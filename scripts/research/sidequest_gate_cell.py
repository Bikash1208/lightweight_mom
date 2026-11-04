"""Side quest step 2: switchable surface-wave GATE cell on the Phase 26 stack. MODEL_ONLY (solver-internal sheet tensor from the MoM).

Stack: RT5880 (eps_r 2.2, tan d 0.0009) 1.575 mm, PEC ground, 6.4 mm lattice, 10.54 GHz. One MA4AGP907 (MEASURED_VENDOR_SPARAM, `src/ris/diodes_measured.py`)
per cell, bridging a gap in a split patch, conduction axis along x (0 deg) or y (90 deg) ONLY (no rotated loads: the Phase 28 rotated-load bug is avoided).
Cell "split patch": two metal rectangles a x h separated by a gap g, the diode footprint (g x wd) sits in the gap (the rest of the gap is empty).
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import json
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.composite_geometry import AxisRect, CompositeGeometry  # noqa: E402
from src.geometry import RectangularLattice, UnitCell  # noqa: E402
from src.hologram.coded_cell_eval import design_loading  # noqa: E402
from src.loads import LumpedLoad  # noqa: E402
from src.materials import Metal, Substrate  # noqa: E402
from src.ris import diodes_measured as dm  # noqa: E402

F0 = 10.54e9
P = 6.4e-3
SUB = Substrate("rt5880", 2.2, 0.0009, 1.575e-3)
STATES = {"OFF": "OFF_0V", "ON": "ON_10mA"}


def diode_load(state, scale, direction, wd, g):
    return dm.load(STATES[state], 0.0, 0.0, g, wd, direction, scale)


def make_cell(a, h, g, wd, state, scale=None, direction="x"):
    """Split patch. direction 'x': halves left/right of a gap g along x (patch a wide along x, h tall); 'y': the same turned by 90 deg."""
    if direction == "x":
        halves = (AxisRect(-(g + (a - g) / 2) / 2 + 0.0 * a, 0.0, (a - g) / 2, h), AxisRect((g + (a - g) / 2) / 2, 0.0, (a - g) / 2, h))
    else:
        halves = (AxisRect(0.0, -(g + (a - g) / 2) / 2, h, (a - g) / 2), AxisRect(0.0, (g + (a - g) / 2) / 2, h, (a - g) / 2))
    bridge = AxisRect(0.0, 0.0, g, wd) if direction == "x" else AxisRect(0.0, 0.0, wd, g)  # the load region is metal carrying the diode impedance
    geo = CompositeGeometry(metal=halves + (bridge,))
    ld = diode_load(state, scale, direction, wd, g)
    return UnitCell(RectangularLattice(P, P), geo, SUB, Metal(), backing="pec", loads=(ld,))


# ------------------------------------------------------------------ sheet model, complex mode, strip transmission (line model)
from scipy.optimize import brentq  # noqa: E402

from src.constants import C0  # noqa: E402,F401
from src.dispersion import sheet_model as sm  # noqa: E402
from src.frequency import angular_frequency, free_space_wavenumber  # noqa: E402
from src.hologram.coded_cell_eval import sheet_tensor  # noqa: E402

K0 = free_space_wavenumber(F0)
W0 = angular_frequency(F0)
EPS0 = 8.8541878128e-12


def _f(beta, y):
    return y + sm._y_air(beta, W0, K0) + sm._y_slab(beta, W0, SUB.eps_complex, K0, SUB.thickness)


def _fg(gam, y):
    beta = np.sqrt(K0**2 + gam**2 + 0j)
    return _f(beta, y), beta


def beta_complex(y, f_unused=None):
    """Complex TM0 constant beta (beta = beta' - j alpha) of a sheet admittance `y` (complex, S): complex Newton in the air decay constant gamma
    (beta^2 = k0^2 + gamma^2, so weakly bound modes beta -> k0 are reached) from several starts. Accepts Re gamma > 0, k0 < Re beta < sqrt(eps_r) k0.
    None = no bound TM0 mode (not hidden)."""
    best = None
    for g0 in K0 * np.array([1e-4, 1e-3, 5e-3, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.8, 1.1]):
        g = g0 + 0j
        ok = False
        for _ in range(80):
            h = 1e-5 * max(abs(g), 1e-3 * K0)
            try:
                d = (_fg(g + h, y)[0] - _fg(g - h, y)[0]) / (2 * h)
                step = _fg(g, y)[0] / d
            except ZeroDivisionError:
                break
            if not np.isfinite(step):
                break
            g = g - 0.7 * step
            if g.real <= 0:
                g = abs(g.real) * 0.5 + 1j * g.imag
            if abs(step) < 1e-11 * max(abs(g), 1.0):
                ok = True
                break
        if not ok:
            continue
        fv, beta = _fg(g, y)
        if beta.imag > 0:
            beta = beta.conjugate()
        if g.real > 0 and K0 < beta.real < np.sqrt(SUB.eps_r) * K0 and abs(fv) < 1e-6:
            if best is None or abs(beta.imag) < abs(best.imag):
                best = beta
    return best


def zm(beta):
    """Modal wave impedance seen from the air side: Z_m = gamma_air / (omega eps0) (reactive for a lossless bound wave)."""
    g = np.sqrt(beta**2 - K0**2 + 0j)
    return g / (W0 * EPS0)


def strip_transmission(y_ref, y_gate, n_rows, pitch=P):
    """|S21|^2 and |S11|^2 (dB) of n_rows gate cells inside the reference sheet (single-mode junction line model; power that the junction
    mismatch would scatter into radiation is NOT removed, so |S21| is an UPPER bound on the transmission = a conservative isolation)."""
    br, bg = beta_complex(y_ref), beta_complex(y_gate)
    if br is None or bg is None:
        return None
    zr, zg = zm(br), zm(bg)
    th = bg * n_rows * pitch
    # ABCD of the gate section (the characteristic impedance zg, propagation constant -j*beta for exp(-j beta x))
    a = np.cos(th)
    b = 1j * zg * np.sin(th)
    c = 1j * np.sin(th) / zg
    d = np.cos(th)
    den = a + b / zr + c * zr + d
    s21 = 2.0 / den
    s11 = (a + b / zr - c * zr - d) / den
    return {"beta_ref_k0": float(br.real / K0), "beta_gate_k0": float(bg.real / K0), "alpha_gate_np_m": float(-bg.imag),
            "S21_db": float(20 * np.log10(abs(s21))), "S11_db": float(20 * np.log10(max(abs(s11), 1e-12)))}


def make_dogbone(e, s, h, g, state, scale=None, direction="x", lt=6.0e-3, wd=None):
    """Dog-bone ('I') cell: two end plates e x h at the cell ends (total length lt along the conduction axis), a centre strip of width s,
    the diode (footprint g x wd, wd = s) in the middle of the strip. Neighbouring cells are lt-to-P apart (gap P - lt)."""
    wd = s if wd is None else wd
    for nm, v in (("e", e), ("s", s), ("h", h), ("g", g), ("lt", lt), ("wd", wd)):  # edges must land on mesh lines (0.2 mm at mesh 32, 0.1 mm at mesh 64)
        if abs(v / 1e-4 - round(v / 1e-4)) > 1e-6:
            raise ValueError(f"{nm}={v*1e3:.4f} mm is not a multiple of 0.1 mm")
    xp = (lt - e) / 2
    xs = (lt / 2 - e)  # strip from -xs..+xs
    if direction == "x":
        rects = (AxisRect(-xp, 0.0, e, h), AxisRect(xp, 0.0, e, h), AxisRect(-(xs + g / 2) / 2 - 0.0, 0.0, xs - g / 2, s), AxisRect((xs + g / 2) / 2, 0.0, xs - g / 2, s), AxisRect(0.0, 0.0, g, wd))
    else:
        rects = (AxisRect(0.0, -xp, h, e), AxisRect(0.0, xp, h, e), AxisRect(0.0, -(xs + g / 2) / 2, s, xs - g / 2), AxisRect(0.0, (xs + g / 2) / 2, s, xs - g / 2), AxisRect(0.0, 0.0, wd, g))
    geo = CompositeGeometry(metal=rects)
    ld = dm.load(STATES[state], 0.0, 0.0, g, wd, direction, scale)
    return UnitCell(RectangularLattice(P, P), geo, SUB, Metal(), backing="pec", loads=(ld,))


def y_xx(cell, mesh):
    y = sheet_tensor(cell, F0, mesh)
    return y, complex(y[0, 0])


class at_freq:
    """Context manager: evaluate the line model at another frequency (module globals K0, W0, F0 swapped, restored on exit)."""

    def __init__(self, f_hz):
        self.f = f_hz

    def __enter__(self):
        g = globals()
        self.old = (g["F0"], g["K0"], g["W0"])
        g["F0"], g["K0"], g["W0"] = self.f, free_space_wavenumber(self.f), angular_frequency(self.f)

    def __exit__(self, *a):
        g = globals()
        g["F0"], g["K0"], g["W0"] = self.old


def no_mode_band(freqs_hz, ys):
    """Frequencies (Hz) at which the sheet admittance has no bound TM0 mode, and the contiguous interval around the middle of them."""
    flags = []
    for f, y in zip(freqs_hz, ys, strict=True):
        with at_freq(f):
            flags.append(beta_complex(y) is None)
    return np.array(flags)


def roots_all(y, re_max=1.7):
    """All decaying (Re gamma > 0) complex roots beta = beta' - j alpha of F(beta) = y + Y_air + Y_slab with 1 < beta' < re_max k0 and 0 <= alpha <= 0.6 k0
    (multi-start least squares). Includes the weakly bound TM0-like mode AND the damped parallel-plate-like (quasi-TEM) mode of a conductive sheet over the
    ground plane, which the single-root `beta_complex` can miss or reject (its upper limit sqrt(eps_r) k0 is wrong for a lossy sheet)."""
    from scipy.optimize import root

    out = []
    for br in np.concatenate([[1.00005, 1.0005, 1.005], np.linspace(1.02, re_max, 36)]):
        for bi in (0.0, -0.002, -0.01, -0.04, -0.15, -0.4):
            def fun(v):
                z = _f((v[0] - 1j * v[1]) * K0, y)
                return [z.real, z.imag]
            try:
                s = root(fun, [br, -bi], method="hybr", tol=1e-13)
            except Exception:  # noqa: BLE001
                continue
            if not s.success or abs(fun(s.x)[0]) + abs(fun(s.x)[1]) > 1e-9:
                continue
            b = (s.x[0] - 1j * s.x[1]) * K0
            if b.imag > 1e-12:
                b = b.conjugate()
            g = np.sqrt(b**2 - K0**2 + 0j)
            if g.real > 1e-6 * K0 and K0 * (1 + 1e-9) < b.real < re_max * K0 and -0.6 * K0 <= b.imag <= 1e-9:
                if not any(abs(b - o) < 1e-6 * K0 for o in out):
                    out.append(b)
    return out


def strip_transmission_modes(y_ref, y_gate, n_rows, pitch=P):
    """Single-mode line model for EVERY root of the gate sheet; returns the list and the root of the largest |S21| (conservative: transmission is carried by the least
    attenuated mode; mode coupling is not modelled)."""
    br = beta_complex(y_ref)
    zr = zm(br)
    res = []
    for bg in roots_all(y_gate):
        zg = zm(bg)
        th = bg * n_rows * pitch
        a, b, c, d = np.cos(th), 1j * zg * np.sin(th), 1j * np.sin(th) / zg, np.cos(th)
        den = a + b / zr + c * zr + d
        res.append({"beta_k0": [float(bg.real / K0), float(-bg.imag / K0)], "alpha_np_m": float(-bg.imag), "S21_db": float(20 * np.log10(abs(2 / den))),
                    "S11_db": float(20 * np.log10(max(abs((a + b / zr - c * zr - d) / den), 1e-12)))})
    return res
