"""Side quest: switchable via-to-ground (mushroom-type) gate cell. LC / transverse-resonance estimate. MODEL_ONLY (no vias in the MoM).

Cell (RT5880 1.575 mm, 6.4 mm lattice, PEC ground): square patch w x w with a square hole (pad 0.8 mm + 0.4 mm gap all round); the pad carries the via (r 0.2 mm) to
the ground; the MA4AGP907 bridges the 0.4 mm gap between pad and patch in plane (conduction axis x). ON = diode forward (R) -> via connected; OFF = diode C -> via open.
Patch sheet admittance Y_sheet(f): MoM (normal incidence) of the patch+pad cell with the pad floating (no via, no diode). The via enters analytically:
* wire-medium (dense via) slab: E_z shorted, the grounded slab seen from the sheet is the TEM line  Z_slab = j eta_s tan(k_s h)  (beta-independent), in series with the diode:
    ON : Y_down = 1 / (Z_slab + Z_on(f))
    OFF: Y_down = Y_plain(beta) + 1 / (Z_slab + Z_off(f))     (plain dielectric slab in parallel with the capacitively-loaded via; approximation)
* Y_up: air half space (design: lid_mm = None) or PEC lid at lid_mm (HFSS eigenmode prediction).
Transverse resonance Y_sheet + Y_down + Y_up = 0 -> beta(f) (lossless parts, real root); a Bloch phase phi along x is phi = beta a.
References (bibliographic-only, NOT verified from this session): D. Sievenpiper et al., IEEE Trans. Microwave Theory Tech. 47(11), 2059-2074, 1999, doi 10.1109/22.798001
(L = mu0 h, C = w eps0 (1+eps_r)/pi acosh((2w+g)/g), f0 = 1/(2 pi sqrt(LC)), BW = sqrt(L/C)/eta0).
"""
import json
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np
from scipy.optimize import brentq

from src.composite_geometry import AxisRect, CompositeGeometry
from src.frequency import angular_frequency, free_space_wavenumber
from src.geometry import RectangularLattice, UnitCell
from src.hologram.coded_cell_eval import sheet_tensor
from src.materials import Metal, Substrate
from src.ris import diodes_measured as dm

P, H = 6.4e-3, 1.575e-3
SUB = Substrate("rt5880", 2.2, 0.0009, H)
EPS0, MU0, C0 = 8.8541878128e-12, 4e-7 * np.pi, 299792458.0
PAD, GD = 0.8e-3, 0.4e-3  # pad side, pad-to-patch gap (diode gap)
ETA0 = 376.730313


def cell(w):
    hh = PAD / 2 + GD
    ring = (AxisRect(-(w / 2 + hh) / 2, 0.0, w / 2 - hh, w), AxisRect((w / 2 + hh) / 2, 0.0, w / 2 - hh, w),
            AxisRect(0.0, (w / 2 + hh) / 2, 2 * hh, w / 2 - hh), AxisRect(0.0, -(w / 2 + hh) / 2, 2 * hh, w / 2 - hh), AxisRect(0.0, 0.0, PAD, PAD))
    return UnitCell(RectangularLattice(P, P), CompositeGeometry(metal=ring), SUB, Metal(), backing="pec")


def sheet_table(w, freqs_ghz, mesh=32):
    out = []
    for f in freqs_ghz:
        y = sheet_tensor(cell(w), f * 1e9, mesh)
        out.append([f, complex(y[0, 0]).imag * 1e3, complex(y[1, 1]).imag * 1e3])
    return out


class Model:
    def __init__(self, table):
        t = np.array(table)
        self.f, self.bxx = t[:, 0] * 1e9, t[:, 1] * 1e-3

    def b_sheet(self, f):
        return float(np.interp(f, self.f, self.bxx))

    # --- admittances at frequency f, Bloch constant beta
    @staticmethod
    def _kz(eps, f, beta):
        k0 = free_space_wavenumber(f)
        return np.sqrt(eps * k0**2 - beta**2 + 0j)

    def y_layer(self, eps, thick, f, beta):
        kz = self._kz(eps, f, beta)
        w = angular_frequency(f)
        return -1j * w * EPS0 * eps / kz / np.tan(kz * thick)

    def y_up(self, f, beta, lid):
        if lid is None:
            return 1j * angular_frequency(f) * EPS0 / np.sqrt(beta**2 - free_space_wavenumber(f) ** 2 + 0j)  # decaying half space
        return self.y_layer(1.0, lid, f, beta)

    def z_slab_tem(self, f):
        k_s = free_space_wavenumber(f) * np.sqrt(SUB.eps_r)
        eta_s = ETA0 / np.sqrt(SUB.eps_r)
        return 1j * eta_s * np.tan(k_s * H)

    def y_down(self, state, f, beta):
        if state == "ON":
            return 1.0 / (self.z_slab_tem(f) + dm.z_of_f("ON_10mA")(f))
        return self.y_layer(SUB.eps_r, H, f, beta) + 1.0 / (self.z_slab_tem(f) + dm.z_of_f("OFF_0V")(f))

    def resid(self, state, f, beta, lid):
        return (1j * self.b_sheet(f) + self.y_down(state, f, beta) + self.y_up(f, beta, lid)).imag

    def beta_at(self, state, f, lid, lo=None, hi=None, n=3000):
        """Real roots beta of the transverse resonance at frequency f (largest first)."""
        k0 = free_space_wavenumber(f)
        betas = np.linspace(0.02 * k0, 4.0 * k0, n) if lid is not None else np.linspace(1.0005 * k0, 4.0 * k0, n)
        v = np.array([self.resid(state, f, b, lid) for b in betas])
        idx = np.flatnonzero(np.diff(np.sign(v)) != 0)
        roots = []
        for i in idx:
            if abs(v[i]) < 1e3 and abs(v[i + 1]) < 1e3:  # skip pole sign flips
                roots.append(brentq(lambda b: self.resid(state, f, b, lid), betas[i], betas[i + 1]))
        return sorted(roots, reverse=True)

    def f_of_phase(self, state, phase_deg, lid, fmin=0.5e9, fmax=20e9, n=400):
        """Frequencies at which a TM-like mode has Bloch phase phase_deg along x (beta = phase / a); sign changes of the transverse resonance at fixed beta."""
        beta = np.radians(phase_deg) / P
        fs = np.linspace(fmin, fmax, n)
        out = []
        prev = None
        for f in fs:
            k0 = free_space_wavenumber(f)
            if lid is None and beta <= k0 * 1.0005:
                prev = None
                continue
            v = self.resid(state, f, beta, lid)
            if prev is not None and np.sign(v) != np.sign(prev[1]) and abs(v) < 1e2 and abs(prev[1]) < 1e2:
                out.append(brentq(lambda ff: self.resid(state, ff, beta, lid), prev[0], f))
            prev = (f, v)
        return out

    def f0(self, state="ON"):
        """Via-patch resonance: Y_sheet + Y_down = 0 (beta -> infinity: Y_up -> 0)."""
        fs = np.linspace(6e9, 18e9, 1200)
        v = [(1j * self.b_sheet(f) + self.y_down(state, f, 5 * free_space_wavenumber(f))).imag for f in fs]
        idx = [i for i in range(len(fs) - 1) if np.sign(v[i]) != np.sign(v[i + 1]) and abs(v[i]) < 0.5]
        return [brentq(lambda f: (1j * self.b_sheet(f) + self.y_down(state, f, 5 * free_space_wavenumber(f))).imag, fs[i], fs[i + 1]) for i in idx]


def lc_estimate(w, g):
    L = MU0 * H
    C = w * EPS0 * (1 + SUB.eps_r) / np.pi * np.arccosh((2 * w + g) / g)
    return {"L_nH": L * 1e9, "C_pF": C * 1e12, "f0_GHz": 1 / (2 * np.pi * np.sqrt(L * C)) / 1e9, "BW_frac": np.sqrt(L / C) / ETA0}


if __name__ == "__main__":
    freqs = [4, 6, 8, 9.5, 10.54, 11.5, 12.5, 14, 16]
    res = {}
    for w in (5.2e-3, 5.6e-3, 6.0e-3):
        tab = sheet_table(w, freqs)
        m = Model(tab)
        g = P - w
        res[f"{w*1e3:.1f}"] = {"table": tab, "lc": lc_estimate(w, g), "f0_ON": [x / 1e9 for x in m.f0("ON")], "b_sheet_10p54_mS": m.b_sheet(10.54e9) * 1e3,
                               "f180_ON_halfspace": [x / 1e9 for x in m.f_of_phase("ON", 180, None)], "beta_OFF_10p54": [b / free_space_wavenumber(10.54e9) for b in m.beta_at("OFF", 10.54e9, None)][:2]}
        print(f"w={w*1e3:.1f}", {k: v for k, v in res[f"{w*1e3:.1f}"].items() if k != "table"}, flush=True)
    Path("results/sidequest_hybrid/mushroom_theory_raw.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline="\n")
