"""Item 2 model: transverse resonance of (patch-grid MoM Bloch Y) + (lumped via/diode branch) on a grounded slab.  MODEL_ONLY.
Branch per square cell (Z_sheet = Z_cell for a square cell):  ON  Z = R_on + j w L_b ;  OFF  Z = j w L_b + 1/(j w (C_off + C_ap)),  L_b = kappa mu0 h + L_ser (kappa = 1: dense-lattice mushroom value).
Radiation: first-order leakage of the -1 harmonic at theta0 = 25 deg (formula of sidequest_a0_direct, with the MoM patch Y at kx = k0 sin theta0 instead of the SW sheet value).
Loss: Re F / (dIm F/dbeta) at the real root, F = Y_air + Y_slab(eps (1 - j tand)) + Y_patch + Y_branch."""
import json
import math
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
C0 = 299792458.0
EPS0 = 8.8541878128e-12
MU0 = 4e-7 * math.pi
F0 = 10.5e9
W0 = 2 * math.pi * F0
K0 = W0 / C0
SIN25 = math.sin(math.radians(25.0))

DIODE = {  # label: (R_on ohm, C_off F)
    "MA4AGP907_measured": (2.44, 28.16e-15),
    "MA4AGFCP910_datasheet_typ": (5.2, 18e-15),  # MACOM V4: CT 0.018 pF typ (-10 V), RS 5.2 ohm typ (+10 mA, 10 GHz)
    "MA4AGFCP910_datasheet_max": (6.0, 21e-15),
}


class Patch:
    """MoM patch Y_uu(kr) tables for one (lam, h, gap, mesh, psi) from a ypatch json."""

    def __init__(self, rows):
        rows = sorted(rows, key=lambda r: r["kr"])
        self.kr = np.array([r["kr"] for r in rows])
        y = np.array([complex(*r["Yuu_S"]) for r in rows])
        self.rad = [r for r in rows if r["kr"] < 0.9][0]
        sw = [r for r in rows if r["kr"] >= 0.9]
        self.krs = np.array([r["kr"] for r in sw])
        ys = np.array([complex(*r["Yuu_S"]) for r in sw])
        self.im = CubicSpline(self.krs, ys.imag)
        self.re = CubicSpline(self.krs, ys.real)
        self.y_rad = complex(*self.rad["Yuu_S"])
        self.yuv_max = max(abs(complex(*r["Yuv_S"])) / max(abs(complex(*r["Yuu_S"])), 1e-12) for r in rows)

    def y_sw(self, kr):
        return self.re(kr) + 1j * self.im(kr)


def y_air(kr):
    return 1j * W0 * EPS0 / (K0 * math.sqrt(kr * kr - 1.0))


def y_slab(kr, er, tand, h):
    ec = er * (1 - 1j * tand)
    kz = K0 * np.sqrt(ec - kr * kr + 0j)
    return -1j * (W0 * EPS0 * ec / kz) / np.tan(kz * h)


def branch(state, h, diode, L_ser=0.0, kappa=1.0, C_ap=0.0, dR=1.0, dC=1.0, dL=1.0):
    R, C = DIODE[diode]
    L = (kappa * MU0 * h + L_ser) * dL
    if state == "ON":
        return 1.0 / (R * dR + 1j * W0 * L)
    return 1.0 / (1j * W0 * L + 1.0 / (1j * W0 * (C * dC + C_ap)))


def solve_beta(patch, er, tand, h, yb, lossless=False):
    """beta/k0 and alpha_loss [Np/m] of the TM wave for branch admittance yb."""
    def F(kr):
        return y_air(kr) + y_slab(kr, er, 0.0 if lossless else tand, h) + patch.y_sw(kr) + yb

    lo, hi = 1.0 + 1e-3, 1.48
    g = lambda kr: F(kr).imag
    ks = np.linspace(lo, hi, 600)
    v = np.array([g(k) for k in ks])
    ix = np.where(np.diff(np.sign(v)) != 0)[0]
    if len(ix) == 0:
        return None
    kr = brentq(g, ks[ix[0]], ks[ix[0] + 1], xtol=1e-9)
    d = (F(kr * (1 + 1e-6)).imag - F(kr * (1 - 1e-6)).imag) / (2e-6 * kr * K0)
    alpha = F(kr).real / d if d != 0 else float("nan")
    return kr, abs(alpha)


def tm_power(kr, er, h):
    """TM0 SW power per |E0|^2 at the sheet plane, W/m (same as sidequest_a0_direct)."""
    b = kr * K0
    aa = math.sqrt(b * b - K0 * K0)
    kz = math.sqrt(er * K0 ** 2 - b * b)
    return b * W0 * EPS0 * (1 / (4 * aa ** 3) + er / (2 * kz ** 2 * math.sin(kz * h) ** 2) * (h / 2 + math.sin(2 * kz * h) / (4 * kz)))


def a0_rad(patch, er, h, yb_on, yb_off, kr_mean):
    """first-order radiative alpha [Np/m] for a 1-bit square modulation (50/50) of the branch; E, beta from the mean sheet."""
    ybar = 0.5 * (yb_on + yb_off)
    dy = abs(yb_on - yb_off)
    kx = SIN25
    yup = W0 * EPS0 / (K0 * math.sqrt(1 - kx * kx))
    ytot = patch.y_rad + ybar + yup + y_slab(kx, er, 0.0, h)
    prad = 0.5 * yup * (dy / math.pi) ** 2 / abs(ytot) ** 2  # (2 h0/pi)^2, h0 = dy/2
    return prad / (2 * tm_power(kr_mean, er, h)), abs(ytot) * 1e3


def evaluate(patch, lam, er, tand, h, diode="MA4AGP907_measured", **kw):
    yon, yoff = branch("ON", h, diode, **kw), branch("OFF", h, diode, **kw)
    ymean = 0.5 * (yon + yoff)
    on, off, mean = (solve_beta(patch, er, tand, h, y) for y in (yon, yoff, ymean))
    if None in (on, off, mean):
        return {"ok": False, "ON": on, "OFF": off, "mean": mean}
    a0, ytot = a0_rad(patch, er, h, yon, yoff, mean[0])
    return {"ok": True, "beta_ON": on[0], "beta_OFF": off[0], "beta_mean": mean[0], "alpha_loss_ON": on[1], "alpha_loss_OFF": off[1], "A0_rad": a0, "dY_mS": abs(yon - yoff) * 1e3,
            "Ypatch_mS": float(patch.y_sw(mean[0]).imag) * 1e3, "ytot_rad_mS": ytot, "Y_ON_mS": yon * 1e3, "Y_OFF_mS": yoff * 1e3}


def load_patches(path, mesh=48):
    rows = json.loads(Path(path).read_text())
    out = {}
    for r in rows:
        if r["mesh"] == mesh:
            out.setdefault((r["lam"], r["h"], r["gap"], r["psi"]), []).append(r)
    return {k: Patch(v) for k, v in out.items()}
