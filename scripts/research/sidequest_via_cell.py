"""Via-switched GB cell on RT5880 1.575 (patch + via through the ground anti-pad to a back diode, thin bias layer, back-side series L / C knob).  MODEL_ONLY.
Patch grid: MoM Bloch sheet tensor Y_p(kr) (grounded RT5880, no via).  Branch per square cell (lumped, Z_sheet = Z_cell):
    ON  : Z = R_on + j w L_tot - j / (w C_s)
    OFF : Z = j w L_tot - j / (w C_s) - j / (w (C_off + C_pad)),   L_tot = kappa mu0 h + L_ser
Knobs: C_s (series MIM in the bias layer) and L_ser (back trace); C_pad (pad / island capacitance in parallel with the diode) is a geometry estimate (25 fF, from the Route B pad estimate), NOT in the measured diode model.
Targets: mean sheet (average of the two sheet levels) gives beta/k0 = 1.062 (tensor cells), switched step s = Y_OFF - Y_ON.   python sidequest_via_cell.py"""
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from scipy.optimize import brentq, fsolve

import sidequest_cellB_model as M

RES = ROOT / "results/sidequest_hybrid"
ER, TD, H = 2.2, 0.0009, 1.575e-3
F0 = 10.5e9
W0 = 2 * math.pi * F0
R_ON, C_OFF = 2.44, 28.16e-15
C_PAD = 25e-15
MU0 = M.MU0


def branch(state, C_s, L_ser, f=F0, kappa=1.0, C_pad=C_PAD, dR=1.0, dC=1.0, dL=1.0, dCs=1.0, dLs=1.0, h=H):
    w = 2 * math.pi * f
    L = kappa * MU0 * h * dL + L_ser * dLs
    x = w * L - 1.0 / (w * C_s * dCs)
    if state == "ON":
        return 1.0 / (R_ON * dR + 1j * x)
    return 1.0 / (1j * (x - 1.0 / (w * (C_OFF * dC + C_pad))))


def y_air(kr, f=F0):
    w = 2 * math.pi * f
    k0 = w / M.C0
    return 1j * w * M.EPS0 / (k0 * np.sqrt(kr * kr - 1.0 + 0j))


def y_slab(kr, f=F0, er=ER, h=H):
    w = 2 * math.pi * f
    k0 = w / M.C0
    ec = er * (1 - 1j * TD)
    kz = k0 * np.sqrt(ec - kr * kr + 0j)
    return -1j * (w * M.EPS0 * ec / kz) / np.tan(kz * h)


def root(patch, y_br, f=F0, er=ER, h=H, scale=1.0, lo=1.0005, hi=1.44):
    """First (lowest-beta) TM root beta/k0 and alpha_loss for a sheet = scale*Y_p + y_br; None if absent below `hi` (the second, inductive branch starts above the slab pole)."""
    w = 2 * math.pi * f
    k0 = w / M.C0

    def F(kr):
        return y_air(kr, f) + y_slab(kr, f, er, h) + scale * patch.y_sw(min(kr, 1.5)) + y_br

    ks = np.linspace(lo, hi, 900)
    v = np.array([F(k).imag for k in ks])
    ix = np.where(np.diff(np.sign(v)) != 0)[0]
    if len(ix) == 0:
        return None
    kr = brentq(lambda k: F(k).imag, ks[ix[0]], ks[ix[0] + 1], xtol=1e-10)
    d = (F(kr * (1 + 1e-6)).imag - F(kr * (1 - 1e-6)).imag) / (2e-6 * kr * k0)
    return kr, abs(F(kr).real / d)


def levels(patch, C_s, L_ser, **kw):
    yp = patch.y_sw(1.06)
    yon, yoff = branch("ON", C_s, L_ser, **kw), branch("OFF", C_s, L_ser, **kw)
    return yp + yon, yp + yoff, yon, yoff


def solve_knobs(patch, step_mS, beta_mean=1.062, Cs0=1.0e-12, Ls0=0.0):
    """(C_s, L_ser) such that the mean sheet gives beta_mean and the branch step is step_mS."""
    target_mean = None

    def res(x):
        C_s, L_ser = max(abs(x[0]), 0.01) * 1e-12, max(x[1], 0.0) * 1e-9
        yon, yoff = branch("ON", C_s, L_ser), branch("OFF", C_s, L_ser)
        r = root(patch, 0.5 * (yon + yoff), lo=1.0005, hi=1.44)
        bm = 1.0 if r is None else r[0]
        return [(bm - beta_mean) * 100.0, ((yoff - yon).imag * 1e3 - step_mS)]

    sol = fsolve(res, [Cs0 * 1e12, Ls0 * 1e9], xtol=1e-10)
    return max(abs(sol[0]), 0.01) * 1e-12, max(sol[1], 0.0) * 1e-9, res(sol)


def main():
    rows = json.loads((RES / "cellB_rt5880_explore.json").read_text())
    P = M.Patch([r for r in rows if abs(r["gap"] - 0.4) < 1e-6])
    for s in (12.0, 14.0, 15.0, 16.0, 20.0):
        C_s, L_ser, rr = solve_knobs(P, s)
        ys_on, ys_off, yon, yoff = levels(P, C_s, L_ser)
        ro, rf = root(P, yon), root(P, yoff)
        rm = root(P, 0.5 * (yon + yoff))
        print("s %4.1f mS: C_s %.3f pF L_ser %.2f nH | branch Y_ON %.2f%+.2fj Y_OFF %.2f mS | sheet ON %.2f OFF %.2f mS | beta_mean %.3f beta_OFF %s beta_ON %s | resid %s" % (
            s, C_s * 1e12, L_ser * 1e9, yon.real * 1e3, yon.imag * 1e3, yoff.imag * 1e3, ys_on.imag * 1e3, ys_off.imag * 1e3, rm[0] if rm else float("nan"), None if rf is None else round(rf[0], 3), None if ro is None else round(ro[0], 3), np.round(rr, 4)))


if __name__ == "__main__":
    main()
