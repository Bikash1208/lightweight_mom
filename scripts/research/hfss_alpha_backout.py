"""Back out the sheet admittances, the surface-wave constant and alpha_beam of the RT5880 coded cell from HFSS eigenmode results (phase 90 deg per cell, one axis per run).
Input: CSV with header `state,axis,f_GHz[,Q]`, state in OFF_OFF / ON_OFF / OFF_ON / ON_ON (x diode, y diode), axis x|y, f_GHz = Re(eigenfrequency) of the bound (TM0) mode, optional Q = Re f / (2 Im f).
    python hfss_alpha_backout.py results.csv            # your HFSS numbers
    python hfss_alpha_backout.py --selftest <predictions.json>   # round trip of the model predictions (should reproduce the model alpha)
Method (transverse resonance, lossless): beta = phase / a = (pi/2) / 6.4 mm for every run. At the eigenfrequency f the sheet susceptance along the propagation axis is
    B_s = - Im[ Y_air(beta, f) + Y_slab(beta, f) ]      (Y_air = TM admittance of the 30 mm air layer under the PEC lid, Y_slab the grounded RT5880 slab, repo `dispersion.sheet_model`)
    beta/k0 = beta / (2 pi f / c);  implied tensor entry Y_xx (x runs) or Y_yy (y runs) = j B_s.
Per axis dY = B_s(ON on that axis) - B_s(OFF/OFF) (other diode OFF); the mean sheet is the average of all runs' B_s (board mean, isotropic). Then
    beta_m from  Im[ j Y_m + Y_air + Y_slab ] = 0 at f0 = 10.54 GHz,   c1 = |dY| / pi,
    alpha_beam = leakage_alpha_first_order(beta_m, j Y_m, c1, theta = 25 deg, zeta = 0, RT5880, f0)       (first-order leakage; alpha ~ c1^2)
    alpha_abs  = (d beta_m / d B) * G_ON / 2,   G_ON = Re Y of the ON state, here from Q:  G ~ (f / (2 Q)) |dB_tot/df| (lossless-sheet slope; order of magnitude).
Caveat: the sheet values are taken at the HFSS eigenfrequency (~10.5-10.9 GHz), not exactly at f0; kt is fixed at 90 deg/cell (beta/k0 ~ 1.05-1.1), the operating beta."""
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np

from src.dispersion import sheet_model as sm
from src.dispersion.leakage_model import leakage_alpha_first_order
from src.materials import Substrate

SUB = Substrate("rt5880", 2.2, 0.0009, 1.575e-3)
A = 6.4e-3
F0 = 10.54e9
BETA = np.radians(90.0) / A
C0 = 299792458.0


LID = 30e-3  # PEC lid above the patch plane in the HFSS model [m]; 0 = open half-space
EPS0 = 8.8541878128e-12


def y_air(beta, om, k0, lid=LID):
    """TM admittance of the air layer looking up from the sheet: open half-space (lid = 0) or a PEC lid at height `lid`: -j om eps0 / (kz tan(kz lid)), kz = -j kappa (bound mode)."""
    if not lid:
        return sm._y_air(beta, om, k0)
    kz = np.sqrt(k0**2 - beta**2 + 0j)
    if kz.imag > 0:
        kz = -kz
    return -1j * om * EPS0 / (kz * np.tan(kz * lid))


def stack_y(beta, f):
    om, k0 = 2 * np.pi * f, 2 * np.pi * f / C0
    return y_air(beta, om, k0) + sm._y_slab(beta, om, SUB.eps_complex, k0, SUB.thickness)


def implied(rows):
    out = {}
    for r in rows:
        f = r["f_GHz"] * 1e9
        b_s = -stack_y(BETA, f).imag
        g = None
        if r.get("Q"):
            df = 1e6
            slope = abs((-stack_y(BETA, f + df).imag - -stack_y(BETA, f - df).imag) / (2 * df)) if False else abs((stack_y(BETA, f + df).imag - stack_y(BETA, f - df).imag) / (2 * df))
            g = f / (2 * r["Q"]) * slope
        out[(r["state"], r["axis"])] = {"f_GHz": r["f_GHz"], "beta_over_k0": BETA / (2 * np.pi * f / C0), "B_s_mS": b_s * 1e3, "G_s_mS": None if g is None else g * 1e3}
    return out


def alpha_from(imp):
    def b(state, axis):
        return imp[(state, axis)]["B_s_mS"] * 1e-3

    dx = b("ON_OFF", "x") - b("OFF_OFF", "x")
    dy = b("OFF_ON", "y") - b("OFF_OFF", "y")
    ym = float(np.mean([v["B_s_mS"] for v in imp.values()])) * 1e-3
    res = {"dY_x_mS": dx * 1e3, "dY_y_mS": dy * 1e3, "Y_mean_mS": ym * 1e3}
    if ("ON_ON", "x") in imp:
        res["ON_ON_additivity_mS"] = (b("ON_ON", "x") - b("OFF_ON", "x") - dx) * 1e3  # extra change of Y_xx when the y diode is also ON
    bm = sm.find_beta_from_sheet_admittance(1j * ym, SUB, F0).beta_rad_per_m  # operating beta of the mean sheet (open half-space, as in the MoM leakage model)
    k0 = 2 * np.pi * F0 / C0
    db = 1e-6
    dbdb = (sm.find_beta_from_sheet_admittance(1j * (ym + db), SUB, F0).beta_rad_per_m - bm) / db
    for nm, d in (("x", dx), ("y", dy)):
        res[f"alpha_beam_{nm}"] = leakage_alpha_first_order(bm, 1j * ym, abs(d) / np.pi, np.radians(25.0), 0.0, SUB, F0)
    res["beta_m_over_k0"] = bm / k0
    gs = [v["G_s_mS"] for k, v in imp.items() if k[0] != "OFF_OFF" and v["G_s_mS"] is not None]
    res["alpha_abs_est"] = None if not gs else dbdb * (np.mean(gs) * 1e-3) / 2
    return res


def read(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return [{"state": r["state"].replace("x", "").replace("y", "").strip("_") if r["state"].startswith("x") else r["state"], "axis": r["axis"], "f_GHz": float(r["f_GHz"]), "Q": float(r["Q"]) if r.get("Q") else None} for r in csv.DictReader(fh)]


if __name__ == "__main__":
    if sys.argv[1] == "--selftest":
        pred = json.loads(Path(sys.argv[2]).read_text(encoding="utf-8"))
        rows = [{"state": "_".join(r["state"].split("/")), "axis": r["axis"], "f_GHz": r["f_GHz"], "Q": None} for r in pred["rows"] if r["variant"] == "A_rim" and r["f_GHz"]]
    else:
        rows = read(sys.argv[1])
    imp = implied(rows)
    for k, v in imp.items():
        print(k, {a: (round(b, 4) if isinstance(b, float) else b) for a, b in v.items()})
    print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in alpha_from(imp).items()}, indent=1))
