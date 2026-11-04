"""Predictions for the four ON variants (calibrated M2 + M1 spread). MODEL_ONLY.  Writes results/sidequest_hybrid/mushroom_variant_predictions.json"""
import json
import math
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np

import sidequest_mushroom_calibrate as cal
import sidequest_mushroom_cell as mc
import sidequest_mushroom_m2 as m2

KAPPA = 4.772  # calibrated on HFSS ON f(180) = 4.92 GHz (sidequest_mushroom_m2.py)
VARS = {"V1": dict(w=2.4, gd=0.4, r=0.2), "V2": dict(w=2.8, gd=0.4, r=0.2), "V3": dict(w=3.2, gd=0.4, r=0.2), "V4": dict(w=2.8, gd=0.8, r=0.5)}
H_MM = 1.575


def job(item):
    k, v = item
    return k, mc.sheet_table(v["w"] * 1e-3, m2.FREQS, 32, v["gd"] * 1e-3, 0.8e-3)


def l_scale(v, scn):
    """Relative change of the series inductance w.r.t. the HFSS baseline (gd 0.4, r 0.2). S0: none. S1: via self-inductance (r) and strap (gd) share ~16 % + 10 % of L (physical estimate, 0.63 + 0.99 nH of 6 nH)."""
    if scn == "S0":
        return 1.0
    via = 0.63 * (math.log(2 * H_MM / v["r"]) - 0.75) / 2.01
    strap = 0.99 * (v["gd"] / 0.4)
    return (6.0 - 0.63 - 0.99 + via + strap) / 6.0


if __name__ == "__main__":
    with ProcessPoolExecutor(4) as ex:
        tabs = dict(ex.map(job, VARS.items()))
    out = {"kappa": KAPPA, "calibration": "HFSS ON f(180 deg) = 4.92 GHz; HFSS ON flat 4.9-5.1 GHz for 90-180 deg; HFSS OFF 3.86/7.66/10.98/12.54/13.17/13.36 GHz", "variants": {}}
    for k, v in VARS.items():
        row = {"geometry_mm": v | {"pad": 0.8, "patch_to_patch_gap": 6.4 - v["w"]}}
        for scn in ("S0", "S1"):
            s = l_scale(v, scn)
            mM2 = m2.ModelM2(tabs[k], KAPPA * s)
            mM1 = cal.ModelPB(tabs[k], 6.0e-9 * s)
            row[scn] = {"L_scale": round(s, 3),
                        "M2": {"ON_90": m2.lowest(mM2, "ON", 90), "ON_180": m2.lowest(mM2, "ON", 180)},
                        "M1": {"ON_90": cal.surface_band(mM1, "ON", 90)[0] / 1e9 if cal.surface_band(mM1, "ON", 90) else None, "ON_180": m2.lowest(mM1, "ON", 180),
                               "ON_second_band_90": sorted(x / 1e9 for x in cal.surface_band(mM1, "ON", 90))[1:3]}}
        mM2 = m2.ModelM2(tabs[k], KAPPA)
        mM1 = cal.ModelPB(tabs[k], 6.0e-9)
        row["OFF_beta_over_k0_10p54"] = {"M2": [b / mc.free_space_wavenumber(10.54e9) for b in mM2.beta_at("OFF", 10.54e9, None) if b / mc.free_space_wavenumber(10.54e9) < 1.4][:1],
                                         "M1": [b / mc.free_space_wavenumber(10.54e9) for b in mM1.beta_at("OFF", 10.54e9, None) if b / mc.free_space_wavenumber(10.54e9) < 1.4][:1]}
        row["OFF_bands_GHz_M2"] = {ph: [round(x / 1e9, 2) for x in cal.surface_band(mM2, "OFF", ph)][-1:] for ph in (90, 180)}
        out["variants"][k] = row
        print(k, v, {s: {m: {a: (round(b, 2) if isinstance(b, float) else b) for a, b in d.items()} for m, d in row[s].items() if isinstance(d, dict)} for s in ("S0", "S1")}, row["OFF_beta_over_k0_10p54"], flush=True)
    Path("results/sidequest_hybrid/mushroom_variant_predictions.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8", newline="\n")
