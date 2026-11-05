"""Evaluate the family scan raw files: direct A0, beta/k0, X (lid) OFF/ON for each (stack, r).  MODEL_ONLY.  -> results/sidequest_hybrid/k2_family_scan.json"""
import json, sys, math
sys.path.insert(0, "scripts/research")
import sidequest_a0_direct as A
import sidequest_leaf_cell as L1
from src.materials import Substrate
STK = {"RT5880_1.575": (2.2, 0.0009, 1.575), "RT5880_3.175": (2.2, 0.0009, 3.175), "RO3006_1.27": (6.15, 0.002, 1.27)}
f = 10.5e9
raw = json.load(open("results/sidequest_hybrid/k2_family_scan_raw.json")) + json.load(open("results/sidequest_hybrid/k2_family_scan_raw_fine.json"))
out = []
for n, (e, t, h) in STK.items():
    for r in sorted({x["r"] for x in raw if x["stack"] == n}):
        d = {x["mode"]: x["Y_mS"] for x in raw if x["stack"] == n and x["r"] == r}
        yo, yn = complex(*d["open"]) * 1e-3, complex(*d["ON"]) * 1e-3
        row = {"stack": n, "r": r, "Y_off_mS": round(yo.imag * 1e3, 2), "Y_on_mS": [round(yn.real * 1e3, 2), round(yn.imag * 1e3, 2)], "dY_mS": round(abs(yo - yn) * 1e3, 2)}
        try:
            yb = (yo.imag + yn.imag) / 2
            a0, bk = A.a0_direct(yb, abs(yo - yn), e, h * 1e-3, f, tand=t); row.update(A0=round(float(a0), 3), beta_k0=round(float(bk), 3))
            sub = Substrate(n, e, t, h * 1e-3); H = 3e8 / f / 4
            row.update(X_off=round(L1.xs_from_y(yo, sub, f, H)["X_lid_ohm"], 1), X_on=round(L1.xs_from_y(yn, sub, f, H)["X_lid_ohm"], 1))
        except Exception as ex:
            row["err"] = str(ex)[:40]
        out.append(row)
json.dump(out, open("results/sidequest_hybrid/k2_family_scan.json", "w"), indent=1)
for o in out:
    if o["r"] >= 2.0 and "A0" in o and o["A0"] > 0.1: print(o)
