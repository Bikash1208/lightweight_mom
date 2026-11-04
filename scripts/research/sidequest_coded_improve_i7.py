"""I7: joint re-optimisation (state reference phase eps_B, state weight w_B, target scale g, lambda) on top of I1 + I2 + I4 + I5 (I3 taper and I6 edge ring lose gain in the model: not adopted), diodes <= 766.
Selection: best mean gain with SLL <= 0 dB-limit of the baseline + 1 dB; also reports the static-maximal alternative. Output coded_improve_i7.json."""
import itertools
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_coded_improve as I
import sidequest_coded_improve_run as R

kw = dict(R.BASE)
for it in ("I1", "I2"):
    kw = R.ITEMS[it](kw)
kw = R.ITEMS["I4"](kw)
kw = R.ITEMS["I5"](kw)
# I3 (taper) and I6 (edge ring) lose gain in the model and are not adopted
rows = []
best = None
for eb, wb, g in itertools.product(np.radians([0, 45, 90, 135, 180, 225, 270, 315]), (1.0,), (1.5, 2.0, 2.7)):
    r, aux = I.full_run(kw["vs"], kw["pairs"], [0.0, eb], kw["mu"], g, kw["launch"], (1.0, wb), kw["scale"], 766, R.bd, kw["ring"], kw["g_ring"])
    gm = float(np.mean([x["gain_dbi"] for x in r["states"]]))
    sll = max(x["SLL_db"] for x in r["states"])
    row = {"eps_B_deg": float(np.degrees(eb)), "w_B": wb, "g": g, "gain": [round(x["gain_dbi"], 2) for x in r["states"]], "SLL": [round(x["SLL_db"], 1) for x in r["states"]], "AR": [round(x["AR_db"], 2) for x in r["states"]],
           "static_pct": round(100 * r["static_frac"], 1), "diodes": r["diodes"], "lines": r["lines_with_inverters"], "SW_to_beam_pct": round(100 * r["SW_to_beam"], 1), "edge_pct": round(100 * r["edge"], 1), "ap_eff_pct": round(100 * float(r["ap_eff"][0]), 2),
           "XPD": [round(x["XPD_db"], 1) for x in r["states"]]}
    rows.append(row)
    print(row, flush=True)
    ok = sll <= R.LIMIT_SLL
    if best is None or (ok, gm) > (best[0], best[1]):
        best = (ok, gm, row)
out = {"rows": rows, "selected": best[2]}
Path("results/sidequest_hybrid/coded_improve_i7.json").write_text(json.dumps(out, indent=1), encoding="utf-8", newline="\n")
print("SELECTED", best[2])
