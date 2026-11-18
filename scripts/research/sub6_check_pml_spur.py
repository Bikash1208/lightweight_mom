"""Compute matched-branch PML spurious attenuation from own v252 unit-cell solves."""

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
R = ROOT / "results" / "sub6_meander_model" / "supercell_hfss"
P_M = 0.0116


def read(cell, top, direction, phase):
    if top == "PEC" and direction == "across" and phase == 85:
        name = f"{cell}_PEC_modes.json"
    else:
        air = "40.035" if top == "PML" else "20.035"
        name = f"{cell}_{top}_{direction}_{phase}_air{air}_modes.json"
    return json.loads((R / name).read_text(encoding="utf-8"))["modes"]


out = {"method": "own 2025.2 unit-cell PEC and PML40 solves; 85/95 degree local slope",
       "formula": "pi*f*(1/Q_PML - 1/Q_PEC)/v_g; v_g=360*p*df/dphase",
       "cells": {}}
for cell in ("ANCHOR", "B2p2_D3p0", "B2p2_D3p6"):
    rows = []
    for direction in ("across", "along"):
        pec85 = [x for x in read(cell, "PEC", direction, 85)
                 if x["f_ghz"] <= 6.5 and x["Q"] < 100000]
        pec95 = [x for x in read(cell, "PEC", direction, 95)
                 if x["f_ghz"] <= 7 and x["Q"] < 100000]
        pml85 = [x for x in read(cell, "PML", direction, 85)
                 if x["f_ghz"] <= 6.5 and x["Q"] < 100000]
        for index, ref in enumerate(pec85):
            # Dispersion rises with phase; preserve order among nearby branches.
            if index >= len(pec95) or index >= len(pml85):
                rows.append({"direction": direction, "branch_index": index + 1,
                             "status": "UNMATCHED"})
                continue
            later, opened = pec95[index], pml85[index]
            slope = (later["f_ghz"] - ref["f_ghz"]) / 10
            vg = 360 * P_M * slope * 1e9
            spur = math.pi * opened["f_ghz"] * 1e9 * (1 / opened["Q"] - 1 / ref["Q"]) / vg
            row = {"direction": direction, "branch_index": index + 1,
                   "pec_85_f_ghz": ref["f_ghz"], "pec_85_Q": ref["Q"],
                   "pec_95_f_ghz": later["f_ghz"],
                   "pml_85_f_ghz": opened["f_ghz"], "pml_85_Q": opened["Q"],
                   "frequency_shift_pct": 100 * (opened["f_ghz"] / ref["f_ghz"] - 1),
                   "vg_m_s": vg, "alpha_spur_np_m": spur,
                   "pass_le_0p1": abs(spur) <= 0.1 and slope > 0}
            rows.append(row)
    out["cells"][cell] = rows

(R / "pml_spur_check.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
print(json.dumps(out, indent=2))
