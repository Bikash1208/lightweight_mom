"""Real beta/k0 of the r 2.30 board from the board-mean cell susceptance: MoM raw, and with the HFSS offsets found at r 2.20 (OFF open-slit branch lower by 0.34-0.85 mS, ON by 0.10 mS).  MODEL_ONLY.  -> results/sidequest_hybrid/k2_beta_real_r2.3.json"""
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import sidequest_a0_direct as A

ER, TD, H, F = 6.15, 0.002, 1.27, 10.5e9
lib = json.loads(Path("results/sidequest_hybrid/ro3006_r2.3_mesh96_lib.json").read_text())["rows"]
cells = json.loads(Path("results/sidequest_hybrid/ro3006_board_design_r2.3.json").read_text())["cells"]


def B(kind, fr):
    return next(x for x in lib if x["kind"] == kind and abs(x["fr"] - fr) < 1e-9)["Y_mS"][1]


def cell_B(c, state):
    ix, iy, kd, tk, fr = c
    if kd == "SW":
        return (B("open", fr), "open") if state == 1 else (B("ON", fr), "ON")
    return (B(tk if tk != "none" else "none", fr if tk != "none" else 0.0), tk)


out = {"scenarios": []}
for d_open, d_on, label in ((0.0, 0.0, "MoM raw"), (0.34, 0.10, "HFSS offset, mesh-80 size (open -0.34 mS, ON -0.10)"), (0.82, 0.10, "HFSS offset, mesh-96 size (open -0.82 mS, ON -0.10)"), (0.85, 0.10, "HFSS offset, mesh-64 size (open -0.85 mS, ON -0.10)")):
    row = {"label": label}
    for state in (1, 2):
        bs = []
        for c in cells:
            b, tk = cell_B(c, state)
            if tk == "ON":
                b -= d_on
            elif tk == "open":
                b -= d_open
            bs.append(b)
        ybar = float(np.mean(bs))
        beta = A.beta_of(ybar * 1e-3, ER, TD, H * 1e-3, F)
        k0 = 2 * np.pi * F / 299792458.0
        row[f"state{state}"] = {"Ybar_mS": ybar, "beta_over_k0": float(beta / k0)}
    out["scenarios"].append(row)
    print(label, row["state1"], row["state2"], flush=True)
Path("results/sidequest_hybrid/k2_beta_real_r2.3.json").write_text(json.dumps(out, indent=1), encoding="utf-8", newline=chr(10))
