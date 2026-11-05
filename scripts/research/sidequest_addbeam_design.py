"""Add-a-beam (b): w2 1.25, eps_B 0, <= 514 diodes, state-dependent beta -> design JSON with the cells and the fields the package script reads.  MODEL_ONLY."""
import json
import math
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sidequest_ro3006_v2 as V

if __name__ == "__main__":
    cfg = {"mode": "add", "decouple": True, "g": 2.0, "mu": 30.0, "eps_B": 0.0, "w2": 1.25, "keep": True, "betas": [1.33, 1.26]}
    e = V.run(("add-a-beam (b)", cfg))
    print(V.summarize(e), flush=True)
    st = []
    for k, s in enumerate(e["states"]):
        b1, b2 = s["beam_gain_dbi"]["b1_(25,0)"], s["beam_gain_dbi"]["b2_(25,180)"]
        st.append({"peak_dir": s["peak_dir"], "gain_design_dir": b1, "gain_peak": max(b1, b2) if k == 1 else b1, "SLL_aperture_db": s["SLL_firstnull_db"], "SLL_floor_db": s["floor_db"], "SLL_total_db": max(s["SLL_firstnull_db"], s["floor_db"]),
                   "other_beam_dbi": b2, "isolation_db": b1 - b2, "XPD_db": float("nan"), "raw": s})
    out = {"cells": e["cells"], "states": st, "A0": e["A0"], "beta_k0": sum(cfg["betas"]) / 2, "cfg": cfg, "lam": e["lam"], "conv": e["conv"], "edge": e["edge"], "NB": e["NB"], "addbeam_states": e["states"], "contrast_db": e["beam2_on_off_contrast_db"], "diodes": e["diodes"]}
    Path("results/sidequest_hybrid/ro3006_board_design_r2.3_addbeam_b.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8", newline=chr(10))
