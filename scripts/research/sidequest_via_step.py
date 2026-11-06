"""Via-switched GB cell on RT5880: (1) board-level required step (levels M +- s/2 about the tensor mean sheet), both signs, 180 / 360 cells.  MODEL_ONLY.   python sidequest_via_step.py"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_a0_direct as A0
import sidequest_cp_hybrid as H

RES = ROOT / "results/sidequest_hybrid"
M_MEAN = A0.ybar_for_beta(1.062, 2.2, 1.575e-3, 10.54e9)  # S: mean sheet susceptance of the tensor cells


def job(a):
    s_mS, rel, ncell, mf, lb = a
    import warnings
    warnings.simplefilter("ignore")
    s = s_mS * 1e-3
    H.GB_PAIR = (0.0002 + 1j * (M_MEAN + s / 2), 0.0004 + 1j * (M_MEAN - s / 2))
    jo = H.joint_opt2(rel, ncell=ncell, mu=4.0, loss_budget_db=lb, eta_design=0.9, xpd_min_db=16.0, max_flips=mf)
    rr, _ = H.run(rel, gb_mask_override=jo["gb_mask"], flip=jo["flip"], cache=jo["cache"], verbose=False)
    sm = H.summarize_run(rr)
    sm.update({"step_mS": s_mS, "rel_sign": rel, "ncell": int(jo["gb_mask"].sum()), "diodes_1_per_cell": int(jo["gb_mask"].sum()), "max_flips": mf, "loss_budget_db": lb, "flips": int(np.sum(jo["flip"] < 0)),
               "ar_worst_on": rr["GB_ON"]["beam1"]["ar_worst_in_3dB_beam_db"], "eta_conv": rr["GB_OFF"]["eta_conv"]})
    return sm


if __name__ == "__main__":
    print("tensor mean sheet susceptance %.2f mS" % (M_MEAN * 1e3), flush=True)
    jobs = [(s, rel, nc, mf, lb) for s in (10.0, 15.0, 20.0, 25.0) for rel in (-1.0, +1.0) for nc in (180, 360) for mf, lb in ((80, 2.0), (160, 3.0))]
    with ProcessPoolExecutor(10) as ex:
        res = list(ex.map(job, jobs))
    (RES / "via_step_required.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline="\n")
    for r in res:
        print("step %4.1f mS rel %+d cells %3d flips %3d | b1 OFF/ON %.1f/%.1f AR %.2f/%.2f worst %.2f | b2 OFF/ON %.1f/%.1f contrast %.1f balance %.1f | SLL %.1f/%.1f" % (
            r["step_mS"], r["rel_sign"], r["ncell"], r["flips"], r["b1_off_dbi"], r["b1_on_dbi"], r["ar_off"], r["ar_on"], r["ar_worst_3dB_beam_off"], r["b2_off_dbi"], r["b2_on_dbi"], r["contrast_db"], r["balance_b1_minus_b2_on_db"], r["sll_off_with_offset"], r["sll_on_with_offset"]), flush=True)
