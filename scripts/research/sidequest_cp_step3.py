"""Step 3 (model only): stronger switched GB cell (larger disc radius = free geometry of the same crossed-slit family, MoM pairs from k2_family_scan_raw*.json), both sign hypotheses.  MODEL_ONLY.
    python sidequest_cp_step3.py"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_cp_hybrid as H

RES = ROOT / "results/sidequest_hybrid"
PAIRS = {3.0: (0.011 + 8.69j, 0.061 + 6.062j), 3.1: (0.019 + 12.118j, 0.088 + 7.52j), 3.12: (0.02 + 12.725j, 0.093 + 7.821j), 3.15: (0.036 + 16.861j, 0.116 + 8.549j)}  # mS


def job(a):
    r, rel, eta, ncell, lb, mf = a
    import warnings
    warnings.simplefilter("ignore")
    yo, yn = PAIRS[r]
    H.GB_PAIR = (yo * 1e-3, yn * 1e-3)
    jo = H.joint_opt2(rel, ncell=ncell, mu=4.0, loss_budget_db=lb, eta_design=eta, xpd_min_db=16.0, max_flips=mf)
    rr, _ = H.run(rel, gb_mask_override=jo["gb_mask"], flip=jo["flip"], cache=jo["cache"], verbose=False)
    sm = H.summarize_run(rr)
    sm.update({"r_disc": r, "rel_sign": rel, "eta_design": eta, "ncell": int(jo["gb_mask"].sum()), "diodes": int(2 * jo["gb_mask"].sum()), "loss_budget_db": lb, "max_flips": mf, "flips": int(np.sum(jo["flip"] < 0)),
               "dY_step_mS": float(abs(yo.imag - yn.imag)), "gb_sites_changed": int(np.sum(jo["gb_mask"] & (jo["grp"] != "GB"))), "eta_conv": rr["GB_OFF"]["eta_conv"],
               "b1_xpd_off": rr["GB_OFF"]["beam1"]["xpd_db"], "ar_worst_on": rr["GB_ON"]["beam1"]["ar_worst_in_3dB_beam_db"]})
    return sm


if __name__ == "__main__":
    jobs = []
    for rel in (-1.0, +1.0):
        for r in (3.12, 3.15):
            for ncell in (180, 270):
                for mf, lb in ((160, 3.0), (320, 4.0)):
                    jobs.append((r, rel, 0.9, ncell, lb, mf))
    with ProcessPoolExecutor(10) as ex:
        res = list(ex.map(job, jobs))
    (RES / "cp_step3_free.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline="\n")
    for s in res:
        print("rel %+d r %.2f (dY %.1f mS) eta_d %.1f cells %d diodes %d | b1 OFF/ON %.1f/%.1f AR %.2f/%.2f worst %.2f | b2 OFF/ON %.1f/%.1f contrast %.1f balance %.1f | SLL %.1f/%.1f | flips %d" % (
            s["rel_sign"], s["r_disc"], s["dY_step_mS"], s["eta_design"], s["ncell"], s["diodes"], s["b1_off_dbi"], s["b1_on_dbi"], s["ar_off"], s["ar_on"], s["ar_worst_3dB_beam_off"], s["b2_off_dbi"], s["b2_on_dbi"], s["contrast_db"], s["balance_b1_minus_b2_on_db"],
            s["sll_off_with_offset"], s["sll_on_with_offset"], s["flips"]), flush=True)
