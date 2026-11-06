"""Calibrated board predictions with the real-hardware via cell (resistive bias feed, cathode via), sign -1, one diode per cell.  MODEL_ONLY.   python sidequest_via_hwboard.py"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_cp_hybrid as H
import sidequest_via_board as B
import sidequest_via_hw as W
import sidequest_via_cell as V

RES = ROOT / "results/sidequest_hybrid"
RB = dict(W.RB, R_b=4700.0, C_p=8e-15, C_A=10e-15, cathode="via")
C_S = 0.216e-12


def pair(**kw):
    P = W.patch()
    yp = P.y_sw(1.062).imag
    yon, yoff = W.branch2("ON", V.F0, C_S, RB, **kw), W.branch2("OFF", V.F0, C_S, RB, **kw)
    return complex(yon.real, yon.imag + yp), complex(yoff.real, yoff.imag + yp)


def job(a):
    ncell, mf, lb = a
    import warnings
    warnings.simplefilter("ignore")
    yon, yoff = pair()
    H.GB_PAIR = (yoff, yon)
    rel = -1.0
    jo = H.joint_opt2(rel, ncell=ncell, mu=4.0, loss_budget_db=lb, eta_design=0.9, xpd_min_db=16.0, max_flips=mf)
    rr, _ = H.run(rel, gb_mask_override=jo["gb_mask"], flip=jo["flip"], cache=jo["cache"], verbose=False)
    sm = H.summarize_run(rr)
    out = {"ncell": int(jo["gb_mask"].sum()), "flips": int(np.sum(jo["flip"] < 0)), "step_mS": (yoff - yon).imag * 1e3, "model": sm, "cal": B.calibrated(sm), "ar_worst_off": rr["GB_OFF"]["beam1"]["ar_worst_in_3dB_beam_db"] * B.AR_F}
    tol = {}
    for nm, kw in (("diode C +20%", dict(dC=1.2)), ("diode C -20%", dict(dC=0.8)), ("drill via L +8%", dict(dL=1.08)), ("drill via L -8%", dict(dL=0.92)), ("C_s +20%", dict(dCs=1.2)), ("C_s -20%", dict(dCs=0.8)),
                   ("island C_A +30%", dict(dCA=1.3)), ("island C_A -30%", dict(dCA=0.7)), ("resistor C_p +30%", dict(dCp=1.3)), ("resistor C_p -30%", dict(dCp=0.7)), ("via-land C_V +30%", dict(dCV=1.3)), ("via-land C_V -30%", dict(dCV=0.7))):
        a2, b2 = pair(**kw)
        H.GB_PAIR = (b2, a2)
        r2, _ = H.run(rel, gb_mask_override=jo["gb_mask"], flip=jo["flip"], cache=jo["cache"], verbose=False)
        tol[nm] = B.calibrated(H.summarize_run(r2)) | {"step_mS": (b2 - a2).imag * 1e3}
    out["tolerance"] = tol
    return out


if __name__ == "__main__":
    jobs = [(nc, mf, lb) for nc in (90, 150, 180) for mf, lb in ((80, 2.0),)]
    with ProcessPoolExecutor(3) as ex:
        res = list(ex.map(job, jobs))
    (RES / "via_hwboard.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline="\n")
    for r in res:
        c = r["cal"]
        bal = [v["balance"] for v in r["tolerance"].values()]
        print("HW cell: diodes %3d flips %3d step %.1f | CAL b1 %.1f/%.1f change %.1f AR %.1f/%.1f (worst %.1f) | b2 %.1f/%.1f contrast %.1f balance %.1f | SLL %.1f/%.1f | tolerance balance %.1f..%.1f, SLL ON worst %.1f" % (
            r["ncell"], r["flips"], r["step_mS"], c["b1_off"], c["b1_on"], c["beam1_change_db"], c["ar_off"], c["ar_on"], r["ar_worst_off"], c["b2_off"], c["b2_on"], c["contrast"], c["balance"], c["sll_off"], c["sll_on"], min(bal), max(bal), max(v["sll_on"] for v in r["tolerance"].values())), flush=True)
