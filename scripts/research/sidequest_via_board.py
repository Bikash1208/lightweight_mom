"""Board predictions with the via-switched GB cell (MA4AGP907 chip PIN, one diode per cell), sign -1 (HFSS-confirmed), diodes vs balance, with the CP/LP-hybrid HFSS calibration applied.  MODEL_ONLY.
Calibration from the HFSS run of the hybrid package (model 22.4 / 22.2 dBi beam 1, AR 2.1, beam 2 7.7 / -3.5 dBi, balance 14.5, SLL -11.7 with the +2.5 offset  vs  HFSS 20.7 / 20.3, AR 1.2 / 1.0, 9.4 / -2.2, 10.9, -13.2 / -12.7):
beam 1 -1.8 dB, beam 2 +1.5 dB, balance -3.3 dB, AR x0.5, SLL: raw model + 1.0 dB (not +2.5).   python sidequest_via_board.py"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_cp_hybrid as H
import sidequest_via_design as D

RES = ROOT / "results/sidequest_hybrid"
OFF_B1, OFF_B2, OFF_BAL, AR_F, SLL_RAW_OFF = -1.8, +1.5, -3.3, 0.5, +1.0


def pair(C_s_pF, L_ser_nH, gap=0.4, diode="MA4AGP907 (measured)", **kw):
    gaps, yp = D.patch_table()
    yon = D.br("ON", diode, C_s_pF * 1e-12, L_ser_nH * 1e-9, **kw)
    yoff = D.br("OFF", diode, C_s_pF * 1e-12, L_ser_nH * 1e-9, **kw)
    return complex(yon.real, yon.imag + yp[gap]), complex(yoff.real, yoff.imag + yp[gap])


def calibrated(sm):
    return {"b1_off": sm["b1_off_dbi"] + OFF_B1, "b1_on": sm["b1_on_dbi"] + OFF_B1, "b2_off": sm["b2_off_dbi"] + OFF_B2, "b2_on": sm["b2_on_dbi"] + OFF_B2, "contrast": sm["contrast_db"], "balance": sm["balance_b1_minus_b2_on_db"] + OFF_BAL,
            "beam1_change_db": sm["b1_on_dbi"] - sm["b1_off_dbi"], "ar_off": sm["ar_off"] * AR_F, "ar_on": sm["ar_on"] * AR_F, "sll_off": sm["sll_off_with_offset"] - 2.5 + SLL_RAW_OFF, "sll_on": sm["sll_on_with_offset"] - 2.5 + SLL_RAW_OFF}


def job(a):
    C_s, L_ser, ncell, mf, lb, perturb = a
    import warnings
    warnings.simplefilter("ignore")
    yon, yoff = pair(C_s, L_ser)
    H.GB_PAIR = (yoff, yon)
    rel = -1.0
    jo = H.joint_opt2(rel, ncell=ncell, mu=4.0, loss_budget_db=lb, eta_design=0.9, xpd_min_db=16.0, max_flips=mf)
    rr, _ = H.run(rel, gb_mask_override=jo["gb_mask"], flip=jo["flip"], cache=jo["cache"], verbose=False)
    sm = H.summarize_run(rr)
    out = {"C_s_pF": C_s, "L_ser_nH": L_ser, "ncell": int(jo["gb_mask"].sum()), "diodes": int(jo["gb_mask"].sum()), "max_flips": mf, "flips": int(np.sum(jo["flip"] < 0)), "loss_budget_db": lb, "model": sm, "cal": calibrated(sm),
           "ar_worst_off": rr["GB_OFF"]["beam1"]["ar_worst_in_3dB_beam_db"] * AR_F, "step_mS": (yoff - yon).imag * 1e3}
    # tolerance: evaluate the SAME design with perturbed cell levels
    tol = {}
    for nm, kw in perturb.items():
        yon2, yoff2 = pair(C_s, L_ser, **kw)
        H.GB_PAIR = (yoff2, yon2)
        r2, _ = H.run(rel, gb_mask_override=jo["gb_mask"], flip=jo["flip"], cache=jo["cache"], verbose=False)
        tol[nm] = calibrated(H.summarize_run(r2)) | {"step_mS": (yoff2 - yon2).imag * 1e3}
    out["tolerance"] = tol
    return out


if __name__ == "__main__":
    C_s, L_ser = float(sys.argv[1]) if len(sys.argv) > 1 else 0.518, float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    perturb = {"diode C +20%": dict(dC=1.2), "diode C -20%": dict(dC=0.8), "via L +8%": dict(dL=1.08), "via L -8%": dict(dL=0.92), "C_s +10%": dict(dCs=1.1), "C_s -10%": dict(dCs=0.9), "pad C +5fF": dict(C_pad=30e-15), "pad C -5fF": dict(C_pad=20e-15), "R_on +20%": dict(dR=1.2)}
    jobs = [(C_s, L_ser, nc, mf, lb, perturb) for nc in (60, 90, 120, 150, 180, 270, 360) for mf, lb in ((80, 2.0), (160, 3.0))]
    with ProcessPoolExecutor(10) as ex:
        res = list(ex.map(job, jobs))
    (RES / "via_board.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline="\n")
    for r in res:
        c = r["cal"]
        print("cells/diodes %3d flips %3d | step %.1f mS | CAL b1 %.1f/%.1f (change %.1f) AR %.1f/%.1f (worst %.1f) | b2 %.1f/%.1f contrast %.1f balance %.1f | SLL %.1f/%.1f" % (
            r["diodes"], r["flips"], r["step_mS"], c["b1_off"], c["b1_on"], c["beam1_change_db"], c["ar_off"], c["ar_on"], r["ar_worst_off"], c["b2_off"], c["b2_on"], c["contrast"], c["balance"], c["sll_off"], c["sll_on"]), flush=True)
