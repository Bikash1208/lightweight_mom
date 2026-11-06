"""Second board scan: p = 4.5 mm cell (open question a) and a retry of the failed swap R = 150 case.  MODEL_ONLY.   python sidequest_swhb_scan2.py"""
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sidequest_swhb_design as DS
import sidequest_swhb_design2 as D2
import sidequest_swhb_run as RUN


def one(a):
    R, variant, mu, seed, pfile, mesh, gap, p, Lser, tag = a
    rows = RUN.patch_rows(RUN.RES / pfile, mesh, 0.0, gap, p=p * 1e3)
    lib = DS.CellLib(rows, 3.38, 0.0027, 1.524e-3, p, L_ser=Lser)
    t0 = time.time()
    des = D2.design2(lib, R, variant, mu_sl=mu, mu_w=2.0, outer=3, sweeps=6, sq=True, seed=seed)
    if not des["ok"]:
        return {"R": R, "variant": variant, "tag": tag, "ok": False, "why": des.get("why")}
    ms = RUN.metrics(des, DS.evaluate(des))
    return {"R": R, "variant": variant, "tag": tag, "mu": mu, "seed": seed, "p": p, "ok": True, "n_cells": des["n"], "diodes": des["nd"], "time_s": time.time() - t0, "states": ms}


if __name__ == "__main__":
    f64 = "cellB_ypatch_final_m32_m64.json"
    f45 = "cellB_ypatch_p45_m60.json"
    jobs = [(0.15, "swap", 3.0, 1, f64, 64, 0.8, 6.4e-3, 1.6e-9, "p6.4_retry_seed1"), (0.15, "swap", 2.0, 2, f64, 64, 0.8, 6.4e-3, 1.6e-9, "p6.4_retry_seed2_mu2")]
    for R in (0.10, 0.125, 0.15):
        for variant in ("swap", "add"):
            jobs.append((R, variant, 3.0, 0, f45, 60, 0.3, 4.5e-3, 1.6e-9, "p4.5"))
    with ProcessPoolExecutor(8) as ex:
        res = list(ex.map(one, jobs))
    (RUN.RES / "swhb_scan_p45_retry.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline="\n")
    for r in res:
        if not r["ok"]:
            print(r)
            continue
        s = r["states"]
        print("%s R %.0f %s diodes %d/%d | S1 G %.1f SLL %.1f | S2 G %.1f SLL %.1f | edge %.2f/%.2f" % (r["tag"], r["R"] * 1e3, r["variant"], r["diodes"], r["n_cells"], s[0]["G_peak_dBi"], s[0]["SLL_db"], s[1]["G_peak_dBi"], s[1]["SLL_db"], s[0]["P_edge"], s[1]["P_edge"]), flush=True)
