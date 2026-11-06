"""Board scan for candidate B: diameter x variant x full/partial.  MODEL_ONLY.   python sidequest_swhb_scan.py [tag]  -> results/sidequest_hybrid/swhb_scan_<tag>.json"""
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_swhb_design as DS
import sidequest_swhb_design2 as D2
import sidequest_swhb_run as RUN


def one(a):
    R, variant, full, mu, patch_file, mesh, gap, p, h, Lser, tag, maxd = a
    rows = RUN.patch_rows(RUN.RES / patch_file, mesh, 0.0, gap, p=p * 1e3)
    lib = DS.CellLib(rows, 3.38, 0.0027, h, p, L_ser=Lser)
    t0 = time.time()
    des = D2.design2(lib, R, variant, full=full, mu_sl=mu, mu_w=2.0, outer=3, sweeps=6, sq=True, max_diodes=maxd)
    if not des["ok"]:
        return {"R": R, "variant": variant, "full": full, "ok": False, "why": des.get("why")}
    ms = RUN.metrics(des, DS.evaluate(des))
    return {"R": R, "variant": variant, "full": full, "mu": mu, "ok": True, "n_cells": des["n"], "diodes": des["nd"], "max_diodes": maxd, "time_s": time.time() - t0, "states": ms, "p": p}


if __name__ == "__main__":
    tag = sys.argv[1] if len(sys.argv) > 1 else "g08"
    jobs = []
    for R in (0.10, 0.125, 0.15, 0.175, 0.20):
        for variant in ("swap", "add"):
            jobs.append((R, variant, False, 3.0, "cellB_ypatch_final_m32_m64.json", 64, 0.8, 6.4e-3, 1.524e-3, 1.6e-9, tag, 514))
    for R in (0.15, 0.20):
        for variant in ("swap", "add"):
            jobs.append((R, variant, False, 3.0, "cellB_ypatch_final_m32_m64.json", 64, 0.8, 6.4e-3, 1.524e-3, 1.6e-9, tag + "_x2", 1028))
    with ProcessPoolExecutor(10) as ex:
        res = list(ex.map(one, jobs))
    (RUN.RES / f"swhb_scan_{tag}.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline="\n")
    for r in res:
        if not r["ok"]:
            print(r); continue
        s = r["states"]
        print("R %.0f %s %s diodes %d/%d | S1 G %.1f SLL %.1f | S2 G %.1f SLL %.1f | edge %.2f/%.2f" % (r["R"] * 1e3, r["variant"], "full" if r["full"] else "partial", r["diodes"], r["n_cells"], s[0]["G_peak_dBi"], s[0]["SLL_db"], s[1]["G_peak_dBi"], s[1]["SLL_db"], s[0]["P_edge"], s[1]["P_edge"]), flush=True)
