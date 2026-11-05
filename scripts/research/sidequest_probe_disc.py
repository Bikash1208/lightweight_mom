"""Probe with a printed top disc (no extra part): pin through the slab ending at the top layer (t = 0) or protruding t mm, PEC disc of radius R_d on the top layer at z = 0, cells from 6.5 mm (bare launch geometry).  Thin-pin FD model; MODEL_ONLY.  -> results/sidequest_hybrid/k2_probe_disc.json"""
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sidequest_launcher_fd as FD
F = 10.5e9
E, TD, H, YS, RIN = 6.15, 0.002, 1.27e-3, 7.4e-3, 6.5e-3


def job(a):
    rd, t = a
    r = FD.solve(F, E, TD, H, t * 1e-3, YS, RIN, disc_r=rd * 1e-3)
    z = complex(r["Z_in"]); g = (z - 50) / (z + 50)
    return {"R_d_mm": rd, "t_mm": t, "R": z.real, "X": z.imag, "S11_db": 20 * math.log10(abs(g)), "launch_rel": r["P_sw"] / r["P_in"]}


if __name__ == "__main__":
    jobs = [(rd, t) for rd in (0.0, 1.5, 2.5, 3.5) for t in (3.0, 4.0, 5.0, 6.0, 7.0)]  # t < 3 mm is not usable in the FD (negative R: the thin-pin FD does not resolve short pins)
    with ProcessPoolExecutor(len(jobs)) as ex:
        res = list(ex.map(job, jobs))
    Path("results/sidequest_hybrid/k2_probe_disc.json").write_text(json.dumps(res, indent=1), encoding="utf-8", newline="\n")
    for r in res:
        print("R_d %.1f t %.0f: R %.1f X %+.1f S11 %.1f dB launch_rel_FD %.2f" % (r["R_d_mm"], r["t_mm"], r["R"], r["X"], r["S11_db"], r["launch_rel"]), flush=True)
