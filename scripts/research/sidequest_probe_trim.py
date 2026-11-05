"""Probe trim scan (thin-pin FD model, `sidequest_launcher_fd`): input impedance of the 50 ohm probe vs pin height t and pin radius a for the RO3006 bare-launch geometry (cells from 6.5 mm) and the RT5880 baseline (calibration point: HFSS bare probe S11 -15 dB).  MODEL_ONLY.
Z_in = R + jX of the prescribed-current pin (no coax opening / feed model): the X offset against HFSS is checked on the baseline.  -> results/sidequest_hybrid/k2_probe_trim.json"""
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sidequest_launcher_fd as FD

F = 10.5e9
CFG = {"RO3006_bare_rho6.5": (6.15, 0.002, 1.27e-3, 7.4e-3, 6.5e-3), "RO3006_pad_rho2.6": (6.15, 0.002, 1.27e-3, 7.4e-3, 2.6e-3), "RT5880_baseline_rho6.5": (2.2, 0.0009, 1.575e-3, 7.38e-3, 6.5e-3)}


def job(a):
    name, t, rad = a
    e, td, h, ys, rin = CFG[name]
    r = FD.solve(F, e, td, h, t * 1e-3, ys, rin, a=rad * 1e-3)
    z = complex(r["Z_in"])
    g = (z - 50) / (z + 50)
    return {"cfg": name, "t_mm": t, "a_mm": rad, "R": z.real, "X": z.imag, "S11_db": 20 * math.log10(abs(g)), "launch_rel_FD": r["P_sw"] / r["P_in"]}


if __name__ == "__main__":
    jobs = [(n, t, 0.635) for n in CFG for t in (3, 4, 5, 6, 7, 8)] + [("RO3006_bare_rho6.5", t, 1.0) for t in (5, 6, 7, 8)]
    with ProcessPoolExecutor(len(jobs)) as ex:
        res = list(ex.map(job, jobs))
    Path("results/sidequest_hybrid/k2_probe_trim.json").write_text(json.dumps(res, indent=1), encoding="utf-8", newline="\n")
    for r in res:
        print(r["cfg"], "t %.0f a %.3f: R %.1f X %+.1f S11 %.1f dB" % (r["t_mm"], r["a_mm"], r["R"], r["X"], r["S11_db"]), flush=True)
