"""SUB6 centre-probe launcher at 5.8 GHz on RT5880 3.175 mm (MODEL_ONLY / NOT_HFSS_VALIDATED).  src/dispersion/pin_launch.py used unmodified.

Model (pin_launch.py): thin filament from the ground (z = -h) to z_top above the substrate surface, optional top-hat disk of radius a; TM0 launch fraction
L = P_sw / (P_sw + P_space); R_rad = 2 (P_sw + P_space) (radiated-power resistance), X = thin-wire transmission-line reactance -Z_w cot(Phi) (ESTIMATE, good to
tens of ohms per the Phase 26 note).  Input impedance estimate  Z(f) = R_rad(f) + j X(f)  (dielectric/conductor loss and the SMA flange / coax aperture are NOT in
the model), reflection against 50 ohm  Gamma = (Z - 50)/(Z + 50), evaluated at 5.725, 5.7625, 5.8, 5.8375, 5.875 GHz.  Pin length from the ground = h + z_top.
Sweep: z_top, disk radius a (<= A_MAX, so the hat stays inside the clear radius around the feed), pin radius r (the SMA centre conductor).
Objective: maximise L(5.8 GHz) subject to max_band |S11| <= -10 dB; also reported: accepted-power share L (1 - |Gamma|^2).
    python sub6_launcher.py run [--workers N]
"""
import argparse
import csv
import json
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

from src.dispersion.pin_launch import Launcher, PinLaunch

RES = ROOT / "results" / "sub6_launcher"
H, EPS_R = 3.175e-3, 2.2
FREQS = (5.725e9, 5.7625e9, 5.8e9, 5.8375e9, 5.875e9)
F0 = 5.8e9
A_MAX_MM = 6.0  # first-ring cell metal starts at 11.6 - 5.44 = 6.2 mm from the feed (pitch 11.6, disc r up to 5.44): hat kept inside 6.0 mm
TOPS_MM = [round(x, 3) for x in np.arange(-2.9, 12.01, 0.25)]
DISKS_MM = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
RADII_MM = [0.25, 0.4, 0.635]
OLD = dict(top_mm=1.219 * (H / 1.575e-3) - H * 1e3, disk_mm=0.0, r_mm=0.305 * (H / 1.575e-3))  # pin used so far (Phase 26 pin scaled by h)


def evalpt(a):
    top, disk, r, nseg = a
    rows = []
    for f in FREQS:
        s = PinLaunch(f, EPS_R, H, Launcher(top * 1e-3, r * 1e-3, disk * 1e-3), n_seg=nseg).split()
        z = complex(s["r_rad_ohm"], s["x_est_ohm"])
        g = abs((z - 50.0) / (z + 50.0))
        rows.append((f, s["launch_fraction"], s["r_rad_ohm"], s["x_est_ohm"], g))
    i0 = FREQS.index(F0)
    g_db = [20 * np.log10(max(x[4], 1e-9)) for x in rows]
    return dict(top_mm=top, disk_mm=disk, r_mm=r, pin_len_mm=top + H * 1e3, L_58=rows[i0][1], L_min=min(x[1] for x in rows), R_58=rows[i0][2], X_58=rows[i0][3],
                S11_58_db=g_db[i0], S11_worst_db=max(g_db), share_58=rows[i0][1] * (1 - rows[i0][4] ** 2),
                share_band_min=min(x[1] * (1 - x[4] ** 2) for x in rows), R_band=[round(x[2], 2) for x in rows], X_band=[round(x[3], 2) for x in rows])


def run(workers):
    RES.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    grid = [(t, a, r, 16) for t in TOPS_MM for a in DISKS_MM for r in RADII_MM if not (a > 0 and t < -0.5) and a <= A_MAX_MM]
    with Pool(workers) as p:
        rows = p.map(evalpt, grid, chunksize=4)
        feas = sorted([r for r in rows if r["S11_worst_db"] <= -10.0], key=lambda r: -r["L_58"])
        # refine around the best feasible candidates (top +-0.25 in 0.05 steps, disk +-0.5 in 0.1, n_seg 24)
        cand = feas[:12] if feas else sorted(rows, key=lambda r: r["S11_worst_db"])[:12]
        ref = set()
        for c in cand:
            for dt in np.arange(-0.25, 0.2501, 0.05):
                for da in np.arange(-0.5, 0.5001, 0.25):
                    a = round(c["disk_mm"] + da, 3)
                    if 0 <= a <= A_MAX_MM and not (a > 0 and c["top_mm"] + dt < -0.5) and (c["disk_mm"] > 0 or da >= 0):
                        ref.add((round(c["top_mm"] + dt, 3), a, c["r_mm"], 24))
        rr = p.map(evalpt, sorted(ref), chunksize=2)
        old = evalpt((OLD["top_mm"], 0.0, OLD["r_mm"], 24))
        old_pts = p.map(evalpt, [(OLD["top_mm"], 0.0, OLD["r_mm"], 24)])
    keys = ["top_mm", "disk_mm", "r_mm", "pin_len_mm", "L_58", "L_min", "R_58", "X_58", "S11_58_db", "S11_worst_db", "share_58", "share_band_min"]
    for name, rs in (("sweep_coarse.csv", rows), ("sweep_refined.csv", rr)):
        with open(RES / name, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(keys)
            for r in rs:
                w.writerow([f"{r[k]:.5g}" for k in keys])
    feas_r = sorted([r for r in rr if r["S11_worst_db"] <= -10.0], key=lambda r: -r["L_58"])
    json.dump(dict(old=old, feasible_coarse_top=feas[:10], feasible_refined_top=feas_r[:10], n_grid=len(grid), n_feasible=len(feas), seconds=time.time() - t0, a_max_mm=A_MAX_MM,
                   best_launch_any=sorted(rows + rr, key=lambda r: -r["L_58"])[:5], best_share=sorted(rows + rr, key=lambda r: -r["share_58"])[:5]),
              open(RES / "launcher_summary.json", "w"), indent=1, default=float)
    print("old", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in old.items() if k not in ("R_band", "X_band")})
    print("grid", len(grid), "feasible", len(feas), "refined", len(rr), f"{time.time() - t0:.0f}s")
    for r in (feas_r or feas)[:8]:
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items() if k not in ("R_band", "X_band")})
    print("best L any:", {k: round(v, 3) for k, v in max(rows + rr, key=lambda r: r["L_58"]).items() if not isinstance(v, list)})


def extend(workers):
    """Second sweep: thicker pins (non-standard SMA centre conductors, flagged), taller tops, small hats."""
    t0 = time.time()
    grid = [(t, a, r, 16) for t in [round(x, 3) for x in np.arange(5.0, 16.01, 0.25)] for a in (0.0, 0.5, 1.0, 2.0) for r in (0.635, 0.8, 1.0, 1.27, 1.6)]
    with Pool(workers) as p:
        rows = p.map(evalpt, grid, chunksize=4)
        feas = sorted([r for r in rows if r["S11_worst_db"] <= -10.0], key=lambda r: -r["L_58"])
        cand = (feas[:10] if feas else sorted(rows, key=lambda r: r["S11_worst_db"])[:10])
        ref = sorted({(round(c["top_mm"] + dt, 3), c["disk_mm"], c["r_mm"], 24) for c in cand for dt in np.arange(-0.25, 0.2501, 0.05)})
        rr = p.map(evalpt, ref, chunksize=2)
    keys = ["top_mm", "disk_mm", "r_mm", "pin_len_mm", "L_58", "L_min", "R_58", "X_58", "S11_58_db", "S11_worst_db", "share_58", "share_band_min"]
    for name, rs in (("extend_coarse.csv", rows), ("extend_refined.csv", rr)):
        with open(RES / name, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(keys)
            for r in rs:
                w.writerow([f"{r[k]:.5g}" for k in keys])
    allr = sorted(rows + rr, key=lambda r: r["S11_worst_db"])
    print("extend", len(grid), "feasible", len(feas), f"{time.time() - t0:.0f}s")
    for r in (sorted([r for r in rows + rr if r["S11_worst_db"] <= -10], key=lambda r: -r["L_58"])[:8] or allr[:8]):
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items() if k not in ("R_band", "X_band")})


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["run", "extend"])
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args()
    run(a.workers) if a.stage == "run" else extend(a.workers)
