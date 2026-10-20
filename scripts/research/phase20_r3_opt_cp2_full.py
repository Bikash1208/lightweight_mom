"""Phase 20 R3 item 3, CP-2 stage 2: Nelder-Mead on the FULL 4-state objective (handedness-preserving r_RR lattice 0/90/180/270 +- 15 deg,
amplitude, opposite-handed <= -15 dB, theta = 0 and 30 deg, tolerance corners inside the objective) from the stage-1 optimum.

    python scripts/research/phase20_r3_opt_cp2_full.py <period_mm> <madp|smp> <out.json> <R> <line> <width> [maxfev] [workers]   (mm)
"""

import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from phase20_r3_opt import corner_kits, lp_responses, lp_score
from phase20_r3_opt_cp2 import BOUNDS, NAMES, make_cell
from scipy.optimize import minimize

from src.ris.cells import madp_kit, smp1345_kit

MESH_N = None  # None: the search pitch (0.25 mm); set from argv for the polish on a finer mesh


def objective(x, period, kits, pool, log):
    lo = np.array([BOUNDS[k][0] * 1e3 for k in NAMES])
    hi = np.array([BOUNDS[k][1] * 1e3 for k in NAMES])
    pen = float(np.sum(np.maximum(0, lo - x) + np.maximum(0, x - hi)))
    xc = np.clip(x, lo, hi)
    resp = lp_responses(make_cell(xc, period), kits, (0.0, 30.0), pool, n=MESH_N, basis="circular")
    val = 60.0 if resp is None else lp_score(resp)[0] + 10 * pen
    log.append((list(map(float, xc)), float(val)))
    return val


def main(period_mm, kit_name, out, x0, maxfev=120, workers=6, mesh_n=None):
    global MESH_N
    MESH_N = mesh_n
    kit = madp_kit() if kit_name == "madp" else smp1345_kit()
    kits = corner_kits(kit)
    period = period_mm * 1e-3
    x0 = np.array(x0, float)
    log = []
    with Pool(workers) as pool:
        sim = np.vstack([x0, x0 + np.diag([0.15, 0.3, 0.12])])
        res = minimize(
            objective,
            x0,
            args=(period, kits, pool, log),
            method="Nelder-Mead",
            options={"maxfev": maxfev, "xatol": 0.01, "fatol": 1e-3, "initial_simplex": sim},
        )
        best = min(log, key=lambda t: t[1])
        resp = lp_responses(make_cell(best[0], period), kits, (0.0, 30.0), pool, n=MESH_N, basis="circular")
        score, metrics = lp_score(resp)
    rec = {
        "period_mm": period_mm,
        "kit": kit.name,
        "names": NAMES,
        "x_mm": best[0],
        "J": best[1],
        "metrics": metrics,
        "nfev": res.nfev,
        "history": log,
    }
    Path(out).write_text(json.dumps(rec, indent=1))
    print(json.dumps({k: v for k, v in rec.items() if k != "history"}), flush=True)


if __name__ == "__main__":
    a = sys.argv
    main(
        float(a[1]),
        a[2],
        a[3],
        [float(v) for v in a[4:7]],
        int(a[7]) if len(a) > 7 else 120,
        int(a[8]) if len(a) > 8 else 6,
        int(a[9]) if len(a) > 9 else None,
    )
