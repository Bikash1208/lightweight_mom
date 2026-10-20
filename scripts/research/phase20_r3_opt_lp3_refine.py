"""Phase 20 R3 item 3: refinement of the best LP-3 stage-1 result: perturbed Nelder-Mead restarts around it.

python scripts/research/phase20_r3_opt_lp3_refine.py <madp|smp> <stage1.json> <out.json> <n_starts> <maxfev> [workers]
"""

import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from phase20_r3_opt_lp3 import SPACE, best_subset, eval_codes, infeasible, make_cell, objective
from scipy.optimize import minimize

from src.ris.cells import madp_kit, smp1345_kit


def refine(args):
    kit_name, x_best, seed, maxfev = args
    kit = madp_kit() if kit_name == "madp" else smp1345_kit()
    lo = np.array([v[0] for v in SPACE.values()])
    hi = np.array([v[1] for v in SPACE.values()])
    rng = np.random.default_rng(3000 + seed)
    x0 = np.array(x_best)
    if seed > 0:
        for _ in range(200):
            xt = np.clip(x0 + rng.normal(0, 0.08, x0.size) * (hi - lo), lo, hi)
            if infeasible(make_cell(xt)) == 0:
                x0 = xt
                break
    log = []
    sim = np.vstack([x0, x0 + np.diag(0.06 * (hi - lo))])
    res = minimize(
        objective,
        x0,
        args=(kit, log),
        method="Nelder-Mead",
        options={"maxfev": maxfev, "xatol": 0.003, "fatol": 1e-4, "initial_simplex": sim},
    )
    best = min(log, key=lambda t: t[1])
    r = eval_codes(make_cell(best[0]), kit)
    score, sub = best_subset(r[0], r[1])
    codes = [make_cell(best[0]).states()[i] for i in sub]
    return {
        "seed": seed,
        "x_mm": best[0],
        "J": best[1],
        "codes": codes,
        "phase": r[0][list(sub)].tolist(),
        "amp_db": (20 * np.log10(r[1][list(sub)])).tolist(),
        "nfev": res.nfev,
    }


if __name__ == "__main__":
    kit_name, src, out, n, maxfev = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
    workers = int(sys.argv[6]) if len(sys.argv) > 6 else 6
    best = min(json.loads(Path(src).read_text())["runs"], key=lambda r: r["J"])
    runs = []
    with Pool(workers) as pool:
        for r in pool.imap_unordered(refine, [(kit_name, best["x_mm"], s, maxfev) for s in range(n)]):
            runs.append(r)
            print(json.dumps(r), flush=True)
            Path(out).write_text(
                json.dumps({"family": "lp3", "kit": kit_name, "names": list(SPACE), "runs": runs}, indent=1)
            )
