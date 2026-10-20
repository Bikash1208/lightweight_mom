"""Phase 20 R3 item 3, LP-1 / LP-2, stage 1: many Nelder-Mead starts in parallel (one process each, serial solves) on the NOMINAL theta = 0
4-state objective with a phase-spreading term (the max-deviation objective is flat when the states coincide); stage 2 (corners, theta = 30) is
`phase20_r3_opt_lp.py` / the finalist script. MODEL_ONLY.

    python scripts/research/phase20_r3_opt_lp_multi.py <lp1|lp2> <madp|smp> <out.json> <n_starts> <maxfev> [workers]
"""

import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from phase20_r3_opt import lp_responses, lp_score
from phase20_r3_opt_lp import CELLS, SPACE, infeasible, make_cell
from scipy.optimize import minimize

from src.ris.cells import madp_kit, smp1345_kit


class SerialPool:
    @staticmethod
    def map(f, xs):
        return [f(x) for x in xs]


def spread_term(ph):
    """0 when the 4 phases are >= 70 deg apart pairwise (circularly), growing as states coincide."""
    ph = np.sort(np.mod(np.asarray(ph), 360.0))
    gaps = np.diff(np.append(ph, ph[0] + 360.0))
    return float(max(0.0, 70.0 - gaps.min()) / 70.0)


def objective(x, fam, kit, log):
    space = SPACE[fam]
    lo = np.array([v[0] for v in space.values()])
    hi = np.array([v[1] for v in space.values()])
    xc = np.clip(x, lo, hi)
    pen = float(np.sum(np.maximum(0, lo - x) + np.maximum(0, x - hi)))
    c = make_cell(fam, xc)
    inf = infeasible(fam, c)
    if inf > 0:
        val = 50.0 + 10 * (pen + inf)
    else:
        resp = lp_responses(c, {"nominal": kit}, (0.0,), SerialPool)
        val = 60.0 if resp is None else lp_score(resp)[0] + spread_term(resp[("nominal", 0.0)][0]) + 10 * pen
    log.append((list(map(float, xc)), float(val)))
    return val


def one_start(args):
    fam, kit_name, seed, maxfev = args
    kit = madp_kit() if kit_name == "madp" else smp1345_kit()
    space = SPACE[fam]
    lo = np.array([v[0] for v in space.values()])
    hi = np.array([v[1] for v in space.values()])
    rng = np.random.default_rng(1000 + seed)
    base = np.array([getattr(CELLS[fam](), k) * 1e3 for k in space])
    for _ in range(500):
        x0 = (
            base
            if seed == 0
            else (np.clip(base * rng.uniform(0.7, 1.3, base.size), lo, hi) if seed % 2 else rng.uniform(lo, hi))
        )
        if infeasible(fam, make_cell(fam, x0)) == 0:
            break
    log = []
    sim = np.vstack([x0, x0 + np.diag(0.12 * (hi - lo))])
    res = minimize(
        objective,
        x0,
        args=(fam, kit, log),
        method="Nelder-Mead",
        options={"maxfev": maxfev, "xatol": 0.005, "fatol": 1e-4, "initial_simplex": sim},
    )
    best = min(log, key=lambda t: t[1])
    resp = lp_responses(make_cell(fam, best[0]), {"nominal": kit}, (0.0,), SerialPool)
    _, m = lp_score(resp)
    return {
        "seed": seed,
        "x_mm": best[0],
        "J": best[1],
        "dev0": m["dev0"],
        "amp_db": m["amp_db"],
        "xpol_db": m["xpol_db"],
        "phase": [float(v) for v in resp[("nominal", 0.0)][0]],
        "nfev": res.nfev,
    }


if __name__ == "__main__":
    fam, kit_name, out, n_starts, maxfev = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
    workers = int(sys.argv[6]) if len(sys.argv) > 6 else 6
    runs = []
    with Pool(workers) as pool:
        for r in pool.imap_unordered(one_start, [(fam, kit_name, s, maxfev) for s in range(n_starts)]):
            runs.append(r)
            print(
                json.dumps(
                    {k: (np.round(v, 2).tolist() if k in ("x_mm", "amp_db", "phase") else v) for k, v in r.items()}
                ),
                flush=True,
            )
            Path(out).write_text(
                json.dumps({"family": fam, "kit": kit_name, "names": list(SPACE[fam]), "runs": runs}, indent=1)
            )
