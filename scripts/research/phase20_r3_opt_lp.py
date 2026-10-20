"""Phase 20 R3 item 3, LP-1 / LP-2: Nelder-Mead over the cell dimensions against the 4-state objective (phase lattice 0/90/180/270 +- 15 deg,
amplitude >= -1.5 dB with spread <= 1 dB, cross-pol <= -15 dB, theta = 0 and 30 deg) with the diode tolerance corners (C_T, L_s one at a time,
step error <= 20 deg) INSIDE the objective. MODEL_ONLY.

    python scripts/research/phase20_r3_opt_lp.py <lp1|lp2> <madp|smp> <out.json> [n_starts] [maxfev]
"""

import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from phase20_r3_opt import corner_kits, lp_responses, lp_score
from scipy.optimize import minimize

from src.ris.cells import Lp1Cell, Lp2Cell, madp_kit, smp1345_kit

# name: (low, high) in mm
SPACE = {
    "lp1": {
        "left": (1.0, 5.0),
        "middle": (1.0, 5.0),
        "right": (0.5, 4.0),
        "width": (2.0, 9.0),
        "strip_width1": (0.4, 3.0),
        "strip_width2": (0.4, 3.0),
        "strip_y1": (-2.5, 2.5),
        "strip_y2": (-2.5, 2.5),
        "period": (10.0, 15.0),
    },
    "lp2": {
        "patch_length": (3.0, 8.0),
        "patch_width": (2.0, 9.0),
        "stub_width": (0.4, 3.0),
        "stub1_length": (0.25, 3.0),
        "stub2_length": (0.25, 3.0),
        "stub1_y": (-2.0, 2.0),
        "stub2_y": (-2.0, 2.0),
        "period": (10.0, 15.0),
    },
}
CELLS = {"lp1": Lp1Cell, "lp2": Lp2Cell}
PERIOD = 12.0e-3


def make_cell(fam, x):
    kw = {k: float(v) * 1e-3 for k, v in zip(SPACE[fam], x, strict=True)}
    return CELLS[fam](**kw)


def infeasible(fam, c):
    """Geometric violation in mm (0 = feasible): everything inside the cell with 0.25 mm clearance to the boundary; strips inside the patch."""
    m = 0.25e-3
    if fam == "lp1":
        total = c.left + c.gap + c.middle + c.gap + c.right
        v = max(0.0, total - (c.period - 2 * m)) + max(0.0, c.width - (c.period - 2 * m))
        v += max(0.0, abs(c.strip_y1) + c.strip_width1 / 2 - c.width / 2) + max(
            0.0, abs(c.strip_y2) + c.strip_width2 / 2 - c.width / 2
        )
        return v * 1e3
    ext = c.patch_length / 2 + c.pad + c.diode_span + c.pad
    v = max(0.0, ext + max(c.stub1_length, c.stub2_length) - (c.period / 2 - m)) + max(
        0.0, c.patch_width - (c.period - 2 * m)
    )
    v += max(
        abs(c.stub1_y) + c.stub_width / 2 - c.patch_width / 2,
        abs(c.stub2_y) + c.stub_width / 2 - c.patch_width / 2,
        0.0,
    )
    return v * 1e3


def objective(x, fam, kits, pool, log, thetas=(0.0, 30.0)):
    lo = np.array([v[0] for v in SPACE[fam].values()])
    hi = np.array([v[1] for v in SPACE[fam].values()])
    xc = np.clip(x, lo, hi)
    pen = float(np.sum(np.maximum(0, lo - x) + np.maximum(0, x - hi)))
    c = make_cell(fam, xc)
    inf = infeasible(fam, c)
    if inf > 0:
        val = 50.0 + 10 * (pen + inf)
    else:
        resp = lp_responses(c, kits, thetas, pool)
        val = 60.0 if resp is None else lp_score(resp)[0] + 10 * pen
    log.append((list(map(float, xc)), float(val)))
    return val


def main(fam, kit_name, out, n_starts=3, maxfev=150):
    kit = madp_kit() if kit_name == "madp" else smp1345_kit()
    kits = corner_kits(kit)
    space = SPACE[fam]
    lo = np.array([v[0] for v in space.values()])
    hi = np.array([v[1] for v in space.values()])
    rng = np.random.default_rng(11)
    runs = []
    with Pool(11) as pool:
        for s in range(n_starts):
            base = np.array([getattr(CELLS[fam](), k) * 1e3 for k in space])  # the item-2 starting cell
            x0 = base
            if s > 0:  # perturbed starts around it
                for _ in range(200):
                    x0 = np.clip(base * rng.uniform(0.75, 1.25, base.size) + rng.normal(0, 0.15, base.size), lo, hi)
                    if infeasible(fam, make_cell(fam, x0)) == 0:
                        break
            log = []
            sim = np.vstack([x0, x0 + np.diag(0.15 * (hi - lo))])
            res = minimize(
                objective,
                x0,
                args=(fam, kits, pool, log),
                method="Nelder-Mead",
                options={"maxfev": maxfev, "xatol": 0.01, "fatol": 1e-3, "initial_simplex": sim},
            )
            best = min(log, key=lambda t: t[1])
            resp = lp_responses(make_cell(fam, best[0]), kits, (0.0, 30.0), pool)
            score, metrics = lp_score(resp) if resp else (None, {})
            runs.append({"x_mm": best[0], "J": best[1], "metrics": metrics, "nfev": res.nfev})
            print(json.dumps(runs[-1]), flush=True)
            Path(out).write_text(
                json.dumps({"family": fam, "kit": kit.name, "names": list(space), "runs": runs}, indent=1)
            )


if __name__ == "__main__":
    main(
        sys.argv[1],
        sys.argv[2],
        sys.argv[3],
        int(sys.argv[4]) if len(sys.argv) > 4 else 3,
        int(sys.argv[5]) if len(sys.argv) > 5 else 150,
    )
