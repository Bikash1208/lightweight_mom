"""Phase 20 R3 item 3, LP-3 (3-diode fallback): many parallel Nelder-Mead starts; per candidate all 8 diode codes are solved (theta = 0, nominal) and the
best 4-code subset (closest to the 0/90/180/270 lattice, amplitude >= -1.5 dB with spread <= 1 dB) is scored. MODEL_ONLY.

    python scripts/research/phase20_r3_opt_lp3.py <madp|smp> <out.json> <n_starts> <maxfev> [workers]
"""

import itertools
import json
import sys
from dataclasses import dataclass
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from phase20_r3_opt import lattice_fit, mesh_for, solve_s11
from phase20_r3_opt_lp_multi import SerialPool, spread_term
from scipy.optimize import minimize

from src.ris.cells import Lp3Cell, madp_kit, smp1345_kit

SPACE = {
    "patch_length": (3.0, 8.0),
    "patch_width": (2.0, 9.0),
    "stub_width": (0.4, 3.0),
    "left_length": (0.25, 3.0),
    "right1_length": (0.25, 3.0),
    "right2_length": (0.25, 3.0),
    "left_y": (-2.0, 2.0),
    "right_y": (-2.0, 2.0),
    "period": (10.0, 15.0),
}


@dataclass(frozen=True)
class SubsetCell:
    """A cell restricted to a chosen list of diode codes (used for scoring / the finalist script)."""

    cell: object
    codes: tuple

    @property
    def period(self):
        return self.cell.period

    def states(self):
        return list(self.codes)

    def unit_cell(self, state, kit):
        return self.cell.unit_cell(state, kit)


def make_cell(x):
    return Lp3Cell(**{k: float(v) * 1e-3 for k, v in zip(SPACE, x, strict=True)})


def infeasible(c):
    m = 0.25e-3
    ext_r = c.patch_length / 2 + 2 * (2 * c.pad + c.diode_span) + c.right1_length + c.right2_length
    ext_l = c.patch_length / 2 + 2 * c.pad + c.diode_span + c.left_length
    v = max(0.0, max(ext_r, ext_l) - (c.period / 2 - m)) + max(0.0, c.patch_width - (c.period - 2 * m))
    v += max(
        0.0, abs(c.left_y) + c.stub_width / 2 - c.patch_width / 2, abs(c.right_y) + c.stub_width / 2 - c.patch_width / 2
    )
    return v * 1e3


def best_subset(phase, amp):
    """(score, code indices) over the 70 four-code subsets."""
    amp_db = 20 * np.log10(np.maximum(amp, 1e-9))
    best = (1e9, None)
    for sub in itertools.combinations(range(len(phase)), 4):
        ph, a = phase[list(sub)], amp_db[list(sub)]
        dev, _ = lattice_fit(ph)
        e = (
            max(0.0, dev - 15) / 15
            + dev / 600
            + max(0.0, -1.5 - a.min()) / 1.5
            + max(0.0, a.max() - a.min() - 1.0)
            + spread_term(ph)
        )
        if e < best[0]:
            best = (e, sub)
    return best


def eval_codes(c, kit):
    n = mesh_for(c.period)
    res = SerialPool.map(solve_s11, [(c, st, kit, 10e9, 0.0, n) for st in c.states()])
    if any(r is None for r in res):
        return None
    return (
        np.degrees(np.angle([r[1, 1] for r in res])),
        np.array([abs(r[1, 1]) for r in res]),
        np.array([abs(r[0, 1]) for r in res]),
    )


def objective(x, kit, log):
    lo = np.array([v[0] for v in SPACE.values()])
    hi = np.array([v[1] for v in SPACE.values()])
    xc = np.clip(x, lo, hi)
    pen = float(np.sum(np.maximum(0, lo - x) + np.maximum(0, x - hi)))
    c = make_cell(xc)
    inf = infeasible(c)
    if inf > 0:
        val = 50.0 + 10 * (pen + inf)
    else:
        r = eval_codes(c, kit)
        val = 60.0 if r is None else best_subset(r[0], r[1])[0] + 10 * pen
    log.append((list(map(float, xc)), float(val)))
    return val


def one_start(args):
    kit_name, seed, maxfev = args
    kit = madp_kit() if kit_name == "madp" else smp1345_kit()
    lo = np.array([v[0] for v in SPACE.values()])
    hi = np.array([v[1] for v in SPACE.values()])
    rng = np.random.default_rng(2000 + seed)
    base = np.array([getattr(Lp3Cell(), k) * 1e3 for k in SPACE])
    for _ in range(500):
        x0 = (
            base
            if seed == 0
            else (np.clip(base * rng.uniform(0.7, 1.3, base.size), lo, hi) if seed % 2 else rng.uniform(lo, hi))
        )
        if infeasible(make_cell(x0)) == 0:
            break
    log = []
    sim = np.vstack([x0, x0 + np.diag(0.12 * (hi - lo))])
    res = minimize(
        objective,
        x0,
        args=(kit, log),
        method="Nelder-Mead",
        options={"maxfev": maxfev, "xatol": 0.005, "fatol": 1e-4, "initial_simplex": sim},
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
    kit_name, out, n_starts, maxfev = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    workers = int(sys.argv[5]) if len(sys.argv) > 5 else 6
    runs = []
    with Pool(workers) as pool:
        for r in pool.imap_unordered(one_start, [(kit_name, s, maxfev) for s in range(n_starts)]):
            runs.append(r)
            print(json.dumps(r), flush=True)
            Path(out).write_text(
                json.dumps({"family": "lp3", "kit": kit_name, "names": list(SPACE), "runs": runs}, indent=1)
            )
