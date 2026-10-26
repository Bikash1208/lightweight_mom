"""Phase 20 R3b optimiser pieces: geometry clamps keep every candidate inside the cell, the cache replay of a Nelder-Mead start is deterministic and covers the committed evaluations."""

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "research"))
import phase20_r3b_opt as opt  # noqa: E402

from src.ris.cells import madp_kit  # noqa: E402


def test_every_candidate_inside_the_bounds_is_a_valid_cell():
    rng = np.random.default_rng(3)
    lo = np.array([b[0] for b in opt.BOUNDS])
    hi = np.array([b[1] for b in opt.BOUNDS])
    xs = [lo, hi, *[lo + rng.random(5) * (hi - lo) for _ in range(40)]]
    for x in xs:
        cell = opt.make_cell(x)
        for st in cell.states():
            cell.unit_cell(st, madp_kit())  # raises ValueError if a stub leaves the cell


def test_cache_replay_is_deterministic_and_matches_the_committed_evaluations():
    cache = opt.cached()
    if not cache:  # a fresh clone without the cache file
        return
    for start in ("baseline", "diag_long", "diag_longer_wide"):
        a = opt.trajectory(start)
        b = opt.trajectory(start)
        assert a == b and len(a) >= 6  # at least the initial simplex
        assert all(k in cache for k in a)
        assert a[0] == tuple(round(v, 4) for v in opt.STARTS[start])
