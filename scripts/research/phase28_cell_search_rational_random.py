"""Phase 28 item 3: random search of `coded_cell_candidates.RationalSectorCell` (rational-slope wedges; mesh-robust) with every dimension a multiple of 0.2 mm and separate
arm / diagonal neck widths. Score = the smaller of the two states' leakage at the beam (alpha_beam); `balance` as in `phase28_cell_search_balanced.py`.

    python scripts/research/phase28_cell_search_rational_random.py <out.json> --period 6.4|8.0 --mode notch|bowtie [--n 1500] [--workers 12] [--seed 5]
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import json  # noqa: E402
import sys  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

from phase28_cell_search_rational import evaluate  # noqa: E402


def sample(rng, half, mode):
    def pick(lo, hi):
        vals = np.round(np.arange(lo, hi + 1e-9, 0.2), 1)
        return float(rng.choice(vals))

    ao = pick(half - 1.0, half - 0.2)
    return {
        "hub_half": pick(1.0, half - 1.6) * 1e-3,
        "neck_gap": pick(0.4, 0.8) * 1e-3,
        "neck_width": pick(0.4, 0.8) * 1e-3,
        "neck_width_diag": pick(0.4, 0.8) * 1e-3,
        "gap_t": float(rng.choice([0.2, 0.4])) * 1e-3,
        "arm_outer": ao * 1e-3,
        "diag_outer": pick(half - 1.4, half - 0.2) * 1e-3,
        "clip": (half - 0.2) * 1e-3,
        "mode": mode,
    }


def main(out, period_mm, mode, n, workers, seed):
    rng = np.random.default_rng(seed)
    half = period_mm / 2
    jobs = [(sample(rng, half, mode), "madp", period_mm) for _ in range(n)]
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i, r in enumerate(pool.map(evaluate, jobs, chunksize=4), 1):
            rows.append(r)
            if i % 100 == 0:
                print(i, flush=True)
                Path(out).write_text(json.dumps(rows), encoding="utf-8")
    Path(out).write_text(json.dumps(rows), encoding="utf-8")


if __name__ == "__main__":
    a = sys.argv[1:]
    main(
        a[0],
        float(a[a.index("--period") + 1]),
        a[a.index("--mode") + 1],
        int(a[a.index("--n") + 1]) if "--n" in a else 1500,
        int(a[a.index("--workers") + 1]) if "--workers" in a else 12,
        int(a[a.index("--seed") + 1]) if "--seed" in a else 5,
    )
