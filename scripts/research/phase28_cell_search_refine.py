"""Phase 28 item 3: refinement search of the hub-and-sector cell at 6.4 mm (or 8.0 mm) around a screening result: wider neck/hub range, narrow gaps, thin edge margin.

    python scripts/research/phase28_cell_search_refine.py <out.json> --period 6.4 [--n 1200] [--workers 8] [--seed 2] [--kit madp|ideal]

Same evaluation and score as `phase28_cell_search_balanced.py`; only the sampling ranges differ (documented in `sample`).
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

from phase28_cell_search_balanced import evaluate  # noqa: E402


def sample(rng, p):
    half = p / 2
    ra = rng.uniform(half - 0.35, half - 0.1)
    return {
        "hub_radius": float(rng.choice([0.4, 0.8, 1.2, 1.6])) * 1e-3,
        "radial_gap": float(rng.choice([0.4, 0.8])) * 1e-3,
        "neck_width": float(rng.choice([0.4, 0.8, 1.2, 1.6])) * 1e-3,
        "sector_gap": 0.2e-3,
        "radius": ra * 1e-3,
        "radius_diag": float(rng.uniform(ra - 0.1, ra + 0.9)) * 1e-3,
        "edge_margin": 0.1e-3,
    }


def main(out, p, n, workers, seed, kit):
    rng = np.random.default_rng(seed)
    jobs = [(sample(rng, p), kit, p) for _ in range(n)]
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i, r in enumerate(pool.map(evaluate, jobs), 1):
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
        int(a[a.index("--n") + 1]) if "--n" in a else 1200,
        int(a[a.index("--workers") + 1]) if "--workers" in a else 8,
        int(a[a.index("--seed") + 1]) if "--seed" in a else 2,
        a[a.index("--kit") + 1] if "--kit" in a else "madp",
    )
