"""Phase 28 item 3: random screening of the octagon + diode-stub cell (`ris.cells.Cp2Cell`, the repo's CP-2 family) scaled to a 6.4 mm cell at 10.0 GHz.

    python scripts/research/phase28_cell_search.py <out.json> [--n 600] [--workers 12] [--kit madp|smp] [--seed 1] [--period 6.4]

Each sample draws the octagon circumradius, stub widths / lengths, pad and diode footprint, evaluates two states (axis stub and diagonal stub) at mesh 32 with
`coded_cell_eval` (normal-incidence tensor sheet) and records b_bar, delta, dX_aniso (the item-2 tensor contrast measure). Screening only; results are not used
before the mesh-convergence study of the survivors.
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import json  # noqa: E402
import sys  # noqa: E402
import warnings  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

F = 10.0e9


def sample(rng, period_mm=6.4):
    r = rng.uniform(0.28 * period_mm, 0.45 * period_mm)
    return {
        "circumradius": r * 1e-3,
        "stub_width": rng.choice([0.3, 0.4, 0.5, 0.6, 0.8, 1.0]) * 1e-3,
        "pad": rng.choice([0.1, 0.15, 0.2]) * 1e-3,
        "diode_span": rng.choice([0.3, 0.4, 0.5, 0.6, 0.8]) * 1e-3,
        "line": rng.uniform(0.05, 0.8) * 1e-3,
        "line_diag": rng.uniform(0.05, 0.8) * 1e-3,
        "stub_width_diag": rng.choice([0.3, 0.4, 0.5, 0.6, 0.8, 1.0]) * 1e-3,
    }


def evaluate(args):
    params, kitname, period_mm = args
    warnings.simplefilter("ignore")
    from src.hologram.coded_cell_eval import evaluate_states
    from src.ris.cells import Cp2Cell, madp_kit, smp1345_kit

    kit = madp_kit() if kitname == "madp" else smp1345_kit()
    try:
        cell = Cp2Cell(period=period_mm * 1e-3, overlap=0.08e-3, **params)
        r = evaluate_states(cell, kit, F, 32, states=cell.states()[:2])
    except Exception as exc:  # noqa: BLE001
        return {"params": params, "error": f"{type(exc).__name__}: {exc}"[:120]}
    return {
        "params": params,
        "kit": kitname,
        "states": {k: {a: v for a, v in m.items() if a != "y_s"} for k, m in r.items()},
    }


def main(out, n, workers, kitname, seed, period_mm=6.4):
    rng = np.random.default_rng(seed)
    jobs = [(sample(rng, period_mm), kitname, period_mm) for _ in range(n)]
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i, r in enumerate(pool.map(evaluate, jobs), 1):
            rows.append(r)
            if i % 50 == 0:
                Path(out).write_text(json.dumps(rows), encoding="utf-8")
                print(i, flush=True)
    Path(out).write_text(json.dumps(rows), encoding="utf-8")


if __name__ == "__main__":
    a = sys.argv[1:]
    main(
        a[0],
        int(a[a.index("--n") + 1]) if "--n" in a else 600,
        int(a[a.index("--workers") + 1]) if "--workers" in a else 12,
        a[a.index("--kit") + 1] if "--kit" in a else "madp",
        int(a[a.index("--seed") + 1]) if "--seed" in a else 1,
        float(a[a.index("--period") + 1]) if "--period" in a else 6.4,
    )
