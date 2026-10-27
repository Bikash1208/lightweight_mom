"""Phase 28 item 3: scan (axis radius, diagonal radius) of the hub-and-sector cell to equalise the axis and diagonal states (rotation equivalence).

python scripts/research/phase28_cell_scan_radii.py <out.json> --period 8.0 [--kit madp|ideal] [--workers 6]
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import itertools  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import warnings  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))


def evaluate(args):
    params, kitname, period_mm = args
    warnings.simplefilter("ignore")
    from src.hologram.coded_cell_candidates import SectorHubCell
    from src.hologram.coded_cell_eval import evaluate_states
    from src.ris.cells import ideal_kit, madp_kit

    kit = ideal_kit() if kitname == "ideal" else madp_kit()
    try:
        cell = SectorHubCell(period=period_mm * 1e-3, **params)
        r = evaluate_states(cell, kit, 10.0e9, 32, states=cell.states()[:2])
    except Exception as exc:  # noqa: BLE001
        return {"params": params, "kit": kitname, "error": f"{type(exc).__name__}: {exc}"[:120]}
    return {
        "params": params,
        "kit": kitname,
        "states": {k: {a: v for a, v in m.items() if a != "y_s"} for k, m in r.items()},
    }


def main(out, period_mm, kit, workers):
    base = {"hub_radius": 0.8e-3, "radial_gap": 0.8e-3, "neck_width": 0.4e-3, "sector_gap": 0.4e-3}
    if period_mm > 7:
        ra, rd = (2.8, 3.0, 3.2, 3.4, 3.6), (3.6, 4.0, 4.4, 4.8, 5.2)
    else:
        ra, rd = (2.4, 2.6, 2.8, 3.0), (2.8, 3.2, 3.6, 4.0, 4.4)
    grid = [
        ({**base, "radius": a * 1e-3, "radius_diag": d * 1e-3}, kit, period_mm) for a, d in itertools.product(ra, rd)
    ]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        rows = list(pool.map(evaluate, grid))
    Path(out).write_text(json.dumps(rows), encoding="utf-8")


if __name__ == "__main__":
    a = sys.argv[1:]
    main(
        a[0],
        float(a[a.index("--period") + 1]),
        a[a.index("--kit") + 1] if "--kit" in a else "madp",
        int(a[a.index("--workers") + 1]) if "--workers" in a else 6,
    )
