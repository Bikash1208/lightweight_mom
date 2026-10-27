"""Phase 28 item 3: grid screening of the hub-and-sector cell (`coded_cell_candidates.SectorHubCell`) at 10.0 GHz.

    python scripts/research/phase28_cell_search_sector.py <out.json> [--period 6.4] [--workers 6]

All lengths are multiples of 0.4 mm (mesh-aligned for mesh 16, 32, 48, 64 on a 6.4 mm cell). Two diode kits: the ideal switch (upper bound of what the
topology can do) and the MADP-000907-14020 model (L_s NOT_AVAILABLE, 0 assumed). Two states per candidate (axis neck and diagonal neck).
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
F = 10.0e9


def evaluate(args):
    params, kitname, period_mm, single = args
    warnings.simplefilter("ignore")
    from src.hologram.coded_cell_candidates import SectorHubCell
    from src.hologram.coded_cell_eval import evaluate_states
    from src.ris.cells import ideal_kit, madp_kit

    kit = ideal_kit() if kitname == "ideal" else madp_kit()
    try:
        cell = SectorHubCell(period=period_mm * 1e-3, single_sided=single, **params)
        r = evaluate_states(cell, kit, F, 32, states=cell.states()[:2])
    except Exception as exc:  # noqa: BLE001
        return {"params": params, "kit": kitname, "single": single, "error": f"{type(exc).__name__}: {exc}"[:120]}
    return {
        "params": params,
        "kit": kitname,
        "single": single,
        "states": {k: {a: v for a, v in m.items() if a != "y_s"} for k, m in r.items()},
    }


def main(out, period_mm, workers):
    radii = (2.8, 3.0) if period_mm < 7 else (3.6, 3.8)
    grid = []
    for hub, gap, neck, sg, rad in itertools.product((0.4, 0.8, 1.2), (0.4, 0.8), (0.4, 0.8, 1.2), (0.4,), radii):
        p = {
            "hub_radius": hub * 1e-3,
            "radial_gap": gap * 1e-3,
            "neck_width": neck * 1e-3,
            "sector_gap": sg * 1e-3,
            "radius": rad * 1e-3,
        }
        for kit in ("ideal", "madp"):
            grid.append((p, kit, period_mm, False))
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i, r in enumerate(pool.map(evaluate, grid), 1):
            rows.append(r)
            if i % 12 == 0:
                print(i, len(grid), flush=True)
    Path(out).write_text(json.dumps(rows), encoding="utf-8")


if __name__ == "__main__":
    a = sys.argv[1:]
    main(
        a[0],
        float(a[a.index("--period") + 1]) if "--period" in a else 6.4,
        int(a[a.index("--workers") + 1]) if "--workers" in a else 6,
    )
