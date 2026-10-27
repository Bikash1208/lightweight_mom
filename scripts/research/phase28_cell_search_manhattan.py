"""Phase 28 item 3: grid screening of the all-axis-aligned (Manhattan) hub cell `coded_cell_candidates.ManhattanHubCell` at 10.0 GHz (every dimension a multiple of 0.4 mm:
exactly on the mesh for mesh 16 k on a 6.4 mm cell, 20 k on an 8.0 mm cell).

    python scripts/research/phase28_cell_search_manhattan.py <out.json> --period 6.4|8.0 [--workers 6] [--kit madp|ideal]
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
    from src.hologram.coded_cell_candidates import ManhattanHubCell
    from src.hologram.coded_cell_eval import evaluate_states
    from src.ris.cells import ideal_kit, madp_kit

    kit = ideal_kit() if kitname == "ideal" else madp_kit()
    mesh = 32 if period_mm < 7 else 40
    try:
        cell = ManhattanHubCell(period=period_mm * 1e-3, **params)
        r = evaluate_states(cell, kit, 10.0e9, mesh, states=cell.states()[:2])
    except Exception as exc:  # noqa: BLE001
        return {"params": params, "kit": kitname, "error": f"{type(exc).__name__}: {exc}"[:120]}
    s0, s1 = r[0], r[1]
    a0, a1 = s0.get("alpha_beam_np_per_m") or 0.0, s1.get("alpha_beam_np_per_m") or 0.0
    bal = max(
        abs(s0["b_bar_ms"] - s1["b_bar_ms"]) / max(abs(s0["b_bar_ms"]), abs(s1["b_bar_ms"]), 1e-9),
        abs(s0["delta_ms"] - s1["delta_ms"]) / max(s0["delta_ms"], s1["delta_ms"], 1e-9),
    )
    return {
        "params": params,
        "kit": kitname,
        "mesh": mesh,
        "score": min(a0, a1),
        "balance": bal,
        "states": {k: {a: v for a, v in m.items() if a != "y_s"} for k, m in r.items()},
    }


def main(out, period_mm, workers, kit, mode="bowtie"):
    half = period_mm / 2
    outers = (half - 0.4, half - 0.8, half - 1.2)
    grid = []
    for h, g, nw, wa, sep, o in itertools.product(
        (0.8, 1.2, 1.6, 2.0, 2.4), (0.4, 0.8), (0.4, 0.8), (0.4, 0.8, 1.2, 1.6, 2.0), (0.4, 0.8), outers
    ):
        if (
            h - (wa + sep) < 0.4 - 1e-9
        ):  # the corner neck needs an overlap of at least 0.4 mm between the hub edge and the corner square
            continue
        if h + g >= o - 0.4:
            continue
        p = {
            "hub_half": h * 1e-3,
            "neck_gap": g * 1e-3,
            "neck_width": nw * 1e-3,
            "arm_half_width": wa * 1e-3,
            "separation": sep * 1e-3,
            "outer": o * 1e-3,
        }
        grid.append(({**p, "mode": mode}, kit, period_mm))
    print(len(grid), "candidates", flush=True)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        rows = list(pool.map(evaluate, grid))
    Path(out).write_text(json.dumps(rows), encoding="utf-8")


if __name__ == "__main__":
    a = sys.argv[1:]
    main(
        a[0],
        float(a[a.index("--period") + 1]),
        int(a[a.index("--workers") + 1]) if "--workers" in a else 6,
        a[a.index("--kit") + 1] if "--kit" in a else "madp",
        a[a.index("--mode") + 1] if "--mode" in a else "bowtie",
    )
