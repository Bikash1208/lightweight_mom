"""Phase 28 item 3: random search of the hub-and-sector cell for BALANCED states (the axis state and the diagonal state must carry the same tensor: rotation equivalence).

    python scripts/research/phase28_cell_search_balanced.py <out.json> --period 6.4|8.0 [--n 400] [--workers 8] [--seed 1] [--kit madp|ideal]

Load dimensions (hub radius, radial gap, neck width) are multiples of 0.4 mm so that they stay on the mesh for the convergence series (6.4 mm: mesh 32, 48, 64;
8.0 mm: mesh 40, 60, 80). Screening mesh: 32 (6.4 mm) / 40 (8.0 mm). Score = the smaller of the two states' leakage `alpha_beam` (beam at 30 deg, Phase 26 closed
form from the cell's own sheet tensor); `balance` = max relative difference of b_bar and delta between the two states.
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


def sample(rng, p):
    half = p / 2
    ra = rng.uniform(half - 0.8, half - 0.1)
    return {
        "hub_radius": float(rng.choice([0.4, 0.8, 1.2])) * 1e-3,
        "radial_gap": float(rng.choice([0.4, 0.8])) * 1e-3,
        "neck_width": float(rng.choice([0.4, 0.8])) * 1e-3,
        "sector_gap": float(rng.choice([0.2, 0.4])) * 1e-3,
        "radius": ra * 1e-3,
        "radius_diag": float(rng.uniform(ra, ra + 1.4)) * 1e-3,
        "edge_margin": float(rng.choice([0.1, 0.2])) * 1e-3,
    }


def evaluate(args):
    params, kitname, p = args
    warnings.simplefilter("ignore")
    from src.hologram.coded_cell_candidates import SectorHubCell
    from src.hologram.coded_cell_eval import evaluate_states
    from src.ris.cells import ideal_kit, madp_kit

    kit = ideal_kit() if kitname == "ideal" else madp_kit()
    mesh = 32 if p < 7 else 40
    try:
        cell = SectorHubCell(period=p * 1e-3, **params)
        r = evaluate_states(cell, kit, 10.0e9, mesh, states=cell.states()[:2])
    except Exception as exc:  # noqa: BLE001
        return {"params": params, "error": f"{type(exc).__name__}: {exc}"[:120]}
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


def main(out, p, n, workers, seed, kit, mode="bowtie"):
    rng = np.random.default_rng(seed)
    jobs = [({**sample(rng, p), "mode": mode}, kit, p) for _ in range(n)]
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i, r in enumerate(pool.map(evaluate, jobs), 1):
            rows.append(r)
            if i % 50 == 0:
                print(i, flush=True)
                Path(out).write_text(json.dumps(rows), encoding="utf-8")
    Path(out).write_text(json.dumps(rows), encoding="utf-8")


if __name__ == "__main__":
    a = sys.argv[1:]
    main(
        a[0],
        float(a[a.index("--period") + 1]),
        int(a[a.index("--n") + 1]) if "--n" in a else 400,
        int(a[a.index("--workers") + 1]) if "--workers" in a else 8,
        int(a[a.index("--seed") + 1]) if "--seed" in a else 1,
        a[a.index("--kit") + 1] if "--kit" in a else "madp",
        a[a.index("--mode") + 1] if "--mode" in a else "bowtie",
    )
