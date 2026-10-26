"""Phase 20 R5-4 (SERVER JOBS `r5-supercell-gradients`, `r5-supercell-random`): array-level coupling check of the designed codings.

Supercell MoM (`src/ris/supercell.py`, Floquet cutoff F = 2n, baseline CP-2 cell, MADP diodes, f = 10.5 GHz, normal incidence) of

* `gradients`: periodic n x 1 (x gradient) and 1 x n (y gradient) windows, n <= 8, whose phase step 2 pi q / n is the best realisable approximation of a designed
  steering angle (sin theta = q lambda / (n P)); the states are the plain design of `src/ris/coding.py` for that gradient;
* `random`: 4 x 4 tiles cut from the designed 20 x 20 codings (dither seed 0 at 30 / 40 / 50 deg; plain at 30 / 50 deg for contrast) and the five R3 random
  codings repeated at F = 2n (they were only run at the minimum cutoff).

Every propagating Floquet order is compared with the isolated-cell prediction used by the design (`predict` of `phase20_r3_locality.py`, which equals the forward
model's order amplitudes). Gate 1 dB / 15 deg on the strongest orders: REPORTED, never tuned. A task = one supercell geometry with several codings that share the
impedance-matrix assembly; each finished task writes its own JSON (checkpoint: finished tasks are skipped on a rerun).

    python scripts/research/phase20_r5_supercell.py <gradients|random> <outdir> [--workers N] [--smoke] [--only NAME]
Run with OPENBLAS_NUM_THREADS=1 per worker. MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET, VIA_OMITTED_APPROXIMATION.
"""

import json
import math
import os
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

for _k in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "1")

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

import phase20_r3_locality as loc  # noqa: E402
from provenance import finish, provenance  # noqa: E402

from src.diagnostics import _ru_maxrss_kb  # noqa: E402
from src.ris import coding as cd  # noqa: E402
from src.ris.calibrated_library import CalibratedLibrary  # noqa: E402
from src.ris.polarization import L, R  # noqa: E402
from src.ris.supercell import build_supercell, default_mesh, solve_orders  # noqa: E402

F_GHZ = 10.5
PITCH = 15e-3
N_PER_CELL = 48
LAMBDA = 299792458.0 / (F_GHZ * 1e9)
loc.F0 = F_GHZ * 1e9  # the R3 helpers read this module global
loc.FLOQUET_FACTOR = 2.0  # converged cutoff (R4 audit)


def realisable_gradients(max_n=8, max_sin=0.8):
    """(q, n) with gcd 1, n <= max_n, sin theta = q lambda / (n P) <= max_sin; the designed angles are matched to the nearest of these."""
    out = []
    for n in range(2, max_n + 1):
        for q in range(1, n):
            if math.gcd(q, n) == 1 and q * LAMBDA / (n * PITCH) <= max_sin:
                out.append((q, n, round(math.degrees(math.asin(q * LAMBDA / (n * PITCH))), 2)))
    return sorted(out, key=lambda t: t[2])


def gradient_tasks(cal):
    table, _ = cd._state_table(cal, F_GHZ, 0.0, 0.0)
    tasks = []
    for q, n, ang in realisable_gradients():
        i = np.arange(n)
        psi = -2 * np.pi * q * i / n
        row = cd.nearest_states(table, psi[None, :]) - 1  # state index 0..3
        tasks.append(
            {
                "name": f"grad_x_q{q}_n{n}",
                "nx": n,
                "ny": 1,
                "angle_deg": ang,
                "axis": "x",
                "codings": [{"label": f"plain q={q}/{n}", "grid": [row[0].tolist()]}],
            }
        )
    # y gradients (1 x n): the cell is not 4-fold symmetric, so the phi_s = 90 plane is a separate physical case
    for q, n in ((1, 4), (2, 5)):
        i = np.arange(n)
        col = cd.nearest_states(table, (-2 * np.pi * q * i / n)[:, None]) - 1
        ang = round(math.degrees(math.asin(q * LAMBDA / (n * PITCH))), 2)
        tasks.append(
            {
                "name": f"grad_y_q{q}_n{n}",
                "nx": 1,
                "ny": n,
                "angle_deg": ang,
                "axis": "y",
                "codings": [{"label": f"plain q={q}/{n} along y", "grid": [[int(v)] for v in col[:, 0]]}],
            }
        )
    return tasks


TILE_ORIGINS = {"corner": (0, 0), "centre": (8, 8), "far_corner": (16, 16), "edge": (0, 16)}


def random_tasks(cal):
    rng = np.random.default_rng(20261020)  # the R3 seed: same five random 4 x 4 codings
    r3 = [[int(v) for v in rng.integers(0, 4, 16)] for _ in range(5)]
    sets = [{"label": f"R3 random #{k}", "grid": [c[4 * j : 4 * j + 4] for j in range(4)]} for k, c in enumerate(r3)]
    for var, th, ph in (
        ("dither", 30, 0),
        ("dither", 40, 0),
        ("dither", 50, 0),
        ("dither", 30, 90),
        ("dither", 50, 90),
        ("plain", 30, 0),
        ("plain", 50, 0),
    ):
        full = cd.design(cal, var, F_GHZ, float(th), float(ph), seed=0).states - 1
        for tname, (r0, c0) in TILE_ORIGINS.items():
            sets.append(
                {
                    "label": f"{var} theta_s={th} phi_s={ph} tile={tname}",
                    "grid": full[r0 : r0 + 4, c0 : c0 + 4].tolist(),
                }
            )
    per = 4  # codings per task share one impedance-matrix assembly
    return [
        {"name": f"tiles4x4_{k // per}", "nx": 4, "ny": 4, "codings": sets[k : k + per]}
        for k in range(0, len(sets), per)
    ]


def run_task(args):
    task, outdir, n_per_cell = args
    out = Path(outdir) / f"{task['name']}.json"
    if out.exists():
        return task["name"], "skipped (done)"
    t0 = time.time()
    cell, kit = loc.cell_and_kit()
    nx, ny = task["nx"], task["ny"]
    iso = loc.isolated_circular(cell, kit, 0.0, n_per_cell)
    states = cell.states()
    grids = [c["grid"] for c in task["codings"]]
    sc = build_supercell(cell, [[states[c] for c in row] for row in grids[0]], kit)
    load_sets = [build_supercell(cell, [[states[c] for c in row] for row in g], kit).loads for g in grids]
    t1 = time.time()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        all_orders = solve_orders(
            sc, F_GHZ * 1e9, default_mesh(n_per_cell, nx, ny, 2.0), loc.PL, 0.0, 0.0, load_sets, max_order=3
        )
    res = []
    for c, g, orders_s in zip(task["codings"], grids, all_orders, strict=True):
        meas = loc.measured(orders_s, 0.0, nx, ny, cell.period)
        keys = list(meas)
        rows = loc.compare(meas, loc.predict(iso, g, keys, R), loc.predict(iso, g, keys, L))
        strong = [r for r in rows if r["above_floor"]]
        res.append(
            {
                "label": c["label"],
                "codes": g,
                "orders": rows,
                "gate_1db_15deg": all(
                    abs(r["err_db"]) <= 1.0 and (r["err_deg"] is None or abs(r["err_deg"]) <= 15.0) for r in strong
                ),
                "max_abs_err_db_above_floor": max((abs(r["err_db"]) for r in strong), default=None),
            }
        )
    rec = {
        "task": {k: v for k, v in task.items() if k != "codings"},
        "n_per_cell": n_per_cell,
        "floquet_factor": 2.0,
        "f_ghz": F_GHZ,
        "results": res,
        "isolated_rr_deg": {str(k): float(np.degrees(np.angle(v[R, R]))) for k, v in iso.items()},
        "seconds_isolated": t1 - t0,
        "seconds_total": time.time() - t0,
        "peak_memory_kb": _ru_maxrss_kb(),
    }
    rec["provenance"] = finish(provenance(extra={"task": task["name"]}))
    out.write_text(json.dumps(rec, indent=1, default=float))
    return task["name"], f"done in {time.time() - t0:.0f} s"


def main(mode, outdir, workers, smoke, only):
    Path(outdir).mkdir(parents=True, exist_ok=True)
    cal = CalibratedLibrary.load()
    tasks = gradient_tasks(cal) if mode == "gradients" else random_tasks(cal)
    n = N_PER_CELL
    if smoke:
        tasks, n = tasks[:1], 8
    if only:
        tasks = [t for t in tasks if only in t["name"]]
    print(f"{len(tasks)} tasks, workers {workers}, n_per_cell {n}", flush=True)
    with Pool(workers) as pool:
        for name, msg in pool.imap_unordered(run_task, [(t, outdir, n) for t in tasks]):
            print(name, msg, flush=True)


if __name__ == "__main__":
    a = sys.argv
    main(
        a[1],
        a[2],
        int(a[a.index("--workers") + 1]) if "--workers" in a else 1,
        "--smoke" in a,
        a[a.index("--only") + 1] if "--only" in a else None,
    )
