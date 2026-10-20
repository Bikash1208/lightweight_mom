"""Phase 20 R3 item 3, CP-2 stage 1: single-state PB condition (stub-0 diode ON, three OFF) at 10 GHz, Nelder-Mead over patch size, stub
line length and stub width, for a fixed period. Writes a JSON of the runs. MODEL_ONLY.

    python scripts/research/phase20_r3_opt_cp2.py <period_mm> <kit: madp|smp> <out.json> [n_starts]
"""

import json
import sys
from dataclasses import replace
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from phase20_r3_opt import ROOT, pb_band, pb_errors, pb_single  # noqa: F401
from scipy.optimize import minimize

from src.ris.cells import Cp2Cell, madp_kit, smp1345_kit

FREQS = np.arange(8.5e9, 11.51e9, 0.25e9)  # band evaluation grid
FREQS_OBJ = np.arange(9.0e9, 11.01e9, 0.5e9)  # objective grid (5 points)
BOUNDS = {"circumradius": (2.5e-3, 5.5e-3), "line": (0.25e-3, 5.0e-3), "stub_width": (0.4e-3, 3.0e-3)}
NAMES = list(BOUNDS)


def make_cell(x, period):
    c = Cp2Cell(period=period, **{k: float(v) * 1e-3 for k, v in zip(NAMES, x, strict=True)})
    room = (
        period / 2 - 0.25e-3 - 1e-6 - (c.stub_extent() - c.line)
    )  # the axis-aligned stubs must stay clear of the cell boundary
    return replace(c, line=max(0.25e-3, min(c.line, room)))


def objective(x, period, kit, pool, log):
    lo = np.array([BOUNDS[k][0] * 1e3 for k in NAMES])
    hi = np.array([BOUNDS[k][1] * 1e3 for k in NAMES])
    pen = float(np.sum(np.maximum(0, lo - x) + np.maximum(0, x - hi)))
    xc = np.clip(x, lo, hi)
    rows = pb_single(make_cell(xc, period), kit, FREQS_OBJ, pool)
    if rows is None:
        val = 1e3 + pen
    else:
        e = pb_errors(rows)
        i10 = int(np.argmin(abs(FREQS_OBJ - 10e9)))
        val = 5.0 * e[i10] + float(np.mean(e)) + 10 * pen
    log.append((list(map(float, xc)), float(val)))
    return val


def main(period_mm, kit_name, out, n_starts=4):
    kit = madp_kit() if kit_name == "madp" else smp1345_kit()
    period = period_mm * 1e-3
    rng = np.random.default_rng(7)
    runs = []
    with Pool(11) as pool:
        for _ in range(n_starts):
            x0 = np.array([rng.uniform(3.0, 5.0), rng.uniform(0.5, 3.0), rng.uniform(0.6, 2.0)])
            log = []
            res = minimize(
                objective,
                x0,
                args=(period, kit, pool, log),
                method="Nelder-Mead",
                options={
                    "maxfev": 90,
                    "xatol": 0.02,
                    "fatol": 1e-3,
                    "initial_simplex": x0 + np.vstack([np.zeros(3), np.eye(3) * [0.5, 0.5, 0.3]]),
                },
            )
            best = min(log, key=lambda t: t[1])
            cell = make_cell(best[0], period)
            rows = pb_single(cell, kit, FREQS, pool)
            runs.append(
                {
                    "x_mm": best[0],
                    "J": best[1],
                    "band_hz": pb_band(FREQS, rows) if rows is not None else 0.0,
                    "nfev": res.nfev,
                }
            )
            print(runs[-1], flush=True)
    Path(out).write_text(json.dumps({"period_mm": period_mm, "kit": kit.name, "names": NAMES, "runs": runs}, indent=1))


if __name__ == "__main__":
    main(float(sys.argv[1]), sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 4)
