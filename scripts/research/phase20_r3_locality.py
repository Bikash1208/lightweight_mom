"""Phase 20 R3 item 4: locality test of the CP-2 cell. Supercell MoM (`src/ris/supercell.py`) with (a) the periodic 0-1-2-3 gradient (4 x 1 supercell,
both directions) and (b) random 2-bit codings on 4 x 4, at theta = 0 and 30 deg, 10 GHz, MADP diodes; every propagating Floquet order is compared
with the isolated-cell prediction (each site's isolated-cell r_RR, summed with the array factor; the project's Phase 13 convention: no element
pattern, no coupling). The opposite-handed term is tracked per order: it carries no PB phase, so it should go to the specular order.

Gate (R3 brief item 4): anomalous-order level within 1 dB and 15 deg of the prediction. Reported, never tuned.
MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET, VIA_OMITTED_APPROXIMATION.

    python scripts/research/phase20_r3_locality.py <gradient|random> <out.json> [n_per_cell] [n_codings]
"""

import json
import math
import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

from phase20_r3_finalist import build  # noqa: E402

from src.excitation import ReferencePlanes  # noqa: E402
from src.ris.cells import madp_kit  # noqa: E402
from src.ris.polarization import L, R, circular_reflection  # noqa: E402
from src.ris.supercell import build_supercell, default_mesh, solve_orders  # noqa: E402
from src.sparameters import compute_s_parameters_oblique  # noqa: E402

PL = ReferencePlanes(0.0, None, 1.575e-3)
F0 = 10e9
FLOQUET_FACTOR = (
    float(os.environ["LOCALITY_FLOQUET_FACTOR"]) if os.environ.get("LOCALITY_FLOQUET_FACTOR") else None
)  # opt-in: 2 = converged cutoff (R4 audit)
PARAMS_JSON = ROOT / "results" / "phase20_r3" / "opt" / "cp2_madp_final.json"


def cell_and_kit():
    d = json.loads(PARAMS_JSON.read_text())
    return build("cp2", d["params_mm"], d["period_mm"]), madp_kit()


def isolated_circular(cell, kit, theta_deg, n):
    """{state index: circular 2x2 [[r_RR, r_RL], [r_LR, r_LL]]} of the isolated cell at (theta, phi = 0), mesh n."""
    out = {}
    from src.basis_functions import MeshSettings

    for k, st in enumerate(cell.states()):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            s = compute_s_parameters_oblique(
                cell.unit_cell(st, kit),
                F0,
                MeshSettings(n, n, n // 2 + 1 if FLOQUET_FACTOR is None else int(FLOQUET_FACTOR * n)),
                PL,
                math.radians(theta_deg),
                0.0,
                warn_grating=False,
            )
        out[k] = circular_reflection(s.s11)
    return out


def predict(iso, code_grid, orders, which):
    """Array-factor prediction b_(m,n) = (1/N) sum_site r_site exp(+j 2 pi (m x + n y)/(N P)) (sites centred as in `build_supercell`),
    with r_site = iso[state][a, R] (a = R: handedness-preserving, a = L: opposite-handed) for R-handed incidence."""
    ny, nx = len(code_grid), len(code_grid[0])
    out = {}
    for m, n in orders:
        acc = 0j
        for j in range(ny):
            for i in range(nx):
                x, y = (i - (nx - 1) / 2) / nx, (j - (ny - 1) / 2) / ny  # in units of the supercell period
                acc += iso[code_grid[j][i]][which, R] * np.exp(1j * 2 * np.pi * (m * x + n * y))
        out[(m, n)] = acc / (nx * ny)
    return out


def measured(orders_s, theta_deg, nx, ny, period):
    """{(m, n): (b_RR, b_LR)}: reflected R- and L-handed amplitudes per unit incident R, from the 2x2 (TE, TM) order matrices. Each order's
    (TM, TE) frame is rotated by that order's azimuth phi_o = atan2(ky, kx), which multiplies the circular amplitudes by exp(+-j phi_o)
    (p - j q = exp(j phi) (x - j y)); the factor is removed so all orders are expressed in the fixed Cartesian frame, the frame of the
    isolated-cell prediction (without it the m < 0 orders would show a spurious 180 deg)."""
    k0 = 2 * math.pi * F0 / 299792458.0
    kx0 = k0 * math.sin(math.radians(theta_deg))
    out = {}
    for (m, n), s in orders_s.items():
        kx, ky = kx0 + 2 * math.pi * m / (nx * period), 2 * math.pi * n / (ny * period)
        phi_o = math.atan2(ky, kx) if math.hypot(kx, ky) > 0 else 0.0
        c = circular_reflection(s)
        out[(m, n)] = (c[R, R] * np.exp(1j * phi_o), c[L, R] * np.exp(-1j * phi_o))
    return out


def compare(meas, pred_rr, pred_lr, floor_db=-30.0):
    rows = []
    for k, (rr, lr) in meas.items():
        p_rr, p_lr = pred_rr[k], pred_lr[k]
        rows.append(
            {
                "order": list(k),
                "meas_rr_db": float(20 * np.log10(max(abs(rr), 1e-9))),
                "pred_rr_db": float(20 * np.log10(max(abs(p_rr), 1e-9))),
                "err_db": float(20 * np.log10(max(abs(rr), 1e-9) / max(abs(p_rr), 1e-9))),
                "err_deg": float(np.degrees(np.angle(rr / p_rr))) if abs(p_rr) > 1e-9 else None,
                "meas_lr_db": float(20 * np.log10(max(abs(lr), 1e-9))),
                "pred_lr_db": float(20 * np.log10(max(abs(p_lr), 1e-9))),
                "above_floor": bool(20 * np.log10(max(abs(p_rr), 1e-9)) > floor_db),
            }
        )
    return rows


def run_cases(name, cell, kit, iso_by_theta, codings, nx, ny, n, out):
    res = {}
    for th in (0.0, 30.0):
        grids = [[[c[j * nx + i] for i in range(nx)] for j in range(ny)] for c in codings]
        states = cell.states()
        sc = build_supercell(cell, [[states[c] for c in row] for row in grids[0]], kit)
        load_sets = []
        for g in grids:
            load_sets.append(build_supercell(cell, [[states[c] for c in row] for row in g], kit).loads)
        t0 = time.time()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            all_orders = solve_orders(
                sc, F0, default_mesh(n, nx, ny, FLOQUET_FACTOR), PL, math.radians(th), 0.0, load_sets, max_order=3
            )
        for g, orders_s in zip(grids, all_orders, strict=True):
            meas = measured(orders_s, th, nx, ny, cell.period)
            keys = list(meas)
            pred_rr = predict(iso_by_theta[th], g, keys, R)
            pred_lr = predict(iso_by_theta[th], g, keys, L)
            rows = compare(meas, pred_rr, pred_lr)
            res.setdefault(f"theta{th:g}", []).append({"codes": g, "orders": rows, "seconds": time.time() - t0})
        Path(out).write_text(json.dumps({"case": name, "n_per_cell": n, "nx": nx, "ny": ny, "results": res}, indent=1))
        print(name, th, "done", time.time() - t0, flush=True)
    return res


def main(mode, out, n=48, n_codings=5):
    cell, kit = cell_and_kit()
    iso = {th: isolated_circular(cell, kit, th, n) for th in (0.0, 30.0)}
    if mode == "gradient":
        # states 0-1-2-3 along x and the reversed gradient 3-2-1-0 (the PB phase ladder rises/falls by ~90 deg per step)
        run_cases("gradient", cell, kit, iso, [(0, 1, 2, 3), (3, 2, 1, 0)], 4, 1, n, out)
    else:
        rng = np.random.default_rng(20261020)
        codings = [tuple(int(v) for v in rng.integers(0, 4, 16)) for _ in range(n_codings)]
        run_cases("random4x4", cell, kit, iso, codings, 4, 4, n, out)


if __name__ == "__main__":
    a = sys.argv
    main(a[1], a[2], int(a[3]) if len(a) > 3 else 48, int(a[4]) if len(a) > 4 else 5)
