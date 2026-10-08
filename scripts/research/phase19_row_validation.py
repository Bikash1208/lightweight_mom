"""Phase 19 D-bis: validate the neighbour-coupling model on REAL cell sequences from the lens layout.

    .venv/bin/python scripts/research/phase19_row_validation.py

Takes 9-cell windows of consecutive cells along x from the committed layout, solves each as an exact periodic supercell
(RCWA, normal incidence, x-polarised), and compares the propagating diffraction orders with
  (a) the isolated-cell prediction  (uniform patch per cell, T_c = isolated response)   and
  (b) the coupling model S          (eps_c = K * mean(d_n), x-neighbours only, K = pool median for N differing neighbours).
Reports the order-spectrum error and the total propagating power ratio.  Writes results/transmitarray_dual30/row_validation.json.
"""

import json
import math
import pathlib
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from phase19_error_bar import D_MIN, k_pools  # noqa: E402
from phase19_local_periodicity import EPS, PLATE, TAND, TOTAL, F, P, isolated_txx  # noqa: E402

from src.dielectric_cell.supercell import CellSpec, solve_grid  # noqa: E402

OUT = ROOT / "results" / "transmitarray_dual30"
K_CELLS = 9
ORDERS = list(range(-4, 5))  # |u| = |m| lambda / (9 p) < 1 for |m| <= 4
HARM = (6 * K_CELLS, 8)


def exact_job(specs):
    grid = [[CellSpec(k, r * 1e-3)] for k, r in specs]
    res = solve_grid(grid, P, TOTAL - PLATE, PLATE, EPS, TAND, True, F, 0.0, 0.0, HARM, orders=[(m, 0) for m in ORDERS])
    return [complex(res["orders"][(m, 0)][0, 1]) for m in ORDERS]  # x-pol in, E_x out


def patch_orders(fields):
    k = len(fields)
    return np.array(
        [
            sum(
                e * np.sinc(m / k) / k * np.exp(1j * 2 * math.pi * m * (c - (k - 1) / 2) / k)
                for c, e in enumerate(fields)
            )
            for m in ORDERS
        ]
    )


def main():
    layout = json.loads((OUT / "layout.json").read_text())["cells"]
    samples = json.loads((OUT / "coupling_samples.json").read_text())["samples"]
    sup = json.loads((OUT / "coupling_superposition.json").read_text())["pairs"]
    pools = k_pools(samples, sup)
    med = {k: complex(np.median(v.real), np.median(v.imag)) for k, v in pools.items()}
    # windows near the lens centre (small incidence angle), different rows and x offsets
    picks = [(2.5, 0), (2.5, 9), (12.5, 3), (22.5, 6)]
    windows = []
    for y, off in picks:
        row = sorted((c for c in layout if abs(c["y_mm"] - y) < 1e-6), key=lambda c: c["x_mm"])
        mid = len(row) // 2 - K_CELLS // 2 + off
        windows.append([(c["kind"], round(c["radius_mm"], 6)) for c in row[mid : mid + K_CELLS]])
    uniq = sorted({s for w in windows for s in w})
    with ProcessPoolExecutor(max_workers=4) as ex:
        t_iso = dict(zip(uniq, ex.map(isolated_txx, [k for k, _ in uniq], [r for _, r in uniq]), strict=True))
        exact = list(ex.map(exact_job, windows))
    out = []
    for w, e_ex in zip(windows, exact, strict=True):
        t = np.array([t_iso[s] for s in w])
        e_iso = patch_orders(list(t))
        eps = np.zeros(len(t), complex)
        for i in range(len(t)):
            diff = [t[j] / t[i] - 1.0 for j in (i - 1, (i + 1) % len(t)) if abs(t[j] / t[i] - 1.0) >= D_MIN]
            if diff:
                n = len(diff)
                eps[i] = med[(n, "x")] * np.mean(diff)
        e_s = patch_orders(list(t * (1.0 + eps)))
        e_ex = np.array(e_ex)
        pw = lambda v: float(np.sum(np.abs(v) ** 2))  # noqa: E731
        out.append(
            dict(
                cells=[list(s) for s in w],
                isolated_phase_steps_deg=[
                    float(abs(np.degrees(np.angle(t[(i + 1) % len(t)] / t[i])))) for i in range(len(t))
                ],
                rel_err_isolated=float(np.linalg.norm(e_ex - e_iso) / np.linalg.norm(e_iso)),
                rel_err_model_S=float(np.linalg.norm(e_ex - e_s) / np.linalg.norm(e_iso)),
                power_ratio_exact_over_isolated=pw(e_ex) / pw(e_iso),
                power_ratio_modelS_over_isolated=pw(e_s) / pw(e_iso),
                orders_exact=[[v.real, v.imag] for v in e_ex],
                orders_isolated=[[v.real, v.imag] for v in e_iso],
                orders_modelS=[[v.real, v.imag] for v in e_s],
            )
        )
        print(
            f"rel err isolated {out[-1]['rel_err_isolated']:.3f}  model S {out[-1]['rel_err_model_S']:.3f}  "
            f"power exact/iso {out[-1]['power_ratio_exact_over_isolated']:.3f} (dB {10 * math.log10(out[-1]['power_ratio_exact_over_isolated']):.2f})  "
            f"model S/iso {out[-1]['power_ratio_modelS_over_isolated']:.3f}",
            flush=True,
        )
    (OUT / "row_validation.json").write_text(
        json.dumps(dict(k_cells=K_CELLS, orders=ORDERS, harmonics=HARM, windows=out), indent=1)
    )


if __name__ == "__main__":
    main()
