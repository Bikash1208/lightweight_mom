"""Phase 19 D-bis: per-cell neighbour-coupling samples from exact 1 x 4 supercell strips [A A B B].

    .venv/bin/python scripts/research/phase19_coupling_samples.py

For each strip the exact diffraction orders m = -1, 0, 1, 2 (a 4-point DFT of the cell-averaged exit field, uniform-patch
representation -- the same representation the lens pattern model uses) are inverted to the four *equivalent per-cell
fields* E_c.  eps_c = E_c / T_c - 1 is then the complex perturbation of cell c relative to its isolated response T_c, and
kappa_c = eps_c / (T_n / T_c - 1) its coupling to the one differing neighbour n.  Both polarisations are recorded:
x-pol (E normal to the A|B interface) and y-pol (E parallel to it).  Writes results/transmitarray_dual30/coupling_samples.json.
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

from phase19_local_periodicity import (  # noqa: E402
    EPS,
    ORDERS,
    PLATE,
    STRIP_HARM,
    TAND,
    TOTAL,
    F,
    P,
    isolated_txx,
    pick_pairs,
)

from src.dielectric_cell.supercell import CellSpec, solve_strip  # noqa: E402

OUT = ROOT / "results" / "transmitarray_dual30"


def strip_job(args):
    specs, harm = args
    res = solve_strip(
        [CellSpec(k, r * 1e-3) for k, r in specs],
        P,
        TOTAL - PLATE,
        PLATE,
        EPS,
        TAND,
        True,
        F,
        0.0,
        0.0,
        harm,
        orders=[(m, 0) for m in ORDERS],
    )
    # [x/y, s/p]: x-pol incidence = p column, co-pol = x row; y-pol incidence = s column, co-pol = y row
    return {m: (complex(res["orders"][(m, 0)][0, 1]), complex(res["orders"][(m, 0)][1, 0])) for m in ORDERS}


def dft_matrix(k=4):
    """A[m, c]: order-m coefficient of a piecewise-constant field with value E_c on cell c (cells centred at c - (k-1)/2)."""
    return np.array(
        [[np.sinc(m / k) / k * np.exp(1j * 2 * math.pi * m * (c - (k - 1) / 2) / k) for c in range(k)] for m in ORDERS]
    )


def main():
    cells = json.loads((OUT / "layout.json").read_text())["cells"]
    steps_all = (20, 40, 60, 80, 100, 120, 140, 160, 175)
    mixed = pick_pairs(cells, True, 0, steps=steps_all)
    pp = pick_pairs(cells, False, 0, kind="pillar", steps=steps_all)
    hh = pick_pairs(cells, False, 0, kind="hole", steps=(40, 60, 100, 140, 175))
    study = [("pillar|hole", p) for p in mixed] + [("pillar|pillar", p) for p in pp] + [("hole|hole", p) for p in hh]
    singles = {}
    for _lab, (a, b) in study:
        singles.setdefault(a, None)
        singles.setdefault(b, None)
    with ProcessPoolExecutor(max_workers=4) as ex:
        for s, t in zip(singles, ex.map(isolated_txx, [k for k, _ in singles], [r for _, r in singles]), strict=True):
            singles[s] = t
        exact = list(ex.map(strip_job, [([a, a, b, b], STRIP_HARM) for _l, (a, b) in study]))
    amat = dft_matrix()
    inv = np.linalg.inv(amat)
    samples, strips = [], []
    for (label, (a, b)), orders in zip(study, exact, strict=True):
        ta, tb = singles[a], singles[b]
        cell_T = [ta, ta, tb, tb]
        nbr_T = [tb, tb, ta, ta]  # the one differing neighbour of each cell
        strips.append(dict(label=label, A=list(a), B=list(b), T_A=[ta.real, ta.imag], T_B=[tb.real, tb.imag]))
        for pol_name, pol_idx in (("x_normal", 0), ("y_parallel", 1)):
            e_orders = np.array([orders[m][pol_idx] for m in ORDERS])
            e_cells = inv @ e_orders
            for c in range(4):
                eps = e_cells[c] / cell_T[c] - 1.0
                d = nbr_T[c] / cell_T[c] - 1.0
                step = float(np.degrees(np.angle(nbr_T[c] / cell_T[c])))
                samples.append(
                    dict(
                        label=label,
                        pol=pol_name,
                        cell=c,
                        kind=(a if c < 2 else b)[0],
                        nbr_kind=(b if c < 2 else a)[0],
                        step_deg=step,
                        amp_ratio=float(abs(nbr_T[c]) / abs(cell_T[c])),
                        eps=[eps.real, eps.imag],
                        kappa=[(eps / d).real, (eps / d).imag],
                        eps_abs=float(abs(eps)),
                    )
                )
    (OUT / "coupling_samples.json").write_text(
        json.dumps(dict(strips=strips, samples=samples, orders=list(ORDERS)), indent=1)
    )
    print(len(samples), "samples written")


if __name__ == "__main__":
    main()
