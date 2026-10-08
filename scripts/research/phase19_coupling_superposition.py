"""Phase 19 D-bis: does neighbour coupling superpose?  Test with 2 and 4 differing neighbours.

    .venv/bin/python scripts/research/phase19_coupling_superposition.py

The 1 x 4 strips ([A A B B]) give the coupling kappa of a cell to ONE differing neighbour.  A lens cell has up to four.
Here the same (A, B) pairs are solved as an alternating strip [A B] (2 differing x-neighbours) and a 2 x 2 checkerboard
(2 x- plus 2 y-neighbours) and the measured per-cell perturbation is compared with the linear-superposition prediction
N * kappa * (T_n / T_c - 1).  Writes results/transmitarray_dual30/coupling_superposition.json.
"""

import json
import math
import pathlib
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.dielectric_cell.supercell import CellSpec, solve_grid  # noqa: E402

OUT = ROOT / "results" / "transmitarray_dual30"
F, P, TOTAL, PLATE, EPS, TAND = 28e9, 5.0e-3, 25.0e-3, 0.6e-3, 2.45, 0.0185


def job(args):
    kind, a, b = args
    sa, sb = CellSpec(*a[:1], a[1] * 1e-3), CellSpec(*b[:1], b[1] * 1e-3)
    if kind == "alt":
        grid, harm, orders = [[sa], [sb]], (16, 8), [(0, 0), (1, 0)]
    else:
        grid, harm, orders = [[sa, sb], [sb, sa]], (16, 16), [(0, 0), (1, 0), (0, 1), (1, 1)]
    res = solve_grid(grid, P, TOTAL - PLATE, PLATE, EPS, TAND, True, F, 0.0, 0.0, harm, orders=orders)
    # x-pol in -> E_x out (p column, x row);  y-pol in -> E_y out (s column, y row)
    return {m: (complex(res["orders"][m][0, 1]), complex(res["orders"][m][1, 0])) for m in orders}


def invert(orders_vals, kind, pol_idx):
    """Per-cell equivalent patch fields from the measured orders (2-point / 2x2 DFT with the sinc patch factor)."""
    if kind == "alt":
        pos = [(-0.5,), (0.5,)]  # cell centres in units of the pitch
        ms = [(0, 0), (1, 0)]
        amat = np.array(
            [[np.sinc(m[0] / 2) / 2 * np.exp(1j * 2 * math.pi * m[0] * pc[0] / 2) for pc in pos] for m in ms]
        )
        cells = ["A", "B"]
    else:
        pos = [(-0.5, -0.5), (-0.5, 0.5), (0.5, -0.5), (0.5, 0.5)]  # (ix, iy) = (0,0)A,(0,1)B,(1,0)B,(1,1)A
        ms = [(0, 0), (1, 0), (0, 1), (1, 1)]
        amat = np.array(
            [
                [
                    np.sinc(m[0] / 2)
                    * np.sinc(m[1] / 2)
                    / 4
                    * np.exp(1j * 2 * math.pi * (m[0] * pc[0] + m[1] * pc[1]) / 2)
                    for pc in pos
                ]
                for m in ms
            ]
        )
        cells = ["A", "B", "B", "A"]
    e = np.linalg.solve(amat, np.array([orders_vals[m][pol_idx] for m in ms]))
    return e, cells


def main():
    data = json.loads((OUT / "coupling_samples.json").read_text())
    strips, samples = data["strips"], data["samples"]
    jobs = []
    for s in strips:
        jobs.append(("alt", tuple(s["A"]), tuple(s["B"])))
        jobs.append(("chk", tuple(s["A"]), tuple(s["B"])))
    with ProcessPoolExecutor(max_workers=4) as ex:
        res = list(ex.map(job, jobs))
    # the sample list is ordered strip-major, then pol (x_normal, y_parallel), then cell 0..3
    per = 8
    out_rows = []
    for i, s in enumerate(strips):
        blk = samples[i * per : (i + 1) * per]
        ta, tb = complex(*s["T_A"]), complex(*s["T_B"])
        kap = {}
        for q in blk:
            kap.setdefault(q["pol"], {"A": [], "B": []})["A" if q["cell"] < 2 else "B"].append(complex(*q["kappa"]))
        kap = {p: {c: complex(np.mean(v)) for c, v in d.items()} for p, d in kap.items()}
        step = float(np.degrees(np.angle(tb / ta)))
        rec = dict(label=s["label"], A=s["A"], B=s["B"], step_deg=step, cases={})
        for kind, r in (("alt", res[2 * i]), ("chk", res[2 * i + 1])):
            for pol_name, pol_idx in (("x_normal", 0), ("y_parallel", 1)):
                e, cells = invert(r, kind, pol_idx)
                for e_c, c in zip(e, cells, strict=True):
                    t_c, t_n = (ta, tb) if c == "A" else (tb, ta)
                    ratio = t_n / t_c - 1.0
                    if kind == "alt":
                        pred = 2 * kap[pol_name][c] * ratio
                        n_diff = 2
                    else:
                        # x-pol checkerboard: x-neighbours across E-normal interfaces, y-neighbours across E-parallel
                        if pol_name != "x_normal":
                            continue
                        pred = 2 * kap["x_normal"][c] * ratio + 2 * kap["y_parallel"][c] * ratio
                        n_diff = 4
                    eps = e_c / t_c - 1.0
                    rec["cases"].setdefault(f"{kind}:{pol_name}", []).append(
                        dict(
                            cell=c,
                            n_diff=n_diff,
                            eps=[eps.real, eps.imag],
                            pred=[pred.real, pred.imag],
                            eps_abs=float(abs(eps)),
                            pred_abs=float(abs(pred)),
                            kappa_eff_per_neighbour=[(eps / (n_diff * ratio)).real, (eps / (n_diff * ratio)).imag],
                            single_kappa=[kap[pol_name][c].real, kap[pol_name][c].imag] if kind == "alt" else None,
                        )
                    )
        out_rows.append(rec)
    (OUT / "coupling_superposition.json").write_text(json.dumps(dict(pairs=out_rows), indent=1))
    # headline: per-neighbour effective kappa vs single-neighbour kappa
    for key in ("alt:x_normal", "alt:y_parallel", "chk:x_normal"):
        ratios, ke, angs = [], [], []
        for r in out_rows:
            for c in r["cases"].get(key, []):
                if abs(r["step_deg"]) < 45:
                    continue
                ratios.append(c["eps_abs"] / max(c["pred_abs"], 1e-12))
                k = complex(*c["kappa_eff_per_neighbour"])
                ke.append(abs(k))
                angs.append(np.degrees(np.angle(k)))
        print(
            key,
            "n",
            len(ratios),
            "median |eps_meas|/|eps_pred|",
            round(float(np.median(ratios)), 3),
            "median |kappa_eff|",
            round(float(np.median(ke)), 3),
            "median arg",
            round(float(np.median(angs)), 1),
        )


if __name__ == "__main__":
    main()
