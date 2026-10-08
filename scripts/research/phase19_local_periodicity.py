"""Phase 19 item 5: how wrong is "each cell behaves as in an infinite array of itself" at pillar|hole boundaries?

    .venv/bin/python scripts/research/phase19_local_periodicity.py

1. Count 4-neighbour adjacencies in results/transmitarray_dual30/layout.json (pillar-pillar, pillar-hole, hole-hole).
2. For representative pairs from the layout, solve 1 x 4 strips [A A B B] (two boundaries per period) by RCWA at normal
   incidence and compare the exact diffraction orders (0, +-1, 2) with the local-periodicity prediction built from the
   two isolated-cell responses (uniform field per cell).  Control strips use same-family pairs with a comparable
   phase step.  Writes results/transmitarray_dual30/local_periodicity.json.
"""

import json
import math
import pathlib
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.dielectric_cell import OctagonPillarCell, solve_cell  # noqa: E402
from src.dielectric_cell.supercell import CellSpec, local_periodicity_prediction, solve_strip  # noqa: E402

OUT = ROOT / "results" / "transmitarray_dual30"
F, P, TOTAL, PLATE, EPS, TAND = 28e9, 5.0e-3, 25.0e-3, 0.6e-3, 2.45, 0.0185
N_HARM_CELL = 8
STRIP_HARM = (24, 8)
ORDERS = (0, 1, -1, 2)


def isolated_txx(kind, r_mm):
    if kind == "hole":
        c = OctagonPillarCell(P, r_mm * 1e-3, TOTAL, EPS, TAND, 0.0, 22.5, "hole")
    else:
        c = OctagonPillarCell(P, r_mm * 1e-3, TOTAL - PLATE, EPS, TAND, PLATE, 22.5, "pillar")
    return complex(solve_cell(c, F, 0.0, 0.0, N_HARM_CELL).t[1, 1])  # normal incidence: p = x


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
    return {m: complex(res["orders"][(m, 0)][0, 1]) for m in ORDERS}  # E_x for x(=p)-polarised incidence


def adjacency(cells):
    p = 5.0
    idx = {(math.floor(c["x_mm"] / p), math.floor(c["y_mm"] / p)): c for c in cells}
    kinds = {k: c["kind"] for k, c in idx.items()}
    counts = {"pillar-pillar": 0, "hole-pillar": 0, "hole-hole": 0}
    steps = []
    boundary = 0
    for (ix, iy), kd in kinds.items():
        for dx, dy in ((1, 0), (0, 1)):
            nb = (ix + dx, iy + dy)
            if nb not in kinds:
                boundary += 1
                continue
            key = "-".join(sorted((kd, kinds[nb])))
            counts[key] += 1
            if key == "hole-pillar":
                a, b = idx[(ix, iy)], idx[nb]
                steps.append(abs(((a["target_deg"] - b["target_deg"]) + 180) % 360 - 180))
    total = sum(counts.values())
    return dict(
        counts=counts,
        total_internal_edges=total,
        outline_edges=boundary,
        pillar_hole_fraction=counts["hole-pillar"] / total,
        pillar_hole_target_phase_step_deg=dict(
            median=float(np.median(steps)), p90=float(np.percentile(steps, 90)), max=float(np.max(steps))
        )
        if steps
        else None,
    )


def pick_pairs(cells, want_mixed, n, kind=None, steps=None):
    """Distinct x-adjacent pairs ((kind, R_mm), (kind, R_mm)).  ``steps``: pick the pair whose layout target-phase
    step is closest to each requested |step| [deg] instead of spreading over R."""
    p = 5.0
    idx = {(math.floor(c["x_mm"] / p), math.floor(c["y_mm"] / p)): c for c in cells}
    found = {}
    for (ix, iy), a in idx.items():
        b = idx.get((ix + 1, iy))
        if b is None:
            continue
        mixed = a["kind"] != b["kind"]
        if mixed != want_mixed or (kind is not None and a["kind"] != kind):
            continue
        step = abs(((a["target_deg"] - b["target_deg"]) + 180) % 360 - 180)
        found[((a["kind"], round(a["radius_mm"], 4)), (b["kind"], round(b["radius_mm"], 4)))] = step
    if steps is not None:
        out = []
        for st in steps:
            best = min((k for k in found if k not in out), key=lambda k: abs(found[k] - st), default=None)
            if best is not None:
                out.append(best)
        return out
    uniq = sorted(found, key=lambda t: (t[0][1], t[1][1]))
    return [uniq[int(i)] for i in np.linspace(0, len(uniq) - 1, min(n, len(uniq)))]


def main():
    cells = json.loads((OUT / "layout.json").read_text())["cells"]
    result = {"adjacency": adjacency(cells)}
    print(json.dumps(result["adjacency"], indent=1), flush=True)
    mixed = pick_pairs(cells, True, 0, steps=(30, 60, 90, 120, 150, 175))
    control = pick_pairs(cells, False, 0, kind="pillar", steps=(30, 60, 90, 120, 150, 175)) + pick_pairs(
        cells, False, 0, kind="hole", steps=(60, 120, 170)
    )
    study = [("pillar|hole", pr) for pr in mixed] + [("same-family", pr) for pr in control]
    singles = {}
    for _kind, (a, b) in study:
        for s in (a, b):
            singles.setdefault(s, None)
    with ProcessPoolExecutor(max_workers=4) as ex:
        for s, t in zip(singles, ex.map(isolated_txx, [k for k, _ in singles], [r for _, r in singles]), strict=True):
            singles[s] = t
        strips = [([a, a, b, b], STRIP_HARM) for _k, (a, b) in study]
        exact = list(ex.map(strip_job, strips))
    rows = []
    for (label, (a, b)), ex_orders in zip(study, exact, strict=True):
        pred = local_periodicity_prediction([singles[a], singles[a], singles[b], singles[b]], ORDERS)
        mean_amp = 0.5 * (abs(singles[a]) + abs(singles[b]))
        row = dict(
            label=label,
            A=list(a),
            B=list(b),
            isolated_phase_step_deg=float(np.degrees(np.angle(singles[b] / singles[a]))),
            isolated_abs=[abs(singles[a]), abs(singles[b])],
            order0_error_over_cell_amplitude=float(abs(ex_orders[0] - pred[0]) / mean_amp),
            orders={},
        )
        for m in ORDERS:
            e, pr = ex_orders[m], pred[m]
            row["orders"][str(m)] = dict(exact_abs=abs(e), pred_abs=abs(pr), abs_error=float(abs(e - pr)))
        rows.append(row)
        print(
            label,
            a,
            b,
            f"step {row['isolated_phase_step_deg']:.0f} deg",
            f"err/amp {row['order0_error_over_cell_amplitude']:.3f}",
            flush=True,
        )
    result["strips"] = rows
    # truncation sanity for the first mixed strip
    a, b = mixed[0]
    conv = {}
    with ProcessPoolExecutor(max_workers=3) as ex:
        hs = [(16, 8), (24, 8), (32, 8)]
        for h, r in zip(hs, ex.map(strip_job, [([a, a, b, b], h) for h in hs]), strict=True):
            conv[str(h)] = {str(m): [abs(r[m]), float(np.degrees(np.angle(r[m])))] for m in ORDERS}
    result["truncation_check_first_mixed_strip"] = conv
    mixed_rows = [r for r in rows if r["label"] == "pillar|hole"]
    ctrl_rows = [r for r in rows if r["label"] == "same-family"]

    # error vs |isolated phase step|, mixed vs same-family: least-squares line through the origin in step^1
    def slope(rs):
        st = np.abs([r["isolated_phase_step_deg"] for r in rs])
        er = np.array([r["order0_error_over_cell_amplitude"] for r in rs])
        return float(np.sum(st * er) / np.sum(st * st)) if len(rs) else None

    result["error_vs_step"] = dict(
        mixed_slope_per_deg=slope(mixed_rows),
        same_family_slope_per_deg=slope(ctrl_rows),
        note="order-0 |exact - local-periodicity| / mean isolated |T|, fitted as slope * |isolated phase step|",
    )
    result["summary"] = dict(
        mixed_order0_error_over_cell_amp_max=max(r["order0_error_over_cell_amplitude"] for r in mixed_rows),
        same_family_order0_error_over_cell_amp_max=max(r["order0_error_over_cell_amplitude"] for r in ctrl_rows),
        n_mixed=len(mixed_rows),
        n_same_family=len(ctrl_rows),
    )
    (OUT / "local_periodicity.json").write_text(json.dumps(result, indent=1))
    print(json.dumps(result["summary"], indent=1))


if __name__ == "__main__":
    main()
