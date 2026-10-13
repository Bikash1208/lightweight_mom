"""Build the 28 GHz nylon octagon cell library (Phase 19 Part B).

    .venv/bin/python examples/build_nylon_cell_library_28ghz.py

Writes results/nylon_cells_28ghz/library.csv (complex t, r for every cell/angle) and summary.json.
Model-only: no HFSS comparison exists yet (NOT_AVAILABLE_FROM_CURRENT_DATA).
"""

import csv
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from src.dielectric_cell.library import (  # noqa: E402
    LibraryParams,
    best_window,
    normal_incidence_table,
    order_by_fill,
    sweep,
)

OUT = pathlib.Path(__file__).resolve().parents[1] / "results" / "nylon_cells_28ghz"

PARAMS = LibraryParams(
    freq_hz=28e9,
    period_m=5.0e-3,
    total_thickness_m=25.0e-3,
    base_plate_m=1.0e-3,
    eps_r=2.45,
    tan_delta=0.0185,
    min_feature_m=0.6e-3,
    n_per_family=28,
    thetas_deg=(0.0, 10.0, 20.0, 25.0),
    phis_deg=(0.0, 45.0),
    n_harm=8,
)


def main(argv=None):
    """No arguments: the Part B library (1 mm plate, blind pockets) -> library.csv.
    ``--plate-mm X [--through]``: the same sweep for another plate / through-hole pockets -> library_<tag>.csv."""
    import argparse
    import dataclasses

    ap = argparse.ArgumentParser()
    ap.add_argument("--plate-mm", type=float, default=None)
    ap.add_argument("--through", action="store_true")
    args = ap.parse_args(argv)
    params, tag = PARAMS, ""
    if args.plate_mm is not None or args.through:
        plate = args.plate_mm if args.plate_mm is not None else 1.0
        params = dataclasses.replace(PARAMS, base_plate_m=plate * 1e-3, through_hole=args.through)
        tag = f"_plate{plate:g}mm_{'through' if args.through else 'blind'}"
    OUT.mkdir(parents=True, exist_ok=True)
    rows = sweep(params)
    names = ["ss", "sp", "ps", "pp"]
    with open(OUT / f"library{tag}.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        head = ["kind", "radius_mm", "fill", "theta_deg", "phi_deg"]
        head += [f"t_{n}_{part}" for n in names for part in ("re", "im")]
        head += [f"r_{n}_{part}" for n in names for part in ("re", "im")]
        w.writerow(head)
        for kind, r, fill, th, ph, t, rr in rows:
            vals = [kind, f"{r * 1e3:.6f}", f"{fill:.6f}", th, ph]
            for m in (t, rr):
                for z in m.reshape(-1):
                    vals += [f"{z.real:.9e}", f"{z.imag:.9e}"]
            w.writerow(vals)
    ordered = order_by_fill(normal_incidence_table(rows))
    win = best_window(ordered)
    summary = {
        "params": {k: (list(v) if isinstance(v, tuple) else v) for k, v in params.__dict__.items()},
        "n_cells": len(ordered),
        "full_library_phase_span_deg": abs(ordered[-1][3] - ordered[0][3]),
        "full_library_loss_db_range": [min(x[4] for x in ordered), max(x[4] for x in ordered)],
        "best_360_window": None
        if win is None
        else {
            "fill_range": [ordered[win[0]][2], ordered[win[1]][2]],
            "span_deg": win[3],
            "worst_case_loss_db": win[2],
            "n_cells": win[1] - win[0] + 1,
        },
        "hfss_comparison": "NOT_AVAILABLE_FROM_CURRENT_DATA",
    }
    (OUT / f"summary{tag}.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
