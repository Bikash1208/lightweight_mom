"""Phase 23 D (revision 2): square vs hexagon vs hexagon + optimiser, with
  * Phase 16 pointing correction on every starting design,
  * the two cell states RE-SOLVED at 10.33 / 10.54 / 10.75 GHz (Phase 22 auto_states, same square side),
  * absolute pointing constrained to <= 0.3 deg of the requested direction,
  * a robust option: worst case over the three frequencies (cells re-solved at each) and, for all-ports
    mode, over 20 seeded random port-phase errors of +-10 deg.
Every number is on the held-out grid (theta 0.5 deg x phi 1 deg) the optimiser never saw.

    python scripts/research/phase23_compare_v2.py <out_dir> [--quick]
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

import phase23_compare as v1  # noqa: E402
from src.app.routers.hmsa import HmsaSynthesizePortsRequest, _run_synthesize_ports_job  # noqa: E402
from src.app.storage import app_data_dir  # noqa: E402
from src.hologram.sll_board import (  # noqa: E402
    HeldOutGrid,
    apply_levels,
    build_problem,
    cells_changed,
    evaluate,
    evaluate_phase_error_draws,
    optimise,
)
from src.hologram.sll_optimizer import OptimiserConfig  # noqa: E402

F0, PERIOD, EPS_R, TAND, H_MM = v1.F0, v1.PERIOD, v1.EPS_R, v1.TAND, v1.H_MM
FREQS = (10.33e9, 10.54e9, 10.75e9)
EXTRA = (10.33e9, 10.75e9)
BEAMS = v1.BEAMS
N_DRAWS, SIGMA = 20, 10.0


class _Job(v1._Job):
    pass


def synthesize(name, layout):
    body = {
        "formulation": "coded",
        "frequency_hz": F0,
        "period_mm": PERIOD,
        "layout": layout,
        "ports": [{"beams": [v1._beam(t, p)]} for t, p in BEAMS],
        "port_eps_r": EPS_R,
        "port_thickness_mm": H_MM,
        "port_tan_delta": TAND,
        "two_state": {"source": "auto"},
        "beta_mode": "model",
        "pointing_correction": True,
        "pointing_iterations": 2,
    }
    return _run_synthesize_ports_job(_Job(f"p23v2_{name}"), HmsaSynthesizePortsRequest(**body).model_dump())


def board_rows(prob, levels, alpha, grid, with_draws=False):
    """Held-out table at the three frequencies (cells re-solved at each)."""
    out = {}
    for f in FREQS:
        key = f"{f / 1e9:.2f}GHz"
        ev = evaluate(prob, levels, alpha, frequency_hz=f, grid=grid)
        out[key] = v1.summarise(ev)
        if with_draws and len(prob.ports) > 1:
            ds = evaluate_phase_error_draws(prob, levels, alpha, frequency_hz=f, grid=grid)
            worst = [max(b["sll_global_db"] for b in d["beams"]) for d in ds]
            out[key]["all_phase_error_draws"] = {
                "n": len(ds),
                "sigma_deg": SIGMA,
                "worst_beam_sll_db": {
                    "median": round(float(np.median(worst)), 2),
                    "worst": round(float(np.max(worst)), 2),
                    "best": round(float(np.min(worst)), 2),
                },
                "min_directivity_dbi": round(float(min(b["directivity_dbi"] for d in ds for b in d["beams"])), 2),
            }
    return out


def main(out_dir, quick=False):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cache = app_data_dir() / "auto_states_cache"
    grid = (
        HeldOutGrid(theta_step_deg=2.0, phi_step_deg=4.0, cut_step_u=0.01, patch_step_deg=0.25)
        if quick
        else HeldOutGrid()
    )
    base_cfg = dict(
        max_flips=120 if quick else 800,
        q_schedule=(6.0, 12.0, 24.0, 48.0),
        accepts_per_round=12,
        seed=0,
        alpha_step_deg=15.0,
        max_pointing_error_deg=0.3,
    )
    boards = {
        "square_300.8_gap0": {
            "preset": "grid",
            "rows": 2,
            "cols": 2,
            "cell_side_mm": 300.8,
            "gap_mm": 0.0,
            "snap_side_to_cells": True,
        },
        "square_243.2_gap0": {
            "preset": "grid",
            "rows": 2,
            "cols": 2,
            "cell_side_mm": 241.8,
            "gap_mm": 0.0,
            "snap_side_to_cells": True,
        },
        "hex4_side150_gap0": {"preset": "hex4", "side_mm": 150.0, "gap_mm": 0.0},
    }
    table, keep = {}, {}
    for name, layout in boards.items():
        t0 = time.time()
        res = synthesize(name, layout)
        prob = build_problem(
            res,
            F0,
            EPS_R,
            TAND,
            H_MM,
            PERIOD,
            [[b] for b in BEAMS],
            extra_frequencies_hz=EXTRA,
            states_cache_dir=cache,
            n_phase_draws=N_DRAWS,
            phase_error_deg=SIGMA,
        )
        t1 = time.time()
        row = {
            "area_mm2_per_port": round(prob.ports[0].region.area_m2 * 1e6, 1),
            "cells_per_port": [int(bp.owned.sum()) for bp in prob.ports],
            "patches": [p["patch_count"] for p in res["ports"]],
            "state_phase_deg": {
                f"{f / 1e9:.2f}GHz": [round(float(np.degrees(v)), 2) for v in prob.freq_psi[f]] for f in FREQS
            },
            "frequencies": board_rows(prob, None, None, grid, with_draws=True),
        }
        table[name] = row
        keep[name] = (res, prob)
        print(f"[{name}] built {t1 - t0:.0f}s evaluated {time.time() - t1:.0f}s", flush=True)
        (out / "comparison_v2.json").write_text(json.dumps(table, indent=1))

    res, prob = keep["hex4_side150_gap0"]
    for mode in ("single", "all"):
        for robust in (False, True):
            label = f"hex4+opt_{mode}" + ("_robust" if robust else "")
            t0 = time.time()
            cfg = OptimiserConfig(**base_cfg, optimise_alpha=(mode == "all"))
            levels, alpha, runs = optimise(prob, mode, cfg, robust=robust)
            t1 = time.time()
            ev0 = evaluate(prob, levels, alpha, grid=grid)
            new = apply_levels(res, prob, levels, alpha, evaluation=ev0)
            row = {
                "cells_changed": cells_changed(prob, levels),
                "alpha_deg": {int(k): round(float(v), 1) for k, v in alpha.items()},
                "optimiser_instances": len(runs[0].before) if runs else 0,
                "frequencies": board_rows(prob, levels, alpha, grid, with_draws=(mode == "all")),
            }
            table[label] = row
            keep[label] = (new, levels, alpha)
            print(f"[{label}] optimised {t1 - t0:.0f}s evaluated {time.time() - t1:.0f}s", flush=True)
            (out / "comparison_v2.json").write_text(json.dumps(table, indent=1))
    if not quick:
        notes = {
            "single": "Cell states optimised for each port excited on its own (one port at a time in HFSS).",
            "all": "Cell states optimised for ALL FOUR ports excited together.",
        }
        tags = {
            "square_243.2_gap0": "four_port_square243_gap0",
            "hex4_side150_gap0": "four_port_hex150_gap0",
            "hex4+opt_single": "four_port_hex150_opt_single",
            "hex4+opt_all": "four_port_hex150_opt_all",
            "hex4+opt_single_robust": "four_port_hex150_opt_single_robust",
            "hex4+opt_all_robust": "four_port_hex150_opt_all_robust",
        }
        for label, tag in tags.items():
            result = keep[label][0]
            alpha = keep[label][2] if len(keep[label]) > 2 else None
            is_all = "opt_all" in label
            extra = (
                notes["all" if is_all else "single"]
                if label.startswith("hex4+opt")
                else "Starting design (Phase 16 pointing-corrected, no sidelobe optimisation)."
            )
            if "robust" in label:
                extra += " ROBUST: the worst case over 10.33 / 10.54 / 10.75 GHz (cells re-solved at each)" + (
                    " and 20 seeded random port-phase errors of +-10 deg." if is_all else "."
                )
            if is_all:
                extra += (
                    "\n\n**All-ports SLL needs ONE phase-locked source split four ways** (the cell states and the port phases "
                    "above were optimised for the four ports radiating coherently with fixed relative phases). If the "
                    "four ports are driven by independent signals the all-ports sidelobe figures do not apply: use the "
                    "single-port metric (each port on its own) instead."
                )
            row = table.get(label) or table[label]
            v1.export_package(ROOT / "hfss_exports" / tag, tag, result, _flat(row), alpha if is_all else None, extra)
            table[label]["step"] = str(ROOT / "hfss_exports" / tag / "board.step")
        (out / "comparison_v2.json").write_text(json.dumps(table, indent=1))
    print(json.dumps(table, indent=1))


def _flat(row):
    """The design-frequency rows in `export_package`'s expected shape."""
    return row["frequencies"]["10.54GHz"] if "frequencies" in row else row


if __name__ == "__main__":
    main(sys.argv[1], quick="--quick" in sys.argv)
