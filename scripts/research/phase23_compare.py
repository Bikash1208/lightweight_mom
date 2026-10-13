"""Phase 23 D: square vs hexagon vs hexagon + sidelobe optimiser, 4-port 1-bit boards, 10.54 GHz.

Every SLL / directivity number is evaluated on the held-out fine grid (theta 0.5 deg x phi 1 deg) that the
optimiser never saw; frequency checks re-evaluate (never re-optimise) the same states at +-2 %.

    python scripts/research/phase23_compare.py <out_dir> [--quick]
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.app.routers.hmsa import HmsaSynthesizePortsRequest, _run_synthesize_ports_job  # noqa: E402
from src.hologram.sll_board import (  # noqa: E402
    HeldOutGrid,
    apply_levels,
    build_problem,
    cells_changed,
    evaluate,
    optimise,
    optimise_alpha_only,
)
from src.hologram.sll_optimizer import OptimiserConfig  # noqa: E402

F0 = 10.54e9
PERIOD = 6.4
EPS_R, TAND, H_MM = 2.2, 0.0009, 1.575
BEAMS = [(30.0, 45.0), (30.0, 135.0), (30.0, 225.0), (30.0, 315.0)]


class _Job:
    def __init__(self, job_id):
        self.id = job_id
        self.progress = 0.0
        self.progress_message = ""


def _beam(theta, phi):
    return {
        "theta_deg": theta,
        "phi_deg": phi,
        "weight_real": 1.0,
        "weight_imag": 0.0,
        "initial_phase_deg": 0.0,
        "desired_polarization": "rhcp",
        "oam_charge": 0,
    }


def synthesize(name, layout):
    body = {
        "formulation": "coded",
        "frequency_hz": F0,
        "period_mm": PERIOD,
        "layout": layout,
        "ports": [{"beams": [_beam(t, p)]} for t, p in BEAMS],
        "port_eps_r": EPS_R,
        "port_thickness_mm": H_MM,
        "port_tan_delta": TAND,
        "two_state": {"source": "auto"},
        "beta_mode": "model",
    }
    return _run_synthesize_ports_job(_Job(f"p23cmp_{name}"), HmsaSynthesizePortsRequest(**body).model_dump())


def summarise(ev):
    """Compact per-scenario table rows from `evaluate` output."""
    rows = {}
    for name, sc in ev.items():
        bs = sc["beams"]
        rows[name] = {
            "directivity_dbi": [round(b["directivity_dbi"], 2) for b in bs],
            "pointing_error_deg": [round(b["pointing_error_deg"], 2) for b in bs],
            "theta_phi_deg": [(round(b["theta_deg"], 2), round(b["phi_deg"], 2)) for b in bs],
            "sll_global_db": [None if b["sll_global_db"] is None else round(b["sll_global_db"], 2) for b in bs],
            "sll_xcut_db": [None if b.get("sll_xcut_db") is None else round(b["sll_xcut_db"], 2) for b in bs],
            "sll_ycut_db": [None if b.get("sll_ycut_db") is None else round(b["sll_ycut_db"], 2) for b in bs],
            "xpd_db": [round(b["xpd_db"], 1) for b in bs],
            "level_spread_db": round(sc["level_spread_db"], 2),
        }
    return rows


def export_package(directory, name, result, table_row, alpha=None, notes=""):
    """Writes the exact STEP + layout files + AEDT script + README of one board into `directory`."""
    from src.export.aedt_script import generate_aedt_script
    from src.export.board_readme import build_board_readme
    from src.export.exact_board import (
        exact_board_from_ports_result,
        write_exact_step,
        write_feeds_csv,
        write_layout_csv,
        write_layout_dxf,
        write_layout_json,
    )

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    board = exact_board_from_ports_result(result, PERIOD, H_MM)
    write_exact_step(board, directory / "board.step")
    write_layout_csv(board, directory / "layout.csv")
    write_feeds_csv(board, directory / "feeds.csv")
    write_layout_json(board, directory / "layout.json")
    write_layout_dxf(board, directory / "layout.dxf")
    from src.app.routers.export import _write_board_png

    _write_board_png(board, directory / "layout.png")
    (directory / "hfss_build_board.py").write_text(generate_aedt_script(board, name, EPS_R, TAND, F0))
    files = ["board.step", "layout.csv", "feeds.csv", "layout.json", "layout.dxf", "layout.png", "hfss_build_board.py"]
    readme = build_board_readme(board, result, name, EPS_R, TAND, F0, files)
    extra = ["", "## Phase 23 model predictions on the held-out fine grid (theta 0.5 deg x phi 1 deg)", ""]
    extra += ["Model-only; none of these numbers has been compared with HFSS.", ""]

    def plain(v):
        return ", ".join(str(x) for x in v) if isinstance(v, list) else v

    for scen in ("port1", "port2", "port3", "port4", "all"):
        row = table_row.get(scen)
        if not row:
            continue
        extra.append(
            "**%s**: directivity %s dBi, pointing error %s deg, co-pol SLL %s dB, x-cut %s dB, y-cut %s dB, XPD %s dB"
            % (
                scen,
                row["directivity_dbi"],
                row["pointing_error_deg"],
                row["sll_global_db"],
                row["sll_xcut_db"],
                row["sll_ycut_db"],
                row["xpd_db"],
            )
        )
        extra.append("")
    if alpha is not None:
        extra += [
            "## Port phases for the all-ports case (HFSS > Excitations > Edit Sources)",
            "",
            "| Port | Magnitude | Phase |",
            "|---|---|---|",
        ]
        for pid in sorted(alpha):
            extra.append("| %d | 1 W | %.1f deg |" % (pid + 1, alpha[pid]))
        extra += [
            "",
            "Port 1 is the 0 deg reference. For a one-port-at-a-time run set that port to 1 W and the others to 0 W "
            "(the cell states were optimised for the mode named in this file).",
            "",
        ]
    if notes:
        extra += [notes, ""]
    (directory / "README.md").write_text(readme + "\n".join(extra))
    return str(directory / "board.step")


def main(out_dir, quick=False):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    grid = (
        HeldOutGrid(theta_step_deg=2.0, phi_step_deg=4.0, cut_step_u=0.01, patch_step_deg=0.25)
        if quick
        else HeldOutGrid()
    )
    cfg = OptimiserConfig(
        max_flips=120 if quick else 800,
        q_schedule=(6.0, 12.0, 24.0, 48.0),
        accepts_per_round=12,
        anneal_steps=0,
        seed=0,
        alpha_step_deg=15.0,
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
    table = {}
    results = {}
    for name, layout in boards.items():
        t0 = time.time()
        res = synthesize(name, layout)
        prob = build_problem(res, F0, EPS_R, TAND, H_MM, PERIOD, [[b] for b in BEAMS])
        t1 = time.time()
        ev = evaluate(prob, grid=grid)
        areas_mm2 = [bp.region.area_m2 * 1e6 for bp in prob.ports]
        table[name] = {
            "area_mm2_per_port": round(areas_mm2[0], 1),
            "cells_per_port": [int(bp.owned.sum()) for bp in prob.ports],
            "patches": [p["patch_count"] for p in res["ports"]],
            **summarise(ev),
        }
        print(f"[{name}] built {t1 - t0:.0f}s, evaluated {time.time() - t1:.0f}s", flush=True)
        results[name] = (res, prob)
        (out / "comparison.json").write_text(json.dumps(table, indent=1))
        if not quick and name != "square_300.8_gap0":
            tag = {"square_243.2_gap0": "four_port_square243_gap0", "hex4_side150_gap0": "four_port_hex150_gap0"}[name]
            table[name]["step"] = export_package(
                ROOT / "hfss_exports" / tag, tag, res, table[name], None, "Starting design (no sidelobe optimisation)."
            )

    res, prob = results["hex4_side150_gap0"]
    hex_start = {p.port_id: p.state for p in prob.opt_ports}
    for mode, label in (("single", "hex4+opt_single"), ("all", "hex4+opt_all")):
        t0 = time.time()
        levels, alpha, runs = optimise(prob, mode, cfg)
        t1 = time.time()
        ev = evaluate(prob, levels, alpha, grid=grid)
        table[label] = {
            "cells_changed": cells_changed(prob, levels),
            "alpha_deg": {int(k): round(float(v), 1) for k, v in alpha.items()},
            "optimisation_grid_report": {
                k: {
                    "sll_db": round(v.sll_db, 2),
                    "beam_peak_db": {int(b): round(x, 2) for b, x in v.beam_peak_db.items()},
                }
                for r in runs
                for k, v in r.after.items()
            },
            **summarise(ev),
        }
        # frequency re-evaluation (no re-optimisation)
        table[label]["frequency_check"] = {
            f"{df:+.0%}": summarise(evaluate(prob, levels, alpha, frequency_hz=F0 * (1 + df), grid=grid))
            for df in (-0.02, 0.02)
        }
        print(f"[{label}] optimised {t1 - t0:.0f}s, evaluated {time.time() - t1:.0f}s", flush=True)
        new = apply_levels(res, prob, levels, alpha, evaluation=ev)
        results[label] = (new, levels, alpha)
        if not quick:
            tag = {"hex4+opt_single": "four_port_hex150_opt_single", "hex4+opt_all": "four_port_hex150_opt_all"}[label]
            note = "Cell states optimised for: " + (
                "each port excited on its own."
                if mode == "single"
                else "all four ports excited together (cells and port phases). Set the phases above in HFSS."
            )
            table[label]["step"] = export_package(
                ROOT / "hfss_exports" / tag, tag, new, table[label], alpha if mode == "all" else None, note
            )
        (out / "comparison.json").write_text(json.dumps(table, indent=1))

    # fair baseline for the all-ports column: port phases only, cells untouched
    alpha_only = optimise_alpha_only(prob, cfg)
    ev = evaluate(prob, hex_start, alpha_only, grid=grid, scenarios=["all"])
    table["hex4+alpha_only"] = {
        "alpha_deg": {int(k): round(float(v), 1) for k, v in alpha_only.items()},
        **summarise(ev),
    }
    (out / "comparison.json").write_text(json.dumps(table, indent=1))
    print(json.dumps(table, indent=1))


if __name__ == "__main__":
    main(sys.argv[1], quick="--quick" in sys.argv)
