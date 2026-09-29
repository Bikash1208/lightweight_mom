"""Phase 12D output writers: deterministic JSON/CSV under ``results/phase12d_local_periodicity/``
(no timestamps, matching every earlier phase)."""

import json
from dataclasses import asdict
from pathlib import Path


def _write_json(obj, path):
    Path(path).write_text(json.dumps(obj, indent=2, default=str))
    return path


def write_task_a_json(curve, centers, step_cases, gate_result, path):
    obj = {
        "reference_curve": [asdict(p) for p in curve],
        "centers": {
            "off_resonance": asdict(centers.off_resonance),
            "near_resonance": asdict(centers.near_resonance),
            "at_resonance": asdict(centers.at_resonance),
        },
        "step_cases": [asdict(c) for c in step_cases],
        "gate": {
            "status": gate_result.status,
            "reason": gate_result.reason,
            "max_step_within_design_tolerance_deg": gate_result.max_step_within_design_tolerance_deg,
            "max_step_within_strict_tolerance_deg": gate_result.max_step_within_strict_tolerance_deg,
            "per_case": [
                {
                    "name": n,
                    "achieved_step_deg": s,
                    "phase_error_deg": p,
                    "magnitude_error_db": m,
                    "within_design_tolerance": wd,
                    "within_strict_tolerance": ws,
                }
                for n, s, p, m, wd, ws in gate_result.per_case
            ],
        },
    }
    return _write_json(obj, path)


def write_task_b_json(aperture_results, requested_sizes, gate_result, path):
    obj = {
        "requested_sizes": list(requested_sizes),
        "results": [
            {
                "n_cells": r.n_cells,
                "mesh_cells_per_subcell": r.mesh_cells_per_subcell,
                "floquet_order": r.floquet_order,
                "achieved_span_deg": r.achieved_span_deg,
                "is_full_2pi": r.is_full_2pi,
                "n_unknowns": r.n_unknowns,
                "runtime_s": r.runtime_s,
                "intended_order": list(r.intended_order),
                "fraction_intended": r.fraction_intended,
                "fraction_specular": r.fraction_specular,
                "total_propagating_fraction": r.total_propagating_fraction,
                "spurious_orders": [{"order": list(mn), "fraction": f} for mn, f in r.spurious_orders],
                "intended_order_phase_deg": r.intended_order_phase_deg,
                "per_cell_phase_errors": [asdict(p) for p in r.per_cell_phase_errors],
            }
            for r in aperture_results
        ],
        "gate": {
            "status": gate_result.status,
            "reason": gate_result.reason,
            "per_size": [
                {"n_cells": n, "fraction_intended": f, "max_spurious_fraction": s, "adequate": a, "attempted": at}
                for n, f, s, a, at in gate_result.per_size
            ],
        },
    }
    return _write_json(obj, path)


def write_task_c_json(comparisons, provenance, geometric_identity, gate_result, path):
    obj = {
        "provenance": provenance,
        "comparisons": [
            {
                "case_id": c.case_id,
                "l_mm": c.l_mm,
                "mean_dphase_mode1_deg": c.mean_dphase_deg,
                "max_abs_dphase_mode1_deg": c.max_abs_dphase_deg,
                "max_dmag_mode1_db": c.max_dmag_db,
                "max_crosspol_mag": c.max_crosspol_mag,
                "max_top_to_bottom_mag": c.max_top_to_bottom_mag,
                "meets_tolerance": c.meets_tolerance,
                "frequency_ghz": list(c.frequency_ghz),
                "dphase_mode1_deg": list(c.dphase_mode1_deg),
            }
            for c in comparisons
        ],
        "geometric_identity": [asdict(g) for g in geometric_identity],
        "gate": {
            "status": gate_result.status,
            "reason": gate_result.reason,
            "cases_with_data": list(gate_result.cases_with_data),
            "cases_pending": list(gate_result.cases_pending),
        },
    }
    return _write_json(obj, path)


def write_gate_status_json(task_a_gate, task_b_gate, task_c_gate, path):
    obj = {
        "task_a_step_limit": {"status": task_a_gate.status, "reason": task_a_gate.reason},
        "task_b_aperture": {"status": task_b_gate.status, "reason": task_b_gate.reason},
        "task_c_hfss": {"status": task_c_gate.status, "reason": task_c_gate.reason},
        "hfss_status": "HFSS_NOT_AVAILABLE_FOR_SUPERCELLS",
        "geometry_inversion": "GEOMETRY_INVERSION_LOCKED",
        "step_export": "STEP_EXPORT_LOCKED",
    }
    return _write_json(obj, path)
