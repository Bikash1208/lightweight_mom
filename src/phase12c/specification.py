"""Phase 12C output writers: deterministic JSON/CSV under ``results/phase12c_grounded_locality/``
and ``results/phase12c_hfss_export/`` (no timestamps, matching every earlier phase)."""

import csv
import json
from dataclasses import asdict
from pathlib import Path


def _write_json(obj, path):
    Path(path).write_text(json.dumps(obj, indent=2, default=str))
    return path


def write_task_a_json(correctness, convergence_series, feasibility_points, path):
    obj = {
        "correctness_proof": asdict(correctness),
        "mesh_convergence_series": [asdict(c) for c in convergence_series],
        "largest_feasible_size": [asdict(p) for p in feasibility_points],
    }
    return _write_json(obj, path)


def write_task_b_json(curve, centers, gradient_results, strong_contrast, periodic_results,
                       periodic_grid_mm, periodic_span_deg, periodic_is_full_2pi, gate_result, path):
    obj = {
        "reference_curve": [asdict(p) for p in curve],
        "centers": {
            "off_resonance": asdict(centers.off_resonance),
            "near_resonance": asdict(centers.near_resonance),
            "at_resonance": asdict(centers.at_resonance),
        },
        "gradient_cases": [asdict(r) for r in gradient_results],
        "strong_contrast": asdict(strong_contrast) if strong_contrast is not None else None,
        "periodic_supercell": {
            "cells": [asdict(r) for r in periodic_results],
            "grid_mm": periodic_grid_mm,
            "achieved_span_deg": periodic_span_deg,
            "is_full_2pi": periodic_is_full_2pi,
        },
        "gate": {
            "status": gate_result.status, "reason": gate_result.reason,
            "max_step_within_tolerance_deg": gate_result.max_step_within_tolerance_deg,
            "phase9_classification": gate_result.phase9_classification,
            "per_case": [{"name": n, "achieved_step_deg": s, "phase_error_deg": p,
                          "magnitude_error_db": m, "within_tolerance": w}
                         for n, s, p, m, w in gate_result.per_case],
        },
    }
    return _write_json(obj, path)


def write_task_c_json(reports, path):
    """``reports``: dict name -> ``coupling_study.CouplingSeriesReport``."""
    obj = {}
    for name, r in reports.items():
        obj[name] = {
            "points": [asdict(p) for p in r.points],
            "increments": [asdict(i) for i in r.increments],
            "convergence": asdict(r.convergence),
            "decay_observation": asdict(r.decay_observation),
            "decay_fit": asdict(r.decay_fit),
        }
    return _write_json(obj, path)


def write_task_d_json(c2_c3, phase9_strong, case_d, coupling_mesh, coupling_points, spot_gate, path):
    obj = {
        "phase10_c2_c3": asdict(c2_c3),
        "phase9_strong_contrast": asdict(phase9_strong),
        "phase11_case_d": asdict(case_d),
        "phase11_coupling_radius": {"mesh_cells_per_subcell": coupling_mesh,
                                      "points": [asdict(p) for p in coupling_points]},
        "gate": {"status": spot_gate.status, "reason": spot_gate.reason,
                 "findings": [{"case": n, "verdict": v, "detail": d} for n, v, d in spot_gate.findings]},
    }
    return _write_json(obj, path)


def write_task_d_spot_check_csv(c2_c3, phase9_strong, case_d, path):
    rows = [
        ("phase10_c2_c3", c2_c3.mesh_cells_per_side, "mask_identical", c2_c3.mask_identical,
         "relative_difference", c2_c3.relative_difference),
        ("phase9_strong_contrast", phase9_strong.mesh_cells_per_subcell, "relative_complex_error",
         phase9_strong.relative_complex_error, "earlier_value", phase9_strong.earlier_relative_complex_error),
        ("phase11_case_d", case_d.mesh_cells_per_subcell, "relative_difference", case_d.relative_difference,
         "substantial", case_d.substantial),
    ]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["case", "mesh_cells_per_subcell", "metric_1_name", "metric_1_value", "metric_2_name", "metric_2_value"])
        w.writerows(rows)
    return path


def write_gate_status_json(task_a_status, task_b_gate, task_c_status, task_d_gate, hfss_status, path):
    obj = {
        "task_a_status": task_a_status,
        "task_b_local_periodicity": {"status": task_b_gate.status, "reason": task_b_gate.reason},
        "task_c_status": task_c_status,
        "task_d_spot_check": {"status": task_d_gate.status, "reason": task_d_gate.reason},
        "hfss_export": hfss_status,
        "hfss_status": "HFSS_NOT_AVAILABLE",
        "geometry_inversion": "GEOMETRY_INVERSION_LOCKED",
        "step_export": "STEP_EXPORT_LOCKED",
    }
    return _write_json(obj, path)
