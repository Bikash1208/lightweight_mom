"""Sec. 27: the required JSON output writers (Phase 11's required-outputs list is JSON-only, no
CSVs mandated)."""

import json
from pathlib import Path


def _write_json(doc, path):
    Path(path).write_text(json.dumps(doc, indent=1) + "\n")
    return path


def _c(z):
    return None if z is None else [z.real, z.imag]


def write_coupling_radius_results_json(results_by_center, path):
    doc = {
        name: {"points": [{"n": p.n, "response": _c(p.response), "n_unknowns": p.n_unknowns} for p in points]}
        for name, points in results_by_center.items()
    }
    return _write_json(doc, path)


def write_incremental_shell_results_json(increments_by_center, path):
    doc = {
        name: [
            {
                "n_from": inc.n_from,
                "n_to": inc.n_to,
                "delta_complex": inc.delta_complex,
                "delta_relative": inc.delta_relative,
            }
            for inc in incs
        ]
        for name, incs in increments_by_center.items()
    }
    return _write_json(doc, path)


def write_directional_coupling_json(baseline, points, path):
    doc = {
        "baseline": _c(baseline),
        "points": [
            {
                "direction": p.direction,
                "distance": p.distance,
                "delta_abs": p.delta_abs,
                "delta_complex": _c(p.delta_complex),
            }
            for p in points
        ],
    }
    return _write_json(doc, path)


def write_diagonal_coupling_json(doc, path):
    return _write_json(doc, path)


def write_pattern_study_json(results, dependence_detected, comparisons, path):
    doc = {
        "results": {name: {"response": _c(r.response), "n_unknowns": r.n_unknowns} for name, r in results.items()},
        "pattern_dependence_detected": dependence_detected,
        "comparisons_relative_to_checkerboard": comparisons,
    }
    return _write_json(doc, path)


def write_frequency_study_json(points, path):
    doc = [
        {
            "frequency_hz": p.frequency_hz,
            "center_name": p.center_name,
            "responses_by_n": {str(n): _c(v) for n, v in p.responses_by_n.items()},
            "max_shell_increment_relative": p.max_shell_increment_relative,
        }
        for p in points
    ]
    return _write_json(doc, path)


def write_polarization_study_json(points, path):
    doc = [
        {
            "polarization": p.polarization,
            "responses_by_n": {str(n): _c(v) for n, v in p.responses_by_n.items()},
            "max_shell_increment_relative": p.max_shell_increment_relative,
        }
        for p in points
    ]
    return _write_json(doc, path)


def write_pairwise_interaction_json(r0, results, path):
    doc = {
        "baseline": _c(r0),
        "results": [
            {
                "pair_label": r.pair_label,
                "delta_a": _c(r.delta_a),
                "delta_b": _c(r.delta_b),
                "delta_ab_measured": _c(r.delta_ab_measured),
                "delta_ab_additive_prediction": _c(r.delta_ab_additive_prediction),
                "interaction_error": r.interaction_error,
                "interaction_error_relative": r.interaction_error_relative,
            }
            for r in results
        ],
    }
    return _write_json(doc, path)


def write_higher_order_interaction_json(result, path):
    doc = {
        "combo_label": result.combo_label,
        "measured": _c(result.measured),
        "pairwise_prediction": _c(result.pairwise_prediction),
        "residual_beyond_pairwise": result.residual_beyond_pairwise,
        "residual_relative": result.residual_relative,
    }
    return _write_json(doc, path)


def write_convergence_json(convergence_results, classification_result, path):
    doc = {
        "per_center": {
            cr.center_name: {
                "status": cr.status,
                "last_increment_relative": cr.last_increment_relative,
                "max_increment_relative": cr.max_increment_relative,
                "converged_at_n": cr.converged_at_n,
            }
            for cr in convergence_results
        },
        "classification": {
            "per_center_case": classification_result.per_center_case,
            "overall_case": classification_result.overall_case,
            "pattern_dependence_detected": classification_result.pattern_dependence_detected,
            "rationale": classification_result.rationale,
        },
        "noise_floor_relative_declared_before_evaluation": True,
    }
    return _write_json(doc, path)


def write_heterogeneous_reference_cases_json(doc, path):
    return _write_json(doc, path)


def write_beam_validation_json(doc, path):
    return _write_json(doc, path)


def write_phase11_summary_json(decision, path):
    doc = {
        "phase": "11",
        "heterogeneous_mom_status": decision.heterogeneous_mom_status,
        "coupling_status": decision.coupling_status,
        "neighborhood_status": decision.neighborhood_status,
        "pattern_status": decision.pattern_status,
        "reduction_status": decision.reduction_status,
        "inverse_status": decision.inverse_status,
        "hfss_status": decision.hfss_status,
        "geometry_inversion_status": decision.geometry_inversion_status,
        "step_export_status": decision.step_export_status,
        "overall_case": decision.overall_case,
    }
    return _write_json(doc, path)


def write_phase11_gate_status_json(decision, path):
    doc = {
        "phase": "11",
        "geometry_inversion_status": decision.geometry_inversion_status,
        "geometry_inversion_conditions": decision.geometry_inversion_conditions,
        "step_export_status": decision.step_export_status,
        "step_export_conditions": decision.step_export_conditions,
    }
    return _write_json(doc, path)
