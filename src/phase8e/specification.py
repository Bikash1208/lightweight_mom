"""Sec. 24: CSV/JSON output writers."""

import csv
import json
from pathlib import Path

import numpy as np


def _write_json(doc, path):
    Path(path).write_text(json.dumps(doc, indent=1, default=_json_default) + "\n")
    return path


def _json_default(o):
    if isinstance(o, complex):
        return {"re": o.real, "im": o.imag}
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.bool_):
        return bool(o)
    return str(o)


def write_formulation_candidates_json(provenance_registry, control_responses, architecture_results, path):
    doc = {}
    for letter, prov in provenance_registry.items():
        cr = control_responses.get(letter)
        arch = architecture_results.get(letter)
        doc[letter] = {
            "formulation_name": prov.formulation_name,
            "mathematical_reference": prov.mathematical_reference,
            "assumptions": list(prov.assumptions),
            "architecture_requirements": list(prov.architecture_requirements),
            "applicability": prov.applicability,
            "incompatibilities": list(prov.incompatibilities),
            "implementation_status": prov.implementation_status,
            "required_quantity": cr.required_quantity if cr else None,
            "value_kind": cr.value_kind if cr else None,
            "architecture_status": arch.status if arch else None,
            "architecture_evidence": arch.evidence if arch else None,
        }
    return _write_json(doc, path)


def write_formulation_comparison_csv(rows, path):
    """A factual comparison table, never a ranked one (mirrors phase8c.compatibility_matrix's
    'no winner' rule) -- columns only, no score column."""
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow([
            "formulation", "required_quantity", "architecture_status", "reference_plane_status",
            "family_exists", "mean_complex_error", "mean_amplitude_error", "mean_phase_error_deg",
            "final_status",
        ])
        for r in rows:
            writer.writerow(r)
    return path


def write_control_variable_analysis_json(decisions, path):
    doc = {}
    for d in decisions:
        doc[d.formulation] = {
            "decision_trail": [list(step) for step in d.decision_trail],
            "gate_checklist": d.gate_checklist,
            "geometry_inversion_permitted": d.geometry_inversion_permitted,
            "final_status": d.final_status,
        }
    return _write_json(doc, path)


def write_reference_spectrum_csv(kx_grid, ky_grid, power, path):
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["kx_per_m", "ky_per_m", "power_linear"])
        for kx, ky, p in zip(kx_grid.ravel(), ky_grid.ravel(), power.ravel()):
            writer.writerow([repr(float(kx)), repr(float(ky)), repr(float(p))])
    return path


def write_object_spectrum_csv(object_spectrum, path):
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["theta_deg", "phi_deg", "kx_per_m", "ky_per_m", "kz_real_per_m", "kz_imag_per_m", "is_propagating"])
        for p in object_spectrum.beams:
            writer.writerow([p.theta_deg, p.phi_deg, repr(p.kx_per_m), repr(p.ky_per_m),
                              repr(p.kz_per_m.real), repr(p.kz_per_m.imag), int(p.is_propagating)])
    return path


def write_required_response_map_csv(entries, path):
    """``entries``: iterable of ``(formulation, x_m, y_m, values)``, one per formulation, combined
    into a single file distinguished by the ``formulation`` column."""
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["x_mm", "y_mm", "formulation", "value_real_or_phase_deg", "value_imag_or_na"])
        for formulation, x_m, y_m, values in entries:
            x_flat, y_flat, v_flat = np.asarray(x_m).ravel(), np.asarray(y_m).ravel(), np.asarray(values).ravel()
            is_complex = np.iscomplexobj(v_flat)
            for x, y, v in zip(x_flat, y_flat, v_flat):
                if is_complex:
                    writer.writerow([repr(x * 1e3), repr(y * 1e3), formulation, repr(float(v.real)), repr(float(v.imag))])
                else:
                    writer.writerow([repr(x * 1e3), repr(y * 1e3), formulation, repr(float(v)), ""])
    return path


def write_library_response_map_csv(entries, path):
    """``entries``: iterable of ``(rows, component_name)``, combined into a single file."""
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["L_mm", "W_mm", "component", "real", "imag", "mag", "phase_deg"])
        for rows, component_name in entries:
            for r in rows:
                writer.writerow([repr(r.l_mm), repr(r.w_mm), component_name, repr(r.real), repr(r.imag), repr(r.mag), repr(r.phase_deg)])
    return path


def write_reconstruction_metrics_json(per_formulation_metrics, path):
    doc = {}
    for name, m in per_formulation_metrics.items():
        doc[name] = m
    return _write_json(doc, path)


def write_required_library_specification_json(period_assessment, path):
    doc = {
        "status": period_assessment.status,
        "registered_periods_mm": list(period_assessment.registered_periods_mm),
        "candidate_periods_mm": list(period_assessment.candidate_periods_mm),
        "registered_period_matches_candidates": period_assessment.registered_period_matches_candidates,
        "required_library_specification": period_assessment.required_library_specification,
    }
    return _write_json(doc, path)
