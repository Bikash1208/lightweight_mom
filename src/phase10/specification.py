"""Sec. 45: the required CSV/JSON output writers, following ``phase9.specification``'s own
pattern (``repr()`` for every float so values round-trip exactly; ``NOT_AVAILABLE``/
``OUT_OF_VALIDATED_DOMAIN`` written literally rather than a blank field)."""

import csv
import json
from pathlib import Path

from .dataset import ALL_POSITIONS, RANDOM_SEED


def _write_json(doc, path):
    Path(path).write_text(json.dumps(doc, indent=1) + "\n")
    return path


# --------------------------------------------------------------------------- neighborhood_study/

def write_five_cell_cases_csv(study_results, path):
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["center_name", "contrast_label", "case_name", "description",
                    "px_real", "px_imag", "n_unknowns"])
        for name, r in study_results.items():
            for label, cases in r.five_cell_cases.items():
                for case_name, c in cases.items():
                    w.writerow([name, label, case_name, c.description, repr(c.px.real),
                                repr(c.px.imag), c.n_unknowns])
    return path


def write_diagonal_cases_csv(study_results, path):
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["center_name", "contrast_label", "case_name", "description",
                    "px_real", "px_imag", "n_unknowns"])
        for name, r in study_results.items():
            for label, cases in r.diagonal_cases.items():
                for case_name, c in cases.items():
                    w.writerow([name, label, case_name, c.description, repr(c.px.real),
                                repr(c.px.imag), c.n_unknowns])
    return path


def write_long_range_cases_csv(study_results, path):
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["center_name", "direction", "px_real", "px_imag", "delta_s_abs",
                    "relative_to_baseline", "n_unknowns"])
        for name, r in study_results.items():
            base = abs(r.long_range_baseline_px)
            for direction, c in r.long_range_cases.items():
                rel = c.delta_s_abs / base if base > 0 else float("inf")
                w.writerow([name, direction, repr(c.px.real), repr(c.px.imag),
                            repr(c.delta_s_abs), repr(rel), c.n_unknowns])
    return path


def write_neighborhood_decision_json(study_results, path):
    doc = {}
    for name, r in study_results.items():
        doc[name] = {
            "center_lw_mm": list(r.center_lw_mm),
            "final_decision": r.decision.final_decision,
            "rationale": r.decision.rationale,
            "per_contrast_classification": r.decision.per_contrast_classification,
            "per_contrast_relative_prediction_error": {
                k: v.relative_to_s0 for k, v in r.decision.per_contrast.items()
            },
            "long_range_relative": r.decision.long_range_relative,
            "long_range_significant": r.decision.long_range_significant,
            "discretization_floor_flags": r.discretization_floor_flags,
        }
    return _write_json(doc, path)


# --------------------------------------------------------------------------- environment_dataset/

def _sample_row(sample, split):
    row = [sample.sample_id, split, sample.source, sample.category,
           repr(sample.center_lw_mm[0]), repr(sample.center_lw_mm[1]), repr(sample.frequency_hz)]
    for pos in ALL_POSITIONS:
        l, w = sample.neighbors_lw_mm[pos]
        row.extend([repr(l), repr(w)])
    row.extend([repr(sample.s_complex.real), repr(sample.s_complex.imag)])
    return row


_DATASET_HEADER = (
    ["sample_id", "split", "source", "category", "center_L_mm", "center_W_mm", "frequency_hz"]
    + [f"{pos}_{axis}_mm" for pos in ALL_POSITIONS for axis in ("L", "W")]
    + ["S_real", "S_imag"]
)


def write_dataset_split_csv(samples, splits, split_name, path):
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(_DATASET_HEADER)
        for s in samples:
            if splits[s.sample_id] == split_name:
                w.writerow(_sample_row(s, split_name))
    return path


def write_dataset_metadata_json(samples, splits, path):
    from collections import Counter
    counts = Counter(splits.values())
    doc = {
        "n_total": len(samples),
        "n_train": counts.get("train", 0), "n_validation": counts.get("validation", 0),
        "n_blind_test": counts.get("blind_test", 0),
        "source_counts": dict(Counter(s.source for s in samples)),
        "category_counts": dict(Counter(s.category for s in samples)),
        "center_counts": dict(Counter(s.center_name for s in samples)),
        "geometry_based_holdout": "center C5 (4.2, 4.2) mm reserved entirely for blind_test (Sec. 21)",
        "symmetry_leakage_policy": (
            "samples related by a y-reflection, x-reflection or 180-degree rotation of the 3x3 "
            "environment are grouped and assigned to the same split (Rule 4); see "
            "dataset.symmetry_canonical_key"
        ),
        "random_seed": RANDOM_SEED,
    }
    return _write_json(doc, path)


# --------------------------------------------------------------------------------------- surrogate/

def write_model_metadata_json(surrogate, best_lambda, val_error, path):
    doc = surrogate.to_dict()
    doc["selected_ridge_lambda"] = best_lambda
    doc["validation_mean_complex_error"] = val_error
    return _write_json(doc, path)


def write_predictions_csv(samples, model_predictions, path):
    """``model_predictions``: {"baseline1": [...], "baseline2": [...], "surrogate": [...]}
    (each aligned index-for-index with ``samples``)."""
    path = Path(path)
    model_names = list(model_predictions.keys())
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        header = ["sample_id", "S_ref_real", "S_ref_imag"]
        for m in model_names:
            header.extend([f"{m}_real", f"{m}_imag"])
        w.writerow(header)
        for i, s in enumerate(samples):
            row = [s.sample_id, repr(s.s_complex.real), repr(s.s_complex.imag)]
            for m in model_names:
                p = model_predictions[m][i]
                if isinstance(p, complex):
                    row.extend([repr(p.real), repr(p.imag)])
                else:
                    row.extend([str(p), str(p)])
            w.writerow(row)
    return path


def write_error_metrics_csv(reports_by_model, path):
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["model", "n", "mean_complex_error", "median_complex_error", "rms_complex_error",
                    "p95_complex_error", "max_complex_error", "mean_relative_error",
                    "p95_relative_error", "max_relative_error", "mean_magnitude_error",
                    "p95_magnitude_error", "max_magnitude_error", "mean_phase_error_deg",
                    "p95_phase_error_deg", "max_phase_error_deg"])
        for model, r in reports_by_model.items():
            w.writerow([model, r.n, repr(r.mean_complex_error), repr(r.median_complex_error),
                        repr(r.rms_complex_error), repr(r.p95_complex_error), repr(r.max_complex_error),
                        repr(r.mean_relative_error), repr(r.p95_relative_error), repr(r.max_relative_error),
                        repr(r.mean_magnitude_error), repr(r.p95_magnitude_error), repr(r.max_magnitude_error),
                        repr(r.mean_phase_error_deg), repr(r.p95_phase_error_deg), repr(r.max_phase_error_deg)])
    return path


def write_generalization_csv(cases, label, path):
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["scenario", "i", "j", "center_L_mm", "center_W_mm", "px_mom_real", "px_mom_imag",
                    "px_surrogate_real", "px_surrogate_imag", "complex_error"])
        for c in cases:
            if isinstance(c.px_surrogate, complex):
                sr, si = repr(c.px_surrogate.real), repr(c.px_surrogate.imag)
                err = repr(abs(c.px_mom - c.px_surrogate))
            else:
                sr = si = err = str(c.px_surrogate)
            w.writerow([label, c.i, c.j, repr(c.center_lw_mm[0]), repr(c.center_lw_mm[1]),
                        repr(c.px_mom.real), repr(c.px_mom.imag), sr, si, err])
    return path


# ---------------------------------------------------------------------------- inverse_validation/

def write_target_responses_csv(results, path):
    """Sec. 45: ``target_responses.csv`` -- the target complex response and fixed environment
    each inversion was asked to reproduce."""
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["case_index", "target_real", "target_imag"]
                         + [f"{pos}_{axis}_mm" for pos in ALL_POSITIONS for axis in ("L", "W")])
        for idx, r in enumerate(results):
            row = [idx, repr(r.target.real), repr(r.target.imag)]
            for pos in ALL_POSITIONS:
                l_mm, w_mm = r.neighbors_lw_mm[pos]
                row.extend([repr(l_mm), repr(w_mm)])
            writer.writerow(row)
    return path


def write_inverse_results_csv(results, path):
    """Sec. 45: ``inverse_results.csv`` -- what the surrogate-driven search proposed."""
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["case_index", "best_candidate_L_mm", "best_candidate_W_mm",
                    "surrogate_predicted_real", "surrogate_predicted_imag"])
        for idx, r in enumerate(results):
            cand = r.best_candidate_lw_mm or (None, None)
            pred = r.surrogate_predicted_at_candidate
            pr, pi = (pred.real, pred.imag) if isinstance(pred, complex) else (None, None)
            w.writerow([idx, cand[0], cand[1], pr, pi])
    return path


def write_mom_verification_csv(results, path):
    """Sec. 45, 39: ``mom_verification.csv`` -- the SAME candidate re-checked against the real
    heterogeneous MoM solver, never accepted on the surrogate's own say-so."""
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["case_index", "mom_achieved_real", "mom_achieved_imag", "mom_relative_error", "status"])
        for idx, r in enumerate(results):
            achieved = r.mom_achieved_at_candidate
            ar, ai = (achieved.real, achieved.imag) if isinstance(achieved, complex) else (None, None)
            w.writerow([idx, ar, ai,
                        repr(r.mom_relative_error) if r.mom_relative_error != float("inf") else "inf",
                        r.status])
    return path


# ------------------------------------------------------------------------------- beam_validation/

def write_beam_validation_json(beam_result, path):
    m_mom, m_sur, comp = beam_result["metrics_mom"], beam_result["metrics_surrogate"], beam_result["comparison"]
    doc = {
        "n_surrogate_excluded": beam_result["n_surrogate_excluded"],
        "mom_pointing_error_deg": list(m_mom.pointing_error_deg),
        "mom_matched_peak_power_db": list(m_mom.matched_peak_power_db),
        "surrogate_pointing_error_deg": list(m_sur.pointing_error_deg) if m_sur else None,
        "surrogate_matched_peak_power_db": list(m_sur.matched_peak_power_db) if m_sur else None,
        "far_field_comparison": (
            {"complex_error_mean": comp.complex_error_mean, "complex_error_max": comp.complex_error_max,
             "amplitude_error_mean": comp.amplitude_error_mean, "phase_error_deg_mean": comp.phase_error_deg_mean}
            if comp else None
        ),
    }
    return _write_json(doc, path)


# --------------------------------------------------------------------------------- top-level

def write_phase10_summary_json(decision, path):
    doc = {
        "phase": "10",
        "final_status": decision.final_status,
        "gates": {
            "10A_control_variable": decision.gate_10a.status,
            "10B_environment_model": decision.gate_10b.status,
            "10C_surrogate_forward": decision.gate_10c.status,
            "10D_inverse_validation": decision.gate_10d.status,
            "10E_beam_level": decision.gate_10e.status,
        },
        "surrogate_beats_baseline1": decision.surrogate_beats_baseline1,
        "inverse_mostly_reachable": decision.inverse_mostly_reachable,
        "beam_level_acceptable": decision.beam_level_acceptable,
        "geometry_inversion": {"status": "UNLOCKED" if decision.geometry_inversion_unlocked else "LOCKED"},
        "step_generation": {"status": "UNLOCKED" if decision.geometry_inversion_unlocked else "LOCKED"},
    }
    return _write_json(doc, path)


def write_phase10_gate_status_json(decision, path):
    doc = {
        "phase": "10",
        "final_status": decision.final_status,
        "geometry_inversion": {"status": "LOCKED" if not decision.geometry_inversion_unlocked else "UNLOCKED"},
        "step_generation": {"status": "LOCKED" if not decision.geometry_inversion_unlocked else "UNLOCKED"},
        "gate_10a_checklist": decision.gate_10a.checklist,
        "gate_10b_checklist": decision.gate_10b.checklist,
        "gate_10c_checklist": decision.gate_10c.checklist,
        "gate_10d_checklist": decision.gate_10d.checklist,
        "gate_10e_checklist": decision.gate_10e.checklist,
    }
    return _write_json(doc, path)
