"""Sec. 27-28: required CSV/JSON output writers."""

import csv
import json
from pathlib import Path


def _write_json(doc, path):
    Path(path).write_text(json.dumps(doc, indent=1) + "\n")
    return path


# ------------------------------------------------------------------------------------- Phase 9A

def write_hfss_reference_cells_csv(solver_rows, hfss_records, path):
    hfss_by_key = {(r.cell_id, round(r.frequency_ghz, 6)): r for r in hfss_records} if hfss_records else {}
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["cell_id", "L_mm", "W_mm", "frequency_GHz",
                          "solver_S11_real", "solver_S11_imag", "solver_S21_real", "solver_S21_imag",
                          "hfss_S11_real", "hfss_S11_imag", "hfss_S21_real", "hfss_S21_imag"])
        for row in solver_rows:
            key = (row.cell_id, round(row.frequency_ghz, 6))
            hfss = hfss_by_key.get(key)
            hfss_vals = ([repr(hfss.s11.real), repr(hfss.s11.imag), repr(hfss.s21.real), repr(hfss.s21.imag)]
                         if hfss is not None else ["NOT_AVAILABLE"] * 4)
            writer.writerow([
                row.cell_id, repr(row.l_mm), repr(row.w_mm), repr(row.frequency_ghz),
                repr(row.s11.real), repr(row.s11.imag), repr(row.s21.real), repr(row.s21.imag),
                *hfss_vals,
            ])
    return path


def write_s_parameter_comparison_csv(solver_rows, hfss_records, component, path):
    """Sec. 9-10: ``s11_comparison.csv`` / ``s21_comparison.csv``. Every row is written even when
    no HFSS data exists -- the comparison columns are `NOT_AVAILABLE`, never omitted or fabricated."""
    from .error_metrics import complex_error as complex_error_fn
    from .error_metrics import magnitude_error as magnitude_error_fn
    from .error_metrics import phase_error_deg as phase_error_fn

    hfss_by_key = {(r.cell_id, round(r.frequency_ghz, 6)): r for r in hfss_records} if hfss_records else {}
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["cell_id", "L_mm", "W_mm", "frequency_GHz",
                          f"solver_{component}_real", f"solver_{component}_imag",
                          f"hfss_{component}_real", f"hfss_{component}_imag",
                          "complex_error", "magnitude_error", "phase_error_deg"])
        for row in solver_rows:
            key = (row.cell_id, round(row.frequency_ghz, 6))
            hfss = hfss_by_key.get(key)
            solver_val = getattr(row, component)
            if hfss is None:
                writer.writerow([row.cell_id, repr(row.l_mm), repr(row.w_mm), repr(row.frequency_ghz),
                                  repr(solver_val.real), repr(solver_val.imag),
                                  "NOT_AVAILABLE", "NOT_AVAILABLE", "NOT_AVAILABLE", "NOT_AVAILABLE", "NOT_AVAILABLE"])
                continue
            hfss_val = getattr(hfss, component)
            writer.writerow([
                row.cell_id, repr(row.l_mm), repr(row.w_mm), repr(row.frequency_ghz),
                repr(solver_val.real), repr(solver_val.imag), repr(hfss_val.real), repr(hfss_val.imag),
                repr(float(complex_error_fn(solver_val, hfss_val))),
                repr(float(magnitude_error_fn(solver_val, hfss_val))),
                repr(float(phase_error_fn(solver_val, hfss_val))),
            ])
    return path


def write_architecture_audit_json(audit_dict, path):
    return _write_json(audit_dict, path)


# ------------------------------------------------------------------------------------- Phase 9B

def write_supercell_definitions_csv(cases, path):
    """``cases``: iterable of (name, n, grid_m) with ``grid_m[i][j] = (L_m, W_m)``."""
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["case_name", "n", "row_i", "col_j", "L_mm", "W_mm"])
        for name, n, grid in cases:
            for i in range(n):
                for j in range(n):
                    l_m, w_m = grid[i][j]
                    writer.writerow([name, n, i, j, repr(l_m * 1e3), repr(w_m * 1e3)])
    return path


def write_supercell_results_csv(results, path):
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["name", "n", "isolated_px_real", "isolated_px_imag", "embedded_px_real",
                          "embedded_px_imag", "complex_error", "magnitude_error", "phase_error_deg",
                          "relative_complex_error", "n_unknowns_isolated", "n_unknowns_embedded"])
        for r in results:
            writer.writerow([
                r.name, r.n, repr(r.isolated_px.real), repr(r.isolated_px.imag),
                repr(r.embedded_px.real), repr(r.embedded_px.imag), repr(r.complex_error),
                repr(r.magnitude_error), repr(r.phase_error_deg), repr(r.relative_complex_error),
                r.n_unknowns_isolated, r.n_unknowns_embedded,
            ])
    return path


def write_lpa_error_metrics_csv(results, path):
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["name", "n", "complex_error", "magnitude_error", "phase_error_deg", "relative_complex_error"])
        for r in results:
            writer.writerow([r.name, r.n, repr(r.complex_error), repr(r.magnitude_error),
                              repr(r.phase_error_deg), repr(r.relative_complex_error)])
    return path


def write_locality_analysis_json(results, path):
    by_name = {r.name: r for r in results}
    pairs = [("moderate", "3x3_moderate", "5x5_moderate"), ("strong", "3x3_strong", "5x5_strong")]
    doc = {}
    for label, name_3x3, name_5x5 in pairs:
        if name_3x3 not in by_name or name_5x5 not in by_name:
            continue
        r3, r5 = by_name[name_3x3], by_name[name_5x5]
        doc[label] = {
            "relative_complex_error_3x3": r3.relative_complex_error,
            "relative_complex_error_5x5": r5.relative_complex_error,
            "change_3x3_to_5x5": r5.relative_complex_error - r3.relative_complex_error,
            "note": (
                "influence appears localized (5x5 close to 3x3)" if abs(r5.relative_complex_error - r3.relative_complex_error) < 0.5 * r3.relative_complex_error
                else "influence extends beyond the immediate 3x3 neighborhood (5x5 differs substantially from 3x3)"
            ),
        }
    return _write_json(doc, path)


# --------------------------------------------------------------------------------- top-level

def write_phase9_summary_json(decision, path):
    doc = {
        "phase": "9",
        "final_status": decision.final_status,
        "gate_a_status": decision.gate_a.status,
        "gate_b_status": decision.gate_b.status,
        "control_variable": {
            "s11": decision.s11_status, "s21": decision.s21_status,
            "selected_control": decision.selected_control,
        },
        "local_periodicity": {"status": decision.local_periodicity_status},
        "geometry_inversion": {"status": "UNLOCKED" if decision.geometry_inversion_unlocked else "LOCKED"},
        "step_generation": {"status": "UNLOCKED" if decision.geometry_inversion_unlocked else "LOCKED"},
    }
    return _write_json(doc, path)


def write_phase9_gate_status_json(decision, path):
    doc = {
        "phase": "9",
        "status": "PASS" if decision.geometry_inversion_unlocked else (
            "PARTIAL" if decision.final_status == "PARTIALLY_VALIDATED" else "BLOCKED"
        ),
        "control_variable": {
            "s11": decision.s11_status, "s21": decision.s21_status,
            "selected_control": decision.selected_control,
        },
        "local_periodicity": {"status": decision.local_periodicity_status},
        "geometry_inversion": {"status": "LOCKED" if not decision.geometry_inversion_unlocked else "UNLOCKED"},
        "step_generation": {"status": "LOCKED" if not decision.geometry_inversion_unlocked else "UNLOCKED"},
        "gate_a_checklist": decision.gate_a.checklist,
        "gate_b_checklist": decision.gate_b.checklist,
    }
    return _write_json(doc, path)
