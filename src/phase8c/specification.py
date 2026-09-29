"""Phase 8C Sec. 47-48: CSV outputs and the machine-readable `control_quantity_validation.json`.
Every field is populated from an already-computed object passed in; nothing is invented here."""

import csv
import json
from pathlib import Path


def write_scattering_csv(s11_rows, s21_rows, path):
    """unit_cell_scattering_12p5GHz.csv (Sec. 47)."""
    by_lw_s21 = {(r.l_mm, r.w_mm): r for r in s21_rows}
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(
            [
                "L_mm",
                "W_mm",
                "S11xx_real",
                "S11xx_imag",
                "S11xx_mag",
                "S11xx_phase_deg",
                "S21xx_real",
                "S21xx_imag",
                "S21xx_mag",
                "S21xx_phase_deg",
            ]
        )
        for r11 in s11_rows:
            r21 = by_lw_s21.get((r11.l_mm, r11.w_mm))
            writer.writerow(
                [
                    repr(r11.l_mm),
                    repr(r11.w_mm),
                    repr(r11.real),
                    repr(r11.imag),
                    repr(r11.mag),
                    repr(r11.phase_deg),
                    repr(r21.real) if r21 else "",
                    repr(r21.imag) if r21 else "",
                    repr(r21.mag) if r21 else "",
                    repr(r21.phase_deg) if r21 else "",
                ]
            )
    return path


def write_control_quantity_candidates_csv(coverages, path):
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(
            [
                "component",
                "n_points",
                "coverage_deg",
                "largest_gap_deg",
                "magnitude_min",
                "magnitude_max",
                "magnitude_median",
                "magnitude_std",
                "phase_magnitude_correlation",
            ]
        )
        for c in coverages:
            writer.writerow(
                [
                    c.component,
                    c.n_points,
                    repr(c.coverage_deg),
                    repr(c.largest_gap_deg),
                    repr(c.magnitude_min),
                    repr(c.magnitude_max),
                    repr(c.magnitude_median),
                    repr(c.magnitude_std),
                    repr(c.phase_magnitude_correlation),
                ]
            )
    return path


def write_phase_state_analysis_csv(analyses, path):
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["component", "n_geometries", "n_distinct_complex_states", "n_distinct_phase_states"])
        for component, a in analyses:
            writer.writerow([component, a.n_geometries, a.n_distinct_complex_states, a.n_distinct_phase_states])
    return path


def build_control_quantity_validation(status, audit, s11_coverage, s21_coverage, s11_investigation, s21_investigation):
    """Sec. 48's schema, populated entirely from already-computed objects."""
    return {
        "status": status,
        "architecture": {
            "classification": audit.architecture_classification,
            "ground_plane": audit.ground_plane,
            "transmission_path": "present (ungrounded, two-port)",
            "reflection_path": "present (ungrounded, two-port)",
        },
        "candidates": {
            "S11xx": {
                "physical_interpretation": s11_investigation.physical_definition,
                "holographic_compatibility": "PARTIAL (Formulation B candidate; unresolved which architecture)",
                "phase_coverage_deg": s11_coverage.coverage_deg,
                "largest_circular_gap_deg": s11_coverage.largest_gap_deg,
                "magnitude_range": [s11_coverage.magnitude_min, s11_coverage.magnitude_max],
                "status": s11_investigation.validated,
            },
            "S21xx": {
                "physical_interpretation": s21_investigation.physical_definition,
                "holographic_compatibility": "PARTIAL (Formulation B candidate; unresolved which architecture)",
                "phase_coverage_deg": s21_coverage.coverage_deg,
                "largest_circular_gap_deg": s21_coverage.largest_gap_deg,
                "magnitude_range": [s21_coverage.magnitude_min, s21_coverage.magnitude_max],
                "status": s21_investigation.validated,
            },
        },
        "hfss_validation": {"performed": False},
    }


def write_json(doc, path):
    path = Path(path)
    path.write_text(json.dumps(doc, indent=1) + "\n")
    return path
