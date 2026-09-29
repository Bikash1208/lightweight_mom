"""Phase 12B required outputs: deterministic CSV/JSON, one writer per file (project convention)."""

import csv
import json
from dataclasses import asdict


def _write_json(doc, path):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(json.dumps(doc, indent=1, allow_nan=False) + "\n")
    return path


def write_alignment_report_json(family, aligned_n, duplicate_report, path):
    doc = {
        "l_values_mm": list(family),
        "aligned_cells": aligned_n,
        "phase7_duplicate_report": {
            "n_geometries": duplicate_report.n_geometries,
            "n_distinct": duplicate_report.n_distinct,
            "n_duplicate_groups": duplicate_report.n_duplicate_groups,
            "groups": [list(g) for g in duplicate_report.groups],
        },
    }
    return _write_json(doc, path)


def write_reassessment_json(
    aligned_points,
    aligned_diffs,
    floquet_points,
    floquet_diffs,
    sweep_points,
    coverage_by_frequency,
    path,
):
    doc = {
        "aligned_mesh_refinement_study": {
            "points": [asdict(p) for p in aligned_points],
            "diffs": list(aligned_diffs),
        },
        "floquet_only_convergence": {
            "points": [asdict(p) for p in floquet_points],
            "diffs": list(floquet_diffs),
        },
        "aligned_sweep_points": [asdict(p) for p in sweep_points],
        "coverage_by_frequency": {
            str(f): {"coverage_deg": c.coverage_deg, "largest_gap_deg": c.largest_gap_deg, "n_points": c.n_points}
            for f, c in coverage_by_frequency.items()
        },
    }
    return _write_json(doc, path)


def write_duplicate_audit_json(
    phase7_report,
    phase7_lookup_cross_check,
    phase9_report,
    phase10_report,
    phase10_leakage,
    phase10_perturbation_report,
    phase11_check,
    path,
):
    doc = {
        "phase7_library": {
            "n_geometries": phase7_report.n_geometries,
            "n_distinct": phase7_report.n_distinct,
            "n_duplicate_groups": phase7_report.n_duplicate_groups,
            "groups": [list(g) for g in phase7_report.groups],
        },
        "phase7_delivered_lookup_cross_check": phase7_lookup_cross_check,
        "phase9_strong_contrast_groups": [list(g) for g in phase9_report],
        "phase10_center_duplicate_groups": [list(g) for g in phase10_report],
        "phase10_holdout_leakage": {
            "holdout_center_names": list(phase10_leakage["holdout_center_names"]),
            "leaks": [list(g) for g in phase10_leakage["leaks"]],
        },
        "phase10_cross_center_perturbation_groups": [
            list(g) for g in phase10_perturbation_report["cross_center_groups"]
        ],
        "phase11_pattern_endpoints_identical": phase11_check,
    }
    return _write_json(doc, path)


def write_period_study_json(period_results, path):
    """``period_results``: {period_mm: {"cells": int, "points": [...], "coverage": {...},
    "resonant_length_mm": float, "grating_surface_wave": GratingSurfaceWaveReport}}."""
    doc = {}
    for period_mm, r in period_results.items():
        gsw = r["grating_surface_wave"]
        doc[str(period_mm)] = {
            "cells": r["cells"],
            "points": [asdict(p) for p in r["points"]],
            "coverage_by_frequency": {
                str(f): {"coverage_deg": c.coverage_deg, "largest_gap_deg": c.largest_gap_deg}
                for f, c in r["coverage"].items()
            },
            "resonant_length_mm": r["resonant_length_mm"],
            "grating_onset_normal_hz": gsw.grating_onset_normal_hz,
            "grating_onset_oblique_hz": {str(k): v for k, v in gsw.grating_onset_oblique_hz.items()},
            "specular_only_across_band": gsw.specular_only_across_band,
            "surface_wave_margins": {str(k): v for k, v in gsw.surface_wave_margins.items()},
            "convergence": {
                "points": [asdict(p) for p in r["convergence_points"]],
                "diffs": list(r["convergence_diffs"]),
            },
        }
    return _write_json(doc, path)


def write_gate_status_json(mesh_gate, audit_gate, period_gate, path):
    doc = {
        "mesh_alignment": asdict(mesh_gate),
        "audit": {
            "status": audit_gate.status,
            "reason": audit_gate.reason,
            "at_risk_findings": list(audit_gate.at_risk_findings),
        },
        "period_study": {
            "status": period_gate.status,
            "reason": period_gate.reason,
            "coverage_by_period_deg": period_gate.coverage_by_period_deg,
        },
    }
    return _write_json(doc, path)


def write_period_sweep_csv(period_results, path):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["period_mm", "length_width_mm", "frequency_hz", "s11_xx_mag", "s11_xx_phase_deg"])
        for period_mm, r in period_results.items():
            for p in sorted(r["points"], key=lambda p: (p.frequency_hz, p.length_width_mm)):
                writer.writerow(
                    [
                        repr(p.period_mm),
                        repr(p.length_width_mm),
                        repr(p.frequency_hz),
                        repr(p.s11_xx_mag),
                        repr(p.s11_xx_phase_deg),
                    ]
                )
    return path
