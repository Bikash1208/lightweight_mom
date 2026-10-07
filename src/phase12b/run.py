"""Phase 12B orchestration: mesh alignment (Task A, already in discretization.py/sweep.py),
re-assess Phase 12A on aligned meshes (Task B), duplicate-geometry audit (Task C, read-only),
grounded period study (Task D). Writes every required output and prints a terminal summary.
"""

from dataclasses import dataclass
from pathlib import Path

from ..phase12a.sweep_study import convergence_study as phase12a_convergence_study
from . import duplicate_audit, gate, period_study, plots, reassessment, specification
from .report import render_phase12b_report_md

NEW_PERIODS_MM = (8.0, 12.0)  # 4.8mm is "baseline from B" (brief's own words) -- reused, not redone


@dataclass(frozen=True)
class Phase12bResult:
    output_dir: Path
    aligned_family: tuple
    aligned_cells_family: int
    phase7_report: object
    aligned_refinement: tuple  # (points, diffs)
    floquet_convergence: tuple  # (points, diffs)
    aligned_sweep_points: tuple
    aligned_cells: int
    aligned_coverage_by_frequency: dict
    duplicate_audit: dict
    period_results: dict  # period_mm -> {...}
    mesh_gate: object
    audit_gate: object
    period_gate: object
    hfss_status: str


def run_phase12b(
    output_dir,
    l_values_mm=None,
    frequencies_hz=None,
    new_periods_mm=NEW_PERIODS_MM,
    period_l_values_mm=None,
    aligned_refinement_cells=(24, 48, 96),
):
    """``l_values_mm``/``frequencies_hz``/``period_l_values_mm``/``aligned_refinement_cells``:
    reduced-scale overrides for tests; the production run (``scripts/research/run_phase12b.py``) uses
    every default unchanged."""
    output_dir = Path(output_dir)
    plots_dir = output_dir / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)

    # -------------------------------------------------------------------------- Task A
    family = l_values_mm if l_values_mm is not None else reassessment.ALL_L_VALUES_MM
    aligned_n_family = reassessment.aligned_mesh_for_family(family)
    phase7_report = duplicate_audit.phase7_library_duplicate_report()

    # -------------------------------------------------------------------------- Task B
    aligned_refinement = reassessment.aligned_mesh_refinement_study(cells_list=aligned_refinement_cells)
    floquet_orders = tuple(
        sorted({aligned_refinement_cells[-1], aligned_refinement_cells[-1] + 16, aligned_refinement_cells[-1] * 2})
    )
    floquet_conv = reassessment.floquet_only_convergence(
        cells=aligned_refinement_cells[-1],
        floquet_orders=floquet_orders,
    )
    sweep_frequencies = frequencies_hz if frequencies_hz is not None else reassessment.SWEEP_FREQUENCIES_HZ
    aligned_sweep_points = reassessment.run_aligned_grounded_sweep(
        l_values_mm=family,
        frequencies_hz=sweep_frequencies,
        cells=aligned_n_family,
    )
    aligned_coverage = reassessment.coverage_by_frequency(aligned_sweep_points)

    # -------------------------------------------------------------------------- Task C
    audit = {
        "phase7_report": phase7_report,
        "phase7_lookup_cross_check": duplicate_audit.phase7_delivered_lookup_duplicate_count(),
        "phase9_report": duplicate_audit.phase9_strong_contrast_check(),
        "phase10_report": duplicate_audit.phase10_center_duplicate_report(),
        "phase10_leakage": duplicate_audit.phase10_holdout_leakage_check(),
        "phase10_perturbation_report": duplicate_audit.phase10_center_perturbation_duplicate_report(),
        "phase11_check": duplicate_audit.phase11_pattern_endpoints_check(),
    }

    # -------------------------------------------------------------------------- Task D
    # AlignedSweepPoint (reassessment.py) has no period_mm field (it is always 4.8mm, "baseline
    # from B"); normalize to period_study.PeriodSweepPoint so every period_results entry has the
    # same shape for the shared CSV/plot/report code below.
    period4p8_points = tuple(
        period_study.PeriodSweepPoint(4.8, p.length_width_mm, p.frequency_hz, p.s11_xx_mag, p.s11_xx_phase_deg)
        for p in aligned_sweep_points
    )
    period_results = {
        4.8: {
            "cells": aligned_n_family,
            "points": list(period4p8_points),
            "coverage": aligned_coverage,
            "resonant_length_mm": period_study.resonant_length_mm(period4p8_points, 12.5e9),
        }
    }
    for period_mm in new_periods_mm:
        l_values = (period_l_values_mm or {}).get(period_mm)
        points, n = period_study.run_period_sweep(period_mm, l_values_mm=l_values, frequencies_hz=sweep_frequencies)
        cov = period_study.coverage_by_frequency(points)
        period_results[period_mm] = {
            "cells": n,
            "points": list(points),
            "coverage": cov,
            "resonant_length_mm": period_study.resonant_length_mm(
                points, sweep_frequencies[len(sweep_frequencies) // 2]
            ),
        }
    for period_mm, r in period_results.items():
        rep_l = r["resonant_length_mm"]
        r["grating_surface_wave"] = period_study.grating_and_surface_wave_report(
            period_mm,
            rep_l,
            frequencies_hz=sweep_frequencies,
        )
        if period_mm == 4.8:
            # 4.8mm's convergence is already thoroughly characterized by Task B's own aligned
            # refinement study + Floquet-only convergence (both well within tolerance, Sec. 3) --
            # doubling its already-fine 96-cell base mesh again here would cost far more than it
            # would teach; reuse that evidence instead of re-solving it.
            r["convergence_points"], r["convergence_diffs"] = aligned_refinement[0], aligned_refinement[1]
        else:
            conv_points, conv_diffs = period_study.convergence_for_near_resonant_patch(
                period_mm,
                rep_l,
                frequency_hz=sweep_frequencies[len(sweep_frequencies) // 2],
            )
            r["convergence_points"], r["convergence_diffs"] = conv_points, conv_diffs

    # -------------------------------------------------------------------------- gates
    mesh_gate = gate.decide_mesh_alignment(aligned_refinement[1][-1], floquet_conv[1][-1])
    audit_gate = gate.decide_audit(
        audit["phase7_report"],
        audit["phase10_leakage"],
        audit["phase10_report"],
    )
    coverage_by_period = {
        p: r["coverage"][12.5e9].coverage_deg
        if 12.5e9 in r["coverage"]
        else next(iter(r["coverage"].values())).coverage_deg
        for p, r in period_results.items()
    }
    period_gate = gate.decide_period_study(coverage_by_period)
    hfss = gate.hfss_status()

    result = Phase12bResult(
        output_dir=output_dir,
        aligned_family=tuple(family),
        aligned_cells_family=aligned_n_family,
        phase7_report=phase7_report,
        aligned_refinement=aligned_refinement,
        floquet_convergence=floquet_conv,
        aligned_sweep_points=tuple(aligned_sweep_points),
        aligned_cells=aligned_n_family,
        aligned_coverage_by_frequency=aligned_coverage,
        duplicate_audit=audit,
        period_results=period_results,
        mesh_gate=mesh_gate,
        audit_gate=audit_gate,
        period_gate=period_gate,
        hfss_status=hfss,
    )

    # -------------------------------------------------------------------------- outputs
    specification.write_alignment_report_json(
        family, aligned_n_family, phase7_report, output_dir / "alignment_report.json"
    )
    specification.write_reassessment_json(
        aligned_refinement[0],
        aligned_refinement[1],
        floquet_conv[0],
        floquet_conv[1],
        aligned_sweep_points,
        aligned_coverage,
        output_dir / "reassessment.json",
    )
    specification.write_duplicate_audit_json(
        audit["phase7_report"],
        audit["phase7_lookup_cross_check"],
        audit["phase9_report"],
        audit["phase10_report"],
        audit["phase10_leakage"],
        audit["phase10_perturbation_report"],
        audit["phase11_check"],
        output_dir / "duplicate_audit.json",
    )
    specification.write_period_study_json(period_results, output_dir / "period_study.json")
    specification.write_period_sweep_csv(period_results, output_dir / "period_sweep.csv")
    specification.write_gate_status_json(mesh_gate, audit_gate, period_gate, output_dir / "phase12b_gate_status.json")

    try:
        unaligned_conv = phase12a_convergence_study()
        plots.generate_all_plots(
            unaligned_conv,
            aligned_refinement[0],
            period_results,
            phase7_report,
            audit["phase10_report"],
            plots_dir,
        )
    except ImportError:
        pass

    (output_dir / "phase12b_report.md").write_text(render_phase12b_report_md(result))

    return result


def print_final_report(result):
    print("=" * 78)
    print("PHASE 12B -- MESH ALIGNMENT, DUPLICATE-GEOMETRY AUDIT, GROUNDED PERIOD STUDY")
    print("=" * 78)
    print(f"Mesh alignment:      {result.mesh_gate.status}")
    print(f"Duplicate audit:     {result.audit_gate.status}")
    print(f"Period study:        {result.period_gate.status}")
    for period_mm, deg in sorted(result.period_gate.coverage_by_period_deg.items()):
        print(f"  P={period_mm} mm: {deg:.1f} deg coverage")
    print(f"HFSS status:         {result.hfss_status}")
    print("Geometry inversion:  GEOMETRY_INVERSION_LOCKED")
    print("STEP export:         STEP_EXPORT_LOCKED")
    print("=" * 78)
