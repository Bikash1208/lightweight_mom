"""Phase 12C orchestration: grounded supercells (Task A), local-periodicity test (Task B),
coupling extent (Task C), fine-aligned-mesh spot-check (Task D, read-only for earlier results),
HFSS-ready exports (Task E). Writes every required output and prints a terminal summary.
"""

from dataclasses import dataclass
from pathlib import Path

from . import coupling_study, exports, gate, locality_study, plots, specification, spot_check
from . import design_point as dp
from .grounded_supercell import homogeneous_correctness_check, supercell_mesh_convergence_check
from .report import render_phase12c_report_md

CORRECTNESS_L_MM = 4.0  # smallest-aligned-mesh representative point (module docstring, locality_study.py)
CONVERGENCE_MESH_SEQUENCE = (4, 8, 16, 32)
FEASIBILITY_SIZES = (1, 3, 5)


@dataclass(frozen=True)
class Phase12cResult:
    output_dir: Path
    task_a_status: str
    correctness: object
    convergence_series: tuple
    feasibility_points: tuple
    curve: tuple
    centers: object
    gradient_results: tuple
    strong_contrast: object
    periodic_results: tuple
    periodic_grid_mm: tuple
    periodic_span_deg: float
    periodic_is_full_2pi: bool
    task_b_gate: object
    task_c_status: str
    coupling_reports: dict
    task_d_gate: object
    single_cell_exports: tuple
    supercell_exports: tuple
    hfss_export_status: str
    recommendation: str


RECOMMENDATION = (
    "Local periodicity, at this grounded 8 mm design point, holds within the stated tolerance only "
    "for small neighbour phase steps and degrades sharply approaching resonance (Sec. 3); a "
    "reflectarray-style lookup-table forward model is therefore defensible ONLY if the aperture's "
    "own phase gradient per cell is kept at or below the measured max-step-within-tolerance "
    "everywhere (a real, checkable design constraint, not an assumption), and the periodic "
    "phase-gradient supercell (Sec. 3) shows this constraint is hard to satisfy near resonance "
    "specifically -- an aperture whose gradient must pass through that region cannot rely on a pure "
    "lookup table there. A full-aperture integral-equation solver (a la Budhu & Grbic 2021) remains "
    "the more broadly defensible choice unless the design is constrained to stay in the "
    "small-step-safe region measured here; this is offered as a recommendation, not a decision."
)


def run_phase12c(
    output_dir,
    correctness_mesh_sequence=CONVERGENCE_MESH_SEQUENCE,
    feasibility_sizes=FEASIBILITY_SIZES,
    feasibility_mesh_cells_per_subcell=20,
    periodic_n_cells=6,
    coupling_sizes=None,
    spot_check_mesh_cells_per_subcell=None,
    export_mesh_cells_per_side=None,
):
    """``correctness_mesh_sequence``/``feasibility_sizes``/``periodic_n_cells``/``coupling_sizes``/
    ``spot_check_mesh_cells_per_subcell``/``export_mesh_cells_per_side``: reduced-scale overrides
    for tests; the production run (``examples/run_phase12c.py``) uses every default unchanged."""
    output_dir = Path(output_dir)
    plots_dir = output_dir / "plots"
    export_dir = Path(str(output_dir).replace("phase12c_grounded_locality", "phase12c_hfss_export"))
    output_dir.mkdir(parents=True, exist_ok=True)
    coupling_sizes = coupling_sizes if coupling_sizes is not None else FEASIBILITY_SIZES

    # ---------------------------------------------------------------------------------- Task A
    correctness = homogeneous_correctness_check(
        CORRECTNESS_L_MM * 1e-3,
        3,
        correctness_mesh_sequence[0],
        dp.PERIOD_M,
        dp.SUBSTRATE,
        dp.METAL,
        dp.DESIGN_FREQUENCY_HZ,
    )
    convergence_series = []
    for base, finer in zip(correctness_mesh_sequence, correctness_mesh_sequence[1:]):
        convergence_series.append(
            supercell_mesh_convergence_check(
                CORRECTNESS_L_MM * 1e-3,
                3,
                base,
                finer,
                dp.PERIOD_M,
                dp.SUBSTRATE,
                dp.METAL,
                dp.DESIGN_FREQUENCY_HZ,
            )
        )
    feasibility_points = [
        homogeneous_correctness_check(
            CORRECTNESS_L_MM * 1e-3,
            n,
            feasibility_mesh_cells_per_subcell,
            dp.PERIOD_M,
            dp.SUBSTRATE,
            dp.METAL,
            dp.DESIGN_FREQUENCY_HZ,
        )
        for n in feasibility_sizes
    ]
    task_a_status = (
        "GROUNDED_SUPERCELLS_VALIDATED" if correctness.relative_error < 1e-9 else "GROUNDED_SUPERCELLS_UNVALIDATED"
    )

    # ---------------------------------------------------------------------------------- Task B
    curve = locality_study.build_reference_curve()
    centers = locality_study.classify_centers(curve)
    gradient_results = []
    for label, center in (
        ("off_resonance", centers.off_resonance),
        ("near_resonance", centers.near_resonance),
        ("at_resonance", centers.at_resonance),
    ):
        steps = locality_study.immediate_neighbor_steps(curve, center)
        west = next((s for s in steps if s.direction == "lower_l"), None)
        east = next((s for s in steps if s.direction == "higher_l"), None)
        w_lw = (west.neighbor.l_mm, west.neighbor.l_mm) if west else (center.l_mm, center.l_mm)
        e_lw = (east.neighbor.l_mm, east.neighbor.l_mm) if east else (center.l_mm, center.l_mm)
        grid_mm = locality_study.build_1d_row_grid((center.l_mm, center.l_mm), w_lw, e_lw)
        step_mags = [abs(s.phase_step_deg) for s in steps]
        achieved = sum(step_mags) / len(step_mags)
        gradient_results.append(
            locality_study.run_gradient_case(f"{label}_1d_row", (center.l_mm, center.l_mm), grid_mm, achieved)
        )
        if label == "near_resonance" and east is not None:
            grid2d, clamped = locality_study.build_2d_gradient_grid(
                (center.l_mm, center.l_mm), curve, center, east.phase_step_deg
            )
            gradient_results.append(
                locality_study.run_gradient_case(
                    f"{label}_2d_gradient", (center.l_mm, center.l_mm), grid2d, abs(east.phase_step_deg)
                )
            )

    lo, hi = locality_study.L_FAMILY_MM[0], locality_study.L_FAMILY_MM[-1]
    center_mid = centers.near_resonance
    # The opposite extreme is whichever domain endpoint (lo/hi) is FARTHER in raw geometry from the
    # centre (Phase 9/11's own "strong_isolated_defect" convention, reused in spirit), not always lo.
    strong_opposite_mm = lo if abs(center_mid.l_mm - hi) < abs(center_mid.l_mm - lo) else hi
    strong_grid_mm = coupling_study.build_checkerboard_pattern_mm(center_mid.l_mm, strong_opposite_mm, 3)
    strong_opposite_point = next(p for p in curve if p.l_mm == strong_opposite_mm)
    strong_step = abs(
        locality_study.wrapped_difference_deg(strong_opposite_point.s11_phase_deg, center_mid.s11_phase_deg)
    )
    strong_contrast = locality_study.run_gradient_case(
        "strong_contrast", (center_mid.l_mm, center_mid.l_mm), strong_grid_mm, strong_step
    )

    periodic_results, periodic_grid_mm, periodic_span_deg, periodic_is_full_2pi = locality_study.run_periodic_supercell(
        periodic_n_cells,
        curve=curve,
    )

    task_b_gate = gate.decide_local_periodicity(gradient_results)

    # ---------------------------------------------------------------------------------- Task C
    off_l, near_l = centers.off_resonance.l_mm, centers.near_resonance.l_mm
    coupling_reports = {}
    for label, center_l in (("off_resonance", off_l), ("near_resonance", near_l)):
        opposite = hi if center_l == lo else lo
        gradient_pattern = coupling_study.build_gradient_pattern_mm(center_l, lo, hi, max(coupling_sizes))
        points = coupling_study.run_coupling_radius_series(gradient_pattern, sizes=coupling_sizes)
        coupling_reports[f"{label}_gradient"] = coupling_study.analyze_coupling_series(f"{label}_gradient", points)

        checker_pattern = coupling_study.build_checkerboard_pattern_mm(center_l, opposite, max(coupling_sizes))
        points_ck = coupling_study.run_coupling_radius_series(checker_pattern, sizes=coupling_sizes)
        coupling_reports[f"{label}_checkerboard"] = coupling_study.analyze_coupling_series(
            f"{label}_checkerboard", points_ck
        )

    any_converged = any(r.convergence.status == "CONVERGED" for r in coupling_reports.values())
    task_c_status = (
        "GROUNDED_COUPLING_LOCAL_FOR_TESTED_SIZES" if any_converged else "GROUNDED_COUPLING_NOT_LOCAL_FOR_TESTED_SIZES"
    )

    # ---------------------------------------------------------------------------------- Task D
    c2_c3 = spot_check.recheck_phase10_c2_c3()
    phase9_strong = spot_check.recheck_phase9_strong_contrast(mesh_cells_per_subcell=spot_check_mesh_cells_per_subcell)
    case_d = spot_check.recheck_phase11_case_d(mesh_cells_per_subcell=spot_check_mesh_cells_per_subcell)
    coupling_mesh, coupling_points = spot_check.recheck_phase11_coupling_radius(
        sizes=coupling_sizes[:2] if len(coupling_sizes) >= 2 else coupling_sizes,
        mesh_cells_per_subcell=spot_check_mesh_cells_per_subcell,
    )
    task_d_gate = gate.decide_spot_check(c2_c3, phase9_strong, case_d, coupling_points)

    # ---------------------------------------------------------------------------------- Task E
    export_kwargs = {"mesh_cells_per_side": export_mesh_cells_per_side} if export_mesh_cells_per_side else {}
    single_cell_exports = exports.export_single_cells(export_dir, **export_kwargs)
    export_sc_kwargs = {"mesh_cells_per_subcell": export_mesh_cells_per_side} if export_mesh_cells_per_side else {}
    supercell_exports = exports.export_supercells(export_dir, **export_sc_kwargs)
    synthetic_path, synthetic_records = exports.write_synthetic_hfss_reference_test(
        export_dir / "synthetic_hfss_reference_test.csv"
    )
    exports.write_export_manifest(single_cell_exports, supercell_exports, synthetic_path, export_dir / "manifest.json")
    hfss_export_status = gate.hfss_export_status(single_cell_exports, supercell_exports, synthetic_records)

    result = Phase12cResult(
        output_dir=output_dir,
        task_a_status=task_a_status,
        correctness=correctness,
        convergence_series=tuple(convergence_series),
        feasibility_points=tuple(feasibility_points),
        curve=curve,
        centers=centers,
        gradient_results=tuple(gradient_results),
        strong_contrast=strong_contrast,
        periodic_results=periodic_results,
        periodic_grid_mm=periodic_grid_mm,
        periodic_span_deg=periodic_span_deg,
        periodic_is_full_2pi=periodic_is_full_2pi,
        task_b_gate=task_b_gate,
        task_c_status=task_c_status,
        coupling_reports=coupling_reports,
        task_d_gate=task_d_gate,
        single_cell_exports=single_cell_exports,
        supercell_exports=supercell_exports,
        hfss_export_status=hfss_export_status,
        recommendation=RECOMMENDATION,
    )

    # ------------------------------------------------------------------------------- outputs
    specification.write_task_a_json(correctness, convergence_series, feasibility_points, output_dir / "task_a.json")
    specification.write_task_b_json(
        curve,
        centers,
        gradient_results,
        strong_contrast,
        periodic_results,
        periodic_grid_mm,
        periodic_span_deg,
        periodic_is_full_2pi,
        task_b_gate,
        output_dir / "task_b.json",
    )
    specification.write_task_c_json(coupling_reports, output_dir / "task_c.json")
    specification.write_task_d_json(
        c2_c3, phase9_strong, case_d, coupling_mesh, coupling_points, task_d_gate, output_dir / "task_d.json"
    )
    specification.write_task_d_spot_check_csv(c2_c3, phase9_strong, case_d, output_dir / "task_d_spot_check.csv")
    specification.write_gate_status_json(
        task_a_status,
        task_b_gate,
        task_c_status,
        task_d_gate,
        hfss_export_status,
        output_dir / "phase12c_gate_status.json",
    )

    try:
        plots.generate_all_plots(gradient_results, periodic_results, coupling_reports, {}, task_d_gate, plots_dir)
    except ImportError:
        pass

    (output_dir / "phase12c_report.md").write_text(render_phase12c_report_md(result))

    return result


def print_final_report(result):
    print("=" * 78)
    print("PHASE 12C -- DOES LOCAL PERIODICITY HOLD FOR THE GROUNDED 8 MM REFLECTARRAY CELL?")
    print("=" * 78)
    print(f"Task A (grounded supercells):  {result.task_a_status}")
    print(f"Task B (local periodicity):    {result.task_b_gate.status}")
    print(f"Task C (coupling extent):      {result.task_c_status}")
    print(f"Task D (fine-mesh spot-check): {result.task_d_gate.status}")
    print(f"Task E (HFSS-ready exports):   {result.hfss_export_status}")
    print("HFSS status:         HFSS_NOT_AVAILABLE")
    print("Geometry inversion:  GEOMETRY_INVERSION_LOCKED")
    print("STEP export:         STEP_EXPORT_LOCKED")
    print("=" * 78)
