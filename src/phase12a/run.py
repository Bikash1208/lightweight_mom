"""Phase 12A orchestration: run the physics-sanity sweep, the analytic/surface-wave/convergence
validation, decide the gate status, write every required output, and print a terminal summary.
"""

from dataclasses import dataclass
from pathlib import Path

from ..discretization import refinement_study
from . import analytic_checks, gate, plots, specification
from .report import render_phase12a_report_md
from .sweep_study import (
    DESIGN_FREQUENCY_HZ,
    SWEEP_FREQUENCIES_HZ,
    _mesh,
    _planes,
    convergence_study,
    coverage_by_frequency,
    grounded_library_cell,
    run_grounded_patch_sweep,
    ungrounded_phase7_coverage,
    ungrounded_phase7_diagonal_table,
)

CONVERGENCE_CENTER_LW_M = 2.8e-3  # C3, this project's own already-used center geometry


@dataclass(frozen=True)
class Phase12aResult:
    output_dir: Path
    sweep_points: tuple
    coverage_by_frequency: dict
    design_frequency_hz: float
    ungrounded_comparison: dict
    convergence_points: tuple
    floquet_error: float
    mesh_error: float
    mesh_floquet_tolerance: float
    slab_check_points: tuple
    max_analytic_error: float
    thin_limit: complex
    surface_wave_points: tuple
    decision: object


def run_phase12a(
    output_dir,
    l_values_mm=None,
    frequencies_hz=SWEEP_FREQUENCIES_HZ,
    convergence_sizes=(8, 16, 32, 64),
    skip_ungrounded_comparison=False,
):
    """``convergence_sizes``/``l_values_mm``/``frequencies_hz``/``skip_ungrounded_comparison``:
    reduced-scale overrides for tests (module docstring pattern of every later phase's own
    ``run_phaseN``); the production run (``examples/run_phase12a_ground_plane.py``) uses every
    default unchanged."""
    output_dir = Path(output_dir)
    plots_dir = output_dir / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)

    kwargs = {"frequencies_hz": frequencies_hz}
    if l_values_mm is not None:
        kwargs["l_values_mm"] = l_values_mm
    sweep_points = run_grounded_patch_sweep(**kwargs)
    coverage = coverage_by_frequency(sweep_points)

    if skip_ungrounded_comparison:
        ungrounded_comparison = {"status": "NOT_AVAILABLE_FROM_CURRENT_DATA"}
        ungrounded_diagonal = []
    else:
        ungrounded_comparison = ungrounded_phase7_coverage(frequency_hz=DESIGN_FREQUENCY_HZ)
        ungrounded_diagonal = ungrounded_phase7_diagonal_table(frequency_hz=DESIGN_FREQUENCY_HZ)

    convergence_points = convergence_study(CONVERGENCE_CENTER_LW_M, sizes=convergence_sizes)

    cell = grounded_library_cell(CONVERGENCE_CENTER_LW_M)
    planes = _planes()
    mesh = _mesh()
    from ..sparameters import compute_s_parameters

    base = compute_s_parameters(cell, DESIGN_FREQUENCY_HZ, mesh, planes, warn_grating=False)
    floquet_error, mesh_error = refinement_study(cell, planes, mesh, DESIGN_FREQUENCY_HZ, base.s11, 2)

    slab_points = analytic_checks.bare_grounded_slab_checks()
    max_analytic_error = max(p.max_error for p in slab_points)
    thin_limit = analytic_checks.vanishing_thickness_limit()

    surface_wave_points = analytic_checks.surface_wave_margin_table()
    any_sw_in_band = any(p.any_in_band for p in surface_wave_points)

    decision = gate.decide(
        max_analytic_error,
        floquet_error,
        mesh_error,
        any_sw_in_band,
        coverage[DESIGN_FREQUENCY_HZ].coverage_deg,
    )

    result = Phase12aResult(
        output_dir=output_dir,
        sweep_points=sweep_points,
        coverage_by_frequency=coverage,
        design_frequency_hz=DESIGN_FREQUENCY_HZ,
        ungrounded_comparison=ungrounded_comparison,
        convergence_points=convergence_points,
        floquet_error=floquet_error,
        mesh_error=mesh_error,
        mesh_floquet_tolerance=gate.MESH_FLOQUET_TOLERANCE,
        slab_check_points=slab_points,
        max_analytic_error=max_analytic_error,
        thin_limit=thin_limit,
        surface_wave_points=surface_wave_points,
        decision=decision,
    )

    specification.write_sweep_csv(sweep_points, output_dir / "sweep.csv")
    specification.write_sweep_json(sweep_points, coverage, output_dir / "sweep.json")
    specification.write_ungrounded_comparison_json(ungrounded_comparison, output_dir / "ungrounded_comparison.json")
    specification.write_analytic_checks_json(slab_points, thin_limit, output_dir / "analytic_checks.json")
    specification.write_surface_wave_json(surface_wave_points, output_dir / "surface_wave_pole.json")
    specification.write_convergence_json(convergence_points, floquet_error, mesh_error, output_dir / "convergence.json")
    specification.write_gate_status_json(decision, output_dir / "phase12a_gate_status.json")

    try:
        plots.generate_all_plots(sweep_points, ungrounded_diagonal, DESIGN_FREQUENCY_HZ, convergence_points, plots_dir)
    except ImportError:
        pass  # matplotlib optional, like every other phase's plotting step

    (output_dir / "phase12a_report.md").write_text(render_phase12a_report_md(result))

    return result


def print_final_report(result):
    d = result.decision
    cov = result.coverage_by_frequency[result.design_frequency_hz]
    print("=" * 78)
    print("PHASE 12A -- PEC GROUND PLANE (GROUNDED-SUBSTRATE STACK)")
    print("=" * 78)
    print(f"Grounded-stack status:   {d.status}")
    print(f"Reason:                  {d.reason}")
    print(f"S11_xx phase coverage:   {cov.coverage_deg:.2f} deg (design frequency)")
    ug = result.ungrounded_comparison
    if ug.get("status") == "AVAILABLE":
        print(f"Ungrounded (Phase 7):    {ug['S11_xx']['coverage_deg']:.2f} deg")
    else:
        print("Ungrounded (Phase 7):    NOT_AVAILABLE_FROM_CURRENT_DATA")
    print(f"HFSS status:             {d.hfss_status}")
    print(f"Geometry inversion:      {d.geometry_inversion_status}")
    print(f"STEP export:             {d.step_export_status}")
    print("-" * 78)
    print(f"Max analytic (bare-slab) error: {result.max_analytic_error:.3g}")
    print(f"Refinement study: floquet_error={result.floquet_error:.4g}  mesh_error={result.mesh_error:.4g}")
    print(f"Surface-wave: any_in_band={any(p.any_in_band for p in result.surface_wave_points)}")
    print("=" * 78)
