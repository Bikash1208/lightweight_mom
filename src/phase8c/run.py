"""Phase 8C orchestrator: runs every analysis module against the real Phase 7 data (read via
`data_source.load_full_database`, never re-solved) and writes `results/phase8c_scattering/`."""

from pathlib import Path

from ..database import CONVENTIONS
from ..phase8b import control_variable as cv
from ..phase8b import library_compatibility as lib_compat
from . import architecture_test as arch
from . import compatibility_matrix as cm
from . import data_source
from . import phase_only_test as pot
from . import plots as p8c_plots
from . import power_balance as pb
from . import reciprocity_symmetry as rs
from . import reference_plane_sensitivity as rps
from . import scattering_analysis as sa
from . import specification as spec
from . import topology_audit as ta
from .report import render_control_quantity_analysis_md

DESIGN_FREQUENCY_GHZ = 12.5
REPRESENTATIVE_CELLS = ((1.0, 1.0), (2.0, 2.0), (3.0, 3.0), (4.0, 4.0))


def run_phase8c(output_dir, checkpoint_path=data_source.DEFAULT_CHECKPOINT_JSON, progress=None, database=None):
    """``database``: an already-loaded `database.Database` (e.g. the session-scoped
    `phase8c_full_database` test fixture) to skip the ~20 s checkpoint parse; if omitted, loads it
    from ``checkpoint_path`` as usual."""
    output_dir = Path(output_dir)
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    if database is None:
        if progress:
            progress("loading full checkpoint database...")
        database = data_source.load_full_database(checkpoint_path)
    check = data_source.cross_check_against_delivered_library(database)
    if not check["ok"]:
        raise RuntimeError(f"Phase 8C: checkpoint database disagrees with the delivered library: {check['mismatches']}")

    freqs = database.frequencies
    f_idx = next(i for i, f in enumerate(freqs) if abs(f / 1e9 - DESIGN_FREQUENCY_GHZ) < 1e-6)
    l_values_mm = sorted({round(c.parameters["length_m"] * 1e3, 9) for c in database.cells})
    w_values_mm = sorted({round(c.parameters["width_m"] * 1e3, 9) for c in database.cells})

    if progress:
        progress("power balance + topology audit...")
    power_summary = pb.library_power_summary(database, f_idx)
    audit = ta.audit_topology(database, power_summary)
    unit_cell = database.unit_cell(next(c.index for c in database.cells if c.status == "ok"))
    mode_check = pb.propagating_mode_check(unit_cell, freqs)

    if progress:
        progress("full-library S11xx/S21xx scattering analysis...")
    s11_rows = sa.extract_component_table(database, "S11_xx", f_idx)
    s21_rows = sa.extract_component_table(database, "S21_xx", f_idx)
    s11_coverage = sa.candidate_coverage(s11_rows, "S11xx")
    s21_coverage = sa.candidate_coverage(s21_rows, "S21xx")
    s11_states = sa.distinct_complex_and_phase_states(s11_rows)
    s21_states = sa.distinct_complex_and_phase_states(s21_rows)
    s11_phase_only = pot.phase_only_justification_test(s11_coverage)
    s21_phase_only = pot.phase_only_justification_test(s21_coverage)

    if progress:
        progress("reciprocity/symmetry checks...")
    reciprocity = rs.check_reciprocity(database, f_idx)
    symmetry = rs.check_lw_swap_symmetry(database, f_idx)

    s11_investigation = cv.investigate_s11xx()
    s21_investigation = cv.investigate_s21xx()
    reflection_test = arch.reflection_architecture_test()
    transmission_test = arch.transmission_architecture_test()
    formulation_a = arch.formulation_a_reevaluation()
    formulation_b = arch.formulation_b_reevaluation()

    architecture_resolved = audit.architecture_classification != "UNRESOLVED"
    status = cm.final_status(
        architecture_resolved=architecture_resolved,
        s11_investigation=s11_investigation, s21_investigation=s21_investigation,
        formulation_a_status=formulation_a["status"],
        reciprocity_ok=reciprocity.reciprocal, propagating_modes_ok=mode_check["only_fundamental_mode_propagates"],
    )

    matrix = cm.build_compatibility_matrix(
        s11_investigation, s21_investigation, s11_coverage, s21_coverage,
        reflection_test, transmission_test, s11_phase_only, s21_phase_only,
    )

    missing_evidence = [
        "HFSS validation data (not performed for any cell)",
        "A quantitative magnitude tolerance for a phase-only approximation, from either Phase 8B formulation",
        "A project decision on reflective-vs-transmissive intended architecture",
    ]

    if progress:
        progress("writing reports...")
    (output_dir / "unit_cell_topology_audit.md").write_text(ta.render_topology_audit_md(audit))
    (output_dir / "s_parameter_definition.md").write_text(ta.render_s_parameter_definition_md())
    report_md = render_control_quantity_analysis_md(
        audit, power_summary, mode_check, s11_investigation, s21_investigation,
        s11_coverage, s21_coverage, s11_states, s21_states, s11_phase_only, s21_phase_only,
        reflection_test, transmission_test, formulation_a, formulation_b, reciprocity, symmetry,
        matrix, status, missing_evidence,
    )
    (output_dir / "control_quantity_analysis.md").write_text(report_md)

    validation_doc = spec.build_control_quantity_validation(status, audit, s11_coverage, s21_coverage, s11_investigation, s21_investigation)
    spec.write_json(validation_doc, output_dir / "control_quantity_validation.json")
    spec.write_scattering_csv(s11_rows, s21_rows, output_dir / "unit_cell_scattering_12p5GHz.csv")
    spec.write_control_quantity_candidates_csv([s11_coverage, s21_coverage], output_dir / "control_quantity_candidates.csv")
    spec.write_phase_state_analysis_csv([("S11xx", s11_states), ("S21xx", s21_states)], output_dir / "phase_state_analysis.csv")

    try:
        if progress:
            progress("plotting...")
        p8c_plots.plot_topology_summary(audit, plots_dir / "01_unit_cell_topology.png")
        p8c_plots.plot_representative_curve(database, REPRESENTATIVE_CELLS, "S11_xx", "mag", plots_dir / "02_representative_S11_magnitude.png")
        p8c_plots.plot_representative_curve(database, REPRESENTATIVE_CELLS, "S11_xx", "phase", plots_dir / "03_representative_S11_phase.png")
        p8c_plots.plot_representative_curve(database, REPRESENTATIVE_CELLS, "S21_xx", "mag", plots_dir / "04_representative_S21_magnitude.png")
        p8c_plots.plot_representative_curve(database, REPRESENTATIVE_CELLS, "S21_xx", "phase", plots_dir / "05_representative_S21_phase.png")
        p8c_plots.plot_library_map(s11_rows, l_values_mm, w_values_mm, lambda r: r.phase_deg, "phase(S11xx) [deg]",
                                    "Library S11xx phase map", plots_dir / "06_library_S11_phase_map.png", -180, 180, "twilight")
        p8c_plots.plot_library_map(s21_rows, l_values_mm, w_values_mm, lambda r: r.phase_deg, "phase(S21xx) [deg]",
                                    "Library S21xx phase map", plots_dir / "07_library_S21_phase_map.png", -180, 180, "twilight")
        p8c_plots.plot_library_map(s11_rows, l_values_mm, w_values_mm, lambda r: r.mag_db, "|S11xx| [dB]",
                                    "Library S11xx magnitude map", plots_dir / "08_library_S11_magnitude_map.png")
        p8c_plots.plot_library_map(s21_rows, l_values_mm, w_values_mm, lambda r: r.mag_db, "|S21xx| [dB]",
                                    "Library S21xx magnitude map", plots_dir / "09_library_S21_magnitude_map.png")
        p8c_plots.plot_complex_plane(s11_rows, "S11xx", plots_dir / "10_S11_complex_plane.png")
        p8c_plots.plot_complex_plane(s21_rows, "S21xx", plots_dir / "11_S21_complex_plane.png")
        p8c_plots.plot_phase_vs_magnitude(s11_rows, "S11xx", plots_dir / "12_S11_phase_vs_magnitude.png")
        p8c_plots.plot_phase_vs_magnitude(s21_rows, "S21xx", plots_dir / "13_S21_phase_vs_magnitude.png")
        from ..phase_extraction import phase_coverage as _cov
        p8c_plots.plot_phase_coverage(s11_rows, _cov([{"phase_deg": r.phase_deg} for r in s11_rows]), "S11xx", plots_dir / "14_S11_phase_coverage.png")
        p8c_plots.plot_phase_coverage(s21_rows, _cov([{"phase_deg": r.phase_deg} for r in s21_rows]), "S21xx", plots_dir / "15_S21_phase_coverage.png")
        p8c_plots.plot_power_balance(power_summary, plots_dir / "16_power_balance.png")
        p8c_plots.plot_phase_state_distribution(s11_states, plots_dir / "17_phase_state_distribution_S11.png", "S11xx")
        p8c_plots.plot_phase_state_distribution(s21_states, plots_dir / "18_phase_state_distribution_S21.png", "S21xx")
    except ImportError:
        pass

    return {
        "status": status, "audit": audit, "mode_check": mode_check, "power_summary": power_summary,
        "s11_investigation": s11_investigation, "s21_investigation": s21_investigation,
        "s11_coverage": s11_coverage, "s21_coverage": s21_coverage,
        "s11_states": s11_states, "s21_states": s21_states,
        "s11_phase_only": s11_phase_only, "s21_phase_only": s21_phase_only,
        "reflection_test": reflection_test, "transmission_test": transmission_test,
        "formulation_a": formulation_a, "formulation_b": formulation_b,
        "reciprocity": reciprocity, "symmetry": symmetry, "matrix": matrix,
        "output_dir": output_dir,
    }


def print_final_report(result):
    a, m = result["audit"], result["mode_check"]
    print("=" * 60)
    print("PHASE 8C -- UNIT-CELL SCATTERING VALIDATION")
    print("=" * 60)
    print()
    print(f"Unit-cell topology: {a.metal_layers}")
    print(f"Ground plane: {a.ground_plane}")
    print(f"Physical architecture: {a.architecture_classification} ({a.architecture_evidence})")
    print("Period: 4.8 x 4.8 mm")
    print(f"Design frequency: {12.5} GHz")
    print("Library geometries: 324")
    print()
    print("-" * 60)
    print("S-PARAMETER DEFINITIONS")
    print("-" * 60)
    print(f"S11: top-port co-pol reflection; S21: transmission port1->port2 ({CONVENTIONS['s_matrix']})")
    print(f"Polarization: {CONVENTIONS['linear_basis']}")
    print(f"Higher-order modes: {m['note']}")
    print()
    print("-" * 60)
    print("REFLECTION CONTROL (S11xx)")
    print("-" * 60)
    print(f"Physical interpretation: {result['s11_investigation'].physical_definition}")
    print(f"Holographic compatibility: PARTIAL/UNVALIDATED ({result['reflection_test'].status})")
    print(f"Phase coverage: {result['s11_coverage'].coverage_deg:.2f} deg")
    print(f"Largest circular gap: {result['s11_coverage'].largest_gap_deg:.2f} deg")
    print(f"Unique phase states: {result['s11_states'].n_distinct_phase_states}")
    print(f"Unique complex states: {result['s11_states'].n_distinct_complex_states}")
    print(f"Magnitude range: [{result['s11_coverage'].magnitude_min:.4f}, {result['s11_coverage'].magnitude_max:.4f}]")
    print(f"Phase-only control: {result['s11_phase_only'].status}")
    print(f"Status: {result['s11_investigation'].validated}")
    print()
    print("-" * 60)
    print("TRANSMISSION CONTROL (S21xx)")
    print("-" * 60)
    print(f"Physical interpretation: {result['s21_investigation'].physical_definition}")
    print(f"Holographic compatibility: PARTIAL/UNVALIDATED ({result['transmission_test'].status})")
    print(f"Phase coverage: {result['s21_coverage'].coverage_deg:.2f} deg")
    print(f"Largest circular gap: {result['s21_coverage'].largest_gap_deg:.2f} deg")
    print(f"Unique phase states: {result['s21_states'].n_distinct_phase_states}")
    print(f"Unique complex states: {result['s21_states'].n_distinct_complex_states}")
    print(f"Magnitude range: [{result['s21_coverage'].magnitude_min:.4f}, {result['s21_coverage'].magnitude_max:.4f}]")
    print(f"Phase-only control: {result['s21_phase_only'].status}")
    print(f"Status: {result['s21_investigation'].validated}")
    print()
    print("-" * 60)
    print("POWER / MODE ANALYSIS")
    print("-" * 60)
    print(f"Power conservation: mean(R+T)={result['power_summary']['mean_R_plus_T']:.4f}")
    print(f"Higher-order modes: {'none in band' if m['only_fundamental_mode_propagates'] else 'PRESENT in band'}")
    print(f"Cross-polarization: R_cross={result['power_summary']['mean_R_cross']:.4f}, T_cross={result['power_summary']['mean_T_cross']:.4f}")
    print(f"Reciprocity: {'OK' if result['reciprocity'].reciprocal else 'FAILED'} ({result['reciprocity'].note})")
    print()
    print("-" * 60)
    print("PHASE 8B COMPATIBILITY")
    print("-" * 60)
    print(f"Formulation A: {result['formulation_a']['status']}")
    print(f"Formulation B: reflection->{result['formulation_b']['corresponds_to_reflection_phase']}, "
          f"transmission->{result['formulation_b']['corresponds_to_transmission_phase']}")
    print()
    print("-" * 60)
    print("LIBRARY ADEQUACY")
    print("-" * 60)
    print(f"S11xx: {cm_adequacy(result['s11_investigation'].validated)}")
    print(f"S21xx: {cm_adequacy(result['s21_investigation'].validated)}")
    print()
    print("-" * 60)
    print("VALIDATION STATUS")
    print("-" * 60)
    print("HFSS validation: NOT PERFORMED")
    print("Solver modified: NO")
    print("Phase 7 modified: NO")
    print("Interpolation: NO")
    print("Extrapolation: NO")
    print("Artificial phase offset: NO")
    print("Artificial phase expansion: NO")
    print()
    print("-" * 60)
    print("FINAL SCIENTIFIC STATUS")
    print("-" * 60)
    print()
    print(result["status"])
    print("=" * 60)


def cm_adequacy(validated):
    from .compatibility_matrix import library_adequacy

    return library_adequacy(validated, None)
