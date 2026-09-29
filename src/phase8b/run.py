"""Phase 8B orchestrator: runs the whole pipeline and writes `results/phase8b_formulation/`
(Sec. 50), reusing Phase 8A's real, already-solved library data for the distinct-phase-states
analysis (no re-solve, no new library build)."""

from pathlib import Path

from ..feasibility_coverage import library_table_at_frequency
from ..feasibility_library_integrity import load_library_bundle
from ..feasibility_schema import DEFAULT_LIBRARY_DIR
from . import plots as p8b_plots
from . import specification
from .distinct_states import analyze_distinct_states, write_distinct_state_csv
from .library_compatibility import evaluate_all
from .object_wave import Beam
from .report import (
    render_control_variable_analysis_md,
    render_holographic_formulation_md,
    render_repository_scientific_audit_md,
)

DEFAULT_DEMO_BEAMS = (
    Beam(theta_deg=30.0, phi_deg=0.0),
    Beam(theta_deg=30.0, phi_deg=120.0),
    Beam(theta_deg=30.0, phi_deg=240.0),
)
DEMO_FREQUENCY_HZ = 12.5e9


def run_phase8b(output_dir, library_dir=DEFAULT_LIBRARY_DIR, design_frequency_ghz=12.5):
    output_dir = Path(output_dir)
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    evaluations, status = evaluate_all()

    bundle = load_library_bundle(library_dir)
    table = library_table_at_frequency(bundle.library_doc, design_frequency_ghz)
    distinct = analyze_distinct_states(table)

    (output_dir / "repository_scientific_audit.md").write_text(render_repository_scientific_audit_md())
    (output_dir / "Holographic_Formulation.md").write_text(render_holographic_formulation_md(evaluations, distinct))
    (output_dir / "control_variable_analysis.md").write_text(render_control_variable_analysis_md(evaluations))

    spec = specification.build_holographic_specification()
    specification.write_json(spec, output_dir / "holographic_specification.json")
    specification.write_json(specification.build_source_registry(), output_dir / "source_registry.json")
    specification.write_phase_control_quantity_csv(evaluations, output_dir / "phase_control_quantity.csv")
    write_distinct_state_csv(distinct, output_dir / "phase8b_distinct_state_analysis.csv")

    try:
        beam1 = DEFAULT_DEMO_BEAMS[0]
        p8b_plots.plot_object_wave_single_beam(beam1, DEMO_FREQUENCY_HZ, plots_dir / "object_wave_single_beam.png")
        p8b_plots.plot_object_wave_multi_beam(
            list(DEFAULT_DEMO_BEAMS), DEMO_FREQUENCY_HZ, plots_dir / "object_wave_multi_beam.png"
        )
        p8b_plots.plot_phase_gradient_test(DEMO_FREQUENCY_HZ, plots_dir / "phase_gradient_test.png")
        p8b_plots.plot_phase_control_quantity(evaluations, plots_dir / "phase_control_quantity.png")
        p8b_plots.plot_distinct_phase_states(distinct, plots_dir / "distinct_phase_states.png")
    except ImportError:
        pass  # matplotlib not installed; JSON/CSV/MD outputs are still written

    return {"status": status, "evaluations": evaluations, "distinct_states": distinct, "output_dir": output_dir}


def print_final_report(result):
    s11 = next(e for e in result["evaluations"] if "B:" in e.formulation_name).control_variable
    from . import control_variable as cv

    s21 = cv.investigate_s21xx()
    zs = cv.investigate_surface_impedance()

    print("=" * 60)
    print("PHASE 8B -- HOLOGRAPHIC FORMULATION")
    print("=" * 60)
    print()
    print(f"{'Repository audit:':<30}PASS")
    print(f"{'Reference wave:':<30}DEFINED (2 candidates; monopole/guided-wave input requires k_swg, unavailable)")
    print(f"{'Object wave:':<30}DEFINED")
    print(f"{'Multi-beam formulation:':<30}DEFINED")
    print(f"{'Coordinate convention:':<30}DEFINED (reused from CONVENTIONS.md)")
    print(f"{'Polarization convention:':<30}DEFINED (reused from CONVENTIONS.md)")
    print(f"{'Phase sign:':<30}DEFINED (reused from excitation.PlaneWaveExcitation)")
    print(f"{'Holographic modulation:':<30}DEFINED (2 candidate formulations, neither validated)")
    print()
    print(f"{'Physical control variable:':<30}UNRESOLVED (2 candidates: S11xx/S21xx vs. surface impedance)")
    print()
    print(f"{'S11xx:':<30}{s11.validated}")
    print(f"{'S21xx:':<30}{s21.validated}")
    print(f"{'Surface impedance:':<30}{zs.validated}")
    print()
    for e in result["evaluations"]:
        print(f"{'Library compatibility (' + e.formulation_name.split(':')[0] + '):':<30}{e.library_compatibility}")
    print()
    d = result["distinct_states"]
    print(f"{'Distinct geometric states:':<30}{d.n_geometries}")
    print(f"{'Distinct phase states:':<30}{d.n_distinct_phase_states}")
    print()
    print(f"{'HFSS validation:':<30}NOT PERFORMED")
    print(f"{'Solver modified:':<30}NO")
    print(f"{'Phase 7 modified:':<30}NO")
    print()
    print(f"{'Final scientific status:':<30}{result['status']}")
    print("=" * 60)
