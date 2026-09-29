"""Phase 8E orchestration: reference field -> object wave -> candidate formulations -> required
control variable -> architecture compatibility -> reference-plane compatibility -> complex
controllability -> universal-database query -> continuous/discrete beam reconstruction ->
geometry-inversion gate (Sec. 29's fixed order). Nothing here modifies Phase 7 or Phase 8D-REV
outputs; no geometry inversion or STEP generation is performed (Sec. 31).
"""

import json
from pathlib import Path

import numpy as np

from ..phase8b.object_wave import combine_beams, combine_beams_scalar
from ..phase8b.run import DEFAULT_DEMO_BEAMS
from ..phase8c import data_source, power_balance
from ..phase8c import topology_audit as topo_mod
from ..phase8c.scattering_analysis import extract_component_table
from ..phase8d_rev import library_entry as le
from ..phase8d_rev.query import UnitCellQuery, search_unit_cells
from ..phase8d_rev.response_matching import nearest_response_match
from ..phase8d_rev.universal_database import UniversalDatabase
from ..phase_extraction import design_frequency_index
from . import architecture_compatibility as arch_mod
from . import beam_reconstruction as br_mod
from . import candidate_library_period as clp_mod
from . import complex_controllability as cc_mod
from . import control_response as cr_mod
from . import decision_tree as dt_mod
from . import local_periodicity as lp_mod
from . import plots as p8e_plots
from . import provenance as prov_mod
from . import reference_field_input as rfi
from . import reference_plane as refplane
from . import report as report_mod
from . import specification as spec_out
from . import spectral as spec_mod

DESIGN_FREQUENCY_HZ = 12.5e9
DEMO_FEED_HEIGHT_M = 0.060  # 60 mm, synthetic demo feed height; not real hardware data
APERTURE_N_CELLS = 15  # per side, odd, centred on the aperture
SPECTRAL_GRID_N = 64
SPECTRAL_GRID_EXTENT_M = 0.032
DEMO_MISMATCH_Z_M = 1.7e-3  # Sec. 3's historical HFSS grid height, demo only
FAR_FIELD_THETA_DEG = np.arange(0.0, 61.0, 2.0)
FAR_FIELD_PHI_DEG = np.arange(0.0, 360.0, 4.0)


def _discrete_aperture_field(e_ref_ap, required_phase_deg, rows):
    """Sec. 19: for each aperture point, find the nearest achievable library complex value to the
    idealized unit-magnitude target `exp(j*required_phase)`, and build the discrete aperture field
    with that REAL (non-unit-magnitude) response -- the actual source of continuous/discrete
    disagreement (Sec. 19: 'do not hide failures behind averaged phase error')."""
    target = np.exp(1j * np.radians(np.asarray(required_phase_deg, dtype=float)))
    candidates = [(f"L{r.l_mm:g}_W{r.w_mm:g}", r.l_mm, r.w_mm, "lib", complex(r.real, r.imag)) for r in rows]
    matched = np.array([nearest_response_match(t, candidates)[0].complex_value for t in target.ravel()]).reshape(
        target.shape
    )
    return np.asarray(e_ref_ap, dtype=complex) * matched, matched


def run_phase8e(
    output_dir, checkpoint_path=data_source.DEFAULT_CHECKPOINT_JSON, library_dir=le.DEFAULT_LIBRARY_DIR, database=None
):
    output_dir = Path(output_dir)
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    plot_paths = {}

    # ---------------------------------------------------------------- real Phase 7 data (read-only)
    if database is None:
        database = data_source.load_full_database(checkpoint_path)
    manifest_path = Path(library_dir) / "manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else None
    entry = le.from_database(
        database, le.BENCHMARK_LIBRARY_ID, le.BENCHMARK_TOPOLOGY, le.REFERENCE_VALIDATION_LIBRARY, manifest
    )
    udb = UniversalDatabase()
    udb.register(entry)

    freq_idx = design_frequency_index(database, DESIGN_FREQUENCY_HZ)
    mean_power = power_balance.library_power_summary(database, freq_idx)
    audit = topo_mod.audit_topology(database, mean_power)
    s11_rows = extract_component_table(database, "S11_xx", freq_idx)
    s21_rows = extract_component_table(database, "S21_xx", freq_idx)

    # ------------------------------------------------------------------- reference field (Sec. 3-4)
    lin = np.linspace(-SPECTRAL_GRID_EXTENT_M, SPECTRAL_GRID_EXTENT_M, SPECTRAL_GRID_N)
    xg, yg = np.meshgrid(lin, lin, indexing="xy")
    ref_field_fine = rfi.spherical_feed_reference_field(xg, yg, DEMO_FEED_HEIGHT_M, DESIGN_FREQUENCY_HZ)
    ref_spectrum = spec_mod.reference_field_spectrum(ref_field_fine, DESIGN_FREQUENCY_HZ)

    plane_compat_actual = refplane.check_reference_plane_compatibility(ref_field_fine.z_reference_m)
    plane_compat_demo_mismatch = refplane.check_reference_plane_compatibility(DEMO_MISMATCH_Z_M)
    _ = refplane.propagate_angular_spectrum(
        ref_field_fine.Ex,
        ref_field_fine.dx_m,
        ref_field_fine.dy_m,
        DESIGN_FREQUENCY_HZ,
        DEMO_MISMATCH_Z_M,
    )  # Sec. 6 demonstration only; the pipeline's own reference field is already z=0-compatible

    # ---------------------------------------------------------------------- object wave (Sec. 4, 8)
    beams = list(DEFAULT_DEMO_BEAMS)
    obj_spectrum = spec_mod.object_wave_spectrum(beams, DESIGN_FREQUENCY_HZ)

    # ------------------------------------------------------- aperture grid at the unit-cell pitch
    p_x_m, p_y_m = entry.p_x_mm * 1e-3, entry.p_y_mm * 1e-3
    half = (APERTURE_N_CELLS - 1) / 2
    ax = (np.arange(APERTURE_N_CELLS) - half) * p_x_m
    ay = (np.arange(APERTURE_N_CELLS) - half) * p_y_m
    ax_grid, ay_grid = np.meshgrid(ax, ay, indexing="xy")
    e_ref_ap = rfi.spherical_feed_reference_field(ax_grid, ay_grid, DEMO_FEED_HEIGHT_M, DESIGN_FREQUENCY_HZ).Ex
    e_obj_ap = combine_beams_scalar(ax_grid, ay_grid, beams, DESIGN_FREQUENCY_HZ)
    ex_obj_vec, ey_obj_vec = combine_beams(ax_grid, ay_grid, beams, DESIGN_FREQUENCY_HZ)

    # ---------------------------------------------------- provenance & architecture compatibility
    provenance_registry = prov_mod.build_provenance_registry()
    architecture_results = arch_mod.assess_all(audit, mean_power)

    # ---------------------------------------------------------------- control responses (Sec. 9-10)
    cr_b = cr_mod.required_phase_b(e_ref_ap, e_obj_ap)
    cr_c = cr_mod.required_phase_c(e_ref_ap, e_obj_ap)
    cr_d = cr_mod.generalized_transfer_operator(e_ref_ap, np.zeros_like(e_ref_ap), ex_obj_vec, ey_obj_vec)
    control_responses = {"B": cr_b, "C": cr_c, "D": cr_d}

    # ----------------------------------------------------------- complex controllability (Sec. 13-14)
    target_unit = np.exp(1j * np.radians(cr_b.values))  # identical target array for B and C (same kernel)
    cc_b = cc_mod.complex_plane_controllability(target_unit, s11_rows)
    cc_c = cc_mod.complex_plane_controllability(target_unit, s21_rows)
    lib_summary_b = cc_mod.library_complex_plane_summary(s11_rows, "S11_xx")
    lib_summary_c = cc_mod.library_complex_plane_summary(s21_rows, "S21_xx")

    # ------------------------------------------------------------- local periodicity (Sec. 15)
    local_periodicity = lp_mod.assess_local_periodicity(p_x_m, obj_spectrum.k_required_max_per_m)

    # ------------------------------------------------------------------ period requirement (Sec. 16-17)
    period_assessment = clp_mod.assess_period_requirement(
        udb,
        entry.topology,
        obj_spectrum.k_required_max_per_m,
        DESIGN_FREQUENCY_HZ,
        entry.substrate["name"],
    )

    # --------------------------------------------------------------- family existence (Sec. 20's box)
    family_query_result = search_unit_cells(
        udb, UnitCellQuery(topology=entry.topology, frequency_hz=DESIGN_FREQUENCY_HZ)
    )
    family_exists = family_query_result.status == "MATCH"

    # ------------------------------------------------------------ continuous/discrete reconstruction
    theta_grid, phi_grid = np.meshgrid(FAR_FIELD_THETA_DEG, FAR_FIELD_PHI_DEG, indexing="ij")

    a_cont_b = e_ref_ap * np.exp(1j * np.radians(cr_b.values))
    a_disc_b, matched_b = _discrete_aperture_field(e_ref_ap, cr_b.values, s11_rows)
    ff_cont_b = br_mod.far_field_pattern(
        ax_grid, ay_grid, a_cont_b, DESIGN_FREQUENCY_HZ, theta_grid, phi_grid, p_x_m, p_y_m
    )
    ff_disc_b = br_mod.far_field_pattern(
        ax_grid, ay_grid, a_disc_b, DESIGN_FREQUENCY_HZ, theta_grid, phi_grid, p_x_m, p_y_m
    )
    metrics_cont_b = br_mod.evaluate_beam_reconstruction(ff_cont_b, theta_grid, phi_grid, beams)
    metrics_disc_b = br_mod.evaluate_beam_reconstruction(ff_disc_b, theta_grid, phi_grid, beams)
    comparison_b = br_mod.compare_continuous_vs_discrete(ff_cont_b, ff_disc_b)

    a_cont_c = e_ref_ap * np.exp(1j * np.radians(cr_c.values))
    a_disc_c, matched_c = _discrete_aperture_field(e_ref_ap, cr_c.values, s21_rows)
    ff_cont_c = br_mod.far_field_pattern(
        ax_grid, ay_grid, a_cont_c, DESIGN_FREQUENCY_HZ, theta_grid, phi_grid, p_x_m, p_y_m
    )
    ff_disc_c = br_mod.far_field_pattern(
        ax_grid, ay_grid, a_disc_c, DESIGN_FREQUENCY_HZ, theta_grid, phi_grid, p_x_m, p_y_m
    )
    metrics_cont_c = br_mod.evaluate_beam_reconstruction(ff_cont_c, theta_grid, phi_grid, beams)
    metrics_disc_c = br_mod.evaluate_beam_reconstruction(ff_disc_c, theta_grid, phi_grid, beams)
    comparison_c = br_mod.compare_continuous_vs_discrete(ff_cont_c, ff_disc_c)

    def _all_finite(metrics):
        return all(np.isfinite(e) for e in metrics.pointing_error_deg)

    # Sec. 20's remaining two checkboxes, sourced from the actual, already-computed Phase 8B/8E
    # results (never hardcoded inside decision_tree.py itself): Phase 8B never reached
    # FORMULATION_VALIDATED for either candidate formulation it evaluates (A or the phase-
    # difference hologram B/C share); Candidate D has no formulation-level evaluation anywhere.
    from ..phase8b.library_compatibility import evaluate_formulation_a, evaluate_formulation_b

    formulation_validated_a = evaluate_formulation_a().final_status == "FORMULATION_VALIDATED"
    formulation_validated_bc = evaluate_formulation_b().final_status == "FORMULATION_VALIDATED"
    local_periodicity_ok = local_periodicity.status == "LOCAL_PERIODICITY_VALIDATED"

    # ---------------------------------------------------------------------- decision tree (Sec. 20, 29)
    decision_a = dt_mod.evaluate_formulation_gate(
        "A",
        architecture_results["A"].status,
        "UNKNOWN_REFERENCE_PLANE",
        False,
        False,
        False,
        False,
        formulation_validated=formulation_validated_a,
        local_periodicity_acceptable=local_periodicity_ok,
    )
    decision_b = dt_mod.evaluate_formulation_gate(
        "B",
        architecture_results["B"].status,
        plane_compat_actual.status,
        family_exists,
        _all_finite(metrics_cont_b),
        _all_finite(metrics_disc_b),
        False,
        formulation_validated=formulation_validated_bc,
        local_periodicity_acceptable=local_periodicity_ok,
    )
    decision_c = dt_mod.evaluate_formulation_gate(
        "C",
        architecture_results["C"].status,
        plane_compat_actual.status,
        family_exists,
        _all_finite(metrics_cont_c),
        _all_finite(metrics_disc_c),
        False,
        formulation_validated=formulation_validated_bc,
        local_periodicity_acceptable=local_periodicity_ok,
    )
    decision_d = dt_mod.evaluate_formulation_gate(
        "D",
        architecture_results["D"].status,
        "UNKNOWN_REFERENCE_PLANE",
        False,
        False,
        False,
        False,
        local_periodicity_acceptable=local_periodicity_ok,
    )
    decisions = [decision_a, decision_b, decision_c, decision_d]
    overall_status = dt_mod.overall_project_status(decisions)

    # ----------------------------------------------------------------------------------- plots (Sec. 23)
    try:
        p8e_plots.plot_reference_field_amplitude_phase(
            ref_field_fine, plots_dir / "01_reference_field_amplitude_phase.png"
        )
        p8e_plots.plot_reference_field_components(ref_field_fine, plots_dir / "02_reference_field_components.png")
        p8e_plots.plot_reference_field_spectrum(
            ref_field_fine, ref_spectrum, plots_dir / "03_reference_field_spectrum.png"
        )
        p8e_plots.plot_object_wave_amplitude_phase(
            ax_grid, ay_grid, e_obj_ap, plots_dir / "04_object_wave_amplitude_phase.png"
        )
        p8e_plots.plot_object_wave_spectrum(
            obj_spectrum, ref_spectrum.k0_per_m, plots_dir / "05_object_wave_spectrum.png"
        )
        p8e_plots.plot_control_response_amplitude_phase(
            ax_grid,
            ay_grid,
            np.exp(1j * np.radians(cr_b.values)),
            "C_req (B/C)",
            plots_dir / "06_control_response_amplitude_phase.png",
        )
        p8e_plots.plot_control_response_complex_plane(
            np.exp(1j * np.radians(cr_b.values)), "B/C", plots_dir / "07_control_response_complex_plane.png"
        )
        p8e_plots.plot_library_complex_state_scatter(
            s11_rows, "S11_xx", plots_dir / "08_library_complex_state_scatter.png"
        )
        p8e_plots.plot_library_magnitude_vs_geometry(
            s11_rows, "S11_xx", plots_dir / "09_library_magnitude_vs_geometry.png"
        )
        p8e_plots.plot_library_phase_vs_geometry(s11_rows, "S11_xx", plots_dir / "10_library_phase_vs_geometry.png")
        p8e_plots.plot_amplitude_phase_relationship(
            s11_rows, "S11_xx", plots_dir / "11_library_amplitude_phase_relationship.png"
        )
        p8e_plots.plot_aperture_field(
            ax_grid, ay_grid, a_cont_b, "continuous aperture (B)", plots_dir / "12_continuous_aperture.png"
        )
        p8e_plots.plot_aperture_field(
            ax_grid, ay_grid, a_disc_b, "discrete aperture (B)", plots_dir / "13_discrete_aperture.png"
        )
        p8e_plots.plot_far_field(
            theta_grid, phi_grid, ff_cont_b, "continuous, B", plots_dir / "14_continuous_far_field.png"
        )
        p8e_plots.plot_far_field(
            theta_grid, phi_grid, ff_disc_b, "discrete, B", plots_dir / "15_discrete_far_field.png"
        )
        phi0_col = int(np.argmin(np.abs(FAR_FIELD_PHI_DEG - 0.0)))
        cont_cut_db = 10 * np.log10(np.abs(ff_cont_b[:, phi0_col]) ** 2 / np.max(np.abs(ff_cont_b) ** 2) + 1e-12)
        disc_cut_db = 10 * np.log10(np.abs(ff_disc_b[:, phi0_col]) ** 2 / np.max(np.abs(ff_disc_b) ** 2) + 1e-12)
        p8e_plots.plot_beam_cuts(FAR_FIELD_THETA_DEG, cont_cut_db, disc_cut_db, 0.0, plots_dir / "16_beam_cuts.png")
        p8e_plots.plot_beam_pointing_comparison(
            beams, metrics_cont_b, metrics_disc_b, plots_dir / "17_beam_pointing_comparison.png"
        )
        plot_paths = {p.name: p for p in plots_dir.glob("*.png")}
    except ImportError:
        pass

    # --------------------------------------------------------------------------- outputs (Sec. 24)
    spec_out.write_formulation_candidates_json(
        provenance_registry,
        control_responses,
        architecture_results,
        output_dir / "formulation_candidates.json",
    )
    spec_out.write_formulation_comparison_csv(
        [
            (
                "A",
                provenance_registry["A"].applicability,
                architecture_results["A"].status,
                "UNKNOWN_REFERENCE_PLANE",
                False,
                None,
                None,
                None,
                decision_a.final_status,
            ),
            (
                "B",
                "reflection phase",
                architecture_results["B"].status,
                plane_compat_actual.status,
                family_exists,
                cc_b.mean_complex_error,
                cc_b.mean_amplitude_error,
                cc_b.mean_phase_error_deg,
                decision_b.final_status,
            ),
            (
                "C",
                "transmission phase",
                architecture_results["C"].status,
                plane_compat_actual.status,
                family_exists,
                cc_c.mean_complex_error,
                cc_c.mean_amplitude_error,
                cc_c.mean_phase_error_deg,
                decision_c.final_status,
            ),
            (
                "D",
                "generalized transfer operator",
                architecture_results["D"].status,
                "UNKNOWN_REFERENCE_PLANE",
                False,
                None,
                None,
                None,
                decision_d.final_status,
            ),
        ],
        output_dir / "formulation_comparison.csv",
    )
    spec_out.write_control_variable_analysis_json(decisions, output_dir / "control_variable_analysis.json")
    kx = np.fft.fftshift(np.fft.fftfreq(SPECTRAL_GRID_N, d=ref_field_fine.dx_m)) * 2 * np.pi
    ky = np.fft.fftshift(np.fft.fftfreq(SPECTRAL_GRID_N, d=ref_field_fine.dy_m)) * 2 * np.pi
    kx_grid, ky_grid = np.meshgrid(kx, ky, indexing="xy")
    spectrum_power = np.abs(np.fft.fftshift(np.fft.fft2(ref_field_fine.Ex))) ** 2
    spec_out.write_reference_spectrum_csv(kx_grid, ky_grid, spectrum_power, output_dir / "reference_spectrum.csv")
    spec_out.write_object_spectrum_csv(obj_spectrum, output_dir / "object_spectrum.csv")
    spec_out.write_required_response_map_csv(
        [("B", ax_grid, ay_grid, cr_b.values), ("C", ax_grid, ay_grid, cr_c.values)],
        output_dir / "required_response_map.csv",
    )
    spec_out.write_library_response_map_csv(
        [(s11_rows, "S11_xx"), (s21_rows, "S21_xx")],
        output_dir / "library_response_map.csv",
    )
    spec_out.write_reconstruction_metrics_json(
        {
            "B": {
                "continuous_pointing_error_deg": list(metrics_cont_b.pointing_error_deg),
                "discrete_pointing_error_deg": list(metrics_disc_b.pointing_error_deg),
                "continuous_matched_peak_power_db": list(metrics_cont_b.matched_peak_power_db),
                "discrete_matched_peak_power_db": list(metrics_disc_b.matched_peak_power_db),
                "comparison_complex_error_mean": comparison_b.complex_error_mean,
                "comparison_amplitude_error_mean": comparison_b.amplitude_error_mean,
                "comparison_phase_error_deg_mean": comparison_b.phase_error_deg_mean,
            },
            "C": {
                "continuous_pointing_error_deg": list(metrics_cont_c.pointing_error_deg),
                "discrete_pointing_error_deg": list(metrics_disc_c.pointing_error_deg),
                "continuous_matched_peak_power_db": list(metrics_cont_c.matched_peak_power_db),
                "discrete_matched_peak_power_db": list(metrics_disc_c.matched_peak_power_db),
                "comparison_complex_error_mean": comparison_c.complex_error_mean,
                "comparison_amplitude_error_mean": comparison_c.amplitude_error_mean,
                "comparison_phase_error_deg_mean": comparison_c.phase_error_deg_mean,
            },
        },
        output_dir / "reconstruction_metrics.json",
    )
    spec_out.write_required_library_specification_json(
        period_assessment, output_dir / "required_library_specification.json"
    )

    result = {
        "output_dir": output_dir,
        "audit": audit,
        "mean_power": mean_power,
        "ref_spectrum": ref_spectrum,
        "obj_spectrum": obj_spectrum,
        "plane_compat_actual": plane_compat_actual,
        "plane_compat_demo_mismatch": plane_compat_demo_mismatch,
        "provenance_registry": provenance_registry,
        "architecture_results": architecture_results,
        "control_responses": control_responses,
        "cc_b": cc_b,
        "cc_c": cc_c,
        "lib_summary_b": lib_summary_b,
        "lib_summary_c": lib_summary_c,
        "local_periodicity": local_periodicity,
        "period_assessment": period_assessment,
        "family_exists": family_exists,
        "metrics_cont_b": metrics_cont_b,
        "metrics_disc_b": metrics_disc_b,
        "comparison_b": comparison_b,
        "metrics_cont_c": metrics_cont_c,
        "metrics_disc_c": metrics_disc_c,
        "comparison_c": comparison_c,
        "decisions": decisions,
        "overall_status": overall_status,
        "entry": entry,
    }
    (output_dir / "report.md").write_text(report_mod.render_report_md(result))
    return result


def print_final_report(result):
    print("=" * 78)
    print("PHASE 8E -- PHYSICAL CONTROL-VARIABLE AND HOLOGRAPHIC SCATTERING BRIDGE")
    print("=" * 78)
    print("Reference field:          SYNTHETIC_DEMO_SPHERICAL_FEED (no real reference-field data exists)")
    print(f"Reference-plane status:   {result['plane_compat_actual'].status}")
    print(
        "Architecture (A/B/C/D):   " + ", ".join(f"{k}={v.status}" for k, v in result["architecture_results"].items())
    )
    print(f"Local periodicity:        {result['local_periodicity'].status}")
    print(f"Period requirement:       {result['period_assessment'].status}")
    print("HFSS validation:          NOT PERFORMED")
    print("Solver/Phase 7 modified:  NO")
    print("Phase 8D-REV modified:    NO")
    print("-" * 78)
    for d in result["decisions"]:
        print(
            f"  [{d.formulation}] final_status={d.final_status}  geometry_inversion_permitted={d.geometry_inversion_permitted}"
        )
    print("-" * 78)
    print("FINAL SCIENTIFIC STATUS")
    print(f"  {result['overall_status']}")
    print("=" * 78)
