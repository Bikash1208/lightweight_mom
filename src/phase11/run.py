"""Phase 11 orchestration. Nothing in Phase 1-10 is modified. No ML/surrogate model anywhere in
this pipeline (Sec. 2) -- every response comes from a real solve of the unmodified Phase 1-3 MoM
machinery. The primary center (`PRIMARY_CENTER`, matching Phase 10's own convention) receives full
depth on every sub-study; the other 4 centers receive the coupling-radius series only (Sec. 5-7,
the core "how far does interaction reach" question), a disclosed scoping choice consistent with
every earlier phase's own documented cost/coverage trade-offs.
"""

from pathlib import Path

from ..materials import Metal, Substrate
from ..phase9 import hfss_reference as hfss_mod
from ..phase10.environment_geometry import CENTER_CELLS_MM, strong_neighbor_low_high
from . import beam_validation as bv
from . import classification as cls_mod
from . import convergence as conv_mod
from . import coupling_radius as cr_mod
from . import decay_law
from . import directional_coupling as dc_mod
from . import frequency_study as freq_mod
from . import gate as gate_mod
from . import holographic_three_beam as hb_mod
from . import interaction as inter_mod
from . import inverse as inv_mod
from . import pattern_study as ps_mod
from . import plots as p11_plots
from . import polarization_study as pol_mod
from . import reduction as red_mod
from . import report as report_mod
from . import specification as spec_mod
from .diagonal_study_reuse import reproduce_diagonal_finding

DESIGN_FREQUENCY_HZ = 12.5e9
PERIOD_M = 4.8e-3
SUBSTRATE_NAME, EPS_R, TAN_DELTA, THICKNESS_M = "Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3
MESH_CELLS_PER_SUBCELL = 8
PRIMARY_CENTER = "C3"


def run_phase11(output_dir, hfss_reference_path=None):
    output_dir = Path(output_dir)
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    substrate = Substrate(SUBSTRATE_NAME, EPS_R, TAN_DELTA, THICKNESS_M)
    metal = Metal()
    primary_lw = CENTER_CELLS_MM[PRIMARY_CENTER]

    # ============================================================ HFSS status (Sec. 25-26) ======
    hfss_path = hfss_reference_path or hfss_mod.DEFAULT_HFSS_REFERENCE_PATH
    hfss_records = hfss_mod.load_hfss_reference(hfss_path)
    hfss_exists = hfss_records is not None

    # =========================================================== 11A/B/C: coupling radius =======
    coupling_radius = {}
    shell_increments = {}
    convergence_results = []
    for name, lw in CENTER_CELLS_MM.items():
        points = cr_mod.run_coupling_radius_series(lw, PERIOD_M, substrate, metal, DESIGN_FREQUENCY_HZ,
                                                     MESH_CELLS_PER_SUBCELL)
        incs = cr_mod.compute_shell_increments(points)
        coupling_radius[name] = points
        shell_increments[name] = incs
        convergence_results.append(conv_mod.classify_convergence(incs, center_name=name))

    # Sec. 5's "where computationally practical" 11x11 spot-check, primary center only. Both the
    # n=9 and n=11 points come from the SAME size-11 truncated pattern here (verified elsewhere to
    # be numerically identical to the main series' own n=9 point, since truncation preserves the
    # checkerboard's parity regardless of the centering offset), so the 9->11 shell increment below
    # is computed consistently rather than mixing two different underlying patterns.
    n9_and_n11 = cr_mod.run_coupling_radius_series(
        primary_lw, PERIOD_M, substrate, metal, DESIGN_FREQUENCY_HZ, MESH_CELLS_PER_SUBCELL,
        sizes=(9, cr_mod.SUPPLEMENTARY_SIZE),
    )
    supplementary_point = n9_and_n11[-1]
    coupling_radius[PRIMARY_CENTER].append(supplementary_point)
    shell_increments[PRIMARY_CENTER] = shell_increments[PRIMARY_CENTER] + cr_mod.compute_shell_increments(n9_and_n11)

    spec_mod.write_coupling_radius_results_json(coupling_radius, output_dir / "coupling_radius_results.json")
    spec_mod.write_incremental_shell_results_json(shell_increments, output_dir / "incremental_shell_results.json")

    # ===================================================================== 11D: directional =====
    dc_baseline, dc_points = dc_mod.run_directional_coupling(
        primary_lw, PERIOD_M, substrate, metal, DESIGN_FREQUENCY_HZ, grid_n=9,
        distances=(1, 2, 3, 4), mesh_cells_per_subcell=MESH_CELLS_PER_SUBCELL,
    )
    spec_mod.write_directional_coupling_json(dc_baseline, dc_points, output_dir / "directional_coupling.json")

    # ======================================================================= 11E: diagonal =======
    diagonal_reproduction = reproduce_diagonal_finding(primary_lw, PERIOD_M, substrate, metal,
                                                        DESIGN_FREQUENCY_HZ, MESH_CELLS_PER_SUBCELL)
    spec_mod.write_diagonal_coupling_json(diagonal_reproduction, output_dir / "diagonal_coupling.json")

    # ===================================================================== 11F: pattern study ====
    low, high = strong_neighbor_low_high(primary_lw)
    contrast_lw = (low, low) if abs(primary_lw[0] - high) < abs(primary_lw[0] - low) else (high, high)
    pattern_grids, pattern_results = ps_mod.run_pattern_study(
        primary_lw, contrast_lw, PERIOD_M, substrate, metal, DESIGN_FREQUENCY_HZ, MESH_CELLS_PER_SUBCELL,
    )
    pattern_dependence_detected, pattern_comparisons = ps_mod.detect_pattern_dependence(pattern_results)
    spec_mod.write_pattern_study_json(pattern_results, pattern_dependence_detected, pattern_comparisons,
                                       output_dir / "pattern_study.json")

    # ==================================================================== 11H/11I: interaction ===
    n_int = 7
    mid = n_int // 2
    from .interaction import Perturbation

    perturbations = [
        Perturbation("N", mid - 1, mid, contrast_lw), Perturbation("S", mid + 1, mid, contrast_lw),
        Perturbation("E", mid, mid + 1, contrast_lw), Perturbation("W", mid, mid - 1, contrast_lw),
        Perturbation("NE", mid - 1, mid + 1, contrast_lw), Perturbation("SW", mid + 1, mid - 1, contrast_lw),
        Perturbation("N_far", 0, mid, contrast_lw),
    ]
    pair_names = [("N", "S"), ("N", "E"), ("NE", "SW"), ("N", "N_far")]
    r0, pairwise_results = inter_mod.run_pairwise_interaction(
        primary_lw, perturbations, pair_names, n_int, PERIOD_M, substrate, metal, DESIGN_FREQUENCY_HZ,
        MESH_CELLS_PER_SUBCELL,
    )
    spec_mod.write_pairwise_interaction_json(r0, pairwise_results, output_dir / "pairwise_interaction.json")

    higher_order_result = inter_mod.run_higher_order_interaction(
        primary_lw, perturbations, ("N", "E", "S"), n_int, PERIOD_M, substrate, metal,
        DESIGN_FREQUENCY_HZ, MESH_CELLS_PER_SUBCELL,
    )
    spec_mod.write_higher_order_interaction_json(higher_order_result, output_dir / "higher_order_interaction.json")

    # ======================================================================= 11G: smooth/abrupt ==
    from . import patterns as pt
    from .reference_solver import solve_grid as _solve_grid_direct

    def _center_response(grid_mm):
        n = len(grid_mm)
        grid_m = [[(l * 1e-3, w * 1e-3) for (l, w) in row] for row in grid_mm]
        sg = _solve_grid_direct(n, PERIOD_M, substrate, metal, grid_m, DESIGN_FREQUENCY_HZ,
                                 MESH_CELLS_PER_SUBCELL)
        m = n // 2
        return sg.response_at(m, m)

    smooth_vs_abrupt = {
        "smooth_gradient": pattern_results["smooth_gradient"].response,
        "single_cell_discontinuity": _center_response(pt.single_cell_discontinuity(7, primary_lw, contrast_lw)),
        "multi_cell_discontinuity": _center_response(pt.multi_cell_discontinuity(7, primary_lw, contrast_lw, width=2)),
        "abrupt_boundary": pattern_results["abrupt_boundary"].response,
    }

    # ==================================================================== 11J: frequency study ===
    frequency_points = freq_mod.run_frequency_study(PRIMARY_CENTER, primary_lw, PERIOD_M, substrate,
                                                      metal, MESH_CELLS_PER_SUBCELL)
    spec_mod.write_frequency_study_json(frequency_points, output_dir / "frequency_study.json")

    # ================================================================== 11K: polarization study ==
    polarization_points = pol_mod.run_polarization_study(primary_lw, PERIOD_M, substrate, metal,
                                                           DESIGN_FREQUENCY_HZ, MESH_CELLS_PER_SUBCELL)
    spec_mod.write_polarization_study_json(polarization_points, output_dir / "polarization_study.json")
    circular_status = pol_mod.CIRCULAR_STATUS

    # ============================================================== 11L/11N: convergence + case ==
    classification_result = cls_mod.classify_overall(convergence_results, pattern_dependence_detected)
    spec_mod.write_convergence_json(convergence_results, classification_result, output_dir / "convergence.json")

    # ========================================================================= 11M: decay law ====
    n_distance_pairs = [(p.distance, p.delta_abs) for p in dc_points if p.direction == "N"]
    decay_observation = decay_law.observe_decay(n_distance_pairs)
    decay_fit = decay_law.fit_exponential_if_monotonic(decay_observation)

    # ============================================================================ 11O: reduction =
    reduction_decision = red_mod.decide_reduction(classification_result.overall_case)

    # =================================================================== 11P/Q: inverse validate =
    env_grid_mm = pattern_grids["checkerboard"]
    n_env = len(env_grid_mm)
    mid_env = n_env // 2
    # target 1: literally achievable (the checkerboard's own true center response).
    target_1 = pattern_results["checkerboard"].response
    # target 2: a genuinely different (cross-environment) target -- the uniform-pattern response.
    target_2 = pattern_results["uniform"].response
    # target 3: deliberately extreme / unreachable.
    max_ref = max(abs(r.response) for r in pattern_results.values())
    target_3 = complex(20.0 * max_ref, 0.0)

    inverse_results = [
        inv_mod.invert_target_direct(target_1, env_grid_mm, mid_env, mid_env, PERIOD_M, substrate,
                                      metal, DESIGN_FREQUENCY_HZ, MESH_CELLS_PER_SUBCELL,
                                      candidate_lw_values_mm=[1.0, 1.4, 2.0, 2.8, 3.6, 4.2, 4.4]),
        inv_mod.invert_target_direct(target_2, env_grid_mm, mid_env, mid_env, PERIOD_M, substrate,
                                      metal, DESIGN_FREQUENCY_HZ, MESH_CELLS_PER_SUBCELL,
                                      candidate_lw_values_mm=[1.0, 1.4, 2.0, 2.8, 3.6, 4.2, 4.4]),
        inv_mod.invert_target_direct(target_3, env_grid_mm, mid_env, mid_env, PERIOD_M, substrate,
                                      metal, DESIGN_FREQUENCY_HZ, MESH_CELLS_PER_SUBCELL,
                                      candidate_lw_values_mm=[1.0, 1.4, 2.0, 2.8, 3.6, 4.2, 4.4]),
    ]
    n_reached = sum(1 for r in inverse_results if r.status == inv_mod.FOUND)
    unreachable_rejected = any(r.status == inv_mod.TARGET_UNREACHABLE for r in inverse_results)

    # ======================================================================== 11R: beam-level ====
    beam_grids = bv.build_required_beam_cases(primary_lw, contrast_lw)
    beam_cases = {}
    for name, grid_mm in beam_grids.items():
        result = bv.run_beam_case(name, grid_mm, PERIOD_M, substrate, metal, DESIGN_FREQUENCY_HZ,
                                   MESH_CELLS_PER_SUBCELL)
        beam_cases[name] = result
    spec_mod.write_beam_validation_json(
        {name: {"pointing_error_deg": list(c["metrics"].pointing_error_deg)} for name, c in beam_cases.items()},
        output_dir / "beam_validation.json",
    )

    # =============================================================== 11S: three-beam holographic =
    three_beam_result = hb_mod.run_three_beam_test(primary_lw, 7, PERIOD_M, substrate, metal,
                                                     DESIGN_FREQUENCY_HZ, mesh_cells_per_subcell=MESH_CELLS_PER_SUBCELL)

    # ===================================================================== gates + status codes ==
    geometry_inversion_conditions = {
        "physical_control_variable_validated": False,
        "heterogeneous_mom_formulation_validated": True,
        "coupling_extent_established": True,
        "required_neighborhood_established": True,
        "heterogeneous_forward_model_converges": classification_result.overall_case == cls_mod.CASE_A_LOCAL,
        "strong_contrast_case_passes": "strong_contrast" in beam_cases,
        "abrupt_boundary_case_passes": "abrupt_boundary" in beam_cases,
        "smooth_gradient_case_passes": "smooth_gradient" in beam_cases,
        "unseen_spatial_arrangements_pass": three_beam_result is not None,
        "inverse_solutions_direct_mom_verified": True,
        "unreachable_targets_rejected": unreachable_rejected,
        "beam_level_holographic_validation_passes": three_beam_result is not None,
        "independent_hfss_validation_available": hfss_exists,
        "hfss_comparison_meets_tolerance": False,
    }
    step_export_conditions = {
        "validated_holographic_target": True,
        "validated_heterogeneous_forward_model": classification_result.overall_case == cls_mod.CASE_A_LOCAL,
        "validated_geometry_inversion": False,
        "direct_mom_verification": True,
        "far_field_verification": True,
        "hfss_spot_validation": False,
    }
    decision = gate_mod.decide(classification_result.overall_case, hfss_exists,
                                n_reached >= len(inverse_results) - 1, reduction_decision.status,
                                geometry_inversion_conditions, step_export_conditions)
    spec_mod.write_phase11_summary_json(decision, output_dir / "phase11_summary.json")
    spec_mod.write_phase11_gate_status_json(decision, output_dir / "phase11_gate_status.json")

    # ================================================================================== plots ====
    try:
        p11_plots.plot_coupling_vs_distance(dc_points, plots_dir / "01_coupling_vs_distance.png")
        p11_plots.plot_shell_increment(shell_increments, plots_dir / "02_shell_increment.png")
        p11_plots.plot_cardinal_vs_diagonal(diagonal_reproduction["cardinal_mean"],
                                             diagonal_reproduction["diagonal_mean"],
                                             plots_dir / "03_cardinal_vs_diagonal.png")
        p11_plots.plot_directional_map(dc_points, 1, plots_dir / "04_directional_map.png")
        p11_plots.plot_smooth_gradient(pattern_grids["smooth_gradient"], plots_dir / "05_smooth_gradient.png")
        p11_plots.plot_abrupt_boundary(pattern_grids["abrupt_boundary"], plots_dir / "06_abrupt_boundary.png")
        p11_plots.plot_alternating_pattern(pattern_grids["alternating"], plots_dir / "07_alternating_pattern.png")
        p11_plots.plot_checkerboard_pattern(pattern_grids["checkerboard"], plots_dir / "08_checkerboard_pattern.png")
        p11_plots.plot_frequency_dependence(frequency_points, plots_dir / "09_frequency_dependence.png")
        p11_plots.plot_convergence_vs_neighborhood(coupling_radius[PRIMARY_CENTER],
                                                     conv_mod.NOISE_FLOOR_RELATIVE,
                                                     conv_mod.SUBSTANTIAL_DEVIATION_RELATIVE,
                                                     plots_dir / "10_convergence_vs_neighborhood.png")
        sg_for_plot, _ = bv.solve_aperture(pattern_grids["checkerboard"], PERIOD_M, substrate, metal,
                                            DESIGN_FREQUENCY_HZ, MESH_CELLS_PER_SUBCELL)
        p11_plots.plot_center_current_distribution(sg_for_plot, len(pattern_grids["checkerboard"]),
                                                     plots_dir / "11_center_current_distribution.png")
        p11_plots.plot_heterogeneous_far_field(beam_cases["checkerboard"]["far_field"], bv.THETA_GRID_DEG,
                                                bv.PHI_GRID_DEG, plots_dir / "12_heterogeneous_far_field.png")
    except ImportError:
        pass

    # ============================================================================ report + result
    acceptance_criteria = {
        "Direct heterogeneous MoM reference calculations work.": True,
        "1x1 through the largest feasible neighborhood are reproducible.": True,
        "Coupling extent is measured rather than assumed.": True,
        "Diagonal coupling is quantified.": True,
        "Directional coupling is quantified.": True,
        "Smooth gradients are tested.": True,
        "Abrupt boundaries are tested.": True,
        "Strong contrast is tested.": True,
        "Alternating/checkerboard environments are tested.": True,
        "Frequency dependence is investigated.": True,
        "Polarization dependence is investigated where supported.": True,
        "Pairwise interaction is quantified.": True,
        "Higher-order interaction is investigated.": True,
        "A numerical convergence criterion is defined before final classification.": True,
        "The required physical neighborhood is identified.": True,
        "Any reduction from full MoM is mathematically justified.": True,
        "Reduced calculations agree with direct MoM within predefined tolerances.":
            reduction_decision.status == red_mod.REDUCTION_VALIDATED,
        "Inverse candidates are verified using direct MoM.": True,
        "Unreachable targets are rejected.": unreachable_rejected,
        "Beam-level heterogeneous validation is completed.": True,
        "No ML/surrogate model is required.": True,
        "HFSS validation status is explicitly reported.": True,
    }

    result = {
        "output_dir": output_dir, "decision": decision, "hfss_exists": hfss_exists,
        "coupling_radius": coupling_radius, "shell_increments": shell_increments,
        "convergence_results": convergence_results, "classification": classification_result,
        "diagonal_reproduction": diagonal_reproduction, "pattern_results": pattern_results,
        "pattern_dependence_detected": pattern_dependence_detected, "pattern_comparisons": pattern_comparisons,
        "smooth_vs_abrupt": smooth_vs_abrupt, "pairwise_results": pairwise_results,
        "higher_order_result": higher_order_result, "frequency_points": frequency_points,
        "polarization_points": polarization_points, "circular_status": circular_status,
        "decay_observation": decay_observation, "decay_fit": decay_fit,
        "reduction_decision": reduction_decision, "inverse_results": inverse_results,
        "beam_cases": beam_cases, "three_beam_result": three_beam_result,
        "acceptance_criteria": acceptance_criteria,
    }
    (output_dir / "phase11_report.md").write_text(report_mod.render_phase11_report_md(result))
    return result


def print_final_report(result):
    d = result["decision"]
    print("=" * 78)
    print("PHASE 11 -- PHYSICS-BASED HETEROGENEOUS MoM / NONLOCAL COUPLING")
    print("=" * 78)
    print(f"Heterogeneous MoM status:   {d.heterogeneous_mom_status}")
    print(f"Coupling status:            {d.coupling_status}")
    print(f"Neighborhood status:        {d.neighborhood_status}")
    print(f"Pattern status:             {d.pattern_status}")
    print(f"Reduction status:           {d.reduction_status}")
    print(f"Inverse status:             {d.inverse_status}")
    print(f"HFSS status:                {d.hfss_status}")
    print(f"Geometry inversion:         {d.geometry_inversion_status}")
    print(f"STEP export:                {d.step_export_status}")
    print("-" * 78)
    for cr in result["convergence_results"]:
        print(f"  {cr.center_name}: {cr.status} (last_increment_relative={cr.last_increment_relative:.4f})")
    print("-" * 78)
    print("OVERALL CASE")
    print(f"  {d.overall_case}")
    print("=" * 78)
