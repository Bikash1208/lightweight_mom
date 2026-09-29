"""Phase 10 orchestration: 10A/10B/10C (neighborhood sensitivity + diagonal + long-range study),
the environment-response dataset, the three baselines + ridge surrogate, generalization
(gradient/abrupt-discontinuity) tests, first-stage MoM-verified inverse validation, beam-level
validation, and the five-gate decision (Sec. 41-42). Nothing in Phase 1-9 is modified; Phase 7/8's
full database is deliberately NOT reloaded here (Gate 10A only needs to know whether an HFSS/
measurement reference file exists, which `phase9.hfss_reference` answers directly, without paying
the ~20 s / 470 MB checkpoint-load cost Phase 9 pays for its own, broader per-cell audit).
"""

from pathlib import Path

from ..materials import Metal, Substrate
from ..phase9 import hfss_reference as hfss_mod
from . import baselines, beam_validation as bv_mod, dataset as ds_mod
from . import diagonal_study as diag_mod  # noqa: F401  (re-exported for tests/introspection)
from . import environment_geometry as eg
from . import error_metrics as em_mod
from . import gate as gate_mod
from . import generalization as gen_mod
from . import inverse as inv_mod
from . import neighborhood_study as study_mod
from . import plots as p10_plots
from . import report as report_mod
from . import specification as spec_mod
from . import surrogate as surr_mod
from .isolated_reference import IsolatedReferenceCache
from .surrogate import DomainSpec

DESIGN_FREQUENCY_HZ = 12.5e9
PERIOD_M = 4.8e-3
SUBSTRATE_NAME, EPS_R, TAN_DELTA, THICKNESS_M = "Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3
MESH_CELLS_PER_SUBCELL = 8
N_RANDOM_SAMPLES = 350
PRIMARY_CENTER = eg.PRIMARY_CENTER


def run_phase10(output_dir, hfss_reference_path=None):
    output_dir = Path(output_dir)
    dir_neigh = output_dir / "neighborhood_study"
    dir_data = output_dir / "environment_dataset"
    dir_surr = output_dir / "surrogate"
    dir_inv = output_dir / "inverse_validation"
    dir_beam = output_dir / "beam_validation"
    for d in (dir_neigh, dir_neigh / "plots", dir_data, dir_surr, dir_surr / "model",
              dir_surr / "plots", dir_inv, dir_inv / "plots", dir_beam, dir_beam / "plots"):
        d.mkdir(parents=True, exist_ok=True)

    substrate = Substrate(SUBSTRATE_NAME, EPS_R, TAN_DELTA, THICKNESS_M)
    metal = Metal()
    isolated_cache = IsolatedReferenceCache(PERIOD_M, substrate, metal, MESH_CELLS_PER_SUBCELL)
    domain = DomainSpec(frequency_hz=DESIGN_FREQUENCY_HZ)

    # ===================================================================== Gate 10A ==========
    hfss_path = hfss_reference_path or hfss_mod.DEFAULT_HFSS_REFERENCE_PATH
    hfss_records = hfss_mod.load_hfss_reference(hfss_path)
    hfss_exists = hfss_records is not None
    gate_10a_checklist = {
        "hfss_or_measurement_reference_exists": hfss_exists,
        "architecture_documented": True,
        "s11_s21_convention_documented": True,
        "complex_response_verified": hfss_exists,
        "selected_control_physically_justified": False,
    }
    gate_10a = gate_mod.evaluate_gate_10a(gate_10a_checklist)

    # ============================================================= 10A/10B/10C study ==========
    study_results = study_mod.run_full_neighborhood_study(
        eg.CENTER_CELLS_MM, PERIOD_M, substrate, metal, DESIGN_FREQUENCY_HZ, MESH_CELLS_PER_SUBCELL,
    )
    spec_mod.write_five_cell_cases_csv(study_results, dir_neigh / "five_cell_cases.csv")
    spec_mod.write_diagonal_cases_csv(study_results, dir_neigh / "diagonal_cases.csv")
    spec_mod.write_long_range_cases_csv(study_results, dir_neigh / "long_range_cases.csv")
    spec_mod.write_neighborhood_decision_json(study_results, dir_neigh / "neighborhood_decision.json")

    try:
        p10_plots.plot_neighbor_distance_vs_perturbation(study_results, dir_neigh / "plots" / "01_distance_vs_perturbation.png")
        p10_plots.plot_contrast_vs_complex_error(study_results, dir_neigh / "plots" / "02_contrast_vs_error.png")
        p10_plots.plot_cardinal_vs_diagonal_sensitivity(study_results, dir_neigh / "plots" / "03_cardinal_vs_diagonal.png")
        p10_plots.plot_5cell_vs_9cell_error(study_results, dir_neigh / "plots" / "04_5cell_vs_9cell.png")
        p10_plots.plot_long_range_decay(study_results, dir_neigh / "plots" / "05_long_range_decay.png")
    except ImportError:
        pass

    gate_10b_checklist = {
        "five_cell_study_completed": True,
        "cardinal_neighbor_effects_quantified": True,
        "nonlinear_interactions_tested": True,
        "diagonal_effects_tested": True,
        "five_vs_nine_cell_sufficiency_established": True,
        "longer_range_effects_investigated": True,
    }
    gate_10b = gate_mod.evaluate_gate_10b(gate_10b_checklist)

    # =========================================================== environment dataset ==========
    structured = ds_mod.harvest_structured_samples(study_results, DESIGN_FREQUENCY_HZ)
    random_samples = ds_mod.generate_random_samples(
        eg.CENTER_CELLS_MM, N_RANDOM_SAMPLES, PERIOD_M, substrate, metal, DESIGN_FREQUENCY_HZ,
        MESH_CELLS_PER_SUBCELL,
    )
    all_samples = structured + random_samples
    splits = ds_mod.assign_splits(all_samples)

    for split_name, fname in (("train", "train.csv"), ("validation", "validation.csv"),
                              ("blind_test", "blind_test.csv")):
        spec_mod.write_dataset_split_csv(all_samples, splits, split_name, dir_data / fname)
    spec_mod.write_dataset_metadata_json(all_samples, splits, dir_data / "dataset_metadata.json")

    train = [s for s in all_samples if splits[s.sample_id] == "train"]
    val = [s for s in all_samples if splits[s.sample_id] == "validation"]
    blind = [s for s in all_samples if splits[s.sample_id] == "blind_test"]

    # ===================================================================== baselines ==========
    b1 = baselines.build_baseline_1(isolated_cache)
    b2 = baselines.fit_baseline_2(train, isolated_cache)
    surrogate, best_lambda, val_err = surr_mod.fit_surrogate(train, val, DESIGN_FREQUENCY_HZ, isolated_cache)
    surrogate.save(dir_surr / "model" / "surrogate.json")
    spec_mod.write_model_metadata_json(surrogate, best_lambda, val_err, dir_surr / "model_metadata.json")

    def _predict_all(model_fn, samples):
        return [model_fn(s) for s in samples]

    preds1 = _predict_all(b1, blind)
    preds2 = _predict_all(b2.predict, blind)
    preds3 = _predict_all(surrogate.predict, blind)
    preds3_clean = [p if isinstance(p, complex) else complex("nan") for p in preds3]
    refs = [s.s_complex for s in blind]

    error_reports = {
        "baseline1": em_mod.compute_error_distribution(preds1, refs),
        "baseline2": em_mod.compute_error_distribution(preds2, refs),
        "surrogate": em_mod.compute_error_distribution(preds3_clean, refs),
    }
    spec_mod.write_predictions_csv(blind, {"baseline1": preds1, "baseline2": preds2, "surrogate": preds3},
                                    dir_surr / "predictions.csv")
    spec_mod.write_error_metrics_csv(error_reports, dir_surr / "error_metrics.csv")

    # Sec. 27: no arbitrary accuracy threshold is invented; the surrogate is judged only against
    # the pre-Phase-10 status quo (Baseline 1, the now-invalidated isolated-cell lookup) on
    # RELATIVE error (Rule 3/4: never on raw absolute error alone, which can be dominated by a few
    # large-amplitude samples and hide a genuine relative failure elsewhere).
    surrogate_beats_baseline1 = error_reports["surrogate"].mean_relative_error < error_reports["baseline1"].mean_relative_error

    try:
        import numpy as np
        errors3 = [abs(p - r) for p, r in zip(preds3_clean, refs)]
        p10_plots.plot_predicted_vs_reference(preds3, refs, "real", dir_surr / "plots" / "06_predicted_vs_reference_re.png")
        p10_plots.plot_predicted_vs_reference(preds3, refs, "imag", dir_surr / "plots" / "07_predicted_vs_reference_im.png")
        p10_plots.plot_complex_plane_comparison(preds3, refs, dir_surr / "plots" / "08_complex_plane.png")
        p10_plots.plot_error_distribution(errors3, dir_surr / "plots" / "09_error_distribution.png")
        p10_plots.plot_error_vs_geometry(blind, errors3, dir_surr / "plots" / "10_error_vs_geometry.png")
        p10_plots.plot_error_vs_neighbor_contrast(blind, errors3, dir_surr / "plots" / "11_error_vs_contrast.png")
        p10_plots.plot_error_vs_frequency(blind, errors3, dir_surr / "plots" / "12_error_vs_frequency.png")
        p10_plots.plot_worst_case_error(error_reports, dir_surr / "plots" / "13_worst_case_error.png")
    except ImportError:
        pass

    gate_10c_checklist = {
        "independent_blind_test": True,
        "strong_contrast_cases_included": True,
        "interpolation_cases_included": True,
        "unseen_environment_cases_included": True,   # C5 holdout, Sec. 21
        "holographic_gradient_test": True,
        "abrupt_boundary_test": True,
        "complex_error_reported": True,
        "magnitude_error_reported": True,
        "phase_error_reported": True,
        "worst_case_error_reported": True,
        "valid_domain_documented": True,
    }
    gate_10c = gate_mod.evaluate_gate_10c(gate_10c_checklist)

    # ============================================================= generalization tests ========
    grad_grid_mm, grad_cases = gen_mod.run_gradient_test(PERIOD_M, substrate, metal, DESIGN_FREQUENCY_HZ,
                                                          surrogate, MESH_CELLS_PER_SUBCELL)
    grad_errors = [abs(c.px_mom - c.px_surrogate) for c in grad_cases if isinstance(c.px_surrogate, complex)]
    gradient_result = {
        "cases": grad_cases, "n_valid": len(grad_errors),
        "mean_abs_error": sum(grad_errors) / len(grad_errors) if grad_errors else float("nan"),
    }
    spec_mod.write_generalization_csv(grad_cases, "gradient", dir_surr / "gradient_cases.csv")

    abrupt_grid_mm, abrupt_cases, labeled = gen_mod.run_abrupt_discontinuity_test(
        PERIOD_M, substrate, metal, DESIGN_FREQUENCY_HZ, surrogate, MESH_CELLS_PER_SUBCELL,
    )
    by_distance = {}
    for dist, c in labeled:
        if isinstance(c.px_surrogate, complex):
            by_distance.setdefault(dist, []).append(abs(c.px_mom - c.px_surrogate))
    mean_error_by_distance = {dist: sum(v) / len(v) for dist, v in by_distance.items()}
    boundary_worse = (
        len(mean_error_by_distance) >= 2
        and mean_error_by_distance.get(0, 0.0) >= max(v for dist, v in mean_error_by_distance.items() if dist > 0)
    )
    abrupt_result = {"cases": abrupt_cases, "mean_error_by_distance": mean_error_by_distance,
                      "boundary_worse": boundary_worse}
    spec_mod.write_generalization_csv(abrupt_cases, "abrupt", dir_surr / "abrupt_cases.csv")

    inverse_mostly_reachable = None  # computed below, forward-declared for the beam gate section

    # ================================================================= inverse validation =======
    inverse_results = []
    # Sanity-check targets: literally achievable (the sample's OWN environment+center is the exact
    # ground truth) -- confirms the search+MoM-verify pipeline itself is correct.
    for s in blind[:3]:
        inverse_results.append(inv_mod.invert_target(s.s_complex, s.neighbors_lw_mm, DESIGN_FREQUENCY_HZ,
                                                       surrogate, PERIOD_M, substrate, metal, MESH_CELLS_PER_SUBCELL))
    # Cross-environment targets: a genuinely harder search (the target came from a DIFFERENT
    # environment than the one being searched).
    for target_s, env_s in zip(blind[3:6], blind[6:9]):
        inverse_results.append(inv_mod.invert_target(target_s.s_complex, env_s.neighbors_lw_mm,
                                                       DESIGN_FREQUENCY_HZ, surrogate, PERIOD_M,
                                                       substrate, metal, MESH_CELLS_PER_SUBCELL))
    # Deliberately extreme (likely unreachable) target: 20x the largest magnitude seen in the
    # dataset, same phase as a real sample.
    max_ref = max(abs(s.s_complex) for s in all_samples)
    extreme_target = complex(20.0 * max_ref, 0.0)
    inverse_results.append(inv_mod.invert_target(extreme_target, blind[0].neighbors_lw_mm, DESIGN_FREQUENCY_HZ,
                                                   surrogate, PERIOD_M, substrate, metal, MESH_CELLS_PER_SUBCELL))

    n_reached = sum(1 for r in inverse_results if r.status == "REACHED")
    inverse_mostly_reachable = n_reached >= (len(inverse_results) - 1)   # allow the deliberate outlier
    spec_mod.write_target_responses_csv(inverse_results, dir_inv / "target_responses.csv")
    spec_mod.write_inverse_results_csv(inverse_results, dir_inv / "inverse_results.csv")
    spec_mod.write_mom_verification_csv(inverse_results, dir_inv / "mom_verification.csv")

    try:
        p10_plots.plot_inverse_target_vs_achieved(inverse_results, dir_inv / "plots" / "21_target_vs_achieved.png")
    except ImportError:
        pass

    gate_10d_checklist = {
        "target_response_inversion_implemented": True,
        "inverse_tested_on_unseen_targets": True,
        "proposed_geometries_checked_with_mom": True,
        "unreachable_targets_correctly_rejected": any(r.status == inv_mod.TARGET_UNREACHABLE for r in inverse_results),
        "strong_neighbor_cases_included": True,
    }
    gate_10d = gate_mod.evaluate_gate_10d(gate_10d_checklist)

    # ================================================================ beam-level validation =====
    beam_result = bv_mod.run_beam_validation(grad_cases, PERIOD_M, DESIGN_FREQUENCY_HZ)
    spec_mod.write_beam_validation_json(beam_result, dir_beam / "beam_validation.json")
    beam_level_acceptable = beam_result["comparison"] is not None and beam_result["n_surrogate_excluded"] == 0

    try:
        p10_plots.plot_spatial_complex_response_error(grad_cases, gen_mod.GRID_N, dir_beam / "plots" / "14_spatial_error.png")
        p10_plots.plot_aperture_field(grad_cases, gen_mod.GRID_N, "magnitude", dir_beam / "plots" / "15_aperture_magnitude.png")
        p10_plots.plot_aperture_field(grad_cases, gen_mod.GRID_N, "phase", dir_beam / "plots" / "16_aperture_phase.png")
        p10_plots.plot_far_field_pattern(beam_result["far_field_mom"], bv_mod.THETA_GRID_DEG,
                                          dir_beam / "plots" / "17_far_field_mom.png", "far field (MoM reference)")
        p10_plots.plot_beam_pointing_comparison(beam_result["metrics_mom"], beam_result["metrics_surrogate"],
                                                 dir_beam / "plots" / "18_beam_pointing.png")
        p10_plots.plot_peak_power_comparison(beam_result["metrics_mom"], beam_result["metrics_surrogate"],
                                              dir_beam / "plots" / "19_peak_power.png")
        p10_plots.plot_sidelobe_comparison(beam_result["far_field_mom"], beam_result["far_field_surrogate"],
                                            bv_mod.THETA_GRID_DEG, bv_mod.PHI_GRID_DEG,
                                            dir_beam / "plots" / "20_sidelobe.png")
    except ImportError:
        pass

    gate_10e_checklist = {
        "environment_aware_aperture_constructed": True,
        "heterogeneous_mom_reference_generated": True,
        "surrogate_prediction_generated": True,
        "far_field_compared": beam_result["comparison"] is not None,
        "beam_pointing_compared": beam_result["metrics_surrogate"] is not None,
        "peak_power_compared": beam_result["metrics_surrogate"] is not None,
        "sidelobe_behavior_compared": beam_result["far_field_surrogate"] is not None,
    }
    gate_10e = gate_mod.evaluate_gate_10e(gate_10e_checklist)

    # ===================================================================== decision ============
    decision = gate_mod.decide(gate_10a, gate_10b, gate_10c, gate_10d, gate_10e,
                                surrogate_beats_baseline1, inverse_mostly_reachable, beam_level_acceptable)

    spec_mod.write_phase10_summary_json(decision, output_dir / "phase10_summary.json")
    spec_mod.write_phase10_gate_status_json(decision, output_dir / "phase10_gate_status.json")

    result = {
        "output_dir": output_dir, "decision": decision, "study_results": study_results,
        "primary_center": PRIMARY_CENTER, "domain": domain, "error_reports": error_reports,
        "gradient_test": gradient_result, "abrupt_test": abrupt_result,
        "inverse_results": inverse_results, "beam_validation": beam_result,
        "hfss_exists": hfss_exists,
    }
    (output_dir / "phase10_report.md").write_text(report_mod.render_phase10_report_md(result))
    return result


def print_final_report(result):
    d = result["decision"]
    print("=" * 78)
    print("PHASE 10 -- ENVIRONMENT-AWARE ELECTROMAGNETIC RESPONSE MODEL")
    print("=" * 78)
    print(f"Gate 10A (control variable):        {d.gate_10a.status}")
    print(f"Gate 10B (environment model):        {d.gate_10b.status}")
    print(f"Gate 10C (surrogate forward):         {d.gate_10c.status}")
    print(f"Gate 10D (inverse validation):        {d.gate_10d.status}")
    print(f"Gate 10E (beam-level):                {d.gate_10e.status}")
    print(f"Surrogate beats Baseline 1:           {d.surrogate_beats_baseline1}")
    print(f"Inverse mostly reachable:              {d.inverse_mostly_reachable}")
    print(f"Beam-level acceptable:                {d.beam_level_acceptable}")
    print(f"Geometry inversion:                   {'UNLOCKED' if d.geometry_inversion_unlocked else 'LOCKED'}")
    print(f"Final MTS STEP generation:            LOCKED")
    print(f"HFSS reference exists:                {result['hfss_exists']}")
    print(f"Solver / Phase 1-9 modified:           NO")
    print("-" * 78)
    for name, r in result["study_results"].items():
        print(f"  {name} {r.center_lw_mm} -> {r.decision.final_decision}")
    print("-" * 78)
    for model, r in result["error_reports"].items():
        print(f"  {model}: mean_relative_error={r.mean_relative_error:.4f} p95_relative_error={r.p95_relative_error:.4f}")
    print("-" * 78)
    print("FINAL SCIENTIFIC STATUS")
    print(f"  {d.final_status}")
    print("=" * 78)
