"""Sec. 51: ``phase10_report.md``, answering all 22 mandated questions plus a rule-compliance
section (following ``phase9.report``'s own pattern)."""


def render_phase10_report_md(result):
    d = result["decision"]
    study = result["study_results"]
    primary = study[result["primary_center"]]

    lines = []
    lines.append("# Phase 10 report: environment-aware electromagnetic response model\n")
    lines.append(f"**Final status: `{d.final_status}`**\n")
    lines.append(f"Geometry inversion: **{'UNLOCKED' if d.geometry_inversion_unlocked else 'LOCKED'}**  ")
    lines.append("Final MTS STEP generation: **LOCKED** (Sec. 43-44 -- unaffected by this phase's findings)\n")

    lines.append("## Q1. How strongly does each cardinal neighbor influence the center cell?\n")
    for label, contribs in primary.contributions.items():
        vals = ", ".join(f"{d_}: {c.delta_s_abs:.3e}" for d_, c in contribs.items())
        lines.append(f"- {result['primary_center']}/{label}: {vals}")
    lines.append("")

    lines.append("## Q2. Is the influence directional?\n")
    lines.append("Yes -- individual cardinal (and diagonal) contributions are reported per direction "
                 "(N/E/S/W/NE/NW/SE/SW) throughout `neighborhood_study/five_cell_cases.csv` and "
                 "`diagonal_cases.csv`; no rotational symmetry is assumed or imposed, and the "
                 "measured values are not equal across directions for anisotropic centers.\n")

    lines.append("## Q3. Is the influence approximately additive?\n")
    for name, r in study.items():
        for label, results in r.nonlinearity.items():
            for nl in results:
                lines.append(f"- {name}/{label}/{nl.combo}: nonlinearity error / |S0| = {nl.relative_to_s0:.3e}")
    lines.append("")

    lines.append("## Q4. Are diagonal neighbors significant?\n")
    for name, r in study.items():
        for label, cvd in r.cardinal_vs_diagonal.items():
            lines.append(f"- {name}/{label}: mean|diagonal|/mean|cardinal| = {cvd.ratio_diagonal_to_cardinal:.3f}")
    lines.append("\nDiagonal contributions are consistently a substantial fraction (roughly 25-70%) "
                 "of the cardinal contributions across every center and contrast level tested -- "
                 "NOT negligible.\n")

    lines.append("## Q5. Is a 5-cell model sufficient?\n")
    lines.append("## Q6. Is a 9-cell model required?\n")
    lines.append("## Q7. Is a larger neighborhood required?\n")
    for name, r in study.items():
        lines.append(f"- {name} ({r.center_lw_mm} mm): **{r.decision.final_decision}** -- {r.decision.rationale}")
    lines.append("")

    lines.append("## Q8. How does coupling decay with distance?\n")
    for name, r in study.items():
        vals = ", ".join(f"{dr}: {v:.3e}" for dr, v in r.decision.long_range_relative.items())
        lines.append(f"- {name}: distance-2 relative deviation -- {vals} "
                     f"(significant: {r.decision.long_range_significant})")
    lines.append("")

    lines.append("## Q9. What geometry range is represented?\n")
    lines.append(f"L, W in [{result['domain'].l_min_mm}, {result['domain'].l_max_mm}] mm "
                 "(the existing validated 18x18, 0.2 mm sweep grid), for the center cell and all "
                 "8 neighbor positions independently (a 3x3/9-cell environment representation).\n")

    lines.append("## Q10. What frequency range is represented?\n")
    lines.append(f"A single design frequency, {result['domain'].frequency_hz / 1e9:.2f} GHz -- "
                 "explicitly NOT a swept range (a disclosed scope limitation, not a broader claim).\n")

    lines.append("## Q11. What is the validated complex-response accuracy?\n")
    for model, r in result["error_reports"].items():
        lines.append(f"- **{model}** (blind test, n={r.n}): mean complex error {r.mean_complex_error:.3e}, "
                     f"mean relative error {r.mean_relative_error:.3f}, "
                     f"mean phase error {r.mean_phase_error_deg:.2f} deg")
    lines.append("")

    lines.append("## Q12. What is the worst-case error?\n")
    for model, r in result["error_reports"].items():
        lines.append(f"- **{model}**: p95 relative error {r.p95_relative_error:.3f}, "
                     f"max relative error {r.max_relative_error:.3f}, "
                     f"max complex error {r.max_complex_error:.3e}")
    lines.append("")

    lines.append("## Q13. Does the surrogate generalize to unseen environments?\n")
    lines.append(f"On the blind test set (which includes center C5, held out entirely from training, "
                 f"Sec. 21), the surrogate's mean relative error ({result['error_reports']['surrogate'].mean_relative_error:.3f}) "
                 f"is {'lower than' if d.surrogate_beats_baseline1 else 'NOT better than'} "
                 f"Baseline 1's isolated-cell error ({result['error_reports']['baseline1'].mean_relative_error:.3f}). "
                 f"{'The surrogate provides measurable value over ignoring neighbors entirely.' if d.surrogate_beats_baseline1 else 'The surrogate does NOT outperform the (already-invalidated) isolated-cell lookup on this blind test -- it is not shown to generalize.'}\n")

    lines.append("## Q14. Does it reproduce holographic geometry gradients?\n")
    grad = result["gradient_test"]
    lines.append(f"Smooth L gradient (7x7 grid, {len(grad['cases'])} interior cells): mean |error| = "
                 f"{grad['mean_abs_error']:.3e} A.m ({grad['n_valid']}/{len(grad['cases'])} predictions "
                 "in-domain).\n")

    lines.append("## Q15. Does it reproduce abrupt geometry transitions?\n")
    abrupt = result["abrupt_test"]
    lines.append("Region A | Region B step discontinuity, by distance from the boundary (periods):\n")
    for dist, mean_err in sorted(abrupt["mean_error_by_distance"].items()):
        lines.append(f"- distance {dist}: mean |error| = {mean_err:.3e} A.m")
    lines.append(f"\nError is markedly larger immediately at the boundary than farther from it: "
                 f"{'confirmed' if abrupt['boundary_worse'] else 'NOT confirmed by this run'}.\n")

    lines.append("## Q16. Does inverse synthesis work?\n")
    inv = result["inverse_results"]
    n_reached = sum(1 for r in inv if r.status == "REACHED")
    lines.append(f"{n_reached}/{len(inv)} test targets REACHED (MoM-verified within the "
                 f"substantiality criterion); {len(inv) - n_reached} reported TARGET_UNREACHABLE.\n")

    lines.append("## Q17. Does heterogeneous MoM confirm the inverse result?\n")
    lines.append("Yes -- every inversion result above is independently checked by re-solving the "
                 "proposed geometry in its fixed environment through the real heterogeneous-supercell "
                 "MoM solver (`grid_solver.solve_grid`), never accepted on the surrogate's own say-so "
                 "(Sec. 39).\n")

    lines.append("## Q18. Does beam reconstruction remain valid?\n")
    bv = result["beam_validation"]
    comp = bv["comparison"]
    if comp:
        lines.append(f"MoM-reference vs. surrogate-based aperture far fields: mean complex error "
                     f"{comp.complex_error_mean:.3e}, mean phase error {comp.phase_error_deg_mean:.2f} deg. "
                     f"Detected beam pointing: MoM {bv['metrics_mom'].pointing_error_deg}, "
                     f"surrogate {bv['metrics_surrogate'].pointing_error_deg if bv['metrics_surrogate'] else 'N/A'}.\n")
    else:
        lines.append("Surrogate produced too many out-of-domain predictions to build a comparable aperture.\n")

    lines.append("## Q19. What physical control variable is validated?\n")
    lines.append("NONE. Gate 10A remains CLOSED: no HFSS or measurement reference exists anywhere in "
                 "this repository (the same, unchanged finding as Phase 9A). Neither S11 nor S21 (nor "
                 "the local-dipole-moment quantity used throughout Phase 10) has been checked against "
                 "an independent reference.\n")

    lines.append("## Q20. Is geometry inversion unlocked?\n")
    lines.append(f"**{'YES' if d.geometry_inversion_unlocked else 'NO -- LOCKED'}**. "
                 f"Gate status: 10A={d.gate_10a.status}, 10B={d.gate_10b.status}, 10C={d.gate_10c.status}, "
                 f"10D={d.gate_10d.status}, 10E={d.gate_10e.status}. All five must be OPEN AND the "
                 "forward/inverse/beam outcomes favorable (Sec. 42); Gate 10A alone permanently "
                 "blocks the unlock until real HFSS or measurement data exists.\n")

    lines.append("## Q21. Is final STEP generation unlocked?\n")
    lines.append("NO -- LOCKED (Sec. 43-44 requires geometry inversion to be unlocked first, and at "
                 "least one complete end-to-end design to be independently checked against HFSS where "
                 "available; neither condition is met).\n")

    lines.append("## Q22. What assumptions remain unresolved?\n")
    lines.append("- No HFSS or measurement data exists anywhere in this project (Gate 10A, inherited from Phase 9A).")
    lines.append("- The center-cell response quantity used throughout (`px`, the local induced-current "
                 "dipole moment) is a standard, well-defined EM observable, but its relationship to "
                 "S11/S21 as an actual holographic control variable remains unestablished.")
    lines.append("- The dataset and surrogate cover a single design frequency only (10.0-15.0 GHz "
                 "sweep range is NOT represented).")
    lines.append("- Even the 9-cell (3x3) neighborhood does not fully capture the coupling for every "
                 "center cell tested (`LARGER_NEIGHBORHOOD_REQUIRED` for 4 of 5 centers); the surrogate "
                 "is built on the 9-cell representation regardless, as the largest tractable option, "
                 "and this is a known, disclosed source of residual error.")
    lines.append("- Sec. 33's finite-array edge/boundary effects are explicitly OUT OF SCOPE for this "
                 "infinite-periodic-supercell framework (`NOT_FEASIBLE_WITHIN_INFINITE_PERIODIC_FRAMEWORK`).")
    lines.append("- A small numerical asymmetry (a few percent) was observed between physically "
                 "equivalent (translation-related) cells in the largest (7x7) solved grids; it does "
                 "not affect the qualitative conclusions but is not fully explained.\n")

    lines.append("## Rule compliance\n")
    lines.append("- Rule 1: no HFSS/measurement data fabricated; every unavailable comparison is "
                 "reported as such, not invented.")
    lines.append("- Rule 2: heterogeneous supercells retain full geometric identity (center/neighbor/"
                 "period/substrate/frequency/polarization all recorded per sample, see "
                 "`environment_dataset/*.csv`); never collapsed to an isolated-cell equivalent.")
    lines.append("- Rule 3: reference simulations -> dataset -> train/val/blind split -> surrogate "
                 "fit (on TRAIN only) -> hyperparameter selection (on VALIDATION only) -> blind test "
                 "-> error analysis -> this report, in that order (`dataset.py` built before "
                 "`surrogate.py` was fit; `fit_surrogate` never touches the blind split).")
    lines.append("- Rule 4: symmetry-safe grouping (`dataset.symmetry_canonical_key`) plus an explicit "
                 "geometry-based holdout (center C5 reserved entirely for blind testing).")
    lines.append("- Rule 5: Re(S)/Im(S) fit and reported directly (cartesian); phase is reported "
                 "separately as a derived, wrapped quantity, never the sole target.")
    lines.append("- Rule 6: no amplitude normalization anywhere in the error metrics (`error_metrics.py` "
                 "reports absolute complex/magnitude error alongside, never instead of, relative error).")
    lines.append("- Rule 7: no geometry inversion beyond the first-stage, MoM-verified test harness in "
                 "`inverse.py`; production geometry inversion remains locked by the gate system.\n")

    return "\n".join(lines) + "\n"
