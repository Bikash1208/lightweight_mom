"""``phase11_report.md``: Sec. 29's acceptance-criteria checklist, Sec. 35's outcome classification,
and the key measured numbers, following Phase 9/10's own report style."""


def render_phase11_report_md(result):
    d = result["decision"]
    lines = []
    lines.append("# Phase 11 report: physics-based heterogeneous MoM / nonlocal coupling\n")
    lines.append(f"**Overall classification: `{d.overall_case}`**\n")
    lines.append(f"Geometry inversion: **{d.geometry_inversion_status}**  ")
    lines.append(f"STEP export: **{d.step_export_status}**\n")

    lines.append("## No ML / no surrogate\n")
    lines.append("Confirmed by construction: every response in this phase comes from a real solve of "
                 "`Z I = V` (`reference_solver.solve_grid`/`solve_isolated`, the unmodified Phase 1-3 "
                 "MoM machinery); no regression, neural network, random forest or learned mapping "
                 "appears anywhere in `src/phase11/`.\n")

    lines.append("## Coupling-radius series (Sec. 5-7)\n")
    for name, points in result["coupling_radius"].items():
        vals = ", ".join(f"n={p.n}:{abs(p.response):.3e}" for p in points)
        lines.append(f"- {name}: {vals}")
    lines.append("")

    lines.append("## Shell increments (Sec. 7)\n")
    for name, incs in result["shell_increments"].items():
        vals = ", ".join(f"{inc.n_from}->{inc.n_to}:{inc.delta_relative:.3e}" for inc in incs)
        lines.append(f"- {name}: {vals}")
    lines.append("")

    lines.append("## Convergence classification (Sec. 16, tolerances declared before evaluation)\n")
    for cr in result["convergence_results"]:
        lines.append(f"- {cr.center_name}: status={cr.status}, last_increment_relative="
                     f"{cr.last_increment_relative:.3e}, converged_at_n={cr.converged_at_n}")
    lines.append(f"\nPer-center case: {result['classification'].per_center_case}")
    lines.append(f"\nOverall case: **{result['classification'].overall_case}** -- {result['classification'].rationale}\n")

    lines.append("## Diagonal coupling (Sec. 9)\n")
    diag = result["diagonal_reproduction"]
    lines.append(f"mean|diagonal|/mean|cardinal| = {diag['ratio']:.3f} "
                 f"(Phase 10's own measured range: 25-70%)\n")

    lines.append("## Pattern study / Case D check (Sec. 10, 18)\n")
    for name, r in result["pattern_results"].items():
        lines.append(f"- {name}: |response|={abs(r.response):.3e}")
    lines.append(f"\nPattern dependence detected: **{result['pattern_dependence_detected']}** "
                 f"({result['pattern_comparisons']})\n")

    lines.append("## Pairwise / higher-order interaction (Sec. 12-13)\n")
    for r in result["pairwise_results"]:
        lines.append(f"- {r.pair_label}: interaction_error_relative={r.interaction_error_relative:.3e}")
    ho = result["higher_order_result"]
    lines.append(f"- {ho.combo_label} (triple): residual_relative={ho.residual_relative:.3e}\n")

    lines.append("## Frequency dependence (Sec. 14)\n")
    for p in result["frequency_points"]:
        lines.append(f"- {p.frequency_hz/1e9:.2f} GHz: max_shell_increment_relative={p.max_shell_increment_relative:.3e}")
    lines.append("")

    lines.append("## Polarization dependence (Sec. 15)\n")
    for p in result["polarization_points"]:
        lines.append(f"- {p.polarization}: max_shell_increment_relative={p.max_shell_increment_relative:.3e}")
    lines.append(f"\nCircular polarization: {result['circular_status']}\n")

    lines.append("## Decay law (Sec. 17 -- no automatic fitting)\n")
    obs = result["decay_observation"]
    fit = result["decay_fit"]
    lines.append(f"Observed trend: {obs.candidate_trend}. Fit attempted: {fit.law}.")
    if fit.law != "NONE_FITTED":
        lines.append(f"  fitted_params={fit.fitted_params}, r_squared={fit.r_squared:.3f}, "
                     f"residual_rms={fit.residual_rms:.3e}")
    lines.append(f"  note: {fit.note}\n")

    lines.append("## Reduction (Sec. 19)\n")
    lines.append(f"{result['reduction_decision'].status}: {result['reduction_decision'].rationale}\n")

    lines.append("## Inverse validation (Sec. 21-22)\n")
    for r in result["inverse_results"]:
        lines.append(f"- target={r.target:.3e}: status={r.status}, relative_error="
                     f"{r.mom_relative_error if r.mom_relative_error != float('inf') else 'inf'}")
    lines.append("")

    lines.append("## Beam-level validation (Sec. 23)\n")
    for name, case in result["beam_cases"].items():
        lines.append(f"- {name}: pointing_error_deg={case['metrics'].pointing_error_deg}")
    lines.append("")

    lines.append("## Three-beam holographic test (Sec. 24)\n")
    tb = result["three_beam_result"]
    lines.append(f"pointing_error_deg={tb.metrics.pointing_error_deg}, "
                 f"matched_peak_power_db={tb.metrics.matched_peak_power_db}\n")

    lines.append("## HFSS status (Sec. 25-26)\n")
    lines.append(f"{d.hfss_status}\n")

    lines.append("## Sec. 29 acceptance criteria\n")
    for key, val in result["acceptance_criteria"].items():
        lines.append(f"- [{'x' if val else ' '}] {key}")
    lines.append("")

    lines.append("## Sec. 35 final objective\n")
    outcome_map = {
        "CASE_A_LOCAL": "Outcome A: a finite spatial neighborhood is sufficient.",
        "CASE_B_MODERATELY_NONLOCAL": "Outcome B (moderate form): a larger but still finite "
                                       "neighborhood is required; isolated-cell inversion remains insufficient.",
        "CASE_C_STRONGLY_NONLOCAL": "Outcome B (strong form): the metasurface exhibits significant "
                                     "nonlocal interaction; isolated-cell inversion is fundamentally insufficient.",
        "CASE_D_PATTERN_DEPENDENT": "Outcome C: the response depends strongly on spatial pattern, "
                                     "not just local statistics.",
    }
    lines.append(outcome_map.get(d.overall_case, "Outcome D: further solver validation is required.") + "\n")

    return "\n".join(lines) + "\n"
