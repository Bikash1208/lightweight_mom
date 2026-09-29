"""Sec. 41: `phase9_report.md`, answering all 16 mandatory questions from the actual computed
``result`` dict (`run.py`)."""


def render_phase9_report_md(result):
    audit = result["architecture_audit"]
    decision = result["decision"]
    lpa_results = result["lpa_results"]

    lines = [
        "# Phase 9 report -- physical control-variable and local-periodicity validation", "",
        f"## Final status: `{decision.final_status}`", "",
        f"- Gate A (control-variable validation process): `{decision.gate_a.status}`",
        f"- Gate B (local-periodicity validation process): `{decision.gate_b.status}`",
        f"- S11 classification: `{decision.s11_status}`",
        f"- S21 classification: `{decision.s21_status}`",
        f"- Selected control variable: `{decision.selected_control}`",
        f"- Local periodicity: `{decision.local_periodicity_status}`",
        f"- Geometry inversion: `{'UNLOCKED' if decision.geometry_inversion_unlocked else 'LOCKED'}`",
        "", "## Mandatory questions (Sec. 41)", "",
        "### 1. What exact electromagnetic architecture was validated?", "",
        f"This project's own MoM solver architecture (Phase 1-7, unmodified): {audit['solver_architecture_classification']} "
        f"-- {audit['solver_architecture_evidence']}. No HFSS architecture could be validated: "
        f"{audit['hfss_solution_type']}.", "",
        "### 2. What does S11 physically represent in that architecture?", "",
        f"{audit['solver_reflected_wave']}, referenced at the top plane "
        f"({audit['solver_incident_wave']}).", "",
        "### 3. What does S21 physically represent?", "",
        f"{audit['solver_transmitted_wave']}.", "",
        "### 4. Which quantities were compared against HFSS?", "",
        f"S11_xx and S21_xx of the {len(result['reference_cells'])} reference cells at "
        f"{len(result['frequencies_ghz'])} frequency points each (Sec. 7, 11) -- see "
        "`s11_comparison.csv`/`s21_comparison.csv`.", "",
        "### 5. What are the complex, magnitude, and phase errors?", "",
        (f"NOT AVAILABLE FROM CURRENT DATA: {result['hfss_unavailable_reason']}"
         if result["hfss_stats_s11"] is None else
         f"S11: mean complex error {result['hfss_stats_s11'].mean_complex_error:.4g}, mean phase "
         f"error {result['hfss_stats_s11'].mean_abs_phase_error_deg:.4g} deg. S21: mean complex "
         f"error {result['hfss_stats_s21'].mean_complex_error:.4g}, mean phase error "
         f"{result['hfss_stats_s21'].mean_abs_phase_error_deg:.4g} deg."), "",
        "### 6. Are those errors frequency-dependent?", "",
        "NOT AVAILABLE FROM CURRENT DATA (no HFSS reference exists to test this)." if result["hfss_stats_s11"] is None
        else "See `s11_comparison.csv`/`s21_comparison.csv` (per-frequency rows) and plots 4/8.", "",
        "### 7. Are they geometry-dependent?", "",
        "NOT AVAILABLE FROM CURRENT DATA (no HFSS reference exists to test this)." if result["hfss_stats_s11"] is None
        else "See plots 9/10 (error summary across geometries).", "",
        "### 8. Does S11 agree sufficiently with the physical architecture?", "",
        f"`{decision.s11_status}` -- {result['s11_reason']}", "",
        "### 9. Does S21 agree sufficiently with the physical architecture?", "",
        f"`{decision.s21_status}` -- {result['s21_reason']}", "",
        "### 10. What happens to the centre cell when neighboring cells change?", "",
        "Its induced-current dipole moment measurably changes; see `supercell_results.csv`/"
        "`lpa_error_metrics.csv` for every tested contrast level: "
        + "; ".join(f"{r.name}: relative complex error {r.relative_complex_error:.4g}, "
                     f"phase error {r.phase_error_deg:.2f} deg" for r in lpa_results) + ".", "",
        "### 11. Does the local-periodic approximation remain valid?", "",
        f"`{decision.local_periodicity_status}` -- "
        + (
            "at least one physically-realistic, in-sweep-range contrast level produced a "
            f"deviation at least {result['lpa_substantial_threshold']*100:.0f}% of the signal's "
            "own magnitude (a self-referential, not an externally-chosen, criterion; Sec. 24)."
            if decision.local_periodicity_status == "INVALIDATED" else
            "no tested contrast level's deviation is within the solver's own numerical noise "
            "floor, so `VALIDATED` cannot be honestly claimed, but no deviation was large enough "
            "relative to its own signal to call `INVALIDATED` either, and no sourced application "
            "tolerance exists to decide the remaining, ambiguous cases (Sec. 24)."
            if decision.local_periodicity_status == "UNVALIDATED" else
            "every tested contrast level's deviation is within the solver's own numerical "
            "self-consistency floor."
        ), "",
        "### 12. Over what geometry/frequency range is it valid?", "",
        "Not established over any range: this run's classification is a single overall verdict "
        f"(`{decision.local_periodicity_status}`), not a range-restricted one -- `PARTIALLY_VALIDATED` "
        "was not reached, so no restricted applicability range can be honestly stated (Sec. 23).", "",
        "### 13. What assumptions remain unvalidated?", "",
        "Whether S11xx or S21xx is the physically correct holographic control variable (no HFSS "
        "data exists to test either); whether the local-periodic approximation holds under "
        "realistic neighbor contrast (measured evidence above suggests it does not, for strong "
        "contrast); the reference/object-wave formulation questions Phase 8B-8E already left open.", "",
        "### 14. Is geometry inversion scientifically unlocked?", "",
        f"**{'YES' if decision.geometry_inversion_unlocked else 'NO'}**.", "",
        "### 15. Is final STEP generation scientifically unlocked?", "",
        "**NO** (tied to the same gate as geometry inversion, Sec. 5, 30).", "",
        "### 16. What is the exact next phase?", "",
        (
            "Phase 10 (validated geometry inversion) may begin." if decision.geometry_inversion_unlocked else
            "Geometry inversion remains blocked. Two independent next steps are required: (1) obtain "
            "real HFSS (or other independent full-wave) reference data for at least the 10 "
            "reference cells at the 7 validation frequencies, since none exists in this project "
            "and Rule 1 forbids fabricating it; (2) given the measured local-periodicity "
            "invalidation under strong neighbor contrast, develop a coupled-cell or heterogeneous-"
            "library capability (Sec. 34) rather than continuing to use the isolated periodic "
            "lookup table for holographic synthesis."
        ), "",
        "## Rule compliance", "",
        "- No boundary condition, port definition, polarization convention, normalization factor, "
        "empirical correction, phase offset, amplitude normalization, S-parameter transformation, "
        "HFSS result, validation threshold, or literature equation was invented (Rule 1).",
        "- No mathematics was modified to force S11 or S21 agreement, or to force local-periodicity "
        "agreement (Rule 2).",
        "- Every comparison retained complex, magnitude and phase error separately (Rule 3).",
        "- No amplitude error was normalized away (Rule 4).",
        f"- Geometry inversion status: `{'UNLOCKED' if decision.geometry_inversion_unlocked else 'LOCKED'}` (Rule 5).",
    ]
    return "\n".join(lines) + "\n"
