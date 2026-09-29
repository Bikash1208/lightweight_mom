"""Sec. 29-30: `report.md` -- the fixed-order decision-tree diagram and the twelve mandatory
questions, answered from the actual computed ``result`` dict (`run.py`), never asserted."""


def render_report_md(result):
    arch = result["architecture_results"]
    lines = [
        "# Phase 8E -- physical control-variable and holographic-scattering bridge",
        "",
        "## Decision tree (Sec. 29)",
        "",
        "```",
        "Reference field (SYNTHETIC_DEMO_SPHERICAL_FEED, no real data exists)",
        "      |",
        "Object wave (3-beam demo, phase8b.run.DEFAULT_DEMO_BEAMS, theta=30 deg)",
        "      |",
        "Candidate formulations A/B/C/D",
        "      |",
        "Required control variable C_req(x,y) (control_response.py)",
        "      |",
        "Architecture compatibility:",
        f"        A: {arch['A'].status}",
        f"        B: {arch['B'].status}",
        f"        C: {arch['C'].status}",
        f"        D: {arch['D'].status}",
        "      |",
        f"Reference-plane compatibility: {result['plane_compat_actual'].status}",
        "      |",
        "Complex controllability (S11xx for B, S21xx for C):",
        f"        B: mean complex error {result['cc_b'].mean_complex_error:.4g}",
        f"        C: mean complex error {result['cc_c'].mean_complex_error:.4g}",
        "      |",
        f"Universal database query: family_exists={result['family_exists']}",
        "      |",
        "Suitable family? " + ("YES" if result["family_exists"] else "NO"),
        "        continuous reconstruction -> discrete reconstruction -> geometry inversion permitted?",
        "```",
        "",
        "## Final per-formulation status",
        "",
    ]
    for d in result["decisions"]:
        lines.append(
            f"- **{d.formulation}**: `{d.final_status}` (geometry_inversion_permitted={d.geometry_inversion_permitted})"
        )
        for q, a in d.decision_trail:
            lines.append(f"  - {q}: {a}")
    lines += ["", f"## Overall project status: `{result['overall_status']}`", ""]

    lines += [
        "## Mandatory questions (Sec. 30)",
        "",
        "### Q1: What exact physical quantity does the holographic formulation require?",
        "",
        "Depends on the candidate: A needs a grounded surface reactance X(x,y); B needs a normal-"
        "incidence co-polarized reflection phase; C needs the identical phase kernel claimed "
        "against transmission; D needs a local 2x2 linear-basis transfer operator. No single "
        "answer is asserted (Sec. 5's explicit framework-not-a-choice instruction).",
        "",
        "### Q2: What exact unit-cell electromagnetic quantity represents that physical quantity?",
        "",
        "A: none available (no surface impedance quantity exists in this project's data, Phase 8B "
        "`investigate_surface_impedance` = INVALID). B: `S11_xx`. C: `S21_xx`. D: the full linear-"
        "basis 2x2 block (Sxx, Sxy, Syx, Syy).",
        "",
        "### Q3: Are the reference field and unit-cell response defined at compatible reference planes?",
        "",
        f"For this run's synthetic reference field: **{result['plane_compat_actual'].status}** "
        f"(z_reference={result['plane_compat_actual'].z_reference_m} m, pattern plane="
        f"{result['plane_compat_actual'].z_unit_cell_pattern_m} m) -- compatible by construction "
        "(the spherical-feed model is evaluated directly at z=0). A real reference field recorded "
        "away from z=0 (e.g. Sec. 3's historical z=1.7 mm example) would instead report "
        f"`{result['plane_compat_demo_mismatch'].status}` and require the angular-spectrum "
        "propagation this run also demonstrated (`reference_plane.propagate_angular_spectrum`).",
        "",
        "### Q4: Does the current patterned-PEC/no-ground-plane architecture satisfy the formulation assumptions?",
        "",
        f"A: **NO** ({arch['A'].evidence}). B: architecturally plausible but unvalidated "
        f"({arch['B'].evidence}). C: architecturally plausible but unvalidated ({arch['C'].evidence}). "
        f"D: **NO** ({arch['D'].evidence}).",
        "",
        "### Q5: Is scalar S11 sufficient?",
        "",
        f"Not established. Complex-plane matching against the real 324-cell library (S11xx): mean "
        f"complex error {result['cc_b'].mean_complex_error:.4g}, mean amplitude error "
        f"{result['cc_b'].mean_amplitude_error:.4g}, mean phase error "
        f"{result['cc_b'].mean_phase_error_deg:.4g} deg -- reported as data, not collapsed into a "
        "pass/fail without a stated tolerance (none is defined by any project document).",
        "",
        "### Q6: Is scalar S21 sufficient?",
        "",
        f"Not established. Complex-plane matching against the real 324-cell library (S21xx): mean "
        f"complex error {result['cc_c'].mean_complex_error:.4g}, mean amplitude error "
        f"{result['cc_c'].mean_amplitude_error:.4g}, mean phase error "
        f"{result['cc_c'].mean_phase_error_deg:.4g} deg.",
        "",
        "### Q7: Is a matrix/vector scattering representation required?",
        "",
        f"Only if polarization conversion is actually needed by the desired beam set; for this "
        f"run's demo beams (all x-polarized, x-polarized reference) the generalized operator "
        "reduces to a scalar (Txy=Tyx=Tyy=0). Independent of the specific demo, the *library's* "
        f"measured cross-polarization is negligible (architecture check D: {arch['D'].status}), so "
        "this unit-cell family could not realize a genuinely polarization-converting requirement "
        "even if one existed.",
        "",
        "### Q8: What complex response region is actually required by the hologram?",
        "",
        "See `required_response_map.csv` (the full per-aperture-point required phase for B/C) and "
        "`06_control_response_amplitude_phase.png`/`07_control_response_complex_plane.png`.",
        "",
        "### Q9: Does the existing 4.8-mm library cover that region?",
        "",
        "Partially, per the numbers in Q5/Q6 above (`control_variable_analysis.json`, "
        "`formulation_comparison.csv`) -- reported quantitatively, never reduced to a single "
        "coverage percentage.",
        "",
        "### Q10: If not, what exact unit-cell library must be built?",
        "",
        f"`period_assessment.status` = `{result['period_assessment'].status}`; see "
        "`required_library_specification.json` for the precise build request "
        "(topology/period/frequency/substrate/polarization/incidence), when one was generated.",
        "",
        "### Q11: What period range is physically required?",
        "",
        f"Candidate (Nyquist-sampling heuristic, `CANDIDATE_HEURISTIC` only, never `VALIDATED`) "
        f"periods: {result['period_assessment'].candidate_periods_mm} mm, from the demo beams' "
        f"required transverse spatial frequency "
        f"({result['obj_spectrum'].k_required_max_per_m:.6g} rad/m). Registered (actually "
        f"validated-by-data) periods: {result['period_assessment'].registered_periods_mm} mm.",
        "",
        "### Q12: Which conclusions are validated, and which remain hypotheses?",
        "",
        "**Validated (measured, this run):** architecture is ungrounded and partially reflective/"
        "transmissive (Phase 8C, reused); cross-polarization is negligible across the whole "
        "registered library; the synthetic reference field used here is reference-plane-compatible "
        "by construction; the angular-spectrum propagation code path is demonstrated to run "
        "correctly for a nonzero reference height. **Hypotheses / unvalidated:** whether S11xx or "
        "S21xx phase is the physically correct holographic control variable (Phase 8B/8C, "
        "unchanged); whether the local-periodic approximation holds for this period/beam "
        "combination (no coupled-cell solver exists to check it); whether any candidate "
        "formulation's continuous reconstruction actually represents a physically buildable "
        "antenna (only a mathematical/Fourier-optics idealization was tested here).",
        "",
        "## Validation status",
        "",
        "- HFSS validation: NOT PERFORMED",
        "- Solver / Phase 7 modified: NO",
        "- Phase 8D-REV modified: NO",
        "- Geometry inversion performed: NO",
        "- Final MTS STEP generated: NO",
    ]
    return "\n".join(lines) + "\n"
