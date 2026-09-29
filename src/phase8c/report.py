"""Phase 8C Sec. 49: `control_quantity_analysis.md`, assembled entirely from already-computed
objects passed in -- this module only formats text."""


def render_control_quantity_analysis_md(
    audit,
    power_summary,
    mode_check,
    s11_investigation,
    s21_investigation,
    s11_coverage,
    s21_coverage,
    s11_states,
    s21_states,
    s11_phase_only,
    s21_phase_only,
    reflection_test,
    transmission_test,
    formulation_a,
    formulation_b,
    reciprocity,
    symmetry,
    matrix,
    status,
    missing_evidence,
):
    lines = [
        "# control_quantity_analysis.md",
        "",
        "## 1. Physical unit-cell architecture",
        "",
        f"{audit.architecture_classification} -- {audit.architecture_evidence}. Ground plane: "
        f"**{audit.ground_plane}** ({audit.ground_plane_evidence})",
        "",
        "## 2. Geometry/topology",
        "",
        audit.metal_layers,
        "",
        audit.substrate_layers,
        "",
        audit.air_regions,
        "",
        "## 3. Port definition",
        "",
        "See `s_parameter_definition.md` (CONVENTIONS.md/S_PARAMETERS.md, reused verbatim).",
        "",
        "## 4. S-parameter definitions",
        "",
        "See `s_parameter_definition.md`.",
        "",
        "## 5. Candidate control quantities",
        "",
        "S11xx, S21xx (investigated in full below); S11yy/S21yy/cross-polarized terms investigated "
        "via the power-balance analysis (Sec. 6-8 below) using the full linear-basis matrix read "
        "from the internal Phase 7 checkpoint (`phase8c.data_source`).",
        "",
        "## 6. Reflection analysis",
        "",
        reflection_test.hypothesis,
        "",
        " -> ".join(f"[{link}: {st}]" for link, st in reflection_test.chain),
        "",
        f"**{reflection_test.status}**",
        "",
        "## 7. Transmission analysis",
        "",
        transmission_test.hypothesis,
        "",
        " -> ".join(f"[{link}: {st}]" for link, st in transmission_test.chain),
        "",
        f"**{transmission_test.status}**",
        "",
        "## 8. Power balance",
        "",
        f"Library-mean (x-incidence): R_co={power_summary['mean_R_co']:.4f}, "
        f"R_cross={power_summary['mean_R_cross']:.4f}, T_co={power_summary['mean_T_co']:.4f}, "
        f"T_cross={power_summary['mean_T_cross']:.4f}, mean(R+T)={power_summary['mean_R_plus_T']:.4f} "
        f"(range [{power_summary['min_R_plus_T']:.4f}, {power_summary['max_R_plus_T']:.4f}]).",
        "",
        "## 9. Modal analysis",
        "",
        mode_check["note"],
        "",
        "## 10. Complex-response analysis",
        "",
        f"S11xx: {s11_states.n_geometries} geometries, {s11_states.n_distinct_complex_states} "
        f"distinct complex states, {s11_states.n_distinct_phase_states} distinct phase states.",
        f"S21xx: {s21_states.n_geometries} geometries, {s21_states.n_distinct_complex_states} "
        f"distinct complex states, {s21_states.n_distinct_phase_states} distinct phase states.",
        "",
        "## 11. Phase coverage",
        "",
        f"S11xx: {s11_coverage.coverage_deg:.2f} deg of 360 deg (largest gap "
        f"{s11_coverage.largest_gap_deg:.2f} deg). S21xx: {s21_coverage.coverage_deg:.2f} deg of "
        f"360 deg (largest gap {s21_coverage.largest_gap_deg:.2f} deg).",
        "",
        "## 12. Magnitude/phase coupling",
        "",
        f"S11xx: magnitude in [{s11_coverage.magnitude_min:.4f}, {s11_coverage.magnitude_max:.4f}] "
        f"(median {s11_coverage.magnitude_median:.4f}, std {s11_coverage.magnitude_std:.4f}); "
        f"phase-magnitude correlation {s11_coverage.phase_magnitude_correlation}.",
        f"S21xx: magnitude in [{s21_coverage.magnitude_min:.4f}, {s21_coverage.magnitude_max:.4f}] "
        f"(median {s21_coverage.magnitude_median:.4f}, std {s21_coverage.magnitude_std:.4f}); "
        f"phase-magnitude correlation {s21_coverage.phase_magnitude_correlation}.",
        f"Phase-only approximation: S11xx {s11_phase_only.status}; S21xx {s21_phase_only.status}.",
        "",
        "## 13. Polarization",
        "",
        "Both S11xx/S21xx are co-polarized (x-in, x-out); an x-polarized response cannot, by "
        "itself, encode an arbitrary object-wave polarization without an additional assumption "
        "(Phase 8B `control_variable_analysis.md` Sec. 6).",
        "",
        "## 14. Reference-plane considerations",
        "",
        "See `reference_plane_sensitivity.py`: phase is reference-plane dependent by a documented, "
        "reused (not re-derived) formula; not corrected or offset anywhere in this analysis.",
        "",
        "## 15. Phase 8B compatibility",
        "",
        f"Formulation A: {formulation_a['status']} (unchanged from Phase 8B).",
        f"Formulation B: corresponds to reflection phase -> {formulation_b['corresponds_to_reflection_phase']}; "
        f"to transmission phase -> {formulation_b['corresponds_to_transmission_phase']}.",
        "",
        "## 16. Library adequacy",
        "",
        "See `compatibility_matrix.py`'s factual table below.",
        "",
        "| Criterion | S11xx | S21xx |",
        "|---|---|---|",
    ]
    for row in matrix:
        lines.append(f"| {row.criterion} | {row.s11xx} | {row.s21xx} |")
    lines += [
        "",
        "## 17. Missing evidence",
        "",
    ]
    lines += [f"- {e}" for e in missing_evidence] or ["- None beyond what is stated above."]
    lines += [
        "",
        "## 18. Final scientific status",
        "",
        f"**{status}**",
        "",
        "Reciprocity check: " + reciprocity.note,
        "",
        "Symmetry check: " + symmetry.note,
    ]
    return "\n".join(lines) + "\n"
