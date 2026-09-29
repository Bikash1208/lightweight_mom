"""Phase 8C Sec. 37-38, 53, 62: the factual compatibility matrix (never a ranking/score) and the
final status decision, following Sec. 62's fixed order exactly: topology -> scattering meaning ->
reflection/transmission interpretation -> Phase 8B compatibility -> complex-response adequacy ->
library adequacy -> final status. No step is ever reversed (Sec. 62's explicit "do not reverse
this order").
"""

from dataclasses import dataclass

FINAL_STATUS_TAXONOMY = (
    "CONTROL_VARIABLE_VALIDATED_LIBRARY_ADEQUATE",
    "CONTROL_VARIABLE_VALIDATED_LIBRARY_INSUFFICIENT",
    "CONTROL_VARIABLE_VALIDATED_LIBRARY_UNSUITABLE",
    "CONTROL_VARIABLE_UNVALIDATED",
    "ARCHITECTURE_UNRESOLVED",
    "SCATTERING_DATA_INSUFFICIENT",
    "FORMULATION_SCATTERING_INCOMPATIBLE",
)


@dataclass(frozen=True)
class MatrixRow:
    criterion: str
    s11xx: str
    s21xx: str


def build_compatibility_matrix(s11_investigation, s21_investigation, s11_coverage, s21_coverage,
                                s11_reflection_test, s21_transmission_test, s11_phase_only, s21_phase_only):
    """Sec. 38's exact table -- one row per criterion, S11xx and S21xx side by side, no score, no
    winner."""
    return (
        MatrixRow("Physical interpretation", s11_investigation.physical_definition, s21_investigation.physical_definition),
        MatrixRow("Architecture compatible", s11_reflection_test.status, s21_transmission_test.status),
        MatrixRow("Phase coverage (deg)", f"{s11_coverage.coverage_deg:.2f}", f"{s21_coverage.coverage_deg:.2f}"),
        MatrixRow("Largest circular gap (deg)", f"{s11_coverage.largest_gap_deg:.2f}", f"{s21_coverage.largest_gap_deg:.2f}"),
        MatrixRow("Magnitude range", f"[{s11_coverage.magnitude_min:.4f}, {s11_coverage.magnitude_max:.4f}]",
                  f"[{s21_coverage.magnitude_min:.4f}, {s21_coverage.magnitude_max:.4f}]"),
        MatrixRow("Phase-only approximation", s11_phase_only.status, s21_phase_only.status),
        MatrixRow("Polarization", "co-polarized (x-in, x-out)", "co-polarized (x-in, x-out)"),
        MatrixRow("Holographic mapping (Formulation B)", "PARTIAL (candidate, unresolved)", "PARTIAL (candidate, unresolved)"),
        MatrixRow("Validation status", s11_investigation.validated, s21_investigation.validated),
    )


def library_adequacy(control_variable_validated, coverage_deg, tolerance_deg=None):
    """Sec. 37, 53: ADEQUATE / INSUFFICIENT / INCOMPATIBLE / UNVALIDATED for one candidate.
    Requires the control variable to already be VALIDATED before "adequate" is even meaningful
    (Sec. 52's "physical validity precedes numerical convenience") -- if it is not, this returns
    UNVALIDATED regardless of how much phase coverage exists (never "passes on coverage alone",
    Sec. 52's explicit prohibition)."""
    if control_variable_validated == "INVALID":
        return "INCOMPATIBLE"
    if control_variable_validated != "VALIDATED":
        return "UNVALIDATED"
    # Only reached if a future project decision validates S11xx or S21xx; no threshold is
    # invented here either (Sec. 22's same "report UNRESOLVED rather than invent a threshold"
    # spirit) -- coverage adequacy needs its own explicit, sourced numeric criterion.
    return "UNVALIDATED"


def final_status(architecture_resolved, s11_investigation, s21_investigation, formulation_a_status,
                  reciprocity_ok, propagating_modes_ok):
    """Sec. 53/62's decision tree, followed in the fixed order (topology -> scattering meaning ->
    Phase 8B compatibility -> complex-response adequacy -> library adequacy -> final status)."""
    if not architecture_resolved:
        return "ARCHITECTURE_UNRESOLVED"
    if not propagating_modes_ok:
        return "SCATTERING_DATA_INSUFFICIENT"
    if not reciprocity_ok:
        return "SCATTERING_DATA_INSUFFICIENT"
    if formulation_a_status == "LIBRARY_CONTROL_VARIABLE_INCOMPATIBLE" and \
       s11_investigation.validated != "VALIDATED" and s21_investigation.validated != "VALIDATED":
        # Formulation A structurally incompatible; Formulation B's own control variable remains
        # unresolved (neither S11xx nor S21xx independently validated) -- the honest, currently
        # correct answer for this project, not a forced choice.
        return "CONTROL_VARIABLE_UNVALIDATED"
    if s11_investigation.validated == "VALIDATED" or s21_investigation.validated == "VALIDATED":
        adequacy = library_adequacy(
            s11_investigation.validated if s11_investigation.validated == "VALIDATED" else s21_investigation.validated,
            None,
        )
        return {
            "ADEQUATE": "CONTROL_VARIABLE_VALIDATED_LIBRARY_ADEQUATE",
            "INSUFFICIENT": "CONTROL_VARIABLE_VALIDATED_LIBRARY_INSUFFICIENT",
            "INCOMPATIBLE": "CONTROL_VARIABLE_VALIDATED_LIBRARY_UNSUITABLE",
            "UNVALIDATED": "CONTROL_VARIABLE_UNVALIDATED",
        }[adequacy]
    return "CONTROL_VARIABLE_UNVALIDATED"
