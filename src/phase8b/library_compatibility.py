"""Phase 8B Sec. 26-28, 52-53: the required-quantity-vs-library comparison and the final
decision-tree logic, evaluated once per candidate formulation (never merged into one artificial
project-wide "the" formulation, since no project document picks between them -- see
`holographic_modulation.py`'s module docstring).
"""

from dataclasses import dataclass

from . import control_variable as cv

FINAL_STATES = (
    "FORMULATION_VALIDATED",
    "FORMULATION_DEFINED_BUT_NOT_VALIDATED",
    "PHASE_SOURCE_UNVALIDATED",
    "LIBRARY_CONTROL_VARIABLE_INCOMPATIBLE",
    "FORMULATION_BLOCKED",
)

LIBRARY_COMPATIBILITY_STATES = ("FULLY_COMPATIBLE", "PARTIALLY_COMPATIBLE", "INCOMPATIBLE", "UNVALIDATED")


@dataclass(frozen=True)
class FormulationEvaluation:
    formulation_name: str
    required_quantity: str
    control_variable: cv.ControlVariable
    library_compatibility: str  # LIBRARY_COMPATIBILITY_STATES
    decision_trail: tuple  # ordered list of (question, answer) -- Sec. 28's tree, auditable
    final_status: str  # FINAL_STATES


def evaluate_formulation_a():
    """Sec. 28's decision tree for Formulation A (Fong et al. 2010 interference-modulation
    index, control variable = grounded surface reactance)."""
    trail = [("Does a holographic formulation exist?", "YES (Fong et al. 2010, source_registry.json fong2010)")]
    trail.append(
        ("Does it define a physical control variable?", "YES (surface reactance X(x,y) of a grounded structure)")
    )
    zs = cv.investigate_surface_impedance()
    can_provide = zs.validated != "INVALID"
    trail.append(("Can the Phase 7 library provide that quantity?", "NO" if not can_provide else "YES"))
    if not can_provide:
        compatibility = "INCOMPATIBLE"
        status = "LIBRARY_CONTROL_VARIABLE_INCOMPATIBLE"
        trail.append(
            (
                "Result",
                "LIBRARY GAP: no ground plane (invalidates the standard "
                "Zs<->Gamma relation) and no surface-wave dispersion solver "
                "(cannot characterize the reference wave either)",
            )
        )
    else:
        trail.append(("Is the relationship physically validated?", "NO"))
        compatibility = "UNVALIDATED"
        status = "PHASE_SOURCE_UNVALIDATED"
    return FormulationEvaluation(
        formulation_name="A: Fong et al. 2010 interference-modulation index (grounded holographic impedance surface)",
        required_quantity="surface reactance X(x,y) of a grounded, surface-wave-supporting structure",
        control_variable=zs,
        library_compatibility=compatibility,
        decision_trail=tuple(trail),
        final_status=status,
    )


def evaluate_formulation_b():
    """Sec. 28's decision tree for Formulation B (phase-difference hologram, control variable =
    a normal-incidence co-polarized reflection/transmission phase -- the *kind* of quantity
    Phase 7 actually stores)."""
    trail = [
        (
            "Does a holographic formulation exist?",
            "YES (phase-difference/phase-conjugate hologram, standard reflectarray/transmitarray design equation)",
        )
    ]
    trail.append(
        (
            "Does it define a physical control variable?",
            "YES (normal-incidence co-pol reflection OR transmission phase -- S11xx or S21xx)",
        )
    )
    s11 = cv.investigate_s11xx()
    s21 = cv.investigate_s21xx()
    # Both S11xx and S21xx are the *kind* of quantity (complex, phase-bearing, normal-incidence
    # co-pol S-parameter) this formulation needs and the library genuinely stores in full, so the
    # library "can provide" candidates -- but WHICH one is still an open architecture decision
    # (investigate_s11xx/investigate_s21xx both conclude INSUFFICIENT_EVIDENCE, not INVALID).
    can_provide = True
    trail.append(
        (
            "Can the Phase 7 library provide that quantity?",
            "YES (both S11xx and S21xx are stored in full, req. 5) -- but which one is correct is not yet decided",
        )
    )
    trail.append(
        (
            "Is the relationship physically validated?",
            "NO (reflective-vs-transmissive architecture undecided; S11xx/S21xx both INSUFFICIENT_EVIDENCE)",
        )
    )
    return FormulationEvaluation(
        formulation_name="B: phase-difference hologram (reflectarray/transmitarray-style)",
        required_quantity="normal-incidence co-polarized reflection or transmission phase (S11xx or S21xx)",
        control_variable=s11,  # the currently spec-mandated (Phase 7) choice; s21 is the alternative
        library_compatibility="PARTIALLY_COMPATIBLE",
        decision_trail=tuple(trail),
        final_status="PHASE_SOURCE_UNVALIDATED",
    )


def overall_project_status(evaluations):
    """The one top-level status this project reports (Sec. 52), chosen from the per-formulation
    evaluations without hiding either one: if *any* candidate formulation reaches
    `PHASE_SOURCE_UNVALIDATED` (an open-but-resolvable architecture decision) rather than
    `LIBRARY_CONTROL_VARIABLE_INCOMPATIBLE` (a structural impossibility), that is reported as the
    overall status -- it is the more actionable, more specific finding, and matches Phase 8A's
    own conclusion (continuity, not coincidence: both phases found the same open question). Every
    formulation's own status remains fully visible in ``evaluations``/`holographic_specification.json`,
    never collapsed away."""
    statuses = [e.final_status for e in evaluations]
    if "PHASE_SOURCE_UNVALIDATED" in statuses:
        return "PHASE_SOURCE_UNVALIDATED"
    if all(s == "LIBRARY_CONTROL_VARIABLE_INCOMPATIBLE" for s in statuses):
        return "LIBRARY_CONTROL_VARIABLE_INCOMPATIBLE"
    if "FORMULATION_VALIDATED" in statuses and all(s == "FORMULATION_VALIDATED" for s in statuses):
        return "FORMULATION_VALIDATED"
    return "FORMULATION_DEFINED_BUT_NOT_VALIDATED"


def evaluate_all():
    evaluations = [evaluate_formulation_a(), evaluate_formulation_b()]
    return evaluations, overall_project_status(evaluations)
