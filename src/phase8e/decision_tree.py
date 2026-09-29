"""Sec. 20, 22, 25, 29: the fixed-order per-formulation decision tree, the 9-checkbox geometry-
inversion gate, and the 13-value status taxonomy. No positive result is forced (Sec. 22): the
project-wide result may legitimately be `NO_VALIDATED_FORMULATION`.
"""

from dataclasses import dataclass

STATUS_TAXONOMY = (
    "FORMULATION_VALIDATED",
    "CONTROL_VARIABLE_VALIDATED",
    "CURRENT_FAMILY_COMPATIBLE",
    "CURRENT_FAMILY_INSUFFICIENT",
    "REQUIRED_LIBRARY_MISSING",
    "REFERENCE_PLANE_MISMATCH",
    "ARCHITECTURE_INCOMPATIBLE",
    "COMPLEX_RESPONSE_INSUFFICIENT",
    "LOCAL_PERIODICITY_UNVALIDATED",
    "CONTINUOUS_RECONSTRUCTION_FAILED",
    "DISCRETE_RECONSTRUCTION_FAILED",
    "GEOMETRY_INVERSION_BLOCKED",
    "NO_VALIDATED_FORMULATION",
)

GATE_CHECKLIST_KEYS = (
    "formulation_validated",
    "control_variable_validated",
    "planes_compatible",
    "required_response_defined",
    "suitable_family_exists",
    "response_within_validated_library_coverage",
    "local_periodicity_acceptable",
    "continuous_reconstruction_passes",
    "discrete_reconstruction_passes",
)


@dataclass(frozen=True)
class FormulationDecision:
    formulation: str
    decision_trail: tuple  # ordered (question, answer), Sec. 29's tree, auditable
    gate_checklist: dict
    geometry_inversion_permitted: bool
    final_status: str


def evaluate_formulation_gate(
    formulation,
    architecture_status,
    reference_plane_status,
    family_exists,
    continuous_reconstruction_ok,
    discrete_reconstruction_ok,
    control_variable_validated,
    formulation_validated=False,
    local_periodicity_acceptable=False,
):
    """Sec. 29's fixed order, stopping at the first structural block and reporting every later
    checkbox as unmet rather than skipped silently."""
    trail = [("architecture compatibility", architecture_status)]
    if architecture_status == "ARCHITECTURE_INCOMPATIBLE":
        checklist = {k: False for k in GATE_CHECKLIST_KEYS}
        checklist["required_response_defined"] = True  # control_response.py always defines C_req
        return FormulationDecision(formulation, tuple(trail), checklist, False, "ARCHITECTURE_INCOMPATIBLE")

    trail.append(("reference-plane compatibility", reference_plane_status))
    planes_ok = reference_plane_status == "COMPATIBLE"
    trail.append(("suitable unit-cell family registered", "YES" if family_exists else "NO"))
    if not family_exists:
        checklist = {k: False for k in GATE_CHECKLIST_KEYS}
        checklist["required_response_defined"] = True
        checklist["planes_compatible"] = planes_ok
        return FormulationDecision(formulation, tuple(trail), checklist, False, "REQUIRED_LIBRARY_MISSING")

    trail.append(("continuous beam reconstruction", "PASSED" if continuous_reconstruction_ok else "FAILED"))
    if not continuous_reconstruction_ok:
        checklist = {k: False for k in GATE_CHECKLIST_KEYS}
        checklist.update(required_response_defined=True, planes_compatible=planes_ok, suitable_family_exists=True)
        return FormulationDecision(formulation, tuple(trail), checklist, False, "CONTINUOUS_RECONSTRUCTION_FAILED")

    trail.append(("discrete beam reconstruction", "PASSED" if discrete_reconstruction_ok else "FAILED"))
    if not discrete_reconstruction_ok:
        checklist = {k: False for k in GATE_CHECKLIST_KEYS}
        checklist.update(
            required_response_defined=True,
            planes_compatible=planes_ok,
            suitable_family_exists=True,
            continuous_reconstruction_passes=True,
        )
        return FormulationDecision(formulation, tuple(trail), checklist, False, "DISCRETE_RECONSTRUCTION_FAILED")

    trail.append(
        (
            "physical control variable validated",
            "YES" if control_variable_validated else "NO (Phase 8B/8C: CONTROL_VARIABLE_UNVALIDATED, unchanged)",
        )
    )
    trail.append(
        (
            "local periodicity acceptable",
            "YES"
            if local_periodicity_acceptable
            else "NO (local_periodicity.py: no coupled-cell solver exists to validate it)",
        )
    )
    checklist = {
        "formulation_validated": formulation_validated,  # phase8b.library_compatibility's own status
        "control_variable_validated": control_variable_validated,
        "planes_compatible": planes_ok,
        "required_response_defined": True,
        "suitable_family_exists": True,
        "response_within_validated_library_coverage": True,  # data exists to compare against (Sec. 20's 6th box)
        "local_periodicity_acceptable": local_periodicity_acceptable,
        "continuous_reconstruction_passes": True,
        "discrete_reconstruction_passes": True,
    }
    permitted = all(checklist.values())
    status = "CONTROL_VARIABLE_VALIDATED" if permitted else "GEOMETRY_INVERSION_BLOCKED"
    return FormulationDecision(formulation, tuple(trail), checklist, permitted, status)


def overall_project_status(formulation_decisions):
    """Sec. 22, 32-style aggregate: `NO_VALIDATED_FORMULATION` unless at least one formulation's
    gate is fully open (never forced positive, Sec. 22)."""
    if any(d.geometry_inversion_permitted for d in formulation_decisions):
        return "CONTROL_VARIABLE_VALIDATED"
    return "NO_VALIDATED_FORMULATION"
