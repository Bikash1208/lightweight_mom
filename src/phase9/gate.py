"""Sec. 29-30, 33-34, 42: the two independent process gates, the geometry-inversion unlock
condition, and the four-value final status taxonomy. Gate A/B track whether the *validation
process itself* was completed (Sec. 29's checkboxes are process steps, not favorable outcomes);
geometry inversion additionally requires the *outcome* to be favorable (Sec. 31-33) -- a gate can
be fully OPEN (the validation was done properly) while still correctly keeping geometry inversion
LOCKED (because what the validation found was UNVALIDATED/INVALIDATED).
"""

from dataclasses import dataclass

GATE_A_CHECKLIST_KEYS = (
    "hfss_architecture_documented",
    "port_mode_conventions_documented",
    "reference_plane_convention_documented",
    "s11_comparison_performed",
    "s21_comparison_performed",
    "complex_errors_calculated",
    "magnitude_errors_calculated",
    "phase_errors_calculated",
    "physical_architecture_compatibility_established",
    "result_reproducible",
)

GATE_B_CHECKLIST_KEYS = (
    "3x3_heterogeneous_test_completed",
    "multiple_neighbor_contrasts_tested",
    "symmetric_environment_tested",
    "5x5_test_completed",
    "center_cell_response_extracted",
    "periodic_cell_reference_extracted",
    "complex_deviation_calculated",
    "magnitude_deviation_calculated",
    "phase_deviation_calculated",
    "neighbor_distance_contrast_behavior_analyzed",
    "applicability_range_documented",
)

FINAL_STATUS_OPTIONS = (
    "VALIDATED_CONTROL_AND_LPA",
    "PARTIALLY_VALIDATED",
    "NO_VALIDATED_FORMULATION",
    "LOCAL_PERIODICITY_INVALIDATED",
)


@dataclass(frozen=True)
class GateResult:
    checklist: dict
    status: str  # "OPEN" | "CLOSED"


def evaluate_gate_a(checklist):
    missing = [k for k in GATE_A_CHECKLIST_KEYS if k not in checklist]
    if missing:
        raise ValueError(f"gate A checklist missing key(s): {missing}")
    return GateResult(dict(checklist), "OPEN" if all(checklist[k] for k in GATE_A_CHECKLIST_KEYS) else "CLOSED")


def evaluate_gate_b(checklist):
    missing = [k for k in GATE_B_CHECKLIST_KEYS if k not in checklist]
    if missing:
        raise ValueError(f"gate B checklist missing key(s): {missing}")
    return GateResult(dict(checklist), "OPEN" if all(checklist[k] for k in GATE_B_CHECKLIST_KEYS) else "CLOSED")


@dataclass(frozen=True)
class Phase9Decision:
    gate_a: GateResult
    gate_b: GateResult
    s11_status: str
    s21_status: str
    local_periodicity_status: str
    selected_control: str | None
    geometry_inversion_unlocked: bool
    final_status: str


def decide(gate_a: GateResult, gate_b: GateResult, s11_status, s21_status, local_periodicity_status):
    """Sec. 30-33, 42. ``selected_control`` stays `None` unless the evidence actually supports a
    selection (Sec. 28's explicit rule) -- never forced (Sec. 44's prohibitions)."""
    selected_control = None
    if s11_status == "VALIDATED" and s21_status != "VALIDATED":
        selected_control = "S11"
    elif s21_status == "VALIDATED" and s11_status != "VALIDATED":
        selected_control = "S21"
    # if both are VALIDATED, or neither is, no automatic selection is made (Sec. 44's explicit
    # prohibition on selecting by convenience/smaller error -- an ambiguous case where both pass
    # is reported, not resolved, by this function).

    lpa_ok = local_periodicity_status in ("VALIDATED", "PARTIALLY_VALIDATED")
    geometry_inversion_unlocked = (
        gate_a.status == "OPEN" and gate_b.status == "OPEN" and selected_control is not None and lpa_ok
    )

    if local_periodicity_status == "INVALIDATED":
        final_status = "LOCAL_PERIODICITY_INVALIDATED"
    elif geometry_inversion_unlocked and local_periodicity_status == "VALIDATED":
        final_status = "VALIDATED_CONTROL_AND_LPA"
    elif geometry_inversion_unlocked and local_periodicity_status == "PARTIALLY_VALIDATED":
        final_status = "PARTIALLY_VALIDATED"
    else:
        final_status = "NO_VALIDATED_FORMULATION"

    return Phase9Decision(
        gate_a=gate_a,
        gate_b=gate_b,
        s11_status=s11_status,
        s21_status=s21_status,
        local_periodicity_status=local_periodicity_status,
        selected_control=selected_control,
        geometry_inversion_unlocked=geometry_inversion_unlocked,
        final_status=final_status,
    )
