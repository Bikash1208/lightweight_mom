"""Sec. 41-42, 52: the five process gates and the final geometry-inversion unlock condition.
Exactly Phase 9's own gate philosophy (`phase9.gate`, reused conceptually, not by import, since
the checklist keys differ): a gate tracks whether the validation PROCESS was completed, which is
independent of whether the outcome was favorable -- Gate 10B can be fully OPEN (the neighborhood
study was genuinely run) while its own finding is that neither 5 nor 9 cells is fully sufficient.
"""

from dataclasses import dataclass

GATE_10A_KEYS = (
    "hfss_or_measurement_reference_exists",
    "architecture_documented",
    "s11_s21_convention_documented",
    "complex_response_verified",
    "selected_control_physically_justified",
)

GATE_10B_KEYS = (
    "five_cell_study_completed",
    "cardinal_neighbor_effects_quantified",
    "nonlinear_interactions_tested",
    "diagonal_effects_tested",
    "five_vs_nine_cell_sufficiency_established",
    "longer_range_effects_investigated",
)

GATE_10C_KEYS = (
    "independent_blind_test",
    "strong_contrast_cases_included",
    "interpolation_cases_included",
    "unseen_environment_cases_included",
    "holographic_gradient_test",
    "abrupt_boundary_test",
    "complex_error_reported",
    "magnitude_error_reported",
    "phase_error_reported",
    "worst_case_error_reported",
    "valid_domain_documented",
)

GATE_10D_KEYS = (
    "target_response_inversion_implemented",
    "inverse_tested_on_unseen_targets",
    "proposed_geometries_checked_with_mom",
    "unreachable_targets_correctly_rejected",
    "strong_neighbor_cases_included",
)

GATE_10E_KEYS = (
    "environment_aware_aperture_constructed",
    "heterogeneous_mom_reference_generated",
    "surrogate_prediction_generated",
    "far_field_compared",
    "beam_pointing_compared",
    "peak_power_compared",
    "sidelobe_behavior_compared",
)

FINAL_STATUS_OPTIONS = (
    "ENVIRONMENT_MODEL_VALIDATED",
    "ENVIRONMENT_MODEL_PARTIALLY_VALIDATED",
    "SURROGATE_INVALIDATED",
    "INVERSE_MODEL_INVALIDATED",
    "END_TO_END_MODEL_INVALIDATED",
    "NO_VALIDATED_FORMULATION",
)


@dataclass(frozen=True)
class GateResult:
    checklist: dict
    status: str   # "OPEN" | "CLOSED"


def _evaluate(keys, checklist, gate_name):
    missing = [k for k in keys if k not in checklist]
    if missing:
        raise ValueError(f"{gate_name} checklist missing key(s): {missing}")
    return GateResult(dict(checklist), "OPEN" if all(checklist[k] for k in keys) else "CLOSED")


def evaluate_gate_10a(checklist):
    return _evaluate(GATE_10A_KEYS, checklist, "gate 10A")


def evaluate_gate_10b(checklist):
    return _evaluate(GATE_10B_KEYS, checklist, "gate 10B")


def evaluate_gate_10c(checklist):
    return _evaluate(GATE_10C_KEYS, checklist, "gate 10C")


def evaluate_gate_10d(checklist):
    return _evaluate(GATE_10D_KEYS, checklist, "gate 10D")


def evaluate_gate_10e(checklist):
    return _evaluate(GATE_10E_KEYS, checklist, "gate 10E")


@dataclass(frozen=True)
class Phase10Decision:
    gate_10a: GateResult
    gate_10b: GateResult
    gate_10c: GateResult
    gate_10d: GateResult
    gate_10e: GateResult
    surrogate_beats_baseline1: bool
    inverse_mostly_reachable: bool
    beam_level_acceptable: object   # bool or None (not evaluated)
    geometry_inversion_unlocked: bool
    final_status: str


def decide(gate_10a, gate_10b, gate_10c, gate_10d, gate_10e, surrogate_beats_baseline1,
           inverse_mostly_reachable, beam_level_acceptable):
    """Sec. 42: ALL FIVE gates must be OPEN for geometry inversion to unlock. Sec. 52's final
    status is chosen by walking the pipeline in order and reporting the FIRST point of failure,
    which is more informative than the generic `NO_VALIDATED_FORMULATION` catch-all (reserved for
    when the control-variable gate alone is what is missing, exactly Phase 9's own finding)."""
    all_open = all(g.status == "OPEN" for g in (gate_10a, gate_10b, gate_10c, gate_10d, gate_10e))
    geometry_inversion_unlocked = (
        all_open and surrogate_beats_baseline1 and inverse_mostly_reachable and bool(beam_level_acceptable)
    )

    if gate_10c.status == "OPEN" and not surrogate_beats_baseline1:
        final_status = "SURROGATE_INVALIDATED"
    elif gate_10d.status == "OPEN" and surrogate_beats_baseline1 and not inverse_mostly_reachable:
        final_status = "INVERSE_MODEL_INVALIDATED"
    elif (gate_10e.status == "OPEN" and surrogate_beats_baseline1 and inverse_mostly_reachable
          and beam_level_acceptable is False):
        final_status = "END_TO_END_MODEL_INVALIDATED"
    elif geometry_inversion_unlocked:
        final_status = "ENVIRONMENT_MODEL_VALIDATED"
    elif (surrogate_beats_baseline1 and inverse_mostly_reachable and beam_level_acceptable
          and gate_10b.status == "OPEN" and gate_10c.status == "OPEN"
          and gate_10d.status == "OPEN" and gate_10e.status == "OPEN"):
        # Everything about the environment-aware forward/inverse/beam chain checks out; only
        # Gate 10A (the persistent, Phase-9-inherited lack of an HFSS/measurement control
        # variable) blocks the final unlock.
        final_status = "ENVIRONMENT_MODEL_PARTIALLY_VALIDATED"
    else:
        final_status = "NO_VALIDATED_FORMULATION"

    return Phase10Decision(
        gate_10a, gate_10b, gate_10c, gate_10d, gate_10e, surrogate_beats_baseline1,
        inverse_mostly_reachable, beam_level_acceptable, geometry_inversion_unlocked, final_status,
    )
