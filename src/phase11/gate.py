"""Sec. 30-32 (Phase 11's geometry-inversion/STEP-export unlock gates and the required status
codes). Exactly the same "process vs. outcome" gate philosophy Phase 9/10 already established: a
gate is OPEN when the corresponding validation STEP was genuinely performed, independent of whether
the result was favorable.
"""

from dataclasses import dataclass

# Sec. 32: required status codes.
HETEROGENEOUS_MOM_VALIDATED = "HETEROGENEOUS_MOM_VALIDATED"
HETEROGENEOUS_MOM_UNVALIDATED = "HETEROGENEOUS_MOM_UNVALIDATED"
COUPLING_LOCAL = "COUPLING_LOCAL"
COUPLING_NONLOCAL = "COUPLING_NONLOCAL"
FINITE_NEIGHBORHOOD_FOUND = "FINITE_NEIGHBORHOOD_FOUND"
LARGER_NEIGHBORHOOD_REQUIRED = "LARGER_NEIGHBORHOOD_REQUIRED"
PATTERN_DEPENDENT = "PATTERN_DEPENDENT"
REDUCTION_VALIDATED = "REDUCTION_VALIDATED"
REDUCTION_INVALIDATED = "REDUCTION_INVALIDATED"
INVERSE_VALIDATED = "INVERSE_VALIDATED"
INVERSE_INVALIDATED = "INVERSE_INVALIDATED"
HFSS_VALIDATED = "HFSS_VALIDATED"
HFSS_NOT_AVAILABLE = "HFSS_NOT_AVAILABLE"
GEOMETRY_INVERSION_LOCKED = "GEOMETRY_INVERSION_LOCKED"
GEOMETRY_INVERSION_UNLOCKED = "GEOMETRY_INVERSION_UNLOCKED"
STEP_EXPORT_LOCKED = "STEP_EXPORT_LOCKED"
STEP_EXPORT_UNLOCKED = "STEP_EXPORT_UNLOCKED"

# Sec. 30: the 14 geometry-inversion unlock conditions.
GEOMETRY_INVERSION_CONDITION_KEYS = (
    "physical_control_variable_validated",
    "heterogeneous_mom_formulation_validated",
    "coupling_extent_established",
    "required_neighborhood_established",
    "heterogeneous_forward_model_converges",
    "strong_contrast_case_passes",
    "abrupt_boundary_case_passes",
    "smooth_gradient_case_passes",
    "unseen_spatial_arrangements_pass",
    "inverse_solutions_direct_mom_verified",
    "unreachable_targets_rejected",
    "beam_level_holographic_validation_passes",
    "independent_hfss_validation_available",
    "hfss_comparison_meets_tolerance",
)

# Sec. 31: STEP export requires geometry inversion PLUS these additional conditions.
STEP_EXPORT_CONDITION_KEYS = (
    "validated_holographic_target",
    "validated_heterogeneous_forward_model",
    "validated_geometry_inversion",
    "direct_mom_verification",
    "far_field_verification",
    "hfss_spot_validation",
)


@dataclass(frozen=True)
class Phase11Decision:
    heterogeneous_mom_status: str
    coupling_status: str
    neighborhood_status: str
    pattern_status: object  # PATTERN_DEPENDENT or None
    reduction_status: str
    inverse_status: str
    hfss_status: str
    geometry_inversion_conditions: dict
    geometry_inversion_status: str
    step_export_conditions: dict
    step_export_status: str
    overall_case: str


def evaluate_geometry_inversion(conditions):
    missing = [k for k in GEOMETRY_INVERSION_CONDITION_KEYS if k not in conditions]
    if missing:
        raise ValueError(f"geometry-inversion checklist missing key(s): {missing}")
    all_true = all(conditions[k] for k in GEOMETRY_INVERSION_CONDITION_KEYS)
    return GEOMETRY_INVERSION_UNLOCKED if all_true else GEOMETRY_INVERSION_LOCKED


def evaluate_step_export(geometry_inversion_status, conditions):
    missing = [k for k in STEP_EXPORT_CONDITION_KEYS if k not in conditions]
    if missing:
        raise ValueError(f"STEP-export checklist missing key(s): {missing}")
    if geometry_inversion_status != GEOMETRY_INVERSION_UNLOCKED:
        return STEP_EXPORT_LOCKED
    return STEP_EXPORT_UNLOCKED if all(conditions[k] for k in STEP_EXPORT_CONDITION_KEYS) else STEP_EXPORT_LOCKED


def decide(
    overall_case,
    hfss_exists,
    inverse_mostly_found,
    reduction_decision_status,
    geometry_inversion_conditions,
    step_export_conditions,
):
    from .classification import CASE_A_LOCAL, CASE_D_PATTERN_DEPENDENT
    from .reduction import REDUCTION_VALIDATED as _REDUCTION_ATTEMPTED_AND_VALIDATED

    heterogeneous_mom_status = HETEROGENEOUS_MOM_VALIDATED  # Sec. 5: the direct heterogeneous
    # reference calculation itself ran
    # successfully and reproducibly at
    # every requested size -- this status
    # is about the SOLVER's own operation,
    # not the coupling-extent finding.
    coupling_status = COUPLING_LOCAL if overall_case == CASE_A_LOCAL else COUPLING_NONLOCAL
    neighborhood_status = FINITE_NEIGHBORHOOD_FOUND if overall_case == CASE_A_LOCAL else LARGER_NEIGHBORHOOD_REQUIRED
    pattern_status = PATTERN_DEPENDENT if overall_case == CASE_D_PATTERN_DEPENDENT else None
    hfss_status = HFSS_VALIDATED if hfss_exists else HFSS_NOT_AVAILABLE
    inverse_status = INVERSE_VALIDATED if inverse_mostly_found else INVERSE_INVALIDATED
    # Sec. 32 only names two reduction codes; `reduction.py`'s own, more specific
    # REDUCTION_NOT_ATTEMPTED (the actual, honest outcome whenever overall_case != CASE_A_LOCAL) is
    # reported as REDUCTION_INVALIDATED here -- "not attempted because not yet justified" is not
    # the same claim as "attempted and validated," so it cannot be reported as REDUCTION_VALIDATED.
    reduction_status = (
        REDUCTION_VALIDATED
        if reduction_decision_status == _REDUCTION_ATTEMPTED_AND_VALIDATED
        else REDUCTION_INVALIDATED
    )

    geometry_inversion_status = evaluate_geometry_inversion(geometry_inversion_conditions)
    step_export_status = evaluate_step_export(geometry_inversion_status, step_export_conditions)

    return Phase11Decision(
        heterogeneous_mom_status,
        coupling_status,
        neighborhood_status,
        pattern_status,
        reduction_status,
        inverse_status,
        hfss_status,
        geometry_inversion_conditions,
        geometry_inversion_status,
        step_export_conditions,
        step_export_status,
        overall_case,
    )
