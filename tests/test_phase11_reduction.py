"""Sec. 19, 28: reduction is only claimed valid for CASE_A_LOCAL; otherwise honestly reported as
not attempted, never silently upgraded."""

from src.phase11.classification import (CASE_A_LOCAL, CASE_B_MODERATELY_NONLOCAL,
                                          CASE_C_STRONGLY_NONLOCAL)
from src.phase11.reduction import REDUCTION_NOT_ATTEMPTED, REDUCTION_VALIDATED, decide_reduction


def test_reduction_validated_for_case_a():
    decision = decide_reduction(CASE_A_LOCAL)
    assert decision.status == REDUCTION_VALIDATED


def test_reduction_not_attempted_for_nonlocal_cases():
    for case in (CASE_B_MODERATELY_NONLOCAL, CASE_C_STRONGLY_NONLOCAL):
        decision = decide_reduction(case)
        assert decision.status == REDUCTION_NOT_ATTEMPTED
        assert "none is implemented" in decision.rationale
