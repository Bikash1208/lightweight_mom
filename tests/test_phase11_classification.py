"""Sec. 18, 28: CASE A/B/C/D classification from measured convergence results."""

from src.phase11.classification import (
    CASE_A_LOCAL,
    CASE_B_MODERATELY_NONLOCAL,
    CASE_C_STRONGLY_NONLOCAL,
    CASE_D_PATTERN_DEPENDENT,
    classify_overall,
)
from src.phase11.convergence import ConvergenceResult


def test_case_a_when_converged_at_small_n():
    cr = ConvergenceResult("C1", "CONVERGED", 1e-9, 1e-9, 3)
    result = classify_overall([cr], pattern_dependence_detected=False)
    assert result.per_center_case["C1"] == CASE_A_LOCAL
    assert result.overall_case == CASE_A_LOCAL


def test_case_b_when_converged_only_at_large_n():
    cr = ConvergenceResult("C1", "CONVERGED", 1e-9, 1e-9, 9)
    result = classify_overall([cr], pattern_dependence_detected=False)
    assert result.per_center_case["C1"] == CASE_B_MODERATELY_NONLOCAL


def test_case_c_when_still_substantially_changing():
    cr = ConvergenceResult("C1", "STILL_SUBSTANTIALLY_CHANGING", 0.6, 0.9, None)
    result = classify_overall([cr], pattern_dependence_detected=False)
    assert result.per_center_case["C1"] == CASE_C_STRONGLY_NONLOCAL
    assert result.overall_case == CASE_C_STRONGLY_NONLOCAL


def test_overall_case_is_the_worst_across_centers():
    good = ConvergenceResult("C1", "CONVERGED", 1e-9, 1e-9, 3)
    bad = ConvergenceResult("C2", "STILL_SUBSTANTIALLY_CHANGING", 0.6, 0.9, None)
    result = classify_overall([good, bad], pattern_dependence_detected=False)
    assert result.overall_case == CASE_C_STRONGLY_NONLOCAL


def test_pattern_dependence_overrides_to_case_d():
    good = ConvergenceResult("C1", "CONVERGED", 1e-9, 1e-9, 3)
    result = classify_overall([good], pattern_dependence_detected=True)
    assert result.overall_case == CASE_D_PATTERN_DEPENDENT
