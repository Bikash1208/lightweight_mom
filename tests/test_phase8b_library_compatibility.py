"""Phase 8B Sec. 26-28, 52-53: the decision-tree logic and final status selection."""

from src.phase8b.library_compatibility import (
    FINAL_STATES,
    evaluate_all,
    evaluate_formulation_a,
    evaluate_formulation_b,
    overall_project_status,
)


def test_formulation_a_is_library_control_variable_incompatible():
    e = evaluate_formulation_a()
    assert e.library_compatibility == "INCOMPATIBLE"
    assert e.final_status == "LIBRARY_CONTROL_VARIABLE_INCOMPATIBLE"
    assert e.final_status in FINAL_STATES


def test_formulation_a_decision_trail_is_fully_ordered_and_auditable():
    e = evaluate_formulation_a()
    questions = [q for q, _ in e.decision_trail]
    assert questions[0].startswith("Does a holographic formulation exist?")
    assert questions[-1] == "Can the Phase 7 library provide that quantity?" or "Result" in questions[-1]


def test_formulation_b_is_phase_source_unvalidated():
    e = evaluate_formulation_b()
    assert e.library_compatibility == "PARTIALLY_COMPATIBLE"
    assert e.final_status == "PHASE_SOURCE_UNVALIDATED"


def test_overall_status_reflects_the_more_actionable_finding():
    evaluations, status = evaluate_all()
    assert status == "PHASE_SOURCE_UNVALIDATED"
    # Formulation A's incompatibility must still be visible, not hidden by the overall status.
    a = next(e for e in evaluations if e.formulation_name.startswith("A:"))
    assert a.final_status == "LIBRARY_CONTROL_VARIABLE_INCOMPATIBLE"


def test_overall_status_never_returns_a_bare_pass():
    _, status = evaluate_all()
    assert status != "PASS"
    assert status in FINAL_STATES


def test_overall_status_all_incompatible_gives_incompatible():
    class _Fake:
        final_status = "LIBRARY_CONTROL_VARIABLE_INCOMPATIBLE"

    assert overall_project_status([_Fake(), _Fake()]) == "LIBRARY_CONTROL_VARIABLE_INCOMPATIBLE"


def test_overall_status_all_validated_gives_validated():
    class _Fake:
        final_status = "FORMULATION_VALIDATED"

    assert overall_project_status([_Fake(), _Fake()]) == "FORMULATION_VALIDATED"
