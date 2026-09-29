"""Sec. 29-30, 33-34, 42: the two-gate decision logic and the geometry-inversion lock. Includes
the two mandatory negative tests from Sec. 39."""

import pytest

from src.phase9.gate import (
    FINAL_STATUS_OPTIONS,
    GATE_A_CHECKLIST_KEYS,
    GATE_B_CHECKLIST_KEYS,
    decide,
    evaluate_gate_a,
    evaluate_gate_b,
)


def _all_true(keys):
    return {k: True for k in keys}


def test_gate_a_open_only_when_every_checkbox_is_true():
    assert evaluate_gate_a(_all_true(GATE_A_CHECKLIST_KEYS)).status == "OPEN"
    partial = _all_true(GATE_A_CHECKLIST_KEYS)
    partial[GATE_A_CHECKLIST_KEYS[0]] = False
    assert evaluate_gate_a(partial).status == "CLOSED"


def test_gate_a_rejects_an_incomplete_checklist():
    with pytest.raises(ValueError):
        evaluate_gate_a({k: True for k in GATE_A_CHECKLIST_KEYS[:-1]})


def test_gate_b_open_only_when_every_checkbox_is_true():
    assert evaluate_gate_b(_all_true(GATE_B_CHECKLIST_KEYS)).status == "OPEN"
    partial = _all_true(GATE_B_CHECKLIST_KEYS)
    partial[GATE_B_CHECKLIST_KEYS[0]] = False
    assert evaluate_gate_b(partial).status == "CLOSED"


def test_status_taxonomy_has_exactly_four_values():
    assert len(FINAL_STATUS_OPTIONS) == 4
    assert "PASS" not in FINAL_STATUS_OPTIONS


def test_negative_unvalidated_s_parameters_cannot_unlock_geometry_inversion():
    """Sec. 39's first mandatory negative test."""
    gate_a = evaluate_gate_a(_all_true(GATE_A_CHECKLIST_KEYS))
    gate_b = evaluate_gate_b(_all_true(GATE_B_CHECKLIST_KEYS))
    decision = decide(gate_a, gate_b, "UNVALIDATED", "UNVALIDATED", "VALIDATED")
    assert decision.geometry_inversion_unlocked is False
    assert decision.selected_control is None
    assert decision.final_status == "NO_VALIDATED_FORMULATION"


def test_negative_unvalidated_local_periodicity_keeps_geometry_inversion_locked():
    """Sec. 39's second mandatory negative test."""
    gate_a = evaluate_gate_a(_all_true(GATE_A_CHECKLIST_KEYS))
    gate_b = evaluate_gate_b(_all_true(GATE_B_CHECKLIST_KEYS))
    decision = decide(gate_a, gate_b, "VALIDATED", "UNVALIDATED", "UNVALIDATED")
    assert decision.geometry_inversion_unlocked is False
    assert decision.final_status == "NO_VALIDATED_FORMULATION"


def test_invalidated_local_periodicity_produces_its_own_final_status():
    gate_a = evaluate_gate_a(_all_true(GATE_A_CHECKLIST_KEYS))
    gate_b = evaluate_gate_b(_all_true(GATE_B_CHECKLIST_KEYS))
    decision = decide(gate_a, gate_b, "VALIDATED", "UNVALIDATED", "INVALIDATED")
    assert decision.final_status == "LOCAL_PERIODICITY_INVALIDATED"
    assert decision.geometry_inversion_unlocked is False


def test_both_s11_and_s21_validated_does_not_auto_select_either_one():
    gate_a = evaluate_gate_a(_all_true(GATE_A_CHECKLIST_KEYS))
    gate_b = evaluate_gate_b(_all_true(GATE_B_CHECKLIST_KEYS))
    decision = decide(gate_a, gate_b, "VALIDATED", "VALIDATED", "VALIDATED")
    assert decision.selected_control is None   # Sec. 44: never select by convenience/smaller error


def test_fully_validated_s11_and_open_gates_and_validated_lpa_unlocks_geometry_inversion():
    gate_a = evaluate_gate_a(_all_true(GATE_A_CHECKLIST_KEYS))
    gate_b = evaluate_gate_b(_all_true(GATE_B_CHECKLIST_KEYS))
    decision = decide(gate_a, gate_b, "VALIDATED", "ARCHITECTURE_INCOMPATIBLE", "VALIDATED")
    assert decision.selected_control == "S11"
    assert decision.geometry_inversion_unlocked is True
    assert decision.final_status == "VALIDATED_CONTROL_AND_LPA"


def test_partially_validated_lpa_yields_partially_validated_final_status():
    gate_a = evaluate_gate_a(_all_true(GATE_A_CHECKLIST_KEYS))
    gate_b = evaluate_gate_b(_all_true(GATE_B_CHECKLIST_KEYS))
    decision = decide(gate_a, gate_b, "ARCHITECTURE_INCOMPATIBLE", "VALIDATED", "PARTIALLY_VALIDATED")
    assert decision.final_status == "PARTIALLY_VALIDATED"
    assert decision.geometry_inversion_unlocked is True


def test_closed_gate_blocks_geometry_inversion_even_with_a_validated_control_and_lpa():
    gate_a = evaluate_gate_a({**_all_true(GATE_A_CHECKLIST_KEYS), GATE_A_CHECKLIST_KEYS[0]: False})
    gate_b = evaluate_gate_b(_all_true(GATE_B_CHECKLIST_KEYS))
    decision = decide(gate_a, gate_b, "VALIDATED", "UNVALIDATED", "VALIDATED")
    assert gate_a.status == "CLOSED"
    assert decision.geometry_inversion_unlocked is False
