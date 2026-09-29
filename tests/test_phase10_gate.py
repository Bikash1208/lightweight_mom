"""Sec. 41-42, 47, 52: the five process gates, the final unlock condition, and the two mandatory
Sec. 47 negative tests (unvalidated control variable -> locked; surrogate not blind-tested/failing
-> locked)."""

import pytest

from src.phase10 import gate as gate_mod


def _open_gate_a():
    return gate_mod.evaluate_gate_10a({k: True for k in gate_mod.GATE_10A_KEYS})


def _closed_gate_a():
    checklist = {k: True for k in gate_mod.GATE_10A_KEYS}
    checklist["hfss_or_measurement_reference_exists"] = False
    return gate_mod.evaluate_gate_10a(checklist)


def _open_gate(keys, evaluator):
    return evaluator({k: True for k in keys})


def test_gate_raises_on_missing_checklist_key():
    with pytest.raises(ValueError):
        gate_mod.evaluate_gate_10a({k: True for k in gate_mod.GATE_10A_KEYS[:-1]})


def test_gate_closed_if_any_checklist_item_is_false():
    checklist = {k: True for k in gate_mod.GATE_10B_KEYS}
    checklist[gate_mod.GATE_10B_KEYS[0]] = False
    result = gate_mod.evaluate_gate_10b(checklist)
    assert result.status == "CLOSED"


# --------------------------------------------------------------- Sec. 47 mandatory negative tests

def test_negative_unvalidated_control_variable_keeps_geometry_inversion_locked():
    """Even with every OTHER gate open and every outcome favorable, an unvalidated (missing HFSS)
    control variable (Gate 10A CLOSED) must keep geometry inversion LOCKED."""
    gate_a = _closed_gate_a()
    gate_b = _open_gate(gate_mod.GATE_10B_KEYS, gate_mod.evaluate_gate_10b)
    gate_c = _open_gate(gate_mod.GATE_10C_KEYS, gate_mod.evaluate_gate_10c)
    gate_d = _open_gate(gate_mod.GATE_10D_KEYS, gate_mod.evaluate_gate_10d)
    gate_e = _open_gate(gate_mod.GATE_10E_KEYS, gate_mod.evaluate_gate_10e)

    decision = gate_mod.decide(gate_a, gate_b, gate_c, gate_d, gate_e,
                                surrogate_beats_baseline1=True, inverse_mostly_reachable=True,
                                beam_level_acceptable=True)
    assert decision.geometry_inversion_unlocked is False
    assert decision.final_status != "ENVIRONMENT_MODEL_VALIDATED"


def test_negative_surrogate_that_fails_blind_test_keeps_geometry_inversion_locked():
    """Even with Gate 10A open (a hypothetical validated control variable), a surrogate that does
    NOT beat the naive isolated-cell baseline on the blind test must keep geometry inversion
    LOCKED and be reported as SURROGATE_INVALIDATED."""
    gate_a = _open_gate_a()
    gate_b = _open_gate(gate_mod.GATE_10B_KEYS, gate_mod.evaluate_gate_10b)
    gate_c = _open_gate(gate_mod.GATE_10C_KEYS, gate_mod.evaluate_gate_10c)
    gate_d = _open_gate(gate_mod.GATE_10D_KEYS, gate_mod.evaluate_gate_10d)
    gate_e = _open_gate(gate_mod.GATE_10E_KEYS, gate_mod.evaluate_gate_10e)

    decision = gate_mod.decide(gate_a, gate_b, gate_c, gate_d, gate_e,
                                surrogate_beats_baseline1=False, inverse_mostly_reachable=True,
                                beam_level_acceptable=True)
    assert decision.geometry_inversion_unlocked is False
    assert decision.final_status == "SURROGATE_INVALIDATED"


def test_fully_validated_chain_unlocks_geometry_inversion():
    gate_a = _open_gate_a()
    gate_b = _open_gate(gate_mod.GATE_10B_KEYS, gate_mod.evaluate_gate_10b)
    gate_c = _open_gate(gate_mod.GATE_10C_KEYS, gate_mod.evaluate_gate_10c)
    gate_d = _open_gate(gate_mod.GATE_10D_KEYS, gate_mod.evaluate_gate_10d)
    gate_e = _open_gate(gate_mod.GATE_10E_KEYS, gate_mod.evaluate_gate_10e)

    decision = gate_mod.decide(gate_a, gate_b, gate_c, gate_d, gate_e,
                                surrogate_beats_baseline1=True, inverse_mostly_reachable=True,
                                beam_level_acceptable=True)
    assert decision.geometry_inversion_unlocked is True
    assert decision.final_status == "ENVIRONMENT_MODEL_VALIDATED"


def test_real_project_state_reports_no_validated_formulation_or_more_specific_failure():
    """The real Phase 10 production run has Gate 10A CLOSED (no HFSS) -- confirm the decision
    function never reports ENVIRONMENT_MODEL_VALIDATED in that case."""
    gate_a = _closed_gate_a()
    gate_b = _open_gate(gate_mod.GATE_10B_KEYS, gate_mod.evaluate_gate_10b)
    gate_c = _open_gate(gate_mod.GATE_10C_KEYS, gate_mod.evaluate_gate_10c)
    gate_d = _open_gate(gate_mod.GATE_10D_KEYS, gate_mod.evaluate_gate_10d)
    gate_e = _open_gate(gate_mod.GATE_10E_KEYS, gate_mod.evaluate_gate_10e)
    decision = gate_mod.decide(gate_a, gate_b, gate_c, gate_d, gate_e,
                                surrogate_beats_baseline1=False, inverse_mostly_reachable=True,
                                beam_level_acceptable=True)
    assert decision.final_status in gate_mod.FINAL_STATUS_OPTIONS
    assert decision.final_status != "ENVIRONMENT_MODEL_VALIDATED"
    assert decision.geometry_inversion_unlocked is False
