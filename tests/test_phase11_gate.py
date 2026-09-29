"""Sec. 30-32, 28: the 14-condition geometry-inversion gate, the 6-condition STEP-export gate, and
the mandatory negative tests (an unvalidated control variable, and an unvalidated heterogeneous
forward model, must each independently keep geometry inversion LOCKED)."""

import pytest

from src.phase11 import gate as gate_mod
from src.phase11.classification import CASE_A_LOCAL, CASE_C_STRONGLY_NONLOCAL
from src.phase11.reduction import REDUCTION_NOT_ATTEMPTED, REDUCTION_VALIDATED


def _all_true_geometry_conditions():
    return {k: True for k in gate_mod.GEOMETRY_INVERSION_CONDITION_KEYS}


def _all_true_step_conditions():
    return {k: True for k in gate_mod.STEP_EXPORT_CONDITION_KEYS}


def test_evaluate_geometry_inversion_raises_on_missing_key():
    with pytest.raises(ValueError):
        gate_mod.evaluate_geometry_inversion({k: True for k in gate_mod.GEOMETRY_INVERSION_CONDITION_KEYS[:-1]})


def test_geometry_inversion_unlocked_only_when_all_14_conditions_true():
    assert gate_mod.evaluate_geometry_inversion(_all_true_geometry_conditions()) == gate_mod.GEOMETRY_INVERSION_UNLOCKED
    conditions = _all_true_geometry_conditions()
    conditions["hfss_comparison_meets_tolerance"] = False
    assert gate_mod.evaluate_geometry_inversion(conditions) == gate_mod.GEOMETRY_INVERSION_LOCKED


def test_step_export_locked_if_geometry_inversion_is_locked_regardless_of_its_own_conditions():
    status = gate_mod.evaluate_step_export(gate_mod.GEOMETRY_INVERSION_LOCKED, _all_true_step_conditions())
    assert status == gate_mod.STEP_EXPORT_LOCKED


def test_step_export_unlocked_only_when_inversion_unlocked_and_all_conditions_true():
    status = gate_mod.evaluate_step_export(gate_mod.GEOMETRY_INVERSION_UNLOCKED, _all_true_step_conditions())
    assert status == gate_mod.STEP_EXPORT_UNLOCKED


# --------------------------------------------------------------------- Sec. 28 mandatory negatives

def test_negative_unvalidated_control_variable_keeps_geometry_inversion_locked():
    """Even with every other condition true, a missing HFSS/control variable must lock inversion."""
    conditions = _all_true_geometry_conditions()
    conditions["physical_control_variable_validated"] = False
    decision = gate_mod.decide(CASE_A_LOCAL, hfss_exists=False, inverse_mostly_found=True,
                                reduction_decision_status=REDUCTION_VALIDATED,
                                geometry_inversion_conditions=conditions,
                                step_export_conditions=_all_true_step_conditions())
    assert decision.geometry_inversion_status == gate_mod.GEOMETRY_INVERSION_LOCKED
    assert decision.step_export_status == gate_mod.STEP_EXPORT_LOCKED


def test_negative_unvalidated_heterogeneous_forward_model_keeps_geometry_inversion_locked():
    """A nonlocal/unconverged forward model (Sec. 18 Case B/C) must lock inversion even if
    everything else claims to pass."""
    conditions = _all_true_geometry_conditions()
    conditions["heterogeneous_forward_model_converges"] = False
    decision = gate_mod.decide(CASE_C_STRONGLY_NONLOCAL, hfss_exists=True, inverse_mostly_found=True,
                                reduction_decision_status=REDUCTION_NOT_ATTEMPTED,
                                geometry_inversion_conditions=conditions,
                                step_export_conditions=_all_true_step_conditions())
    assert decision.geometry_inversion_status == gate_mod.GEOMETRY_INVERSION_LOCKED
    assert decision.coupling_status == gate_mod.COUPLING_NONLOCAL
    assert decision.reduction_status == gate_mod.REDUCTION_INVALIDATED


def test_fully_validated_chain_unlocks_geometry_inversion_and_step_export():
    decision = gate_mod.decide(CASE_A_LOCAL, hfss_exists=True, inverse_mostly_found=True,
                                reduction_decision_status=REDUCTION_VALIDATED,
                                geometry_inversion_conditions=_all_true_geometry_conditions(),
                                step_export_conditions=_all_true_step_conditions())
    assert decision.geometry_inversion_status == gate_mod.GEOMETRY_INVERSION_UNLOCKED
    assert decision.step_export_status == gate_mod.STEP_EXPORT_UNLOCKED
    assert decision.coupling_status == gate_mod.COUPLING_LOCAL
    assert decision.neighborhood_status == gate_mod.FINITE_NEIGHBORHOOD_FOUND
