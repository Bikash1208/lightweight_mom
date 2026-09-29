"""Sec. 20, 22, 25, 29: the geometry-inversion gate and status taxonomy."""

from src.phase8e.decision_tree import STATUS_TAXONOMY, evaluate_formulation_gate, overall_project_status


def test_status_taxonomy_has_exactly_thirteen_values():
    assert len(STATUS_TAXONOMY) == 13
    assert "NO_VALIDATED_FORMULATION" in STATUS_TAXONOMY
    assert "PASS" not in STATUS_TAXONOMY


def test_architecture_incompatible_short_circuits_the_gate():
    d = evaluate_formulation_gate("A", "ARCHITECTURE_INCOMPATIBLE", "COMPATIBLE", True, True, True, True)
    assert d.final_status == "ARCHITECTURE_INCOMPATIBLE"
    assert d.geometry_inversion_permitted is False


def test_missing_family_blocks_before_reconstruction_is_even_attempted():
    d = evaluate_formulation_gate("B", "ARCHITECTURE_COMPATIBLE_UNVALIDATED", "COMPATIBLE", False, True, True, True)
    assert d.final_status == "REQUIRED_LIBRARY_MISSING"


def test_continuous_reconstruction_failure_is_reported_distinctly():
    d = evaluate_formulation_gate("B", "ARCHITECTURE_COMPATIBLE_UNVALIDATED", "COMPATIBLE", True, False, True, True)
    assert d.final_status == "CONTINUOUS_RECONSTRUCTION_FAILED"


def test_discrete_reconstruction_failure_is_reported_distinctly():
    d = evaluate_formulation_gate("B", "ARCHITECTURE_COMPATIBLE_UNVALIDATED", "COMPATIBLE", True, True, False, True)
    assert d.final_status == "DISCRETE_RECONSTRUCTION_FAILED"


def test_unvalidated_control_variable_blocks_geometry_inversion_even_if_everything_else_passes():
    d = evaluate_formulation_gate("B", "ARCHITECTURE_COMPATIBLE_UNVALIDATED", "COMPATIBLE", True, True, True, False)
    assert d.final_status == "GEOMETRY_INVERSION_BLOCKED"
    assert d.geometry_inversion_permitted is False


def test_fully_validated_chain_opens_the_gate():
    d = evaluate_formulation_gate("B", "ARCHITECTURE_COMPATIBLE_UNVALIDATED", "COMPATIBLE", True, True, True, True,
                                  formulation_validated=True, local_periodicity_acceptable=True)
    assert d.geometry_inversion_permitted is True
    assert d.final_status == "CONTROL_VARIABLE_VALIDATED"


def test_default_formulation_validated_and_local_periodicity_are_false_matching_real_project_state():
    # this project's own real values (phase8b: never FORMULATION_VALIDATED; local_periodicity.py:
    # always LOCAL_PERIODICITY_UNVALIDATED) -- the gate must not silently open without them.
    d = evaluate_formulation_gate("B", "ARCHITECTURE_COMPATIBLE_UNVALIDATED", "COMPATIBLE", True, True, True, True)
    assert d.geometry_inversion_permitted is False
    assert d.gate_checklist["formulation_validated"] is False
    assert d.gate_checklist["local_periodicity_acceptable"] is False


def test_overall_project_status_is_no_validated_formulation_unless_one_gate_is_fully_open():
    blocked = evaluate_formulation_gate("B", "ARCHITECTURE_COMPATIBLE_UNVALIDATED", "COMPATIBLE", True, True, True, False)
    assert overall_project_status([blocked]) == "NO_VALIDATED_FORMULATION"

    open_gate = evaluate_formulation_gate("B", "ARCHITECTURE_COMPATIBLE_UNVALIDATED", "COMPATIBLE", True, True, True, True,
                                           formulation_validated=True, local_periodicity_acceptable=True)
    assert overall_project_status([blocked, open_gate]) == "CONTROL_VARIABLE_VALIDATED"


def test_gate_never_forces_a_positive_result_by_default_arguments():
    # every checkbox missing/false by default reasoning: a formulation that fails architecture
    # never accidentally reports geometry_inversion_permitted=True.
    d = evaluate_formulation_gate("D", "ARCHITECTURE_INCOMPATIBLE", "UNKNOWN_REFERENCE_PLANE", False, False, False, False)
    assert d.geometry_inversion_permitted is False
