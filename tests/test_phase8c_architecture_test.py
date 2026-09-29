"""Phase 8C Sec. 23-28: reflection/transmission architecture tests, Phase 8B re-evaluation."""

from src.phase8c.architecture_test import (
    formulation_a_reevaluation,
    formulation_b_reevaluation,
    reflection_architecture_test,
    transmission_architecture_test,
)


def test_reflection_architecture_is_unvalidated_matching_s11_investigation():
    result = reflection_architecture_test()
    assert result.status == "REFLECTION_ARCHITECTURE_UNVALIDATED"
    assert result.chain[-2][0] == "S11xx interpretation as that quantity"
    assert result.chain[-2][1] == "INSUFFICIENT_EVIDENCE"


def test_transmission_architecture_is_unvalidated_matching_s21_investigation():
    result = transmission_architecture_test()
    assert result.status == "TRANSMISSION_ARCHITECTURE_UNVALIDATED"


def test_formulation_a_reevaluation_is_unchanged_from_phase8b():
    result = formulation_a_reevaluation()
    assert result["status"] == "LIBRARY_CONTROL_VARIABLE_INCOMPATIBLE"
    assert "grounded" in result["ground_plane_requirement"]
    assert "monopole" in result["reference_wave_requirement"]


def test_formulation_b_reevaluation_does_not_validate_by_resemblance_alone():
    result = formulation_b_reevaluation()
    assert result["corresponds_to_reflection_phase"] == "REFLECTION_ARCHITECTURE_UNVALIDATED"
    assert result["corresponds_to_transmission_phase"] == "TRANSMISSION_ARCHITECTURE_UNVALIDATED"
    assert "neither correspondence is validated" in result["note"]
