"""Phase 8A Sec. 4 gate: `validate_holographic_phase_source` must not accept S11xx as
automatically valid merely because Phase 7 extracted it (brief Sec. 44 Test 8)."""

import pytest

from src.feasibility_phase_source import HOLOGRAPHIC_GENERATOR_IMPLEMENTED, validate_holographic_phase_source


def test_holographic_generator_is_not_implemented_in_this_project():
    # The gate's whole conclusion depends on this being an inspected fact, not an assumption;
    # this test pins it so a future change to the project is forced to update it deliberately.
    assert HOLOGRAPHIC_GENERATOR_IMPLEMENTED is False


def test_phase_source_is_reported_unvalidated_not_automatically_accepted():
    result = validate_holographic_phase_source()
    assert result.status == "PHASE_SOURCE_UNVALIDATED"
    assert result.status != "PHYSICALLY_VALIDATED"


def test_phase_source_quantity_required_is_unknown():
    result = validate_holographic_phase_source()
    assert result.quantity_required is None
    assert result.physically_equivalent is None


def test_phase_source_quantity_provided_matches_phase_extraction():
    from src.phase_extraction import PHASE_SOURCE_COMPONENT

    result = validate_holographic_phase_source()
    assert result.quantity_provided == PHASE_SOURCE_COMPONENT == "S11_xx"


def test_phase_source_time_convention_is_reported_consistent():
    result = validate_holographic_phase_source()
    assert result.time_convention_ok is True


def test_phase_source_note_mentions_unimplemented_holographic_synthesis():
    result = validate_holographic_phase_source()
    assert "no implemented holographic-synthesis formulation" in result.note


def test_reject_generator_implemented_without_updating_this_module():
    import src.feasibility_phase_source as mod

    original = mod.HOLOGRAPHIC_GENERATOR_IMPLEMENTED
    mod.HOLOGRAPHIC_GENERATOR_IMPLEMENTED = True
    try:
        with pytest.raises(NotImplementedError):
            mod.validate_holographic_phase_source()
    finally:
        mod.HOLOGRAPHIC_GENERATOR_IMPLEMENTED = original
