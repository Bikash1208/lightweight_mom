"""Phase 8A schema: `FeasibilityConfig` validation and the fixed states/exit-code tables."""

import pytest

from src.feasibility_schema import (
    DEFAULT_LIBRARY_DIR,
    EXIT_CODES,
    STATES,
    FeasibilityConfig,
)


def test_default_library_dir_is_the_phase7_production_library():
    assert DEFAULT_LIBRARY_DIR.name == "unit_cell_library"
    assert DEFAULT_LIBRARY_DIR.parent.name == "results"


def test_every_state_has_an_exit_code():
    # Sec. 47's own list omits three of Sec. 35's states; feasibility_schema.py documents the
    # extension rather than silently reusing 0-8 for them.
    for state in STATES:
        assert state in EXIT_CODES, state


def test_exit_codes_are_unique():
    assert len(set(EXIT_CODES.values())) == len(EXIT_CODES)


def test_design_feasible_is_exit_code_zero():
    assert EXIT_CODES["DESIGN_FEASIBLE"] == 0


def test_config_rejects_unsupported_polarization():
    with pytest.raises(ValueError):
        FeasibilityConfig(polarization="circular")


def test_config_rejects_non_positive_tolerance():
    with pytest.raises(ValueError):
        FeasibilityConfig(phase_match_tolerance_deg=0.0)
    with pytest.raises(ValueError):
        FeasibilityConfig(phase_match_tolerance_deg=-1.0)


def test_config_default_polarization_is_x():
    assert FeasibilityConfig().polarization == "x"
