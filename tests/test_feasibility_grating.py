"""Phase 8A Sec. 28: grating-lobe / period check, cross-checked against
`discretization.grating_onset_frequency` (Phase 4/5, unmodified) on a real solved `UnitCell` so
this module's inlined formula is verified reuse, not an independent re-derivation."""

import warnings

import pytest

from src.constants import C0
from src.discretization import grating_onset_frequency
from src.feasibility_grating import grating_check
from src.library_builder import build_library_database
from src.library_schema import LibraryConfig

warnings.filterwarnings("ignore", category=UserWarning)

_SUBSTRATE = {"name": "Rogers RT/Duroid 5880", "eps_r": 2.2, "tan_delta": 0.0009, "thickness_m": 1.575e-3}


def test_grating_check_rejects_non_positive_inputs():
    with pytest.raises(ValueError):
        grating_check(0.0, 4.8e-3, 12.5e9)
    with pytest.raises(ValueError):
        grating_check(4.8e-3, 4.8e-3, 0.0)


def test_grating_onset_matches_discretization_for_the_locked_period():
    cfg = LibraryConfig(
        period_m=4.8e-3, substrate=_SUBSTRATE,
        f_min_hz=12.0e9, f_max_hz=13.0e9, f_step_hz=1.0e9, f_design_hz=12.5e9,
        l_values_mm=(1.0,), w_values_mm=(1.0,),
    )
    build = build_library_database(cfg, out_dir=None)
    unit_cell = build.database.unit_cell(0)
    expected_onset_hz = grating_onset_frequency(unit_cell)

    result = grating_check(4.8e-3, 4.8e-3, 12.5e9)
    assert result.grating_onset_hz == pytest.approx(expected_onset_hz, rel=1e-12)


def test_locked_library_period_has_no_grating_risk_in_10_to_15_ghz():
    result = grating_check(4.8e-3, 4.8e-3, 12.5e9)
    assert result.is_grating_risk is False
    assert result.grating_onset_hz == pytest.approx(C0 / 4.8e-3, rel=1e-12)
    assert result.grating_onset_hz / 1e9 > 15.0   # onset (~62.5 GHz) is far above the 10-15 GHz band


def test_grating_risk_is_true_above_onset():
    onset_hz = C0 / 30.0e-3
    result = grating_check(30.0e-3, 30.0e-3, onset_hz * 1.01)
    assert result.is_grating_risk is True


def test_period_over_wavelength_is_computed_correctly():
    result = grating_check(4.8e-3, 4.8e-3, 12.5e9)
    wavelength_m = C0 / 12.5e9
    assert result.period_over_wavelength == pytest.approx(4.8e-3 / wavelength_m, rel=1e-12)
