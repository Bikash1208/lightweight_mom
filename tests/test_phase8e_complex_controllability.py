"""Sec. 13-14: complex-plane controllability, on real S11xx/S21xx library data."""

import numpy as np
import pytest

from src.phase8c import data_source
from src.phase8c.scattering_analysis import extract_component_table
from src.phase8e.complex_controllability import (
    complex_plane_controllability,
    library_complex_plane_summary,
    phase_only_controllability,
)
from src.phase_extraction import design_frequency_index

pytestmark = pytest.mark.skipif(not data_source.DEFAULT_CHECKPOINT_JSON.exists(), reason="Phase 7 checkpoint not present")


@pytest.fixture
def s11_rows(phase8c_full_database):
    idx = design_frequency_index(phase8c_full_database, 12.5e9)
    return extract_component_table(phase8c_full_database, "S11_xx", idx)


def test_library_complex_plane_summary_matches_324_geometries(s11_rows):
    summary = library_complex_plane_summary(s11_rows, "S11_xx")
    assert summary.n_geometries == 324
    assert summary.amplitude_min < summary.amplitude_max


def test_complex_plane_controllability_zero_error_when_target_is_a_real_library_value(s11_rows):
    target = complex(s11_rows[0].real, s11_rows[0].imag)
    result = complex_plane_controllability(target, s11_rows)
    assert result.mean_complex_error == pytest.approx(0.0, abs=1e-9)
    assert result.mean_amplitude_error == pytest.approx(0.0, abs=1e-9)
    assert result.mean_phase_error_deg == pytest.approx(0.0, abs=1e-6)


def test_complex_plane_controllability_reports_separate_amplitude_and_phase_error(s11_rows):
    # a target on the unit circle (magnitude 1) is outside the library's achievable magnitude
    # range (Sec. 13's own measured [0.231, 0.831]), so amplitude error must be nonzero.
    target = 1.0 + 0j
    result = complex_plane_controllability(target, s11_rows)
    assert result.mean_amplitude_error > 0.0


def test_phase_only_controllability_reports_uncontrolled_amplitude_range(s11_rows):
    result = phase_only_controllability(np.array([0.0, 90.0, 180.0]), s11_rows)
    assert result.n_required_points == 3
    assert result.uncontrolled_amplitude_min < result.uncontrolled_amplitude_max
    assert result.mean_phase_error_deg >= 0.0
