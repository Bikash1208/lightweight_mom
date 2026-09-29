"""Phase 8C Sec. 29, 55: reference-plane sensitivity, demonstrated only, never applied to stored
data (no arbitrary phase offset introduced anywhere else in this package)."""

import numpy as np
import pytest

from src.phase8c.reference_plane_sensitivity import demonstrate_sensitivity


def test_reference_plane_shift_changes_phase_by_the_documented_formula():
    s11 = 0.3 + 0.1j
    s21 = 0.9 + 0.05j
    result = demonstrate_sensitivity(s11, s21, frequency_hz=12.5e9, shift_m=1.0e-3)
    assert result.delta_deg_S11 != pytest.approx(0.0, abs=1e-6)
    assert result.delta_deg_S21 != pytest.approx(0.0, abs=1e-6)


def test_zero_shift_gives_zero_phase_change():
    s11 = 0.3 + 0.1j
    s21 = 0.9 + 0.05j
    result = demonstrate_sensitivity(s11, s21, frequency_hz=12.5e9, shift_m=0.0)
    assert result.delta_deg_S11 == pytest.approx(0.0, abs=1e-9)
    assert result.delta_deg_S21 == pytest.approx(0.0, abs=1e-9)


def test_sensitivity_magnitude_is_unchanged_by_a_pure_phase_shift():
    s11 = 0.3 + 0.1j
    result = demonstrate_sensitivity(s11, 0.9 + 0.0j, frequency_hz=12.5e9, shift_m=2.0e-3)
    # A pure reference-plane shift must not alter magnitude (only phase).
    from src.frequency import free_space_wavenumber

    k0 = free_space_wavenumber(12.5e9)
    s11_after = s11 * np.exp(-2j * k0 * 2.0e-3)
    assert abs(s11_after) == pytest.approx(abs(s11))
