"""Sec. 3, 9-12: complex/magnitude/phase error metrics."""

import math

import pytest

from src.phase9.error_metrics import complex_error, compute_error_stats, magnitude_error, phase_error_deg


def test_complex_error_is_plain_distance():
    assert complex_error(0 + 0j, 3 + 4j) == pytest.approx(5.0)


def test_magnitude_error_never_normalizes_amplitude():
    # Rule 4's own example: |S|=0.62 vs |S|=0.98
    assert magnitude_error(0.98 + 0j, 0.62 + 0j) == pytest.approx(0.36)


def test_phase_error_deg_uses_wrapped_difference_not_naive_subtraction():
    # 179 deg vs -179 deg: naive subtraction gives 358, wrapped gives 2.
    a = complex(math.cos(math.radians(179)), math.sin(math.radians(179)))
    b = complex(math.cos(math.radians(-179)), math.sin(math.radians(-179)))
    assert phase_error_deg(a, b) == pytest.approx(2.0, abs=1e-6)


def test_compute_error_stats_zero_for_identical_arrays():
    values = [1 + 1j, 2 - 1j, 0.5j]
    stats = compute_error_stats(values, values)
    assert stats.mean_complex_error == pytest.approx(0.0, abs=1e-12)
    assert stats.max_abs_phase_error_deg == pytest.approx(0.0, abs=1e-9)


def test_compute_error_stats_rejects_shape_mismatch():
    with pytest.raises(ValueError):
        compute_error_stats([1 + 1j], [1 + 1j, 2 + 2j])


def test_compute_error_stats_reports_rms_and_max_separately_from_mean():
    a = [1.0 + 0j, 1.0 + 0j, 1.0 + 0j]
    b = [1.1 + 0j, 1.0 + 0j, 2.0 + 0j]
    stats = compute_error_stats(a, b)
    assert stats.max_complex_error > stats.mean_complex_error
    assert stats.rms_complex_error != stats.mean_complex_error
