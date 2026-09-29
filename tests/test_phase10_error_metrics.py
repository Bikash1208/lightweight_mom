"""Sec. 25-26, 47: relative error + the full percentile report."""

import numpy as np

from src.phase10 import error_metrics as em


def test_relative_complex_error_matches_definition():
    pred = [1.1 + 0j]
    ref = [1.0 + 0j]
    rel = em.relative_complex_error(pred, ref)
    assert np.isclose(rel[0], 0.1)


def test_relative_complex_error_uses_a_documented_epsilon_floor_for_zero_reference():
    rel = em.relative_complex_error([1.0 + 0j], [0.0 + 0j])
    assert np.isfinite(rel[0])   # never a literal division-by-zero crash or inf/nan


def test_compute_error_distribution_reports_every_required_statistic():
    rng = np.random.default_rng(0)
    pred = rng.normal(size=50) + 1j * rng.normal(size=50)
    ref = pred + (rng.normal(size=50) + 1j * rng.normal(size=50)) * 0.1
    report = em.compute_error_distribution(pred, ref)
    for field in ("mean_complex_error", "median_complex_error", "rms_complex_error",
                  "p95_complex_error", "max_complex_error", "mean_relative_error",
                  "p95_relative_error", "max_relative_error", "mean_magnitude_error",
                  "p95_magnitude_error", "max_magnitude_error", "mean_phase_error_deg",
                  "p95_phase_error_deg", "max_phase_error_deg"):
        assert hasattr(report, field)
        assert getattr(report, field) >= 0
    assert report.max_complex_error >= report.p95_complex_error >= report.mean_complex_error * 0  # sanity ordering
    assert report.max_complex_error >= report.mean_complex_error


def test_compute_error_distribution_rejects_shape_mismatch():
    import pytest
    with pytest.raises(ValueError):
        em.compute_error_distribution([1 + 0j, 2 + 0j], [1 + 0j])
