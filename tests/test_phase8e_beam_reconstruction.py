"""Sec. 18-19: continuous/discrete aperture -> far-field reconstruction."""

import numpy as np
import pytest

from src.phase8b.object_wave import Beam
from src.phase8e.beam_reconstruction import (
    compare_continuous_vs_discrete,
    detect_local_maxima,
    evaluate_beam_reconstruction,
    far_field_pattern,
)

FREQ_HZ = 12.5e9


def test_single_point_aperture_radiates_a_constant_pattern():
    x, y = np.array([[0.0]]), np.array([[0.0]])
    aperture = np.array([[1.0 + 0j]])
    theta, phi = np.meshgrid(np.linspace(0, 60, 7), np.linspace(0, 360, 9), indexing="ij")
    dx = dy = 1e-3
    field = far_field_pattern(x, y, aperture, FREQ_HZ, theta, phi, dx, dy)
    assert np.allclose(field, dx * dy, atol=1e-15)


def test_far_field_broadside_peak_for_a_uniform_in_phase_aperture():
    lin = np.linspace(-0.02, 0.02, 9)
    x, y = np.meshgrid(lin, lin, indexing="xy")
    aperture = np.ones_like(x, dtype=complex)
    theta = np.array([[0.0, 30.0]])
    phi = np.array([[0.0, 0.0]])
    dx = float(lin[1] - lin[0])
    field = far_field_pattern(x, y, aperture, FREQ_HZ, theta, phi, dx, dx)
    assert abs(field[0, 0]) > abs(field[0, 1])   # broadside stronger than a squinted direction


def test_detect_local_maxima_finds_a_single_known_peak():
    theta_grid, phi_grid = np.meshgrid(np.arange(5), np.arange(5), indexing="ij")
    power = np.zeros((5, 5))
    power[2, 3] = 1.0
    peaks = detect_local_maxima(power, theta_grid, phi_grid, min_relative_db=-40.0)
    assert any(p.theta_deg == 2 and p.phi_deg == 3 for p in peaks)


def test_evaluate_beam_reconstruction_zero_pointing_error_for_an_exact_peak():
    theta_grid, phi_grid = np.meshgrid(np.arange(0, 10, dtype=float), np.arange(0, 10, dtype=float), indexing="ij")
    power = np.full((10, 10), 1e-6)
    power[3, 4] = 1.0
    far_field = np.sqrt(power).astype(complex)
    beam = Beam(theta_deg=3.0, phi_deg=4.0)
    metrics = evaluate_beam_reconstruction(far_field, theta_grid, phi_grid, [beam])
    assert metrics.pointing_error_deg[0] == pytest.approx(0.0, abs=1e-6)
    assert metrics.matched_peak_power_db[0] == pytest.approx(0.0, abs=1e-6)


def test_evaluate_beam_reconstruction_reports_unwanted_peaks_separately():
    theta_grid, phi_grid = np.meshgrid(np.arange(0, 10, dtype=float), np.arange(0, 10, dtype=float), indexing="ij")
    power = np.full((10, 10), 1e-6)
    power[3, 4] = 1.0     # intended
    power[8, 8] = 0.5     # unwanted, dominant sidelobe
    far_field = np.sqrt(power).astype(complex)
    beam = Beam(theta_deg=3.0, phi_deg=4.0)
    metrics = evaluate_beam_reconstruction(far_field, theta_grid, phi_grid, [beam])
    assert metrics.highest_unwanted_relative_db is not None
    assert len(metrics.unwanted_peaks) >= 1


def test_compare_continuous_vs_discrete_is_zero_for_identical_fields():
    field = np.array([1 + 1j, 2 - 1j, 0.5j])
    result = compare_continuous_vs_discrete(field, field)
    assert result.complex_error_mean == pytest.approx(0.0, abs=1e-12)
    assert result.amplitude_error_mean == pytest.approx(0.0, abs=1e-12)
    assert result.phase_error_deg_mean == pytest.approx(0.0, abs=1e-9)
