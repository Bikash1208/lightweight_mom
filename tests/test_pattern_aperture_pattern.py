"""Part 3: aperture far-field pattern prediction, checked against standard array-theory facts
(uniform aperture -> broadside; linear phase gradient -> steered beam at the expected angle)."""

import numpy as np
import pytest

from src.pattern.aperture_pattern import compute_pattern

FREQ = 10e9
C0 = 299792458.0
LAMBDA0 = C0 / FREQ
K0 = 2 * np.pi / LAMBDA0


def _square_aperture_grid(n=21, period_m=LAMBDA0 / 2):
    coords = (np.arange(n) - (n - 1) / 2) * period_m
    x, y = np.meshgrid(coords, coords, indexing="xy")
    return x, y, period_m


def test_uniform_aperture_points_broadside():
    x, y, d = _square_aperture_grid()
    field = np.ones_like(x, dtype=complex)
    result = compute_pattern(x, y, field, FREQ, d, d, n_theta=46, n_phi=37)
    assert result.main_beam_theta_deg < 5.0   # broadside, within grid resolution


def test_linear_phase_gradient_steers_the_beam():
    x, y, d = _square_aperture_grid()
    theta0_deg = 20.0
    kx_target = K0 * np.sin(np.radians(theta0_deg))
    field = np.exp(-1j * kx_target * x)   # standard progressive phase for a beam steered to theta0 along phi=0
    result = compute_pattern(x, y, field, FREQ, d, d, n_theta=46, n_phi=37)
    assert result.main_beam_theta_deg == pytest.approx(theta0_deg, abs=3.0)
    assert result.main_beam_phi_deg in (0.0, pytest.approx(360.0, abs=1e-6)) or result.main_beam_phi_deg < 5.0


def test_directivity_increases_with_larger_uniform_aperture():
    x1, y1, d = _square_aperture_grid(n=11)
    x2, y2, _ = _square_aperture_grid(n=31)
    r1 = compute_pattern(x1, y1, np.ones_like(x1, dtype=complex), FREQ, d, d, n_theta=46, n_phi=37)
    r2 = compute_pattern(x2, y2, np.ones_like(x2, dtype=complex), FREQ, d, d, n_theta=46, n_phi=37)
    assert r2.directivity_estimate_dbi > r1.directivity_estimate_dbi


def test_zero_field_raises():
    x, y, d = _square_aperture_grid(n=5)
    with pytest.raises(ValueError, match="zero power"):
        compute_pattern(x, y, np.zeros_like(x, dtype=complex), FREQ, d, d)
