"""Phase 8B Sec. 18, 40: object-wave polarization is a real vector direction, never silently
collapsed to a scalar phase."""

import numpy as np
import pytest

from src.phase8b.object_wave import Beam, single_beam_vector_field

F_HZ = 12.5e9
X, Y = np.array([0.0, 0.01]), np.array([0.0, -0.02])


def test_x_polarized_beam_has_zero_y_field():
    beam = Beam(theta_deg=15.0, phi_deg=30.0, polarization="x")
    ex, ey = single_beam_vector_field(X, Y, beam, F_HZ)
    assert np.all(ex != 0)
    assert np.all(ey == 0)


def test_y_polarized_beam_has_zero_x_field():
    beam = Beam(theta_deg=15.0, phi_deg=30.0, polarization="y")
    ex, ey = single_beam_vector_field(X, Y, beam, F_HZ)
    assert np.all(ey != 0)
    assert np.all(ex == 0)


def test_x_and_y_polarized_beams_share_the_same_scalar_phase_pattern():
    # Polarization must not change the *spatial* phase pattern, only which component carries it.
    beam_x = Beam(theta_deg=15.0, phi_deg=30.0, polarization="x")
    beam_y = Beam(theta_deg=15.0, phi_deg=30.0, polarization="y")
    ex, _ = single_beam_vector_field(X, Y, beam_x, F_HZ)
    _, ey = single_beam_vector_field(X, Y, beam_y, F_HZ)
    assert np.allclose(np.angle(ex), np.angle(ey))


def test_only_x_and_y_are_supported_no_silent_circular_collapse():
    with pytest.raises(ValueError):
        Beam(theta_deg=10.0, phi_deg=0.0, polarization="R")
