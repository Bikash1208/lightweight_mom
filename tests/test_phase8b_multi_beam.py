"""Phase 8B Sec. 15, 38: multi-beam field-level superposition, verified for 1/2/3 beams by
comparing the underlying complex fields, not only phase maps (Sec. 38's explicit requirement)."""

import numpy as np
import pytest

from src.phase8b.object_wave import Beam, combine_beams, single_beam_object_field, single_beam_vector_field

F_HZ = 12.5e9
X = np.array([0.0, 0.01, -0.02, 0.03])
Y = np.array([0.0, -0.01, 0.02, 0.01])


def test_one_beam_combination_equals_the_single_beam_field():
    beam = Beam(theta_deg=20.0, phi_deg=45.0)
    ex, ey = combine_beams(X, Y, [beam], F_HZ)
    expected_x, expected_y = single_beam_vector_field(X, Y, beam, F_HZ)
    assert np.allclose(ex, expected_x) and np.allclose(ey, expected_y)


def test_two_beams_combine_as_exact_complex_sum_not_phase_average():
    b1 = Beam(theta_deg=20.0, phi_deg=0.0, polarization="x")
    b2 = Beam(theta_deg=20.0, phi_deg=90.0, polarization="x")
    ex, ey = combine_beams(X, Y, [b1, b2], F_HZ)
    f1 = single_beam_object_field(X, Y, b1, F_HZ)
    f2 = single_beam_object_field(X, Y, b2, F_HZ)
    assert np.allclose(ex, f1 + f2)
    # The naive "average the two phases" alternative must NOT match (proves this isn't phase averaging).
    naive_phase_avg = np.angle(np.exp(1j * (np.angle(f1) + np.angle(f2)) / 2))
    assert not np.allclose(np.angle(ex), naive_phase_avg, atol=1e-6)


def test_three_beams_with_different_weights_is_a_weighted_complex_sum():
    beams = [
        Beam(theta_deg=15.0, phi_deg=0.0, weight=1.0),
        Beam(theta_deg=15.0, phi_deg=120.0, weight=0.5 + 0.2j),
        Beam(theta_deg=15.0, phi_deg=240.0, weight=0.3),
    ]
    ex, _ = combine_beams(X, Y, beams, F_HZ)
    expected = sum(single_beam_object_field(X, Y, b, F_HZ) for b in beams)
    assert np.allclose(ex, expected)


def test_combine_beams_rejects_empty_beam_list():
    with pytest.raises(ValueError):
        combine_beams(X, Y, [], F_HZ)


def test_combine_beams_scalar_requires_uniform_polarization():
    from src.phase8b.object_wave import combine_beams_scalar

    with pytest.raises(ValueError):
        combine_beams_scalar(
            X,
            Y,
            [Beam(theta_deg=10.0, phi_deg=0.0, polarization="x"), Beam(theta_deg=10.0, phi_deg=90.0, polarization="y")],
            F_HZ,
        )
