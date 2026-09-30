"""Phase 15 Part A: phi_ext,b (polarization compensation + OAM) and the polarization-compensated
multi-beam combination -- unit tests derived directly from `hologram.polarization_compensation`'s
own module docstring (the brief's eq., algebraically reduced per label)."""

import numpy as np
import pytest

from src.hologram.polarization_compensation import (
    combine_beams_with_polarization,
    polarization_compensation_term,
    polarization_jones_vector,
)
from src.phase8b.object_wave import Beam, combine_beams_scalar
from src.polarization import E_L, E_R, handedness

FREQ = 10e9


def _test_points():
    # four cardinal points about the origin feed: +x, -x, +y, -y (phi' = 0, 180, 90, -90 deg)
    x = np.array([1.0, -1.0, 0.0, 0.0])
    y = np.array([0.0, 0.0, 1.0, -1.0])
    return x, y


# ------------------------------------------------------------------------------------ natural / OAM


def test_natural_polarization_has_unit_compensation_everywhere():
    x, y = np.meshgrid(np.linspace(-0.02, 0.02, 9), np.linspace(-0.02, 0.02, 9))
    phi_ext = polarization_compensation_term(x, y, (0.0, 0.0), "natural")
    assert np.allclose(phi_ext, 1.0 + 0.0j)


def test_natural_plus_oam_charge_is_a_pure_spiral():
    x, y = _test_points()
    l = 2
    phi_ext = polarization_compensation_term(x, y, (0.0, 0.0), "natural", oam_charge=l)
    phi_prime = np.arctan2(y, x)
    assert np.allclose(phi_ext, np.exp(1j * l * phi_prime))


def test_oam_charge_zero_is_pencil_matches_beam_type_default():
    beam = Beam(theta_deg=10.0, phi_deg=0.0)
    assert beam.oam_charge == 0
    assert beam.desired_polarization == "natural"


# --------------------------------------------------------------------------------------- HP / VP


def test_horizontal_polarization_is_a_0_pi_flip_by_sign_of_cos_phi_prime():
    x, y = _test_points()
    phi_ext = polarization_compensation_term(x, y, (0.0, 0.0), "hp")
    phi_prime = np.arctan2(y, x)
    expected = np.sign(np.cos(phi_prime))
    assert np.allclose(phi_ext.real, expected, atol=1e-9)
    assert np.allclose(phi_ext.imag, 0.0, atol=1e-9)


def test_vertical_polarization_is_a_0_pi_flip_by_sign_of_sin_phi_prime():
    # off-axis points only: sin(phi') is exactly 0 at phi'=0/180deg (the two `_test_points()` on
    # the x-axis), where sign() and arg(0) legitimately disagree at that single measure-zero point
    # -- not a real physics disagreement, so this test avoids it rather than assert through it.
    x, y = np.array([1.0, -1.0]), np.array([1.0, -1.0])
    phi_ext = polarization_compensation_term(x, y, (0.0, 0.0), "vp")
    phi_prime = np.arctan2(y, x)
    expected = np.sign(np.sin(phi_prime))
    assert np.allclose(phi_ext.real, expected, atol=1e-9)
    assert np.allclose(phi_ext.imag, 0.0, atol=1e-9)


def test_slant_polarizations_flip_by_sign_of_cos_phi_prime_mp_45deg():
    x, y = _test_points()
    phi_prime = np.arctan2(y, x)
    phi_ext_p45 = polarization_compensation_term(x, y, (0.0, 0.0), "slant_p45")
    phi_ext_m45 = polarization_compensation_term(x, y, (0.0, 0.0), "slant_m45")
    assert np.allclose(phi_ext_p45.real, np.sign(np.cos(phi_prime - np.radians(45.0))), atol=1e-9)
    assert np.allclose(phi_ext_m45.real, np.sign(np.cos(phi_prime + np.radians(45.0))), atol=1e-9)


# ------------------------------------------------------------------------------------------- CP


def test_circular_polarizations_are_spirals_with_opposite_sign():
    x, y = _test_points()
    phi_prime = np.arctan2(y, x)
    rhcp = polarization_compensation_term(x, y, (0.0, 0.0), "rhcp")
    lhcp = polarization_compensation_term(x, y, (0.0, 0.0), "lhcp")
    # spirals of the brief's e^{-+j*phi'} form -- opposite winding for opposite handedness, and each
    # is a genuine spiral (not the HP/VP binary +-1 flip): compare against both candidate windings.
    assert np.allclose(rhcp, np.exp(-1j * phi_prime)) or np.allclose(rhcp, np.exp(1j * phi_prime))
    assert np.allclose(lhcp, np.conj(rhcp))


def test_circular_jones_vector_handedness_is_derived_not_hardcoded():
    """The picked p_hat for "rhcp" must actually BE the E_R/E_L that `polarization.handedness`
    calls RHCP for a wave travelling +z (the object-wave propagation direction) -- this is the
    project's own convention, not the paper's, and must not be hardcoded independently of it."""
    rhcp_p_hat = polarization_jones_vector("rhcp")
    lhcp_p_hat = polarization_jones_vector("lhcp")
    expected_rhcp = E_R if handedness("R", "+z") == "RHCP" else E_L
    expected_lhcp = E_L if expected_rhcp is E_R else E_R
    assert np.allclose(rhcp_p_hat, expected_rhcp)
    assert np.allclose(lhcp_p_hat, expected_lhcp)


def test_cp_phase_offset_deg_adds_a_constant_phase():
    x, y = _test_points()
    base = polarization_compensation_term(x, y, (0.0, 0.0), "rhcp")
    offset = polarization_compensation_term(x, y, (0.0, 0.0), "rhcp", cp_phase_offset_deg=45.0)
    assert np.allclose(offset, base * np.exp(1j * np.radians(45.0)))


# ---------------------------------------------------------------------------------- unit modulus


@pytest.mark.parametrize("label", ["natural", "hp", "vp", "slant_p45", "slant_m45", "rhcp", "lhcp"])
def test_compensation_term_is_always_unit_modulus(label):
    x, y = np.meshgrid(np.linspace(-0.02, 0.02, 11), np.linspace(-0.02, 0.02, 11))
    phi_ext = polarization_compensation_term(x, y, (0.001, -0.0005), label, oam_charge=1)
    assert np.allclose(np.abs(phi_ext), 1.0)


# -------------------------------------------------------------------------- feed position matters


def test_compensation_term_is_about_the_given_feed_not_the_origin():
    x, y = np.array([1.0]), np.array([0.0])
    at_origin = polarization_compensation_term(x, y, (0.0, 0.0), "hp")
    at_offset_feed = polarization_compensation_term(x, y, (2.0, 0.0), "hp")  # point is now "behind" the feed
    assert not np.allclose(at_origin, at_offset_feed)


# --------------------------------------------------------------- combine_beams_with_polarization


def test_natural_pencil_beams_reproduce_combine_beams_scalar_bit_identically():
    x = np.linspace(-0.02, 0.02, 9)
    y = np.linspace(-0.02, 0.02, 9)
    xg, yg = np.meshgrid(x, y, indexing="xy")
    beams = [
        Beam(theta_deg=20.0, phi_deg=0.0, weight=complex(1.2, -0.3)),
        Beam(theta_deg=10.0, phi_deg=180.0, weight=complex(0.5, 0.1)),
    ]
    old = combine_beams_scalar(xg, yg, beams, FREQ)
    new = combine_beams_with_polarization(xg, yg, beams, FREQ)
    assert np.array_equal(old, new)


def test_a_compensated_beam_changes_the_combined_field():
    x = np.linspace(-0.02, 0.02, 9)
    y = np.linspace(-0.02, 0.02, 9)
    xg, yg = np.meshgrid(x, y, indexing="xy")
    natural = [Beam(theta_deg=20.0, phi_deg=0.0, desired_polarization="natural")]
    hp = [Beam(theta_deg=20.0, phi_deg=0.0, desired_polarization="hp")]
    e_natural = combine_beams_with_polarization(xg, yg, natural, FREQ)
    e_hp = combine_beams_with_polarization(xg, yg, hp, FREQ)
    assert not np.allclose(e_natural, e_hp)


def test_sector_owner_index_zeroes_out_non_owning_beams():
    x = np.array([1.0, -1.0])
    y = np.array([0.0, 0.0])
    beams = [Beam(theta_deg=20.0, phi_deg=0.0), Beam(theta_deg=20.0, phi_deg=180.0)]
    owner = np.array([0, 1])  # cell 0 -> beam 0 only, cell 1 -> beam 1 only
    combined = combine_beams_with_polarization(x, y, beams, FREQ, sector_owner_index=owner)
    from src.phase8b.object_wave import single_beam_object_field

    expected0 = single_beam_object_field(x[:1], y[:1], beams[0], FREQ)
    expected1 = single_beam_object_field(x[1:], y[1:], beams[1], FREQ)
    assert np.allclose(combined[0], expected0[0])
    assert np.allclose(combined[1], expected1[0])


def test_per_beam_feed_xy_m_sequence_uses_each_beams_own_feed():
    x, y = np.array([1.0]), np.array([0.0])
    beams = [Beam(theta_deg=20.0, phi_deg=0.0, desired_polarization="hp")]
    at_origin = combine_beams_with_polarization(x, y, beams, FREQ, feed_xy_m=(0.0, 0.0))
    per_beam_feed = combine_beams_with_polarization(x, y, beams, FREQ, feed_xy_m=([2.0], [0.0]))
    assert not np.allclose(at_origin, per_beam_feed)
