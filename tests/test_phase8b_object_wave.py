"""Phase 8B Sec. 37-41: numerical tests of the object-wave equations -- propagation direction,
azimuth, broadside, dimensional consistency, and the phase-sign convention under `exp(+j*omega*t)`."""

import math

import numpy as np
import pytest

from src.excitation import PlaneWaveExcitation
from src.phase8b.object_wave import Beam, combine_beams, single_beam_object_field, transverse_wavenumber_deg

F_HZ = 12.5e9


def test_transverse_wavenumber_matches_plane_wave_excitation():
    # Cross-check against the project's own existing, tested class wherever it IS constructible
    # (theta=0, or a polarization perpendicular to k_inc) -- see object_wave.py's module comment
    # for why general (theta, phi, polarization) triples cannot all go through that class.
    for theta_deg, phi_deg in [(0.0, 0.0), (30.0, 90.0), (30.0, 270.0), (45.0, 0.0)]:
        pol = "y" if phi_deg in (0.0, 180.0) or theta_deg == 0.0 else "x"
        exc = PlaneWaveExcitation(frequency=F_HZ, theta=math.radians(theta_deg), phi=math.radians(phi_deg), polarization=pol)
        expected = exc.transverse_wavenumber
        got = transverse_wavenumber_deg(F_HZ, theta_deg, phi_deg)
        assert got == pytest.approx(expected, abs=1e-9)


# ------------------------------------------------------------------------------- Test 37: broadside

def test_broadside_has_zero_transverse_phase_gradient():
    beam = Beam(theta_deg=0.0, phi_deg=0.0)
    x = np.linspace(-0.05, 0.05, 21)
    y = np.zeros_like(x)
    field = single_beam_object_field(x, y, beam, F_HZ)
    phase = np.degrees(np.angle(field))
    assert np.allclose(phase, phase[0], atol=1e-9)


# --------------------------------------------------------------------------------- Test 37: azimuth

@pytest.mark.parametrize("phi_deg,axis", [(0.0, "x"), (90.0, "y"), (180.0, "x"), (270.0, "y")])
def test_azimuth_orients_the_transverse_wavevector_correctly(phi_deg, axis):
    beam = Beam(theta_deg=30.0, phi_deg=phi_deg)
    kx0, ky0 = transverse_wavenumber_deg(F_HZ, 30.0, phi_deg)
    if axis == "x":
        assert abs(ky0) < 1e-6 and abs(kx0) > 1e-6
    else:
        assert abs(kx0) < 1e-6 and abs(ky0) > 1e-6


def test_opposite_azimuths_reverse_the_phase_gradient_sign():
    kx0_0, _ = transverse_wavenumber_deg(F_HZ, 30.0, 0.0)
    kx0_180, _ = transverse_wavenumber_deg(F_HZ, 30.0, 180.0)
    assert kx0_0 == pytest.approx(-kx0_180, rel=1e-9)


# -------------------------------------------------------------------------------- Test 39: phase sign

def test_phase_sign_matches_exp_minus_jkr_convention():
    # E(x,y) = exp(-j(kx0 x)); at x>0 with kx0>0 the phase should be NEGATIVE (a delay), matching
    # CONVENTIONS.md's "a delay is a negative phase" and the fixed exp(-j k.r) plane-wave rule.
    beam = Beam(theta_deg=30.0, phi_deg=0.0)   # kx0 > 0
    kx0, _ = transverse_wavenumber_deg(F_HZ, 30.0, 0.0)
    assert kx0 > 0.0
    x = np.array([0.0, 0.01])
    y = np.array([0.0, 0.0])
    field = single_beam_object_field(x, y, beam, F_HZ)
    phase = np.angle(field)
    assert phase[1] < phase[0]   # phase decreases (more negative / more delayed) with increasing x


def test_initial_phase_deg_enters_as_a_constant_additive_offset():
    beam0 = Beam(theta_deg=20.0, phi_deg=10.0, initial_phase_deg=0.0)
    beam1 = Beam(theta_deg=20.0, phi_deg=10.0, initial_phase_deg=45.0)
    x, y = np.array([0.01]), np.array([0.02])
    f0 = single_beam_object_field(x, y, beam0, F_HZ)
    f1 = single_beam_object_field(x, y, beam1, F_HZ)
    from src.phase import wrapped_difference_deg

    delta = wrapped_difference_deg(np.degrees(np.angle(f1)), np.degrees(np.angle(f0)))
    assert delta == pytest.approx(45.0, abs=1e-6)


# ------------------------------------------------------------------------------------- dimensional

def test_wavenumber_has_correct_units_and_magnitude():
    from src.constants import C0

    k0 = 2 * math.pi * F_HZ / C0
    kx0, ky0 = transverse_wavenumber_deg(F_HZ, 90.0, 0.0)   # grazing: |k_t| = k0
    assert math.hypot(kx0, ky0) == pytest.approx(k0, rel=1e-9)


def test_beam_rejects_out_of_range_angles():
    with pytest.raises(ValueError):
        Beam(theta_deg=120.0, phi_deg=0.0)
    with pytest.raises(ValueError):
        Beam(theta_deg=10.0, phi_deg=400.0)
    with pytest.raises(ValueError):
        Beam(theta_deg=10.0, phi_deg=0.0, polarization="circular")
