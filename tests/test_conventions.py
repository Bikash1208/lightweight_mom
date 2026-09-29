"""Complex-number / sign conventions: time factor exp(+j*omega*t)."""

import math

import numpy as np
import pytest

from src import constants
from src.excitation import PlaneWaveExcitation
from src.materials import Substrate
from src.periodic_bc import normal_wavenumber


# ---------------------------------------------------------------- constants

def test_speed_of_light_is_exact_si():
    assert constants.C0 == 299_792_458.0


def test_vacuum_constants_are_mutually_consistent():
    assert constants.EPS0 * constants.MU0 * constants.C0**2 == pytest.approx(1.0, rel=1e-15)
    assert constants.ETA0 == pytest.approx(math.sqrt(constants.MU0 / constants.EPS0), rel=1e-14)


def test_free_space_impedance_matches_codata_2018():
    # CODATA 2018: Z0 = 376.730 313 668(57) ohm
    assert constants.ETA0 == pytest.approx(376.730313668, rel=1e-9)
    assert constants.EPS0 == pytest.approx(8.8541878128e-12, rel=1e-10)


def test_time_convention_is_documented():
    assert constants.TIME_CONVENTION == "exp(+jwt)"


# ------------------------------------------------- physical meaning of signs

def test_forward_wave_phase_front_moves_toward_plus_z():
    """With exp(+jwt), the phasor exp(-j*k*z) is a wave travelling toward +z."""
    omega, k = 2 * math.pi * 10e9, 2 * math.pi * 10e9 / constants.C0
    t, dt, z = 0.3e-10, 1e-12, 0.01

    def field(zz, tt):
        return (np.exp(-1j * k * zz) * np.exp(1j * omega * tt)).real

    assert field(z, t) == pytest.approx(math.cos(omega * t - k * z), abs=1e-12)
    # follow a constant-phase point: after dt it has moved by (omega/k)*dt in +z
    assert field(z + (omega / k) * dt, t + dt) == pytest.approx(field(z, t), abs=1e-9)
    # ... and NOT by the same amount in -z
    assert abs(field(z - (omega / k) * dt, t + dt) - field(z, t)) > 1e-3


def test_lossy_wave_decays_along_propagation():
    """Im(k) < 0 => |exp(-j*k*z)| = exp(Im(k)*z) decays for z > 0."""
    k = Substrate("t", 2.2, 0.0009, 1.575e-3).wavenumber(10e9)
    z = np.linspace(0, 0.5, 50)
    amp = np.abs(np.exp(-1j * k * z))
    assert np.all(np.diff(amp) < 0.0)
    assert amp == pytest.approx(np.exp(k.imag * z), rel=1e-12)


def test_dielectric_loss_sign_gives_positive_dissipation():
    """Time-average dissipation density = (omega*eps0/2) * (-Im eps_r) * |E|^2 must be > 0."""
    s = Substrate("t", 2.2, 0.0009, 1.575e-3)
    omega = 2 * math.pi * 10e9
    p_diss = 0.5 * omega * constants.EPS0 * (-s.eps_complex.imag) * 1.0**2
    assert p_diss > 0.0
    # the opposite convention (exp(-iwt)) would have Im(eps) > 0; make sure we did not use it
    assert s.eps_complex.imag < 0.0


def test_evanescent_reflected_harmonic_decays_away_from_surface():
    """Reflected/scattered harmonic ~ exp(-j*kz*z) for z>0 must decay when evanescent."""
    k0 = 2 * math.pi * 10e9 / constants.C0
    kz = normal_wavenumber(k0, 3 * k0, 0.0)
    z = np.array([1e-4, 1e-3, 1e-2])
    amp = np.abs(np.exp(-1j * kz * z))
    assert np.all(np.diff(amp) < 0.0) and np.all(amp < 1.0)
    # transmitted (z<0 side) harmonic ~ exp(+j*kz*z) decays as z -> -inf
    assert np.all(np.abs(np.exp(1j * kz * -z)) < 1.0)


# ------------------------------------------------------ incident wave geometry

def test_normal_incidence_wavevector_points_toward_minus_z():
    e = PlaneWaveExcitation(frequency=10e9)
    kx, ky, kz = e.wavevector
    assert (kx, ky) == (0.0, 0.0)
    assert kz == pytest.approx(-e.k0, rel=1e-15)
    assert e.transverse_wavenumber == (0.0, 0.0)
    assert np.dot(e.e_field_direction, e.wavevector) == 0.0


def test_oblique_wavevector_definition():
    """theta from the normal, phi from +x: kx0 = k0 sin(th) cos(ph), ky0 = k0 sin(th) sin(ph)."""
    th, ph = math.radians(30), math.radians(90)
    e = PlaneWaveExcitation(10e9, theta=th, phi=ph, polarization="x")   # E = x_hat, k in y-z plane
    kx0, ky0 = e.transverse_wavenumber
    assert kx0 == pytest.approx(0.0, abs=1e-9)
    assert ky0 == pytest.approx(e.k0 * 0.5, rel=1e-14)
    assert np.linalg.norm(e.wavevector) == pytest.approx(e.k0, rel=1e-14)
    assert e.wavevector[2] < 0.0
    # kz0 from the dispersion relation agrees with the branch used elsewhere
    assert -normal_wavenumber(e.k0, kx0, ky0).real == pytest.approx(e.wavevector[2], rel=1e-12)


def test_polarization_must_be_transverse():
    with pytest.raises(ValueError, match="perpendicular"):
        PlaneWaveExcitation(10e9, theta=math.radians(30), phi=0.0, polarization="x")
    PlaneWaveExcitation(10e9, theta=math.radians(30), phi=0.0, polarization="y")  # E = y_hat OK


@pytest.mark.parametrize(
    "kw",
    [
        dict(frequency=0.0), dict(frequency=-1.0), dict(frequency=float("nan")),
        dict(frequency=10e9, theta=-0.1), dict(frequency=10e9, theta=math.pi / 2),
        dict(frequency=10e9, phi=float("inf")), dict(frequency=10e9, polarization="z"),
        dict(frequency=10e9, polarization="TE"),
    ],
)
def test_excitation_rejects_invalid(kw):
    with pytest.raises(ValueError):
        PlaneWaveExcitation(**kw)


def test_excitation_derived_quantities():
    e = PlaneWaveExcitation(frequency=10e9)
    assert e.omega == pytest.approx(2 * math.pi * 10e9, rel=1e-15)
    assert e.k0 * e.wavelength == pytest.approx(2 * math.pi, rel=1e-14)
