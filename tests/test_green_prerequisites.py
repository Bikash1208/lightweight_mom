"""Consistency checks that any later Green's function must be built on.

The Green's function itself is NOT implemented in Phase 1. These tests check
only (a) the symmetry of the Floquet-mode kz grid and (b) that the kz branch
chosen in periodic_bc.normal_wavenumber yields a valid outgoing 1-D spectral
kernel

    g(z) = exp(-j*kz*|z|) / (2j*kz),      (d^2/dz^2 + kz^2) g = -delta(z),

which is the building block of every spectral-domain layered Green's function
under exp(+j*omega*t). The kernel is defined *inside this test file only*.
"""

import numpy as np
import pytest

from src.frequency import free_space_wavenumber
from src.materials import Substrate
from src.periodic_bc import RectangularLattice, normal_wavenumber

K0 = free_space_wavenumber(10e9)
K_SUB = Substrate("t", 2.2, 0.0009, 1.575e-3).wavenumber(10e9)


def g(kz, z):
    return np.exp(-1j * kz * np.abs(z)) / (2j * kz)


CASES = {
    "propagating": normal_wavenumber(K0, 0.4 * K0, 0.0),
    "evanescent": normal_wavenumber(K0, 2.5 * K0, 0.0),
    "lossy_propagating": normal_wavenumber(K_SUB, 0.4 * K0, 0.0),
    "lossy_evanescent": normal_wavenumber(K_SUB, 3.0 * K0, 1.0 * K0),
}


# ----------------------------------------------------- kz grid symmetries


@pytest.mark.parametrize("kx0, ky0", [(0.0, 0.0)])
def test_kz_grid_symmetric_for_square_lattice_at_normal_incidence(kx0, ky0):
    lat = RectangularLattice(8e-3, 8e-3)
    m, n = lat.mode_grid(6)
    kx, ky = lat.transverse_wavenumbers(kx0, ky0, m, n)
    kz = normal_wavenumber(K_SUB, kx, ky)
    assert np.allclose(kz, kz[::-1, :])  # m -> -m
    assert np.allclose(kz, kz[:, ::-1])  # n -> -n
    assert np.allclose(kz, kz.T)  # m <-> n (square lattice only)
    assert kz[6, 6] == pytest.approx(K_SUB, rel=1e-14)  # (0,0) mode


def test_kz_depends_only_on_transverse_magnitude():
    rng = np.random.default_rng(11)
    kt = rng.uniform(0, 4 * K0, 50)
    ang = rng.uniform(0, 2 * np.pi, 50)
    a = normal_wavenumber(K_SUB, kt * np.cos(ang), kt * np.sin(ang))
    b = normal_wavenumber(K_SUB, kt, 0.0)
    assert np.allclose(a, b, rtol=1e-12, atol=1e-9)


def test_only_finitely_many_propagating_orders_and_they_are_lowest():
    lat = RectangularLattice(8e-3, 8e-3)  # ~0.27*lambda0 at 10 GHz: only (0,0) radiates in air
    m, n = lat.mode_grid(8)
    kx, ky = lat.transverse_wavenumbers(0.0, 0.0, m, n)
    kz = normal_wavenumber(K0, kx, ky)
    propagating = kz.imag == 0.0
    assert propagating.sum() == 1 and propagating[8, 8]
    # every higher order is strictly evanescent (negative imaginary kz)
    assert np.all(kz[~propagating].imag < 0.0) and np.all(kz[~propagating].real == 0.0)


# ------------------------------------------------ 1-D spectral kernel checks


@pytest.mark.parametrize("name", CASES)
def test_kernel_is_even_in_z(name):
    kz = CASES[name]
    z = np.array([1e-4, 7e-4, 3e-3])
    assert np.allclose(g(kz, z), g(kz, -z), rtol=1e-14)


@pytest.mark.parametrize("name", CASES)
def test_kernel_reciprocity(name):
    kz = CASES[name]
    for z, zp in [(0.3e-3, -1.1e-3), (2e-3, 5e-4)]:
        assert g(kz, z - zp) == pytest.approx(g(kz, zp - z), rel=1e-14)


@pytest.mark.parametrize("name", CASES)
def test_kernel_source_jump_condition(name):
    """g'(0+) - g'(0-) = -1  (unit point source), evaluated by finite differences."""
    kz = CASES[name]
    h = 1e-7 / abs(kz)
    d_plus = (g(kz, h) - g(kz, 0.0)) / h
    d_minus = (g(kz, 0.0) - g(kz, -h)) / h
    assert (d_plus - d_minus) == pytest.approx(-1.0, rel=1e-5)


@pytest.mark.parametrize("name", CASES)
def test_kernel_satisfies_homogeneous_equation_away_from_source(name):
    kz = CASES[name]
    z0 = 0.5 / abs(kz)
    h = 1e-3 / abs(kz)
    d2 = (g(kz, z0 + h) - 2 * g(kz, z0) + g(kz, z0 - h)) / h**2
    assert d2 + kz**2 * g(kz, z0) == pytest.approx(0.0, abs=1e-5 * abs(kz**2 * g(kz, z0)))


@pytest.mark.parametrize("name", ["evanescent", "lossy_propagating", "lossy_evanescent"])
def test_kernel_decays_away_from_source(name):
    kz = CASES[name]
    z = np.linspace(0, 20 / abs(kz.imag), 40)
    assert np.all(np.diff(np.abs(g(kz, z))) < 0.0)


def test_kernel_lossless_propagating_has_constant_magnitude():
    kz = CASES["propagating"]
    z = np.linspace(0, 0.1, 30)
    assert np.allclose(np.abs(g(kz, z)), abs(g(kz, 0.0)), rtol=1e-13)


def test_wrong_branch_would_be_detected():
    """Sanity check on the tests: the opposite kz sign grows, so decay tests are meaningful."""
    kz = -CASES["evanescent"]
    z = np.linspace(0, 20 / abs(kz.imag), 40)
    assert np.all(np.diff(np.abs(g(kz, z))) > 0.0)
