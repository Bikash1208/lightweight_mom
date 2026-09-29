import math

import numpy as np
import pytest

from src.constants import C0
from src.frequency import (
    angular_frequency,
    free_space_wavelength,
    free_space_wavenumber,
    validate_frequency,
)
from src.periodic_bc import normal_wavenumber

# ---------------------------------------------------------------- frequency


def test_known_wavelength_at_10ghz():
    assert free_space_wavelength(10e9) == pytest.approx(29.9792458e-3, rel=1e-15)


@pytest.mark.parametrize("f", [1e6, 2.45e9, 10e9, 94e9, 1e12])
def test_frequency_identities(f):
    k0 = free_space_wavenumber(f)
    assert k0 == pytest.approx(angular_frequency(f) / C0, rel=1e-15)
    assert k0 * free_space_wavelength(f) == pytest.approx(2 * math.pi, rel=1e-14)
    assert free_space_wavelength(f) * f == pytest.approx(C0, rel=1e-14)
    assert angular_frequency(f) == pytest.approx(2 * math.pi * f, rel=1e-15)


def test_frequency_scaling():
    assert free_space_wavenumber(20e9) == pytest.approx(2 * free_space_wavenumber(10e9), rel=1e-15)
    assert free_space_wavelength(20e9) == pytest.approx(free_space_wavelength(10e9) / 2, rel=1e-15)


def test_validate_frequency_accepts_int_and_numpy():
    assert validate_frequency(10_000_000_000) == 1e10
    assert validate_frequency(np.float64(5e9)) == 5e9


@pytest.mark.parametrize(
    "bad",
    [0.0, -1e9, float("nan"), float("inf"), -float("inf"), 1e9 + 1j, True, "10e9", None],
)
def test_validate_frequency_rejects(bad):
    with pytest.raises(ValueError):
        validate_frequency(bad)


# -------------------------------------------------------------- wave numbers

K0 = free_space_wavenumber(10e9)


@pytest.mark.parametrize(
    "k",
    [K0, 1.48 * K0, 1.48 * K0 * (1 - 0.0009j) ** 0.5, K0 * (1 - 0.3j)],
)
@pytest.mark.parametrize(
    "kx, ky", [(0.0, 0.0), (0.3 * K0, 0.4 * K0), (0.9 * K0, 0.0), (1.5 * K0, 0.7 * K0), (5 * K0, 5 * K0)]
)
def test_dispersion_relation_and_branch(k, kx, ky):
    kz = normal_wavenumber(k, kx, ky)
    assert kz * kz + kx**2 + ky**2 == pytest.approx(k * k, rel=1e-12, abs=1e-9)
    assert kz.imag <= 0.0
    assert kz.real >= 0.0


def test_normal_incidence_kz_equals_k():
    assert normal_wavenumber(K0, 0.0, 0.0) == pytest.approx(K0, rel=1e-15)
    k_lossy = 1.5 * K0 * (1 - 0.01j)
    assert normal_wavenumber(k_lossy, 0.0, 0.0) == pytest.approx(k_lossy, rel=1e-14)


def test_lossless_propagating_kz_is_real_positive():
    kz = normal_wavenumber(K0, 0.6 * K0, 0.0)
    assert kz.imag == 0.0
    assert kz.real == pytest.approx(0.8 * K0, rel=1e-14)


def test_lossless_evanescent_kz_is_negative_imaginary():
    kt = 2.0 * K0
    kz = normal_wavenumber(K0, kt, 0.0)
    assert kz.real == 0.0
    assert kz.imag == pytest.approx(-math.sqrt(kt**2 - K0**2), rel=1e-14)


def test_kz_vanishes_at_cutoff_and_is_continuous():
    assert normal_wavenumber(K0, K0, 0.0) == 0.0
    below = normal_wavenumber(K0, K0 * (1 - 1e-8), 0.0)
    above = normal_wavenumber(K0, K0 * (1 + 1e-8), 0.0)
    assert abs(below) < 1e-3 * K0 and abs(above) < 1e-3 * K0
    assert below.imag == 0.0 and below.real > 0.0
    assert above.real == 0.0 and above.imag < 0.0


def test_branch_is_negative_imag_for_dense_kt_sweep():
    kt = np.linspace(0.0, 10 * K0, 2001)
    for k in (K0, 2.2**0.5 * K0 * (1 - 0.0009j) ** 0.5):
        kz = normal_wavenumber(k, kt, 0.0)
        assert np.all(kz.imag <= 0.0)
        assert np.all(kz.real >= 0.0)


def test_normal_wavenumber_preserves_array_shape():
    kx = np.zeros((3, 4))
    assert normal_wavenumber(K0, kx, kx).shape == (3, 4)


def test_normal_wavenumber_rejects_active_medium():
    with pytest.raises(ValueError):
        normal_wavenumber(K0 * (1 + 0.01j), 0.0, 0.0)
