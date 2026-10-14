"""Phase 26 B2.4 correction: Floquet balance of a modulated tensor sheet."""

import numpy as np
import pytest

from src.dispersion import tensor_sheet as ts
from src.dispersion.leakage_model import bound_beta, floquet_complex_beta, two_level_fourier
from src.dispersion.tensor_floquet import floquet_complex_beta_tensor, staircase_fourier_tensor
from src.frequency import free_space_wavenumber
from src.materials import Substrate

F = 10.54e9
SUB = Substrate("rt5880", 2.2, 0.0009, 1.575e-3)
K0 = free_space_wavenumber(F)


def test_fourier_of_a_two_cell_staircase_matches_the_scalar_two_level_series():
    y0, y1 = 3e-3j, 9e-3j
    tens = np.array([y0 * np.eye(2), y1 * np.eye(2)])
    got = staircase_fourier_tensor(tens, 3)[:, 0, 0]
    # the series differ only by the profile origin (a phase per harmonic): compare magnitudes
    assert np.allclose(np.abs(got), np.abs(two_level_fourier(y0, y1, 3, 0.5)), atol=1e-15)


def test_staircase_fourier_reconstructs_the_profile_and_conserves_the_mean():
    rng = np.random.default_rng(3)
    tens = rng.normal(size=(8, 2, 2)) * 1e-3j
    c = staircase_fourier_tensor(tens, 40)
    assert np.allclose(c[80], tens.mean(axis=0))
    x = (np.arange(8) + 0.5) / 8.0
    rec = sum(c[p + 80] * np.exp(-2j * np.pi * p * x[3]) for p in range(-80, 81))
    assert np.allclose(rec, tens[3], atol=0.05 * np.abs(tens).max())


def test_uniform_tensor_reproduces_the_bound_constant_of_the_tensor_dispersion():
    y = np.array([[9e-3j, 2e-3j], [2e-3j, 5e-3j]])
    ref, _ = ts.find_beta_tensor(y, 0.0, SUB, F)
    got, _ = floquet_complex_beta_tensor(np.array([y] * 4), 25.6e-3, SUB, F, ref, n_max=3, steps=1)
    assert got.real == pytest.approx(ref, rel=1e-6)
    assert abs(got.imag) < 1e-3 * ref


def test_scalar_two_level_staircase_equals_the_scalar_floquet_balance():
    y0, y1 = 5e-3j, 9e-3j
    period = 12.8e-3
    b0 = bound_beta(0.5 * (y0 + y1), SUB, F)
    scal, _ = floquet_complex_beta(two_level_fourier(y0, y1, 4, 0.5), period, 0.0, SUB, F, b0, n_max=4)
    tens = np.array([y0 * np.eye(2), y1 * np.eye(2)])
    got, _ = floquet_complex_beta_tensor(tens, period, SUB, F, b0, n_max=4)
    assert got == pytest.approx(scal, rel=1e-8)
