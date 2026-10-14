"""Phase 26 A: the leakage constant of a modulated patch sheet -- closed form against the exact Floquet
harmonic balance, scaling, and the surface-wave power identity."""

import numpy as np
import pytest

from src.constants import ETA0
from src.dispersion import leakage_model as lm
from src.dispersion.analytic_tm0 import tm0_beta
from src.dispersion.sheet_model import _y_air, _y_slab
from src.frequency import angular_frequency, free_space_wavenumber
from src.materials import Substrate

F = 10.54e9
SUB = Substrate("rt5880", 2.2, 0.0, 1.575e-3)
K0 = free_space_wavenumber(F)
Y_MEAN = 0.5j * 2.9 / ETA0  # a mean patch sheet of the Phase 22 pair's order of magnitude


def _hb_alpha(amp_rel, gamma=0.0, n_max=6):
    bm = lm.bound_beta(Y_MEAN, SUB, F)
    period = 2.0 * np.pi / (bm + 0.5 * K0)
    coeffs = lm.sinusoid_fourier(Y_MEAN, 1j * amp_rel * 2.9 / ETA0, n_max)
    beta, _ = lm.floquet_complex_beta(coeffs, period, gamma, SUB, F, bm, n_max=n_max)
    return bm, period, beta


def _closed_alpha(amp_rel, bm, period, gamma=0.0):
    th, ze, _kt = lm.radiated_harmonic_geometry(bm, 2.0 * np.pi / period, gamma, F)
    return lm.leakage_alpha_first_order(bm, Y_MEAN, 0.5 * amp_rel * 2.9 / ETA0, th, ze, SUB, F)


def test_bare_sheet_recovers_analytic_tm0():
    beta = lm.bound_beta(0.0, SUB, F)
    assert beta == pytest.approx(tm0_beta(F, SUB.eps_r, SUB.thickness), rel=1e-6)


def test_surface_wave_power_equals_susceptance_derivative_identity():
    bm = lm.bound_beta(Y_MEAN, SUB, F)
    omega = angular_frequency(F)

    def b_total(b):
        return (_y_air(b, omega, K0) + _y_slab(b, omega, SUB.eps_complex, K0, SUB.thickness)).imag

    h = 1e-4 * K0
    db = (b_total(bm + h) - b_total(bm - h)) / (2 * h)
    e0 = ETA0 * np.sqrt(bm**2 - K0**2) / K0
    assert lm.surface_wave_power_per_h2(bm, SUB.eps_r, SUB.thickness, F) == pytest.approx(abs(db) * e0**2 / 4, rel=2e-3)


def test_closed_form_matches_floquet_balance_weak_modulation():
    bm, period, beta = _hb_alpha(0.05)
    a_hb = -beta.imag
    assert a_hb > 0.0
    assert _closed_alpha(0.05, bm, period) == pytest.approx(a_hb, rel=0.02)
    assert beta.real == pytest.approx(bm, rel=1e-4)


def test_leakage_scales_with_modulation_depth_squared():
    _, _, b1 = _hb_alpha(0.05)
    _, _, b2 = _hb_alpha(0.10)
    assert (-b2.imag) / (-b1.imag) == pytest.approx(4.0, rel=0.03)


def test_oblique_propagation_checks_the_tm_te_split():
    # propagation 40 deg off the modulation vector: the radiated harmonic is a TM/TE mix (zeta != 0)
    bm, period, beta = _hb_alpha(0.05, gamma=np.radians(40.0))
    th, ze, _ = lm.radiated_harmonic_geometry(bm, 2.0 * np.pi / period, np.radians(40.0), F)
    assert np.isfinite(th) and abs(ze) > 0.1
    assert _closed_alpha(0.05, bm, period, np.radians(40.0)) == pytest.approx(-beta.imag, rel=0.04)


def test_no_leakage_when_the_harmonic_does_not_radiate():
    bm = lm.bound_beta(Y_MEAN, SUB, F)
    period = 2.0 * np.pi / (3.0 * K0)  # K = 3 k0: n = -1 has |k_t| = |beta - 3 k0| > k0, evanescent
    coeffs = lm.sinusoid_fourier(Y_MEAN, 0.1j * 2.9 / ETA0, 6)
    beta, _ = lm.floquet_complex_beta(coeffs, period, 0.0, SUB, F, bm, n_max=6)
    assert abs(beta.imag) < 1e-6 * K0
    th, _, _ = lm.radiated_harmonic_geometry(bm, 2.0 * np.pi / period, 0.0, F)
    assert np.isnan(th)


def test_zero_contrast_gives_zero_leakage():
    assert lm.leakage_alpha_first_order(1.03 * K0, Y_MEAN, 0.0, 0.5, 0.0, SUB, F) == 0.0


def test_two_level_fourier_reproduces_the_profile():
    n_max = 40
    c = lm.two_level_fourier(0.0, 1.0, n_max, duty=0.5)
    x = np.array([0.1, 0.3, 0.7, 0.9])
    p = np.arange(-2 * n_max, 2 * n_max + 1)
    y = np.array([np.sum(c * np.exp(-2j * np.pi * p * xi)) for xi in x]).real
    # level 1 is centred on x = 0 and fills `duty` of the period
    assert np.allclose(y, [1.0, 0.0, 0.0, 1.0], atol=0.08)
    d = lm.two_level_fourier(0.0, 1.0, n_max, duty=0.25)
    y2 = np.array([np.sum(d * np.exp(-2j * np.pi * p * xi)) for xi in (0.05, 0.2, 0.5, 0.95)]).real
    assert np.allclose(y2, [1.0, 0.0, 0.0, 1.0], atol=0.1)
    assert lm.two_level_c1(0.0, 1.0, 0.5) == pytest.approx(1.0 / np.pi)
