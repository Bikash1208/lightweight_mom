"""Phase 15 Part E checks: the surface-wave pair rule, the patch-plane reference, the quantization
definition, and the secondary/mirror beam checks. Synthetic inputs throughout (no MoM solve)."""

import numpy as np
import pytest

from src.constants import ETA0
from src.dispersion.analytic_tm0 import tm0_beta
from src.dispersion.sheet_model import Y0, SheetModelMode, _y_slab, find_beta_from_gamma
from src.errors import LwmomError
from src.frequency import angular_frequency, free_space_wavenumber
from src.hologram.coded_synthesis import (
    aperture_coupling_loss_db,
    beta_reference_field,
    textbook_quantization_loss_db,
)
from src.hologram.two_state_cells import (
    LIBRARY_REFERENCE_OFFSET_M,
    TwoStateCell,
    best_valid_pair,
    patch_plane_gamma,
    surface_wave_of_state,
)
from src.labels import quantization_loss_db
from src.materials import Substrate
from src.pattern.vector_aperture_pattern import (
    VectorPatternResult,
    first_null_half_width_deg,
    mirror_beam_level_db,
    strongest_secondary_beam,
)

F_HZ = 10.54e9
SUB = Substrate("user", 2.2, 0.0009, 1.575e-3)


def _bare_slab_gamma():
    k0 = free_space_wavenumber(F_HZ)
    omega = angular_frequency(F_HZ)
    y_slab = _y_slab(0.0, omega, SUB.eps_complex, k0, SUB.thickness)
    return (Y0 - y_slab) / (Y0 + y_slab)


def _mode(beta):
    return SheetModelMode(beta, ETA0 * np.sqrt((beta / free_space_wavenumber(F_HZ)) ** 2 - 1.0), 0j, 0.0, None, False)


def _cell(label, phase_deg, geometry=None):
    return TwoStateCell(
        label, label, complex(np.exp(1j * np.radians(phase_deg))), geometry or {}, reference_offset_m=0.0
    )


def test_sheet_model_from_gamma_recovers_the_analytic_bare_slab_tm0():
    mode = find_beta_from_gamma(_bare_slab_gamma(), SUB, F_HZ)
    assert mode.beta_rad_per_m == pytest.approx(tm0_beta(F_HZ, 2.2, 1.575e-3), rel=1e-6)


def test_patch_plane_gamma_undoes_the_library_reference_offset():
    gamma_patch = _bare_slab_gamma()
    k0 = free_space_wavenumber(F_HZ)
    gamma_top = gamma_patch * np.exp(-2j * k0 * LIBRARY_REFERENCE_OFFSET_M)
    assert patch_plane_gamma(gamma_top, F_HZ, LIBRARY_REFERENCE_OFFSET_M) == pytest.approx(gamma_patch, abs=1e-12)
    cell = TwoStateCell("c", "c", complex(gamma_top), {}, reference_offset_m=LIBRARY_REFERENCE_OFFSET_M)
    assert cell.patch_plane_value(F_HZ) == pytest.approx(gamma_patch, abs=1e-12)


def test_a_short_has_no_surface_wave_and_is_not_solved():
    cell = TwoStateCell("short", "short", complex(-1.0), {})
    assert surface_wave_of_state(cell, F_HZ, SUB, period_m=8e-3) is None


def test_pair_rule_skips_states_without_a_surface_wave():
    a, b, c = _cell("a", 0.0), _cell("b", 180.0), _cell("c", 170.0)
    surface_waves = [_mode(230.0), None, _mode(231.0)]
    choice = best_valid_pair([a, b, c], surface_waves, F_HZ, SUB)
    assert {choice.state_0.label, choice.state_1.label} == {"a", "c"}
    assert choice.achieved_separation_deg == pytest.approx(170.0, abs=1e-6)


def test_pair_rule_rejects_betas_too_far_apart_even_with_a_bigger_separation():
    a, b, c = _cell("a", 0.0), _cell("b", 180.0), _cell("c", 170.0)
    surface_waves = [_mode(230.0), _mode(260.0), _mode(231.0)]
    choice = best_valid_pair([a, b, c], surface_waves, F_HZ, SUB)
    assert {choice.state_0.label, choice.state_1.label} == {"a", "c"}


def test_pair_rule_reports_mean_beta_and_spread():
    a, b = _cell("a", 0.0), _cell("b", 175.0)
    choice = best_valid_pair([a, b], [_mode(230.0), _mode(234.0)], F_HZ, SUB)
    assert choice.beta_bar_rad_per_m == pytest.approx(232.0)
    assert choice.beta_spread_fraction == pytest.approx(4.0 / 232.0)
    assert choice.surface_wave_0.beta_rad_per_m == 230.0


def test_no_valid_pair_is_an_explicit_error():
    a, b = _cell("a", 0.0), _cell("b", 180.0)
    with pytest.raises(LwmomError) as excinfo:
        best_valid_pair([a, b], [_mode(230.0), None], F_HZ, SUB)
    assert excinfo.value.code == "TWO_STATE_NO_VALID_PAIR"


def test_textbook_quantization_loss_values():
    assert textbook_quantization_loss_db(1) == pytest.approx(3.9224, abs=1e-3)
    assert textbook_quantization_loss_db(2) == pytest.approx(0.9121, abs=1e-3)


def test_aperture_coupling_reproduces_textbook_for_uniformly_distributed_phase():
    n = 201
    index = np.arange(n * n, dtype=float)
    continuous = np.mod(index * 0.6180339887 * 360.0, 360.0)
    quantized = np.where(np.cos(np.radians(continuous)) >= 0.0, 0.0, 180.0)
    loss = aperture_coupling_loss_db(np.exp(1j * np.radians(quantized)), np.exp(1j * np.radians(continuous)))
    assert loss == pytest.approx(textbook_quantization_loss_db(1), abs=0.05)


def test_aperture_coupling_is_zero_for_identical_apertures():
    e = np.exp(1j * np.linspace(0, 6, 50))
    assert aperture_coupling_loss_db(e, e) == pytest.approx(0.0, abs=1e-9)


def test_beta_reference_field_has_phase_minus_beta_rho_about_the_feed():
    x = np.linspace(-0.05, 0.05, 5)
    X, Y = np.meshgrid(x, x)

    class _Cell:
        x_m, y_m, dx_m, dy_m = X, Y, 0.02, 0.02

    beta = 231.8
    ref = beta_reference_field(_Cell, (0.0, 0.0), beta)
    rho = np.hypot(X, Y)
    assert np.allclose(np.exp(1j * np.angle(ref.Ex)), np.exp(-1j * beta * rho))


def test_quantization_label_states_the_textbook_value():
    label = quantization_loss_db(3.86, 1, 3.92)
    assert label == "QUANTIZATION_LOSS(3.86dB,1bit,textbook 3.92dB)"


def _vector_result(power, theta_grid, phi_grid, main_theta, main_phi):
    amplitude = np.sqrt(power)
    return VectorPatternResult(
        theta_deg=theta_grid,
        phi_deg=phi_grid,
        e_theta=amplitude.astype(complex),
        e_phi=np.zeros_like(amplitude, dtype=complex),
        e_r=amplitude.astype(complex),
        e_l=np.zeros_like(amplitude, dtype=complex),
        power_db=10 * np.log10(np.maximum(power, 1e-300) / power.max()),
        main_beam_theta_deg=main_theta,
        main_beam_phi_deg=main_phi,
        sll_db=None,
        directivity_estimate_dbi=0.0,
    )


def _gaussian_beam(theta_grid, phi_grid, theta0, phi0, width_deg=3.0):
    dphi = (phi_grid - phi0 + 180.0) % 360.0 - 180.0
    return np.exp(
        -((theta_grid - theta0) ** 2 + (dphi * np.sin(np.radians(theta_grid + 1e-9))) ** 2) / (2 * width_deg**2)
    )


def _two_beam_result(secondary_amplitude=0.5, mirror_amplitude=0.0):
    theta = np.linspace(0.0, 90.0, 91)
    phi = np.linspace(0.0, 360.0, 72, endpoint=False)
    theta_grid, phi_grid = np.meshgrid(theta, phi, indexing="ij")
    power = (
        _gaussian_beam(theta_grid, phi_grid, 20.0, 0.0) ** 2
        + secondary_amplitude**2 * _gaussian_beam(theta_grid, phi_grid, 50.0, 90.0) ** 2
        + mirror_amplitude**2 * _gaussian_beam(theta_grid, phi_grid, 20.0, 180.0) ** 2
    )
    return _vector_result(power, theta_grid, phi_grid, 20.0, 0.0)


def test_strongest_secondary_beam_is_found_outside_the_main_lobe():
    result = _two_beam_result(secondary_amplitude=0.5)
    found = strongest_secondary_beam(result, exclusion_deg=10.0)
    assert found is not None
    theta, phi, level_db = found
    assert theta == pytest.approx(50.0, abs=3.0)
    assert phi == pytest.approx(90.0, abs=3.0)
    assert level_db == pytest.approx(10.0 * np.log10(0.25), abs=0.6)


def test_main_lobe_shoulder_is_excluded_from_the_secondary_search():
    result = _two_beam_result(secondary_amplitude=0.0)
    assert strongest_secondary_beam(result, exclusion_deg=10.0) is None


def test_mirror_beam_level_is_sampled_at_phi_plus_180():
    result = _two_beam_result(secondary_amplitude=0.0, mirror_amplitude=0.3)
    level_db, theta, phi = mirror_beam_level_db(result)
    assert theta == pytest.approx(20.0, abs=3.0)
    assert phi == pytest.approx(180.0, abs=3.0)
    assert level_db == pytest.approx(10.0 * np.log10(0.09), abs=0.3)


def test_first_null_half_width_matches_the_airy_formula():
    assert first_null_half_width_deg(0.3, 0.0284, 0.0) == pytest.approx(np.degrees(np.arcsin(1.22 * 0.0284 / 0.3)))
    assert first_null_half_width_deg(0.3, 0.0284, 60.0) == pytest.approx(
        2.0 * first_null_half_width_deg(0.3, 0.0284, 0.0), rel=1e-6
    )
