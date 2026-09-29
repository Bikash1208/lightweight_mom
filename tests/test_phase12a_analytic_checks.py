"""Phase 12A: analytic bare-slab and surface-wave-pole checks."""

import pytest

from src.phase12a.analytic_checks import (
    bare_grounded_slab_checks,
    surface_wave_margin_table,
    textbook_grounded_slab_gamma,
    vanishing_thickness_limit,
)


def test_bare_grounded_slab_matches_textbook_to_noise_floor():
    points = bare_grounded_slab_checks(frequencies_hz=(10.0e9, 12.5e9, 15.0e9))
    for p in points:
        assert p.max_error < 1e-9


def test_lossless_bare_slab_has_unit_reflection():
    points = bare_grounded_slab_checks(frequencies_hz=(12.5e9,))
    lossless = [p for p in points if p.tan_delta == 0.0]
    assert lossless
    for p in lossless:
        assert abs(p.gamma_mine) == pytest.approx(1.0, abs=1e-9)


def test_lossy_bare_slab_has_sub_unit_reflection():
    points = bare_grounded_slab_checks(frequencies_hz=(12.5e9,))
    lossy = [p for p in points if p.tan_delta > 0.0]
    assert lossy
    for p in lossy:
        assert abs(p.gamma_mine) < 1.0


def test_vanishing_thickness_limit_is_minus_one():
    gamma = vanishing_thickness_limit(thickness_m=1e-9)
    assert abs(gamma - (-1.0)) < 1e-5


def test_textbook_formula_is_symmetric_in_sign_of_tan_delta_zero():
    g1 = textbook_grounded_slab_gamma(12.5e9, 2.2, 0.0, 1.575e-3)
    assert abs(g1) == pytest.approx(1.0, abs=1e-12)


def test_surface_wave_margin_table_reports_no_in_band_harmonic_for_library_stack():
    points = surface_wave_margin_table(periods_mm=(4.8, 8.0), frequencies_hz=(10.0e9, 12.5e9, 15.0e9))
    assert len(points) == 6
    for p in points:
        assert p.any_in_band is False
        assert p.nearest_margin_over_k0 > 0.0
