"""Phase 12B Task D: grounded period study."""

import pytest

from src.phase12b.period_study import (
    aligned_cells_for_period,
    convergence_for_near_resonant_patch,
    coverage_by_frequency,
    grating_and_surface_wave_report,
    l_family_for_period,
    resonant_length_mm,
    run_period_sweep,
)


@pytest.mark.parametrize("period_mm,expected_n", [(4.8, 48), (8.0, 40), (12.0, 60)])
def test_aligned_cells_for_period_matches_measurement(period_mm, expected_n):
    assert aligned_cells_for_period(period_mm) == expected_n


def test_l_family_spans_roughly_0p2_to_0p94_of_the_period():
    fam = l_family_for_period(8.0)
    assert fam[0] == pytest.approx(1.6, abs=0.01)
    assert fam[-1] <= 8.0 * 0.94 + 0.2  # +0.2mm: rounded to the nearest 0.4mm step


def test_run_period_sweep_rejects_an_unaligned_cells_override():
    with pytest.raises(ValueError, match="does not align"):
        run_period_sweep(8.0, l_values_mm=(2.0, 3.0), frequencies_hz=(12.5e9,), cells=17)


def test_run_period_sweep_small_scale_and_resonant_length():
    points, n = run_period_sweep(8.0, l_values_mm=(2.0, 4.0, 6.0), frequencies_hz=(12.5e9,))
    assert n == aligned_cells_for_period(8.0, (2.0, 4.0, 6.0))
    assert len(points) == 3
    l = resonant_length_mm(points, 12.5e9)
    assert l in (2.0, 4.0, 6.0)


def test_coverage_by_frequency_small_scale():
    points, _ = run_period_sweep(8.0, l_values_mm=(2.0, 4.0, 6.0), frequencies_hz=(12.5e9,))
    cov = coverage_by_frequency(points)
    assert 12.5e9 in cov
    assert cov[12.5e9].n_points == 3


def test_grating_and_surface_wave_report_for_12mm_shows_a_shrinking_margin_with_frequency():
    """The brief's own predicted finding: the surface-wave margin shrinks as frequency (hence k0)
    grows, for the larger period."""
    r = grating_and_surface_wave_report(12.0, 6.0, frequencies_hz=(10.0e9, 15.0e9))
    assert r.surface_wave_margins[15.0e9] < r.surface_wave_margins[10.0e9]


def test_grating_onset_oblique_45deg_falls_inside_the_band_for_12mm_period():
    r = grating_and_surface_wave_report(12.0, 6.0, frequencies_hz=(10.0e9, 15.0e9))
    assert r.grating_onset_oblique_hz[45.0] < 15.0e9


def test_grating_report_specular_only_across_band_for_4p8_and_8mm():
    for period_mm in (4.8, 8.0):
        r = grating_and_surface_wave_report(period_mm, period_mm * 0.4, frequencies_hz=(10.0e9, 15.0e9))
        assert r.specular_only_across_band is True


def test_convergence_for_near_resonant_patch_small_scale():
    points, diffs = convergence_for_near_resonant_patch(8.0, 4.0, factors=(1, 2))
    assert len(points) == 2
    assert points[1].cells == 2 * points[0].cells
