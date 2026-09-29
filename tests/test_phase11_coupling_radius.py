"""Sec. 5-7, 28: the growing-pattern truncation series and shell increments."""

import pytest

from src.materials import Metal, Substrate
from src.phase11 import coupling_radius as cr

SUB = Substrate("Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3)
MET = Metal()
P = 4.8e-3
F = 12.5e9
CENTER = (2.8, 2.8)


def test_build_growing_pattern_places_true_center_at_the_center():
    grid = cr.build_growing_pattern_mm(CENTER, 9)
    mid = 4
    assert grid[mid][mid] == CENTER


def test_truncate_centered_preserves_the_center_value():
    full = cr.build_growing_pattern_mm(CENTER, 9)
    small = cr.truncate_centered(full, 3)
    assert len(small) == 3 and len(small[0]) == 3
    assert small[1][1] == CENTER


def test_truncate_centered_rejects_n_larger_than_the_pattern():
    full = cr.build_growing_pattern_mm(CENTER, 5)
    with pytest.raises(ValueError):
        cr.truncate_centered(full, 7)


def test_run_coupling_radius_series_is_monotonic_in_size_and_deterministic():
    points = cr.run_coupling_radius_series(CENTER, P, SUB, MET, F, sizes=(1, 3, 5))
    assert [p.n for p in points] == [1, 3, 5]
    for p in points:
        assert isinstance(p.response, complex)
        assert p.n_unknowns > 0

    points2 = cr.run_coupling_radius_series(CENTER, P, SUB, MET, F, sizes=(1, 3, 5))
    assert [p.response for p in points] == [p.response for p in points2]   # deterministic, no ML


def test_compute_shell_increments_matches_direct_subtraction():
    points = cr.run_coupling_radius_series(CENTER, P, SUB, MET, F, sizes=(1, 3, 5))
    incs = cr.compute_shell_increments(points)
    assert len(incs) == 2
    assert incs[0].n_from == 1 and incs[0].n_to == 3
    expected = abs(points[1].response - points[0].response)
    assert incs[0].delta_complex == pytest.approx(expected)
