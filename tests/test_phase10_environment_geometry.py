"""Sec. 6-9, 12, 17, 47: center-cell definitions, the validated domain, and the 5-/9-/25-cell
environment grid builders."""

import pytest

from src.phase10 import environment_geometry as eg


def test_all_five_center_cells_are_on_the_validated_sweep_grid():
    for name, (l, w) in eg.CENTER_CELLS_MM.items():
        assert eg.in_domain(l, w)
        # exact 0.2 mm grid alignment (matches library_schema.py's L_VALUES_MM/W_VALUES_MM)
        assert round((l - eg.DOMAIN_L_MIN_MM) / eg.DOMAIN_STEP_MM, 6) % 1 == 0
        assert round((w - eg.DOMAIN_W_MIN_MM) / eg.DOMAIN_STEP_MM, 6) % 1 == 0


def test_in_domain_rejects_outside_the_sweep_range():
    assert not eg.in_domain(0.5, 2.0)
    assert not eg.in_domain(2.0, 5.0)
    assert eg.in_domain(1.0, 4.4)


def test_weak_neighbor_prefers_plus_direction_when_valid():
    result = eg.weak_neighbor((2.8, 2.8))
    assert result.direction == "+"
    assert result.lw_mm == pytest.approx((3.0, 3.0))


def test_moderate_neighbor_falls_back_to_minus_near_the_upper_domain_edge():
    result = eg.moderate_neighbor((4.2, 4.2))   # 4.2+0.6=4.8 is out of domain
    assert result.direction == "-"
    assert result.lw_mm == pytest.approx((3.6, 3.6))


def test_perturbed_neighbor_raises_when_neither_direction_is_valid():
    with pytest.raises(ValueError):
        eg.perturbed_neighbor((2.8, 2.8), 5.0)   # +-5mm leaves [1.0, 4.4] entirely


def test_strong_neighbor_extremes_are_center_independent():
    for name, center in eg.CENTER_CELLS_MM.items():
        low, high = eg.strong_neighbor_low_high(center)
        assert (low, high) == (eg.DOMAIN_L_MIN_MM, eg.DOMAIN_L_MAX_MM)


def test_build_5cell_environment_leaves_diagonals_at_center():
    grid = eg.build_5cell_environment((2.8, 2.8), {"N": (3.0, 3.0)})
    assert grid[0][0] == (2.8, 2.8)   # NW diagonal == center
    assert grid[0][2] == (2.8, 2.8)   # NE diagonal == center
    assert grid[0][1] == (3.0, 3.0)   # N cardinal == perturbed value
    assert grid[1][1] == (2.8, 2.8)   # center


def test_build_9cell_environment_allows_independent_diagonals():
    grid = eg.build_9cell_environment((2.8, 2.8), diagonal_mm={"NE": (4.4, 4.4)})
    assert grid[0][2] == (4.4, 4.4)
    assert grid[0][0] == (2.8, 2.8)   # NW untouched, defaults to center


def test_build_5x5_long_range_environment_only_touches_distance_two_cardinals():
    grid = eg.build_5x5_long_range_environment((2.8, 2.8), {"N2": (1.0, 1.0)})
    assert grid[0][2] == (1.0, 1.0)
    for i in range(5):
        for j in range(5):
            if (i, j) != (0, 2):
                assert grid[i][j] == (2.8, 2.8)


def test_grid_to_m_converts_every_entry():
    grid_mm = [[(2.8, 2.8)] * 3 for _ in range(3)]
    grid_m = eg.grid_to_m(grid_mm)
    assert grid_m[1][1] == pytest.approx((2.8e-3, 2.8e-3))
