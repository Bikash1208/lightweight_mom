"""Phase 12C Task C: coupling extent, grounded."""

import pytest

from src.phase12c import coupling_study as cs
from src.phase12c import locality_study as ls


def test_gradient_pattern_snaps_to_the_aligned_family():
    grid = cs.build_gradient_pattern_mm(2.4, 1.6, 7.2, 7)
    assert len(grid) == 7 and len(grid[0]) == 7
    for row in grid:
        for l, w in row:
            assert l in ls.L_FAMILY_MM and l == w
    # every row is identical (a 1-D gradient along x, held fixed along y)
    assert all(row == grid[0] for row in grid)


def test_gradient_pattern_center_is_exact():
    grid = cs.build_gradient_pattern_mm(4.0, 1.6, 7.2, 5)
    mid = 5 // 2
    assert grid[0][mid] == (4.0, 4.0)


def test_checkerboard_pattern_matches_phase11_definition():
    grid = cs.build_checkerboard_pattern_mm(2.0, 6.0, 3)
    assert grid[1][1] == (2.0, 2.0)  # centre keeps the centre value
    assert grid[0][1] == (6.0, 6.0)  # cardinal neighbour is the contrast value
    assert grid[0][0] == (2.0, 2.0)  # diagonal neighbour equals the centre (5-cell pattern)


@pytest.mark.slow
def test_coupling_radius_series_is_monotone_truncation_of_one_pattern():
    grid = cs.build_gradient_pattern_mm(2.4, 1.6, 7.2, 3)
    points = cs.run_coupling_radius_series(grid, sizes=(1, 3))
    assert [p.n for p in points] == [1, 3]
    assert points[0].n_unknowns < points[1].n_unknowns


@pytest.mark.slow
def test_analyze_coupling_series_reuses_phase11_machinery_unmodified():
    grid = cs.build_gradient_pattern_mm(2.4, 1.6, 7.2, 3)
    points = cs.run_coupling_radius_series(grid, sizes=(1, 3))
    report = cs.analyze_coupling_series("test", points)
    assert report.convergence.status in ("CONVERGED", "NOT_CONVERGED", "STILL_SUBSTANTIALLY_CHANGING")
    assert len(report.increments) == 1
    assert report.decay_fit.law in ("NONE_FITTED", "EXPONENTIAL_CANDIDATE_FIT")
