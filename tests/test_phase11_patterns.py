"""Sec. 10, 28: the required controlled spatial-pattern generators."""

from src.phase11 import patterns as pt


def test_uniform_is_all_one_geometry():
    grid = pt.uniform(5, (2.8, 2.8))
    assert all(cell == (2.8, 2.8) for row in grid for cell in row)


def test_single_perturbation_changes_exactly_one_cell():
    grid = pt.single_perturbation(5, (2.8, 2.8), (4.4, 4.4))
    flat = [cell for row in grid for cell in row]
    assert flat.count((4.4, 4.4)) == 1
    assert flat.count((2.8, 2.8)) == 24
    assert grid[2][2] == (4.4, 4.4)   # defaults to center


def test_smooth_gradient_is_monotonic_along_x():
    grid = pt.smooth_gradient(7, 1.4, 4.2, 2.8)
    row = [grid[0][j][0] for j in range(7)]
    assert row == sorted(row)
    assert row[0] == 1.4 and row[-1] == 4.2


def test_abrupt_boundary_splits_into_two_regions():
    grid = pt.abrupt_boundary(6, (1.4, 1.4), (4.2, 4.2), boundary_col=3)
    for row in grid:
        assert row[:3] == [(1.4, 1.4)] * 3
        assert row[3:] == [(4.2, 4.2)] * 3


def test_checkerboard_and_alternating_are_identical():
    a = pt.checkerboard(5, (1.0, 1.0), (4.4, 4.4))
    b = pt.alternating(5, (1.0, 1.0), (4.4, 4.4))
    assert a == b
    # verify the actual parity rule directly.
    for i in range(5):
        for j in range(5):
            expected = (1.0, 1.0) if (i + j) % 2 == 0 else (4.4, 4.4)
            assert a[i][j] == expected


def test_strong_isolated_defect_uses_the_domain_extreme_farthest_from_center():
    grid = pt.strong_isolated_defect(5, (1.4, 1.4))
    assert grid[2][2] == (4.4, 4.4)   # 4.4 is farther from 1.4 than 1.0 is
    grid2 = pt.strong_isolated_defect(5, (4.2, 4.2))
    assert grid2[2][2] == (1.0, 1.0)   # 1.0 is farther from 4.2 than 4.4 is


def test_single_and_multi_cell_discontinuity_widths():
    single = pt.single_cell_discontinuity(7, (2.8, 2.8), (4.4, 4.4))
    multi = pt.multi_cell_discontinuity(7, (2.8, 2.8), (4.4, 4.4), width=3)
    for row in single:
        assert sum(1 for cell in row if cell == (4.4, 4.4)) == 1
    for row in multi:
        assert sum(1 for cell in row if cell == (4.4, 4.4)) == 3


def test_all_pattern_values_stay_within_the_validated_domain():
    from src.phase10.environment_geometry import in_domain

    grids = [
        pt.uniform(5, (2.8, 2.8)), pt.smooth_gradient(7), pt.abrupt_boundary(7),
        pt.checkerboard(7), pt.strong_isolated_defect(5, (2.8, 2.8)),
    ]
    for grid in grids:
        for row in grid:
            for l, w in row:
                assert in_domain(l, w)
