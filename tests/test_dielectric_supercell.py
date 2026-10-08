"""Phase 19 item 5 machinery -- multi-region layers, rectangular truncation, strips of mixed cells."""

import numpy as np
import pytest

from src.dielectric_cell import OctagonPillarCell, solve_cell
from src.dielectric_cell.supercell import CellSpec, local_periodicity_prediction, solve_strip

F, P, TOTAL, PLATE, EPS, TAND = 28e9, 5e-3, 25e-3, 0.6e-3, 2.45, 0.0185
H = TOTAL - PLATE


def test_identical_pillar_cells_in_a_strip_reproduce_the_single_cell():
    c = OctagonPillarCell(P, 1.8e-3, H, EPS, TAND, PLATE, 22.5, "pillar")
    one = solve_cell(c, F, 0, 0, 6).t[1, 1]
    strip = solve_strip([CellSpec("pillar", 1.8e-3)] * 2, P, H, PLATE, EPS, TAND, True, F, 0, 0, (12, 6))["t"][1, 1]
    assert abs(strip / one - 1) < 2e-4


@pytest.mark.parametrize("through", [False, True])
def test_identical_hole_cells_in_a_strip_reproduce_the_single_cell(through):
    kw = dict(height_m=TOTAL, plate=0.0) if through else dict(height_m=H, plate=PLATE)
    c = OctagonPillarCell(P, 1.8e-3, kw["height_m"], EPS, TAND, kw["plate"], 22.5, "hole")
    one = solve_cell(c, F, 0, 0, 6).t[1, 1]
    strip = solve_strip([CellSpec("hole", 1.8e-3)] * 2, P, H, PLATE, EPS, TAND, through, F, 0, 0, (12, 6))["t"][1, 1]
    assert abs(strip / one - 1) < 2e-4


def test_local_periodicity_prediction_of_identical_cells_has_only_the_zeroth_order():
    pred = local_periodicity_prediction([0.7 + 0.2j] * 4, (0, 1, -1, 2))
    assert pred[0] == pytest.approx(0.7 + 0.2j)
    assert all(abs(pred[m]) < 1e-12 for m in (1, -1, 2))


def test_zeroth_order_of_a_mixed_strip_is_energy_bounded_and_orders_are_requested_ones():
    res = solve_strip(
        [CellSpec("pillar", 1.2e-3), CellSpec("pillar", 1.2e-3), CellSpec("hole", 2.0e-3), CellSpec("hole", 2.0e-3)],
        P,
        H,
        PLATE,
        EPS,
        TAND,
        True,
        F,
        0,
        0,
        (16, 4),
        orders=[(0, 0), (1, 0), (-1, 0)],
    )
    assert set(res["orders"]) == {(0, 0), (1, 0), (-1, 0)}
    amp = np.sqrt(sum(abs(v[0, 1]) ** 2 for v in res["orders"].values()))
    assert amp <= 1.0 + 1e-6  # transmitted power in the propagating orders cannot exceed the incident power
    # +1 / -1 orders mirror each other for a normal-incidence strip that is symmetric in |m|
    assert abs(res["orders"][(1, 0)][0, 1]) == pytest.approx(abs(res["orders"][(-1, 0)][0, 1]), rel=1e-6)


def test_checkerboard_of_identical_cells_reproduces_the_single_cell():
    from src.dielectric_cell.supercell import solve_grid

    c = OctagonPillarCell(P, 1.8e-3, TOTAL, EPS, TAND, 0.0, 22.5, "hole")
    one = solve_cell(c, F, 0, 0, 8).t[1, 1]
    a = CellSpec("hole", 1.8e-3)
    chk = solve_grid([[a, a], [a, a]], P, H, PLATE, EPS, TAND, True, F, 0, 0, (10, 10))["t"][1, 1]
    assert abs(chk / one - 1) < 3e-3  # truncation-level agreement (5 harmonics per cell here)


def test_strip_is_exactly_the_one_column_grid():
    from src.dielectric_cell.supercell import solve_grid

    a, b = CellSpec("pillar", 1.2e-3), CellSpec("hole", 2.0e-3)
    s1 = solve_strip([a, b], P, H, PLATE, EPS, TAND, True, F, 0, 0, (8, 4))["t"][1, 1]
    s2 = solve_grid([[a], [b]], P, H, PLATE, EPS, TAND, True, F, 0, 0, (8, 4))["t"][1, 1]
    assert s1 == s2


def test_checkerboard_transmits_x_and_y_polarisation_equally():
    from src.dielectric_cell.supercell import solve_grid

    a, b = CellSpec("pillar", 1.2e-3), CellSpec("hole", 2.0e-3)
    t = solve_grid([[a, b], [b, a]], P, H, PLATE, EPS, TAND, True, F, 0, 0, (8, 8))["t"]
    # the checkerboard is symmetric under x <-> y, so the zeroth-order transmission is the same for x- and y-pol
    assert abs(t[1, 1] - t[0, 0]) < 2e-3 * abs(t[0, 0])
