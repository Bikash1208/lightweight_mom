"""Sec. 3, 28: the Phase-11-owned solve wrapper -- homogeneous equivalence, and that adding a
polarization parameter reproduces Phase 10's own results at the default "x" polarization (i.e. this
new wrapper is not a silent physics change)."""

import pytest

from src.materials import Metal, Substrate
from src.phase9.supercell_solver import solve_isolated_cell
from src.phase10.grid_solver import solve_grid as p10_solve_grid
from src.phase11.reference_solver import solve_grid, solve_isolated

SUB = Substrate("Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3)
MET = Metal()
P = 4.8e-3
F = 12.5e9
CENTER = (2.8e-3, 2.8e-3)


def test_solve_grid_matches_phase10s_own_wrapper_at_default_polarization():
    grid_m = [[CENTER] * 3 for _ in range(3)]
    sg11 = solve_grid(3, P, SUB, MET, grid_m, F, mesh_cells_per_subcell=8)
    sg10 = p10_solve_grid(3, P, SUB, MET, grid_m, F, mesh_cells_per_subcell=8)
    d11 = sg11.dipole_at(1, 1)
    d10 = sg10.dipole_at(1, 1)
    assert d11.px == pytest.approx(d10.px)
    assert sg11.response_at(1, 1) == pytest.approx(d10.px)  # x-pol response is px


def test_solve_grid_homogeneous_matches_isolated_cell():
    grid_m = [[CENTER] * 3 for _ in range(3)]
    sg = solve_grid(3, P, SUB, MET, grid_m, F, mesh_cells_per_subcell=8)
    isolated = solve_isolated_cell(P, P, CENTER[0], CENTER[1], SUB, MET, F, 8, 16)
    rel = abs(sg.response_at(1, 1) - isolated.dipole.px) / abs(isolated.dipole.px)
    assert rel < 1e-8


def test_solve_grid_n1_matches_solve_isolated():
    grid_m = [[CENTER]]
    sg = solve_grid(1, P, SUB, MET, grid_m, F, mesh_cells_per_subcell=8)
    response, _ = solve_isolated(P, CENTER[0], CENTER[1], SUB, MET, F, mesh_cells_per_side=8)
    rel = abs(sg.response_at(0, 0) - response) / abs(response)
    assert rel < 1e-8


def test_y_polarization_reads_py_not_px():
    grid_m = [[(1.4e-3, 4.2e-3)] * 3 for _ in range(3)]  # anisotropic patch: px != py expected
    sg_x = solve_grid(3, P, SUB, MET, grid_m, F, mesh_cells_per_subcell=8, polarization="x")
    sg_y = solve_grid(3, P, SUB, MET, grid_m, F, mesh_cells_per_subcell=8, polarization="y")
    dx = sg_x.dipole_at(1, 1)
    dy = sg_y.dipole_at(1, 1)
    assert sg_x.response_at(1, 1) == pytest.approx(dx.px)
    assert sg_y.response_at(1, 1) == pytest.approx(dy.py)
