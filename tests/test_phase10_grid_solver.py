"""Sec. 28-32, 47: the generic multi-cell-extraction solve wrapper (one solve, many read-outs)."""

from src.materials import Metal, Substrate
from src.phase9.supercell_solver import solve_isolated_cell
from src.phase10.grid_solver import solve_grid

SUB = Substrate("Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3)
MET = Metal()
P = 4.8e-3
F = 12.5e9


def test_solve_grid_homogeneous_matches_isolated_cell():
    grid_m = [[(2.8e-3, 2.8e-3)] * 3 for _ in range(3)]
    sg = solve_grid(3, P, SUB, MET, grid_m, F, mesh_cells_per_subcell=8)
    embedded = sg.dipole_at(1, 1)
    # matches phase9.local_periodicity's own convention: floquet_order = FLOQUET_FACTOR *
    # mesh_cells_per_subcell for the isolated (n=1-equivalent) reference, not scaled by n.
    isolated = solve_isolated_cell(P, P, 2.8e-3, 2.8e-3, SUB, MET, F, 8, 2 * 8)
    rel = abs(embedded.px - isolated.dipole.px) / abs(isolated.dipole.px)
    assert rel < 1e-8


def test_solve_grid_dipole_at_reads_different_subcells_from_one_solve():
    grid_m = [[(1.4e-3, 1.4e-3)] * 3 for _ in range(3)]
    grid_m[0][0] = (4.4e-3, 4.4e-3)  # NW corner is very different
    sg = solve_grid(3, P, SUB, MET, grid_m, F, mesh_cells_per_subcell=8)
    center = sg.dipole_at(1, 1)
    corner = sg.dipole_at(0, 0)
    assert abs(center.px - corner.px) > 0  # genuinely different sub-cell responses from ONE solve
