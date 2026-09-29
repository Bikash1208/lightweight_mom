"""Sec. 20: solving the isolated cell and a heterogeneous supercell through the exact same,
unmodified `solver.assemble_system`/`solve` MoM pipeline. The homogeneous-supercell exact-match
check is the load-bearing correctness proof for the whole Phase 9B methodology."""

import pytest

from src.materials import Metal, Substrate
from src.phase9.supercell_geometry import build_supercell
from src.phase9.supercell_solver import solve_isolated_cell, solve_supercell_center

SUB = Substrate("Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3)
MET = Metal()
P = 4.8e-3
FREQ = 12.5e9
L = W = 2.8e-3


def test_homogeneous_supercell_reproduces_the_isolated_cell_exactly():
    """A 3x3 supercell of 9 IDENTICAL cells, repeated periodically, is mathematically identical to
    one cell repeated periodically -- if the centre-cell dipole moment does not match the isolated
    solve to near machine precision (given matched mesh/Floquet accuracy), something in the
    supercell geometry/solve pipeline is wrong."""
    cells = 8
    iso = solve_isolated_cell(P, P, L, W, SUB, MET, FREQ, cells, 2 * cells)
    grid = [[(L, W)] * 3 for _ in range(3)]
    sc = build_supercell(3, P, P, SUB, MET, grid)
    emb = solve_supercell_center(sc, (L, W), FREQ, cells, 2 * cells * 3)
    rel_error = abs(emb.dipole.px - iso.dipole.px) / abs(iso.dipole.px)
    assert rel_error < 1e-10


def test_solve_supercell_center_rejects_even_n():
    grid = [[(L, W)] * 2 for _ in range(2)]
    sc = build_supercell(2, P, P, SUB, MET, grid)
    with pytest.raises(ValueError):
        solve_supercell_center(sc, (L, W), FREQ, 6, 12)


def test_solve_supercell_center_rejects_mismatched_center_lw():
    grid = [[(L, W)] * 3 for _ in range(3)]
    sc = build_supercell(3, P, P, SUB, MET, grid)
    with pytest.raises(ValueError):
        solve_supercell_center(sc, (1.0e-3, 1.0e-3), FREQ, 6, 36)


def test_heterogeneous_supercell_differs_measurably_from_the_isolated_cell():
    """A real, physically-different neighbor environment must produce a real, nonzero deviation
    (a basic sanity check that the composite geometry is not somehow solving the isolated cell
    over again)."""
    cells = 8
    iso = solve_isolated_cell(P, P, L, W, SUB, MET, FREQ, cells, 2 * cells)
    grid = [
        [(1.0e-3, 1.0e-3), (4.4e-3, 4.4e-3), (1.0e-3, 4.4e-3)],
        [(4.4e-3, 1.0e-3), (L, W), (1.0e-3, 1.0e-3)],
        [(4.4e-3, 4.4e-3), (1.0e-3, 4.4e-3), (4.4e-3, 1.0e-3)],
    ]
    sc = build_supercell(3, P, P, SUB, MET, grid)
    emb = solve_supercell_center(sc, (L, W), FREQ, cells, 2 * cells * 3)
    rel_error = abs(emb.dipole.px - iso.dipole.px) / abs(iso.dipole.px)
    assert rel_error > 0.1   # a real, substantial, physically-driven difference
