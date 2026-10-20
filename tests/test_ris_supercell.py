"""Supercell MoM with per-site diode states and Floquet-order extraction (Phase 20 R3 item 4, src/ris/supercell.py)."""

import math
import warnings
from dataclasses import replace

import numpy as np
import pytest

from src.excitation import ReferencePlanes
from src.materials import Substrate
from src.ris import cells as cells_mod
from src.ris.cells import OFF, ON, Lp2Cell, ideal_kit, madp_kit
from src.ris.supercell import build_supercell, default_mesh, solve_orders, sparse_load_matrix
from src.sparameters import compute_s_parameters_oblique

PL = ReferencePlanes(0.0, None, 1.575e-3)
CELL = Lp2Cell()
N = 24  # cells per period (0.5 mm): fast


def _mesh(nx, ny):
    return default_mesh(N, nx, ny)


def test_one_by_one_supercell_reproduces_the_isolated_cell_s11():
    st = (ON, OFF)
    uc = CELL.unit_cell(st, madp_kit())
    sc = build_supercell(CELL, [[st]], madp_kit())
    th = math.radians(20.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ref = compute_s_parameters_oblique(uc, 10e9, _mesh(1, 1), PL, th, 0.0).s11
        got = solve_orders(sc, 10e9, _mesh(1, 1), PL, th, 0.0, [sc.loads])[0][(0, 0)]
    assert np.allclose(got, ref, atol=1e-9)


def test_sparse_load_matrix_matches_the_dense_one():
    from src.basis_functions import RooftopBasis
    from src.loads import load_matrix

    uc = CELL.unit_cell((ON, OFF), madp_kit())
    basis = RooftopBasis(uc, N, N)
    dense = load_matrix(basis, uc.loads, 2 * math.pi * 10e9)
    d = sparse_load_matrix(basis, uc.loads, 2 * math.pi * 10e9)
    rebuilt = np.zeros_like(dense)
    for (i, j), v in d.items():
        rebuilt[i, j] = v
    assert np.allclose(rebuilt, dense, rtol=0, atol=1e-12 * np.abs(dense).max())


def test_uniform_two_by_one_supercell_has_no_anomalous_orders_and_the_isolated_specular_response():
    st = (ON, OFF)
    sc = build_supercell(CELL, [[st, st]], madp_kit())
    th = math.radians(25.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        iso = compute_s_parameters_oblique(CELL.unit_cell(st, madp_kit()), 14e9, _mesh(1, 1), PL, th, 0.0).s11
        orders = solve_orders(sc, 14e9, _mesh(2, 1), PL, th, 0.0, [sc.loads])[0]
    assert (-1, 0) in orders or (1, 0) in orders  # an anomalous order propagates at 14 GHz, theta 25 deg
    for key, s in orders.items():
        if key != (0, 0):
            assert np.abs(s).max() < 5e-3
    # the Floquet truncation differs (the supercell sums twice as many harmonics along x): agreement to ~1e-2, not machine precision
    assert np.allclose(orders[(0, 0)], iso, atol=1.5e-2)


def test_power_balance_over_all_orders_for_a_lossless_nonuniform_supercell():
    lossless = replace(CELL)  # the cell's substrate is replaced below
    kit = ideal_kit()
    codes = [[(ON, OFF), (OFF, ON)]]
    sc = build_supercell(lossless, codes, kit)
    sc = replace(sc, substrate=Substrate("lossless", 2.2, 0.0, 1.575e-3))
    th = math.radians(20.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        orders = solve_orders(sc, 14e9, _mesh(2, 1), PL, th, 0.0, [sc.loads])[0]
    assert len(orders) >= 2
    for mode in range(2):
        total = sum(float(np.sum(np.abs(s[:, mode]) ** 2)) for s in orders.values())
        assert total == pytest.approx(1.0, abs=2e-3)
    assert (
        max(float(np.abs(s).max()) for k, s in orders.items() if k != (0, 0)) > 0.05
    )  # the coding does scatter into other orders


def test_supercell_geometry_translates_loads_and_metal():
    c = cells_mod.Cp2Cell()
    sc = build_supercell(c, [[c.states()[0], c.states()[1]], [c.states()[2], c.states()[3]]], madp_kit())
    assert sc.lattice.period_x == pytest.approx(2 * c.period) and sc.lattice.period_y == pytest.approx(2 * c.period)
    assert len(sc.loads) == 16
    xs = sorted({round(ld.x, 9) for ld in sc.loads})
    assert xs[0] < -c.period / 4 and xs[-1] > c.period / 4
