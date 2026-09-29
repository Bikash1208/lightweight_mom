"""Sec. 30-33, 47: the gradient and abrupt-discontinuity generalization tests."""

import numpy as np

from src.materials import Metal, Substrate
from src.phase10 import generalization as gen
from src.phase10 import surrogate as surr_mod
from src.phase10.isolated_reference import IsolatedReferenceCache

SUB = Substrate("Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3)
MET = Metal()
P = 4.8e-3
F = 12.5e9


def _dummy_surrogate():
    cache = IsolatedReferenceCache(P, SUB, MET, 8)
    return surr_mod.EnvironmentSurrogate(cache, np.zeros(surr_mod.N_FEATURES), np.zeros(surr_mod.N_FEATURES),
                                          0.0, surr_mod.DomainSpec(frequency_hz=F))


def test_build_gradient_grid_is_monotonic_along_x():
    grid = gen.build_gradient_grid_mm()
    row = [grid[0][j][0] for j in range(gen.GRID_N)]
    assert row == sorted(row)


def test_build_abrupt_grid_has_exactly_two_regions():
    grid = gen.build_abrupt_grid_mm(region_a_mm=(1.0, 1.0), region_b_mm=(4.4, 4.4), boundary_col=3)
    values = {grid[0][j] for j in range(gen.GRID_N)}
    assert values == {(1.0, 1.0), (4.4, 4.4)}
    assert grid[0][0] == (1.0, 1.0)
    assert grid[0][gen.GRID_N - 1] == (4.4, 4.4)


def test_run_gradient_test_produces_one_case_per_interior_cell():
    surrogate = _dummy_surrogate()
    grid_mm, cases = gen.run_gradient_test(P, SUB, MET, F, surrogate, mesh_cells_per_subcell=8)
    n_interior = len(gen.INTERIOR_RANGE)
    assert len(cases) == n_interior * n_interior
    for c in cases:
        assert isinstance(c.px_mom, complex)


def test_run_abrupt_discontinuity_test_labels_distance_from_boundary():
    surrogate = _dummy_surrogate()
    grid_mm, cases, labeled = gen.run_abrupt_discontinuity_test(P, SUB, MET, F, surrogate,
                                                                 mesh_cells_per_subcell=8, boundary_col=3)
    distances = {dist for dist, _ in labeled}
    assert distances == {0, 1, 2}   # interior columns 1,2 | 4,5 at increasing distance from col 3
