"""Sec. 17, 47: the distance-2 (two-period) long-range test."""

from src.materials import Metal, Substrate
from src.phase10 import long_range as lr

SUB = Substrate("Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3)
MET = Metal()
P = 4.8e-3
F = 12.5e9
CENTER = (2.8, 2.8)


def test_solve_5x5_homogeneous_is_physically_reasonable():
    s0 = lr.solve_5x5_homogeneous(CENTER, P, SUB, MET, F, mesh_cells_per_subcell=8)
    assert isinstance(s0, complex)
    assert abs(s0) > 0


def test_run_long_range_test_covers_all_four_directions():
    s0 = lr.solve_5x5_homogeneous(CENTER, P, SUB, MET, F, mesh_cells_per_subcell=8)
    results = lr.run_long_range_test(CENTER, (1.0, 1.0), s0, P, SUB, MET, F, mesh_cells_per_subcell=8)
    assert set(results.keys()) == {"N2", "S2", "W2", "E2"}
    for c in results.values():
        assert c.delta_s_abs == abs(c.delta_s)
        assert c.px != s0   # a distance-2 perturbation must change the solved system somewhat
