"""Sec. 24, 28: the three-beam holographic test -- reference+object wave -> deterministic geometry
assignment -> ONE direct heterogeneous MoM solve -> far field."""

from src.materials import Metal, Substrate
from src.phase11 import holographic_three_beam as hb

SUB = Substrate("Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3)
MET = Metal()
P = 4.8e-3
F = 12.5e9
CENTER = (2.8, 2.8)


def test_build_holographic_target_grid_covers_the_interior_region():
    targets = hb.build_holographic_target_grid(5, P, SUB, MET, F, margin=1, mesh_cells_per_subcell=8)
    assert set(targets.keys()) == {(i, j) for i in range(1, 4) for j in range(1, 4)}
    for v in targets.values():
        assert isinstance(v, complex)


def test_assign_geometry_by_isolated_search_picks_from_the_response_table():
    table = {(1.0, 1.0): 1 + 0j, (2.8, 2.8): 5 + 0j, (4.4, 4.4): 10 + 0j}
    targets = {(0, 0): 4.9 + 0j, (0, 1): 0.9 + 0j}
    assigned = hb.assign_geometry_by_isolated_search(targets, table)
    assert assigned[(0, 0)] == (2.8, 2.8)
    assert assigned[(0, 1)] == (1.0, 1.0)


def test_run_three_beam_test_produces_a_solved_aperture_and_metrics():
    result = hb.run_three_beam_test(CENTER, 5, P, SUB, MET, F, margin=1, mesh_cells_per_subcell=8,
                                     candidate_lw_values_mm=[1.0, 2.8, 4.4])
    assert len(result.assigned) == 9   # 3x3 interior of a 5x5 aperture
    assert result.metrics.n_intended_beams == 3
    for lw in result.assigned.values():
        assert lw in [(1.0, 1.0), (2.8, 2.8), (4.4, 4.4)]
