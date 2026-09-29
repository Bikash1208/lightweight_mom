"""Sec. 8, 28: single-neighbor-at-increasing-distance perturbation, all 8 compass directions."""

from src.materials import Metal, Substrate
from src.phase11 import directional_coupling as dc

SUB = Substrate("Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3)
MET = Metal()
P = 4.8e-3
F = 12.5e9
CENTER = (2.8, 2.8)


def test_run_directional_coupling_covers_all_directions_and_distances():
    baseline, points = dc.run_directional_coupling(
        CENTER, P, SUB, MET, F, grid_n=5, distances=(1, 2), mesh_cells_per_subcell=8
    )
    seen = {(p.direction, p.distance) for p in points}
    expected = {(d, dist) for d in dc.DIRECTIONS for dist in (1, 2)}
    assert seen == expected
    assert isinstance(baseline, complex)
    for p in points:
        assert p.delta_abs == abs(p.delta_complex)


def test_run_directional_coupling_rejects_distance_exceeding_grid():
    import pytest

    with pytest.raises(ValueError):
        dc.run_directional_coupling(CENTER, P, SUB, MET, F, grid_n=3, distances=(1, 2), mesh_cells_per_subcell=8)
