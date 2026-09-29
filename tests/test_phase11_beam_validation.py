"""Sec. 23, 28: beam-level validation generated entirely from a direct heterogeneous MoM solve."""

import pytest

from src.materials import Metal, Substrate
from src.phase11 import beam_validation as bv
from src.phase11 import patterns as pt

SUB = Substrate("Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3)
MET = Metal()
P = 4.8e-3
F = 12.5e9
CENTER = (2.8, 2.8)
CONTRAST = (4.4, 4.4)


@pytest.mark.slow
def test_solve_aperture_returns_only_interior_points():
    grid_mm = pt.checkerboard(5, CENTER, CONTRAST)
    sg, points = bv.solve_aperture(grid_mm, P, SUB, MET, F, mesh_cells_per_subcell=8)
    n = 5
    for p in points:
        assert bv.INTERIOR_MARGIN <= p.i < n - bv.INTERIOR_MARGIN
        assert bv.INTERIOR_MARGIN <= p.j < n - bv.INTERIOR_MARGIN
    assert len(points) == (n - 2 * bv.INTERIOR_MARGIN) ** 2


@pytest.mark.slow
def test_run_beam_case_produces_far_field_and_metrics():
    grid_mm = pt.checkerboard(5, CENTER, CONTRAST)
    result = bv.run_beam_case("checkerboard", grid_mm, P, SUB, MET, F, mesh_cells_per_subcell=8)
    assert result["far_field"].shape == bv.THETA_GRID_DEG.shape
    assert result["metrics"].n_intended_beams == 1


def test_build_required_beam_cases_covers_the_four_named_patterns():
    grids = bv.build_required_beam_cases(CENTER, CONTRAST, n=5)
    assert set(grids.keys()) == {"smooth_gradient", "abrupt_boundary", "strong_contrast", "checkerboard"}
