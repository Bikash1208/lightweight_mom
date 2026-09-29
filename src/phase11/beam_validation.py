"""Sec. 23 (Phase 11R): beam-level validation, generated ENTIRELY from a real, direct heterogeneous
MoM solution -- never inferred from a local lookup table (Sec. 23's own explicit rule) and never
from Phase 10's now-invalidated surrogate. One large grid is solved once per case (the same
one-solve/many-read-out reuse `reference_solver.solve_grid` already provides); every interior
sub-cell's response becomes that cell's aperture-field contribution, fed into Phase 8E's own
unmodified planar-aperture far-field integral.
"""

from dataclasses import dataclass

import numpy as np

from ..phase8b.object_wave import Beam
from ..phase8e.beam_reconstruction import evaluate_beam_reconstruction, far_field_pattern
from . import patterns as pt
from .reference_solver import solve_grid

APERTURE_N = 7                  # Sec. 23 cases 1-4: keeps each solve <= ~4 s (measured).
INTERIOR_MARGIN = 1              # read out every cell except the outer ring (needs full neighbors
                                  # to be a fair "aperture element," matching Phase 10's own choice).
DEFAULT_BEAM = Beam(theta_deg=20.0, phi_deg=0.0)
THETA_GRID_DEG, PHI_GRID_DEG = np.meshgrid(np.linspace(0.0, 70.0, 36), np.linspace(0.0, 359.0, 72))


@dataclass(frozen=True)
class AperturePoint:
    i: int
    j: int
    center_lw_mm: tuple
    response: complex


def solve_aperture(grid_mm, period_m, substrate, metal, frequency_hz, mesh_cells_per_subcell=8):
    n = len(grid_mm)
    grid_m = [[(l * 1e-3, w * 1e-3) for (l, w) in row] for row in grid_mm]
    sg = solve_grid(n, period_m, substrate, metal, grid_m, frequency_hz, mesh_cells_per_subcell)
    points = []
    for i in range(INTERIOR_MARGIN, n - INTERIOR_MARGIN):
        for j in range(INTERIOR_MARGIN, n - INTERIOR_MARGIN):
            points.append(AperturePoint(i, j, grid_mm[i][j], sg.response_at(i, j)))
    return sg, points


def far_field_from_aperture(points, period_m, frequency_hz, n):
    xs = np.array([-n * period_m / 2 + (p.j + 0.5) * period_m for p in points])
    ys = np.array([-n * period_m / 2 + (p.i + 0.5) * period_m for p in points])
    values = np.array([p.response for p in points], dtype=complex)
    return far_field_pattern(xs, ys, values, frequency_hz, THETA_GRID_DEG, PHI_GRID_DEG, period_m, period_m)


def run_beam_case(name, grid_mm, period_m, substrate, metal, frequency_hz, mesh_cells_per_subcell=8,
                   intended_beams=(DEFAULT_BEAM,)):
    """Sec. 23 cases 1-4 (smooth gradient / abrupt boundary / strong contrast / checkerboard):
    solves ``grid_mm`` once, builds its far field, and evaluates the requested intended beam(s)."""
    n = len(grid_mm)
    sg, points = solve_aperture(grid_mm, period_m, substrate, metal, frequency_hz, mesh_cells_per_subcell)
    far_field = far_field_from_aperture(points, period_m, frequency_hz, n)
    metrics = evaluate_beam_reconstruction(far_field, THETA_GRID_DEG, PHI_GRID_DEG, intended_beams)
    return {"name": name, "grid_mm": grid_mm, "points": points, "far_field": far_field, "metrics": metrics}


def build_required_beam_cases(center_lw_mm, contrast_lw_mm, n=APERTURE_N):
    """Sec. 23 cases 1-4's grids (case 5, the three-beam holographic aperture, is built separately
    in ``holographic_three_beam.py`` since it needs an actual per-cell geometry ASSIGNMENT, not just
    a named spatial pattern)."""
    return {
        "smooth_gradient": pt.smooth_gradient(n, center_lw_mm[0], contrast_lw_mm[0], center_lw_mm[1]),
        "abrupt_boundary": pt.abrupt_boundary(n, center_lw_mm, contrast_lw_mm),
        "strong_contrast": pt.strong_isolated_defect(n, center_lw_mm),
        "checkerboard": pt.checkerboard(n, center_lw_mm, contrast_lw_mm),
    }
