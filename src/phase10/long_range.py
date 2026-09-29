"""Sec. 17: does coupling decay by the second ring (two periods)? A 5x5 grid identical to the
center everywhere except one distance-2 cardinal position at a time
(`environment_geometry.build_5x5_long_range_environment`)."""

from dataclasses import dataclass

from . import environment_geometry as eg
from .grid_solver import solve_grid

LONG_RANGE_DIRECTIONS = ("N2", "S2", "W2", "E2")


@dataclass(frozen=True)
class LongRangeCase:
    direction: str
    px: complex
    delta_s: complex
    delta_s_abs: float
    n_unknowns: int


def solve_5x5_homogeneous(center_lw_mm, period_m, substrate, metal, frequency_hz, mesh_cells_per_subcell=8):
    """The n=5 zero-perturbation reference (Sec. 17's own baseline at the SAME mesh/Floquet
    resolution as the long-range test cases -- a 3x3 homogeneous solve is not numerically
    identical to a 5x5 homogeneous solve of the same infinite-periodic structure, so each grid
    size needs its own matched baseline; Phase 9's homogeneous-exact-match proof already
    established both reproduce the isolated cell almost exactly, so either is physically valid --
    this just keeps the *numerical* comparison apples-to-apples)."""
    grid_mm = [[center_lw_mm] * 5 for _ in range(5)]
    grid_m = eg.grid_to_m(grid_mm)
    sg = solve_grid(5, period_m, substrate, metal, grid_m, frequency_hz, mesh_cells_per_subcell)
    return sg.dipole_at(2, 2).px


def run_long_range_test(
    center_lw_mm, contrast_lw_mm, s0_5x5, period_m, substrate, metal, frequency_hz, mesh_cells_per_subcell=8
):
    results = {}
    for d in LONG_RANGE_DIRECTIONS:
        grid_mm = eg.build_5x5_long_range_environment(center_lw_mm, {d: contrast_lw_mm})
        grid_m = eg.grid_to_m(grid_mm)
        sg = solve_grid(5, period_m, substrate, metal, grid_m, frequency_hz, mesh_cells_per_subcell)
        px = sg.dipole_at(2, 2).px
        delta = px - s0_5x5
        results[d] = LongRangeCase(d, px, delta, abs(delta), sg.n_unknowns)
    return results
