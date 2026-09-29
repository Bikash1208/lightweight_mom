"""Sec. 8 (Phase 11D): does one neighboring cell's influence on the center depend on which
direction it sits in, and how does that influence change with distance? A single perturbed cell is
placed at increasing distance from a fixed center, independently for each of the 8 compass
directions, on one grid large enough to hold the largest tested distance (Sec. 8's own worked
example: N, S, E, W, NE, NW, SE, SW, "then repeat for increasing distance").
"""

from dataclasses import dataclass

from ..phase10.environment_geometry import strong_neighbor_low_high
from .patterns import uniform
from .reference_solver import solve_grid

DIRECTIONS = ("N", "S", "E", "W", "NE", "NW", "SE", "SW")
_STEP = {
    "N": (-1, 0), "S": (1, 0), "E": (0, 1), "W": (0, -1),
    "NE": (-1, 1), "NW": (-1, -1), "SE": (1, 1), "SW": (1, -1),
}


def _contrast_for(center_lw_mm):
    low, high = strong_neighbor_low_high(center_lw_mm)
    return (low, low) if abs(center_lw_mm[0] - high) < abs(center_lw_mm[0] - low) else (high, high)


@dataclass(frozen=True)
class DirectionalPoint:
    direction: str
    distance: int
    delta_abs: float
    delta_complex: complex


def run_directional_coupling(center_lw_mm, period_m, substrate, metal, frequency_hz, grid_n,
                              distances, mesh_cells_per_subcell=8, contrast_mm=None):
    """One perturbed-cell solve per (direction, distance) pair, plus one shared homogeneous
    baseline. ``grid_n`` must be large enough that ``mid +/- max(distances)`` stays in-bounds for
    every direction tested."""
    mid = grid_n // 2
    if max(distances) > mid:
        raise ValueError(f"run_directional_coupling: grid_n={grid_n} too small for distance {max(distances)}")
    contrast_mm = contrast_mm or _contrast_for(center_lw_mm)

    def solve(grid_mm):
        grid_m = [[(l * 1e-3, w * 1e-3) for (l, w) in row] for row in grid_mm]
        sg = solve_grid(grid_n, period_m, substrate, metal, grid_m, frequency_hz, mesh_cells_per_subcell)
        return sg.response_at(mid, mid)

    baseline = solve(uniform(grid_n, center_lw_mm))

    points = []
    for direction in DIRECTIONS:
        di, dj = _STEP[direction]
        for d in distances:
            grid_mm = uniform(grid_n, center_lw_mm)
            i, j = mid + di * d, mid + dj * d
            grid_mm[i][j] = contrast_mm
            response = solve(grid_mm)
            delta = response - baseline
            points.append(DirectionalPoint(direction, d, abs(delta), delta))
    return baseline, points
