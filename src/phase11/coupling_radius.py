"""Sec. 5-7 (Phase 11A/B/C): direct heterogeneous MoM reference calculations at growing simulated
extent, and the incremental per-shell contribution derived from them.

The "environment that grows" is one fixed, deterministic, maximally-demanding pattern per center:
a checkerboard between the center's own geometry and the sweep domain's opposite extreme (Sec. 10's
own "strong isolated defect" contrast, repeated -- the most information-dense, least-likely-to-hit-
the-mesh's-own-discretization-floor choice, per Phase 10's measured finding that only strong
contrast reliably clears that floor). The center cell of the LARGEST grid is always the true center
geometry (checkerboard parity places it there automatically); truncating to a smaller centered
n x n sub-block and re-solving it as its own periodic structure is exactly "how much of the true
surrounding pattern is included in the simulated extent," Sec. 6's own framing.
"""

from dataclasses import dataclass

from . import patterns as pt
from .reference_solver import solve_grid

SIZES = (1, 3, 5, 7, 9)          # Sec. 5: the main, full-per-center convergence series, run for
                                  # every center cell.
SUPPLEMENTARY_SIZE = 11          # Sec. 5's own "where computationally practical" qualifier: 11x11
                                  # is run as a spot-check (primary center only, `run.py`) rather
                                  # than in the routine per-center loop -- a single solve at this
                                  # size was measured at ~13-38 s depending on system load, an
                                  # order of magnitude more expensive than the rest of the series
                                  # combined, and a disclosed scoping choice matching Phase 10's own
                                  # documented cost/coverage trade-offs.


def build_growing_pattern_mm(center_lw_mm, max_n):
    """The checkerboard(center, opposite-extreme) pattern at its largest size; smaller sizes are
    centered truncations of this SAME grid (never independently regenerated), so "growing the
    environment" means literally revealing more of one fixed underlying pattern, not solving a
    different structure at each size."""
    from ..phase10.environment_geometry import strong_neighbor_low_high

    low, high = strong_neighbor_low_high(center_lw_mm)
    opposite = (low, low) if abs(center_lw_mm[0] - high) < abs(center_lw_mm[0] - low) else (high, high)
    grid = pt.checkerboard(max_n, a_mm=center_lw_mm, b_mm=opposite)
    mid = max_n // 2
    assert grid[mid][mid] == center_lw_mm, "checkerboard parity must place the true center at the center"
    return grid


def truncate_centered(grid_mm, n):
    """The centered n x n sub-block of ``grid_mm`` (n odd, <= len(grid_mm))."""
    max_n = len(grid_mm)
    if n > max_n:
        raise ValueError(f"truncate_centered: n={n} exceeds the pattern size {max_n}")
    lo = (max_n - n) // 2
    return [row[lo:lo + n] for row in grid_mm[lo:lo + n]]


@dataclass(frozen=True)
class CouplingRadiusPoint:
    n: int
    response: complex
    n_unknowns: int


def run_coupling_radius_series(center_lw_mm, period_m, substrate, metal, frequency_hz,
                                mesh_cells_per_subcell=8, sizes=SIZES, polarization="x"):
    """Sec. 5-6: one `CouplingRadiusPoint` per requested size, all truncations of the SAME
    underlying pattern (``build_growing_pattern_mm``)."""
    max_n = max(sizes)
    full_grid_mm = build_growing_pattern_mm(center_lw_mm, max_n)
    points = []
    for n in sorted(sizes):
        grid_mm = truncate_centered(full_grid_mm, n)
        grid_m = [[(l * 1e-3, w * 1e-3) for (l, w) in row] for row in grid_mm]
        sg = solve_grid(n, period_m, substrate, metal, grid_m, frequency_hz, mesh_cells_per_subcell,
                         polarization=polarization)
        mid = n // 2
        points.append(CouplingRadiusPoint(n, sg.response_at(mid, mid), sg.n_unknowns))
    return points


@dataclass(frozen=True)
class ShellIncrement:
    n_from: int
    n_to: int
    delta_complex: float
    delta_relative: float   # relative to the R_(N-2) magnitude (the "before adding this shell" value)


def compute_shell_increments(points):
    """Sec. 7: DeltaR_N = R_N - R_(N-2), one entry per consecutive pair in ``points`` (already
    sorted by n)."""
    increments = []
    for prev, cur in zip(points, points[1:]):
        delta = cur.response - prev.response
        rel = abs(delta) / abs(prev.response) if abs(prev.response) > 0 else float("inf")
        increments.append(ShellIncrement(prev.n, cur.n, abs(delta), rel))
    return increments
