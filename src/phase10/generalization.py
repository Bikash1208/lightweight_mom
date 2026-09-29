"""Sec. 30-33: does the surrogate generalize to the spatial STRUCTURE a real hologram actually
has, rather than only to randomly sampled environments? One 7x7 heterogeneous supercell is solved
ONCE (`grid_solver.solve_grid`, the same unmodified MoM pipeline) per scenario; every interior
sub-cell's local dipole moment is then read out from that single solve (`SolvedGrid.dipole_at`) --
this is the correct physical procedure (one self-consistent current distribution over the whole
patch, not independently-solved cells), not a shortcut.

Sec. 33's boundary/edge-effect study is explicitly OUT OF SCOPE for this infinite-periodic-
supercell framework and is reported as such (`NOT_FEASIBLE_WITHIN_INFINITE_PERIODIC_FRAMEWORK`),
per Sec. 33's own warning not to claim an infinite-periodic environment-aware model automatically
validates finite-array edges.
"""

from dataclasses import dataclass

from . import environment_geometry as eg
from .grid_solver import solve_grid

GRID_N = 7   # odd, >= 3 interior rings on each side so every interior cell has a full 3x3
             # neighborhood strictly inside the solved grid (Sec. 30-32); n=7 keeps one solve
             # under ~4 s (measured) while still giving a 5x5 = 25 interior test points.
INTERIOR_RANGE = range(1, GRID_N - 1)   # indices 1..5 for n=7

BOUNDARY_STUDY_STATUS = "NOT_FEASIBLE_WITHIN_INFINITE_PERIODIC_FRAMEWORK"


def _interior_neighbors(grid_mm, i, j):
    from .dataset import ALL_POSITIONS
    offsets = {"N": (-1, 0), "S": (1, 0), "E": (0, 1), "W": (0, -1),
               "NE": (-1, 1), "NW": (-1, -1), "SE": (1, 1), "SW": (1, -1)}
    return {pos: grid_mm[i + di][j + dj] for pos, (di, dj) in offsets.items()}


@dataclass(frozen=True)
class SpatialCase:
    i: int
    j: int
    center_lw_mm: tuple
    neighbors_lw_mm: dict
    px_mom: complex
    px_surrogate: object   # complex, or "OUT_OF_VALIDATED_DOMAIN"


def _build_and_solve(grid_mm, period_m, substrate, metal, frequency_hz, mesh_cells_per_subcell, surrogate):
    grid_m = eg.grid_to_m(grid_mm)
    sg = solve_grid(GRID_N, period_m, substrate, metal, grid_m, frequency_hz, mesh_cells_per_subcell)
    cases = []
    for i in INTERIOR_RANGE:
        for j in INTERIOR_RANGE:
            center_lw_mm = grid_mm[i][j]
            neighbors_lw_mm = _interior_neighbors(grid_mm, i, j)
            px_mom = sg.dipole_at(i, j).px
            from .dataset import EnvironmentSample
            fake_sample = EnvironmentSample(f"gen_{i}_{j}", "n/a", center_lw_mm, neighbors_lw_mm,
                                             "generalization", frequency_hz, px_mom, "generalization")
            px_surrogate = surrogate.predict(fake_sample)
            cases.append(SpatialCase(i, j, center_lw_mm, neighbors_lw_mm, px_mom, px_surrogate))
    return cases


def build_gradient_grid_mm(l_min_mm=1.4, l_max_mm=4.2, w_fixed_mm=2.8):
    """Sec. 31: a smooth L gradient (small -> medium -> large) along x (columns), W held fixed."""
    import numpy as np
    values = np.round(np.linspace(l_min_mm, l_max_mm, GRID_N), 2)
    return [[(float(values[j]), w_fixed_mm) for j in range(GRID_N)] for _ in range(GRID_N)]


def build_abrupt_grid_mm(region_a_mm=(1.4, 1.4), region_b_mm=(4.2, 4.2), boundary_col=3):
    """Sec. 32: Region A (columns < boundary_col) | Region B (columns >= boundary_col), a strong
    step discontinuity, so interior cells at increasing distance from the boundary (0, 1, 2
    columns/periods away) can be compared separately."""
    grid = []
    for _ in range(GRID_N):
        row = [region_a_mm if j < boundary_col else region_b_mm for j in range(GRID_N)]
        grid.append(row)
    return grid


def run_gradient_test(period_m, substrate, metal, frequency_hz, surrogate, mesh_cells_per_subcell=8,
                       **grid_kwargs):
    grid_mm = build_gradient_grid_mm(**grid_kwargs)
    cases = _build_and_solve(grid_mm, period_m, substrate, metal, frequency_hz, mesh_cells_per_subcell, surrogate)
    return grid_mm, cases


def run_abrupt_discontinuity_test(period_m, substrate, metal, frequency_hz, surrogate,
                                   mesh_cells_per_subcell=8, **grid_kwargs):
    grid_mm = build_abrupt_grid_mm(**grid_kwargs)
    cases = _build_and_solve(grid_mm, period_m, substrate, metal, frequency_hz, mesh_cells_per_subcell, surrogate)
    boundary_col = grid_kwargs.get("boundary_col", 3)
    # Sec. 32: label each interior case by its column-distance from the A|B boundary.
    labeled = [(abs(c.j - boundary_col) if c.j < boundary_col else c.j - boundary_col, c) for c in cases]
    return grid_mm, cases, labeled
