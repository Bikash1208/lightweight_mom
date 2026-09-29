"""Sec. 6-9, 12: center-cell definitions, the validated geometry domain, and reproducible neighbor
contrast levels, all built on top of the existing 18x18, 1.0-4.4 mm / 0.2 mm sweep grid
(``library_schema.py``, unmodified) so every geometry used here is one the frozen solver/database
has already exercised at production scale.
"""

from dataclasses import dataclass

DOMAIN_L_MIN_MM, DOMAIN_L_MAX_MM = 1.0, 4.4
DOMAIN_W_MIN_MM, DOMAIN_W_MAX_MM = 1.0, 4.4
DOMAIN_STEP_MM = 0.2

# Sec. 6: five center cells spanning near-small, mid-range, near-large and anisotropic geometry,
# all exact points on the existing 0.2 mm sweep grid (verified in tests).
CENTER_CELLS_MM = {
    "C1": (1.4, 1.4),   # near-small
    "C2": (2.0, 2.8),   # anisotropic, mid-small
    "C3": (2.8, 2.8),   # mid-range (Phase 9's own default center)
    "C4": (3.6, 2.8),   # anisotropic, mid-large
    "C5": (4.2, 4.2),   # near-large
}
PRIMARY_CENTER = "C3"

WEAK_DELTA_MM = 0.2     # Sec. 7
MODERATE_DELTA_MM = 0.6  # Sec. 8
STRONG_LW_MM = (1.0, 4.4)  # Sec. 9: the sweep's own domain extremes


def in_domain(l_mm, w_mm, tol=1e-9):
    return (DOMAIN_L_MIN_MM - tol <= l_mm <= DOMAIN_L_MAX_MM + tol
            and DOMAIN_W_MIN_MM - tol <= w_mm <= DOMAIN_W_MAX_MM + tol)


@dataclass(frozen=True)
class PerturbedNeighbor:
    lw_mm: tuple
    direction: str   # "+" or "-", whichever stayed in-domain
    delta_mm: float


def perturbed_neighbor(center_lw_mm, delta_mm):
    """Sec. 7-8: ``center +/- delta`` on both L and W, preferring "+"; falls back to "-" if "+"
    would leave the validated domain, and raises if neither direction is valid (Sec. 7's own
    proviso, "subject to valid geometry limits"). Never silently clamps or extrapolates."""
    l0, w0 = center_lw_mm
    plus = (round(l0 + delta_mm, 6), round(w0 + delta_mm, 6))
    minus = (round(l0 - delta_mm, 6), round(w0 - delta_mm, 6))
    if in_domain(*plus):
        return PerturbedNeighbor(plus, "+", delta_mm)
    if in_domain(*minus):
        return PerturbedNeighbor(minus, "-", delta_mm)
    raise ValueError(f"perturbed_neighbor: neither {center_lw_mm}+-{delta_mm} stays in the "
                      f"validated {DOMAIN_L_MIN_MM}-{DOMAIN_L_MAX_MM} mm domain")


def weak_neighbor(center_lw_mm):
    return perturbed_neighbor(center_lw_mm, WEAK_DELTA_MM)


def moderate_neighbor(center_lw_mm):
    return perturbed_neighbor(center_lw_mm, MODERATE_DELTA_MM)


def strong_neighbor_low_high(center_lw_mm):
    """Sec. 9: the two sweep extremes are always in-domain regardless of center, since they
    define the domain's own boundary."""
    return STRONG_LW_MM[0], STRONG_LW_MM[1]


# --------------------------------------------------------------------- 3x3 grid construction

# Row-major 3x3 index layout used throughout this module (matches ``phase9.local_periodicity``):
#   (0,0) NW  (0,1) N  (0,2) NE
#   (1,0) W   (1,1) C  (1,2) E
#   (2,0) SW  (2,1) S  (2,2) SE
CARDINAL_POSITIONS = {"N": (0, 1), "E": (1, 2), "S": (2, 1), "W": (1, 0)}
DIAGONAL_POSITIONS = {"NE": (0, 2), "NW": (0, 0), "SE": (2, 2), "SW": (2, 0)}


def build_3x3_environment(center_lw_mm, cardinal_mm=None, diagonal_mm=None):
    """Sec. 8, 12: a full 3x3 neighbor grid (mm), ``center_lw_mm`` at (1,1). ``cardinal_mm``/
    ``diagonal_mm`` are dicts keyed by "N"/"E"/"S"/"W"/"NE"/"NW"/"SE"/"SW"; any key not given
    defaults to the center geometry itself. Sec. 12's "5-cell model" is exactly this 3x3 grid
    with every diagonal position left at its default (== center) -- a real, well-defined MoM
    configuration (corner sub-cells identical to the center), not an approximation of a smaller
    grid; it isolates the four cardinal neighbors' effect from any diagonal contribution, which
    is the whole point of comparing it against the full 3x3 ("9-cell") grid in Sec. 16."""
    cardinal_mm = cardinal_mm or {}
    diagonal_mm = diagonal_mm or {}
    grid = [[center_lw_mm] * 3 for _ in range(3)]
    for name, (i, j) in CARDINAL_POSITIONS.items():
        grid[i][j] = cardinal_mm.get(name, center_lw_mm)
    for name, (i, j) in DIAGONAL_POSITIONS.items():
        grid[i][j] = diagonal_mm.get(name, center_lw_mm)
    return grid


def build_5cell_environment(center_lw_mm, cardinal_mm):
    """Sec. 12: the 5-cell (center + 4 cardinal neighbors) model, represented as a 3x3 grid whose
    diagonals equal the center (see ``build_3x3_environment`` docstring)."""
    return build_3x3_environment(center_lw_mm, cardinal_mm=cardinal_mm, diagonal_mm=None)


def build_9cell_environment(center_lw_mm, cardinal_mm=None, diagonal_mm=None):
    """Sec. 13: the full 3x3 ("9-cell") neighborhood, cardinals and diagonals independent."""
    return build_3x3_environment(center_lw_mm, cardinal_mm=cardinal_mm, diagonal_mm=diagonal_mm)


def grid_to_m(grid_mm):
    return [[(l * 1e-3, w * 1e-3) for (l, w) in row] for row in grid_mm]


# ---------------------------------------------------------------- Sec. 17: long-range (5x5) grid

# distance-2 (two periods) cardinal positions in a 0-indexed 5x5 grid, center at (2, 2).
LONG_RANGE_POSITIONS = {"N2": (0, 2), "S2": (4, 2), "W2": (2, 0), "E2": (2, 4)}


def build_5x5_long_range_environment(center_lw_mm, long_range_mm):
    """Sec. 17: a 5x5 grid identical to the center everywhere except the four positions exactly
    two periods from the center (cardinal directions only); isolates whether coupling has decayed
    by the second ring, with the immediate (3x3) neighborhood held homogeneous so it cannot
    contaminate the measurement."""
    grid = [[center_lw_mm] * 5 for _ in range(5)]
    for name, (i, j) in LONG_RANGE_POSITIONS.items():
        grid[i][j] = long_range_mm.get(name, center_lw_mm)
    return grid
