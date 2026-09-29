"""Sec. 10: the required controlled spatial-pattern generators, for an arbitrary odd grid size n,
built on the same validated 1.0-4.4 mm domain ``phase10.environment_geometry`` already defines.
Every pattern is a deterministic function of (n, center geometry) -- no randomness anywhere in this
module (Sec. 10's "controlled environments rather than arbitrary random data").
"""

from ..phase10.environment_geometry import strong_neighbor_low_high

A_DEFAULT_MM = (1.4, 1.4)
B_DEFAULT_MM = (4.2, 4.2)


def uniform(n, geom_mm):
    """Sec. 10: every cell identical (the homogeneous control)."""
    return [[geom_mm] * n for _ in range(n)]


def single_perturbation(n, base_mm, perturb_mm, i=None, j=None):
    """Sec. 10: one cell (default: the center) differs from an otherwise uniform surround."""
    grid = uniform(n, base_mm)
    i = n // 2 if i is None else i
    j = n // 2 if j is None else j
    grid[i][j] = perturb_mm
    return grid


def smooth_gradient(n, l_min_mm=A_DEFAULT_MM[0], l_max_mm=B_DEFAULT_MM[0], w_fixed_mm=2.8, axis="x"):
    """Sec. 10: A A B B C C D-style smooth gradient along one axis, W held fixed."""
    import numpy as np

    values = [float(v) for v in np.round(np.linspace(l_min_mm, l_max_mm, n), 3)]
    if axis == "x":
        return [[(values[j], w_fixed_mm) for j in range(n)] for _ in range(n)]
    return [[(values[i], w_fixed_mm) for _ in range(n)] for i in range(n)]


def abrupt_boundary(n, region_a_mm=A_DEFAULT_MM, region_b_mm=B_DEFAULT_MM, boundary_col=None):
    """Sec. 10: AAAA|BBBB, a strong step discontinuity along columns."""
    boundary_col = n // 2 if boundary_col is None else boundary_col
    return [[region_a_mm if j < boundary_col else region_b_mm for j in range(n)] for _ in range(n)]


def checkerboard(n, a_mm=A_DEFAULT_MM, b_mm=B_DEFAULT_MM):
    """Sec. 10: A B A B / B A B A -- the standard 2-D checkerboard, cell (i,j) = A iff i+j is even."""
    return [[(a_mm if (i + j) % 2 == 0 else b_mm) for j in range(n)] for i in range(n)]


def alternating(n, a_mm=A_DEFAULT_MM, b_mm=B_DEFAULT_MM):
    """Sec. 10: ABABABA / BABABAB row-and-column-alternating stripes -- verified (Sec. 10's ASCII
    diagram) to be the exact same cell(i, j) = A iff (i + j) even pattern as ``checkerboard``, so
    this simply calls it; kept as a separate, explicitly named function because Sec. 10 lists the
    two as distinct required test cases and callers/tests should be able to name either by its own
    section name."""
    return checkerboard(n, a_mm, b_mm)


def strong_isolated_defect(n, base_mm, i=None, j=None):
    """Sec. 10: "the largest permitted geometry contrast" -- the sweep's own opposite extreme from
    ``base_mm`` (whichever of 1.0/4.4 mm is farther from it), as a single-cell defect."""
    low, high = strong_neighbor_low_high(base_mm)
    defect = (low, low) if abs(base_mm[0] - high) < abs(base_mm[0] - low) else (high, high)
    return single_perturbation(n, base_mm, defect, i, j)


def single_cell_discontinuity(n, base_mm, contrast_mm, boundary_col=None):
    """Sec. 11 (Phase 11G): the narrowest possible discontinuity -- one column differs, otherwise
    uniform (distinct from ``abrupt_boundary``'s half-and-half split)."""
    boundary_col = n // 2 if boundary_col is None else boundary_col
    grid = uniform(n, base_mm)
    for i in range(n):
        grid[i][boundary_col] = contrast_mm
    return grid


def multi_cell_discontinuity(n, base_mm, contrast_mm, width=2, boundary_col=None):
    """Sec. 11 (Phase 11G): a ``width``-column band of contrast geometry."""
    boundary_col = n // 2 if boundary_col is None else boundary_col
    grid = uniform(n, base_mm)
    for i in range(n):
        for w in range(width):
            j = boundary_col + w
            if 0 <= j < n:
                grid[i][j] = contrast_mm
    return grid
