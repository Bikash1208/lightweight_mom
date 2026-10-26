"""Automatic mesh / Floquet validation.

Two layers, both reported (never silently "fixed"): the solver keeps the discretization the user
configured and records how adequate it was.

A priori checks (cheap, every cell and frequency; no solve)
    n_unknowns        rooftop unknowns; 0 means the mesh cannot represent the pattern at all
                      (``RooftopBasis`` raises) and the cell is invalid, never "no metal".
    metal_area_error  |A_mesh - A_true| / A_true, A_mesh = (#metal cells) dx dy from the centroid
                      classification. 0 when the pattern edges lie on mesh lines; otherwise the
                      staircase error. Limit ``max_area_error`` (a heuristic, not a theorem).
    cells_per_wavelength  lambda_sub / max(dx, dy) with lambda_sub = lambda0 / Re(sqrt(eps_c)).
                      Limit ``min_cells_per_wavelength`` (a rule of thumb).
    specular_only     k0 < 2 pi / max(Px, Py): no grating order propagates in air, so the
                      specular (0,0) S-matrix is the complete scattering description. Above it S
                      is still passive but leaves out power radiated into higher orders.

Mesh alignment and duplicate-geometry detection (Phase 12B; cheap, no solve)
    A mesh cell is metal iff its centroid is metal (CONVENTIONS.md), so a feature edge that does
    not sit exactly on a mesh line is staircased. ``feature_edges_aligned``/``smallest_aligned_cells``
    answer "does this edge sit on a mesh line" / "what is the smallest mesh that puts every edge
    of this whole family of sizes on a line", derived from the mesh-line positions
    ``-P/2 + k*P/N`` (Sec. below), not guessed. Two geometries that discretize to the *same*
    metal mask on a given mesh solve to bit-identical Z, V and hence S, regardless of how
    different their continuous geometry is; ``metal_mask_key``/``duplicate_geometry_groups``
    detect this before any solve (``sweep.plan_sweep``'s own ``mesh_duplicate_groups``).

Refinement study (a posteriori; extra solves; opt-in)
    With factor q (integer >= 2), base mesh (nx, ny) and Floquet order M, at one frequency:

        A = solve(nx,   ny,   q M)        B = solve(q nx, q ny, q M)

        floquet_error = max |S_A - S_base|      (truncation of the Floquet sum)
        mesh_error    = max |S_B - S_A|         (mesh size, at the finer Floquet order)

    with max over the 16 complex entries of the 4x4 S-matrix (absolute; |S| <= 1, so it is also
    the scale of the phase error in radians for |S| ~ 1). These are *successive differences*, not
    error bars: for a first-order-convergent quantity the error of the base mesh is larger than
    the successive difference (about 2x for the measured rate, MOM_FORMULATION.md). The
    tolerance is a policy the user sets; passing it means "the result moved by less than this
    when the discretization was refined", not "the result is this accurate".
"""

import math

import numpy as np

from .basis_functions import MeshSettings, RooftopBasis
from .constants import C0
from .frequency import free_space_wavelength, free_space_wavenumber
from .sparameters import compute_s_parameters

# Policy defaults. These are numerical-quality heuristics chosen by this project, not physical
# constants; the values actually used are recorded in every database.
DEFAULT_MIN_CELLS_PER_WAVELENGTH = 10.0
DEFAULT_MAX_AREA_ERROR = 0.05
DEFAULT_REFINEMENT_TOLERANCE = 0.02
DEFAULT_REFINEMENT_FACTOR = 2


def check_cell(unit_cell, mesh, max_area_error):
    """Frequency-independent checks. Raises ValueError if the mesh has no unknowns."""
    basis = RooftopBasis(unit_cell, mesh.cells_x, mesh.cells_y)
    mesh_area = float(basis.metal_mask.sum()) * basis.dx * basis.dy
    area_error = abs(mesh_area - unit_cell.metal_area) / unit_cell.metal_area
    return {
        "n_unknowns": int(basis.n_unknowns),
        "metal_area_error": float(area_error),
        "area_ok": bool(area_error <= max_area_error),
    }


def cells_per_wavelength(unit_cell, mesh, frequency):
    """lambda_sub / (largest mesh edge), lambda_sub = lambda0 / Re(refractive index)."""
    lam_sub = free_space_wavelength(frequency) / unit_cell.substrate.refractive_index.real
    lat = unit_cell.lattice
    return float(lam_sub / max(lat.period_x / mesh.cells_x, lat.period_y / mesh.cells_y))


def is_specular_only(unit_cell, frequency):
    """True iff no Floquet order other than (0,0) propagates in air at normal incidence.

    The nearest nonzero orders are (+-1, 0) and (0, +-1), cut off at k0 = 2 pi / P, so the
    condition is k0 < 2 pi / max(Px, Py) (equivalently f < c0 / max(Px, Py))."""
    lat = unit_cell.lattice
    return bool(free_space_wavenumber(frequency) < min(lat.reciprocal_x, lat.reciprocal_y))


def grating_onset_frequency(unit_cell):
    """f = c0 / max(Px, Py) [Hz]: first frequency at which a grating order propagates."""
    lat = unit_cell.lattice
    return C0 / max(lat.period_x, lat.period_y)


def grating_onset_frequency_oblique(unit_cell, theta):
    """f [Hz]: grating onset for a plane wave incident at angle ``theta`` [rad] from the surface
    normal (Phase 12B Task D), in the azimuth most restrictive for the larger period axis: the
    forward diffracted order's transverse wavenumber is ``kx0 + 2*pi/P = k0*sin(theta) + 2*pi/P``;
    it starts to propagate (equals ``-k0``, the worst-case, most restrictive sign) when
    ``P*(1 + sin(theta)) = lambda``, i.e. ``f = c0 / (P*(1 + sin(theta)))``. ``theta = 0`` reduces
    exactly to ``grating_onset_frequency``."""
    lat = unit_cell.lattice
    return C0 / (max(lat.period_x, lat.period_y) * (1.0 + math.sin(theta)))


def check_frequency(unit_cell, mesh, frequency, min_cells_per_wavelength):
    cpw = cells_per_wavelength(unit_cell, mesh, frequency)
    return {
        "cells_per_wavelength": cpw,
        "wavelength_ok": bool(cpw >= min_cells_per_wavelength),
        "specular_only": is_specular_only(unit_cell, frequency),
    }


# ------------------------------------------------------- mesh alignment (Phase 12B)

# Tolerance, in units of one mesh cell (dx or dy), for "this position lies on a mesh line": a
# policy value like the other a priori thresholds above, generous compared to the double-precision
# rounding (~1e-15) an exact-decimal mm value converted through units.parse_length actually has.
ALIGNMENT_RTOL = 1e-9

# Upper bound on the linear search in ``smallest_aligned_cells``: a mesh finer than this is not
# practical (P/4096 is far below any feature size this project ever meshes), so exceeding it means
# "no aligned mesh in any practical range", reported as None, not silently searched forever.
MAX_ALIGNMENT_SEARCH_CELLS = 4096


def is_position_mesh_aligned(position, period, cells, rtol=ALIGNMENT_RTOL):
    """True iff ``position`` sits on a mesh line of a ``cells``-cell mesh over ``period``.

    Mesh lines are at ``-period/2 + k*period/cells`` for integer ``k = 0..cells``
    (``basis_functions.RooftopBasis``'s own cell edges); ``position`` is on one iff
    ``(position + period/2)/(period/cells)`` is an integer, within ``rtol`` of one mesh cell.
    """
    d = period / cells
    idx = (position + period / 2.0) / d
    return abs(idx - round(idx)) <= rtol


def feature_edges_aligned(period, cells, length, center, rtol=ALIGNMENT_RTOL):
    """True iff both edges (``center -+ length/2``) of a feature lie on mesh lines of a
    ``cells``-cell mesh over ``period`` (one axis)."""
    return is_position_mesh_aligned(center - length / 2.0, period, cells, rtol) and is_position_mesh_aligned(
        center + length / 2.0, period, cells, rtol
    )


def geometry_edges_aligned(geometry, lattice, mesh, rtol=ALIGNMENT_RTOL):
    """True iff ``geometry``'s feature edges (x and y) lie on ``mesh``'s lines. Exact only for an
    axis-aligned (``rotation == 0``) rectangle; any rotated feature's edges are not, in general,
    parallel to the mesh at all, so this returns False for nonzero rotation without attempting a
    rotated-edge alignment test. A geometry with its own ``edges_aligned(lattice, mesh, rtol)`` (Phase 20 R1:
    ``composite_geometry.CompositeGeometry``) answers for itself."""
    if hasattr(geometry, "edges_aligned"):
        return geometry.edges_aligned(lattice, mesh, rtol)
    if geometry.rotation != 0.0:
        return False
    return feature_edges_aligned(
        lattice.period_x, mesh.cells_x, geometry.length, geometry.center_x, rtol
    ) and feature_edges_aligned(lattice.period_y, mesh.cells_y, geometry.width, geometry.center_y, rtol)


def smallest_aligned_cells(period, features, max_cells=MAX_ALIGNMENT_SEARCH_CELLS, rtol=ALIGNMENT_RTOL):
    """Smallest integer ``cells`` (>= 2) such that every ``(length, center)`` pair in ``features``
    (one axis) is edge-aligned on a ``cells``-cell mesh over ``period``, or ``None`` if none exists
    up to ``max_cells`` -- an honest negative result (report it; never silently pick an unaligned
    mesh or search past a documented bound). Example: ``period=4.8e-3``, a family of lengths
    1.0-4.4 mm in 0.2 mm steps, all centred at 0: the smallest aligned mesh is 48 cells (0.1 mm).
    """
    if max_cells < 2:
        raise ValueError(f"max_cells must be >= 2, got {max_cells!r}")
    for n in range(2, max_cells + 1):
        if all(feature_edges_aligned(period, n, length, center, rtol) for length, center in features):
            return n
    return None


def metal_mask_key(unit_cell, mesh):
    """A hashable signature of the discretized metal pattern (``RooftopBasis.metal_mask``) on
    ``mesh``. Two unit cells with the same key solve to bit-identical Z, V and hence S on that
    mesh (Phase 12B), regardless of any continuous-geometry difference between them."""
    basis = RooftopBasis(unit_cell, mesh.cells_x, mesh.cells_y)
    return basis.metal_mask.tobytes()


def duplicate_geometry_groups(items):
    """``items``: iterable of ``(key, mask_key)`` (``key`` identifies the item, e.g. a cell index
    or cell_id; ``mask_key`` from ``metal_mask_key``). Returns a tuple of tuples of ``key``s
    sharing the same ``mask_key``, one tuple per group of size >= 2 (a singleton is not a
    duplicate and is omitted), in first-seen order."""
    groups = {}
    for key, mask_key in items:
        groups.setdefault(mask_key, []).append(key)
    return tuple(tuple(v) for v in groups.values() if len(v) > 1)


def max_abs_difference(a, b):
    """max over entries of |a - b| for two complex arrays of equal shape."""
    a, b = np.asarray(a, dtype=complex), np.asarray(b, dtype=complex)
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch {a.shape} vs {b.shape}")
    return float(np.max(np.abs(a - b)))


def refinement_study(unit_cell, planes, mesh, frequency, base_matrix, factor):
    """(floquet_error, mesh_error) at one frequency; see the module docstring.

    ``base_matrix`` (and the comparison) is the full 4x4 S-matrix, or (Phase 12A) just the (2, 2)
    S11 block for a grounded stack (``unit_cell.backing == "pec"``, no port 2 to compare)."""
    if isinstance(factor, bool) or not isinstance(factor, (int, np.integer)) or factor < 2:
        raise ValueError(f"refinement factor must be an integer >= 2, got {factor!r}")
    fine_order = factor * mesh.floquet_order
    a = compute_s_parameters(
        unit_cell,
        frequency,
        MeshSettings(mesh.cells_x, mesh.cells_y, fine_order),
        planes,
        warn_grating=False,
    )
    b = compute_s_parameters(
        unit_cell,
        frequency,
        MeshSettings(factor * mesh.cells_x, factor * mesh.cells_y, fine_order),
        planes,
        warn_grating=False,
    )
    a_m = a.s11 if a.grounded else a.matrix()
    b_m = b.s11 if b.grounded else b.matrix()
    return (
        max_abs_difference(a_m, base_matrix),
        max_abs_difference(b_m, a_m),
    )
