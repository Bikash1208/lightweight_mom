"""Phase 12C Task C: coupling extent, grounded (short Phase 11 series). Reuses Phase 11's own
functions verbatim -- ``phase11.coupling_radius.truncate_centered``/``compute_shell_increments``,
``phase11.convergence.classify_convergence``, ``phase11.decay_law.observe_decay``/
``fit_exponential_if_monotonic`` -- every one of them is a pure function of a grid/series of
responses, backing-agnostic. Only the geometry-building and solve call are new: a grounded grid
solve (:mod:`phase12c.grounded_supercell`) instead of Phase 11's ungrounded ``reference_solver.
solve_grid``, and a GRADIENT pattern option (the brief's "include a gradient pattern, not only a
checkerboard") alongside Phase 11's own checkerboard-of-extremes pattern.

Mesh/domain: the same aligned family and per-subcell resolution as :mod:`phase12c.locality_study`
(``L_FAMILY_MM``, ``MESH_CELLS_PER_SUBCELL``) -- reused, not re-derived, so every Phase 12C module
speaks about the same discretization.
"""

from dataclasses import dataclass

from .. import solver as solver_mod
from ..basis_functions import MeshSettings
from ..excitation import PlaneWaveExcitation
from ..phase9.dipole_extraction import local_dipole_moment
from ..phase11.convergence import classify_convergence
from ..phase11.coupling_radius import CouplingRadiusPoint, compute_shell_increments, truncate_centered
from ..phase11.decay_law import fit_exponential_if_monotonic, observe_decay
from . import design_point as dp
from .grounded_supercell import build_grounded_supercell
from .locality_study import L_FAMILY_MM, MESH_CELLS_PER_SUBCELL, floquet_for

SIZES = (1, 3, 5, 7)  # extended only if a size proves feasible in practice (Sec. 'why', module-level cost note)


def solve_grounded_grid(n, lw_grid_m, mesh_cells_per_subcell=MESH_CELLS_PER_SUBCELL):
    """A grounded analogue of ``phase11.reference_solver.solve_grid``: one MoM solve of the whole
    n x n grounded supercell, returning the basis/coefficients needed to read out any sub-cell's
    dipole moment."""
    supercell = build_grounded_supercell(n, dp.PERIOD_M, dp.PERIOD_M, dp.SUBSTRATE, dp.METAL, lw_grid_m)
    mesh_cells = n * mesh_cells_per_subcell
    floquet = floquet_for(mesh_cells)
    excitation = PlaneWaveExcitation(frequency=dp.DESIGN_FREQUENCY_HZ, polarization="x")
    mesh = MeshSettings(mesh_cells, mesh_cells, floquet)
    system = solver_mod.assemble_system(supercell, excitation, mesh)
    currents = solver_mod.solve(system.z, system.v)
    return supercell, system, currents


def _center_response(supercell, system, currents, n):
    mid = n // 2
    cx, cy = supercell.inner.sub_cell_center(mid, mid)
    dipole = local_dipole_moment(
        system.basis, currents, cx, cy, supercell.inner.sub_period_x / 2, supercell.inner.sub_period_y / 2
    )
    return dipole.px


def build_gradient_pattern_mm(center_l_mm, low_l_mm, high_l_mm, max_n, family_mm=L_FAMILY_MM):
    """A smooth L=W gradient (Sec. C's "include a gradient pattern, not only a checkerboard"):
    row-major linear ramp from ``low_l_mm`` to ``high_l_mm`` along x, held fixed along y, with the
    true centre placed at the exact geometric centre position (``max_n`` odd) -- unlike Phase 11's
    checkerboard, a smooth gradient has no reflection symmetry about the centre column in general,
    so the centre value is set explicitly to ``center_l_mm`` rather than relying on parity.

    Every ramp value is SNAPPED to the nearest member of ``family_mm`` (this module's own aligned
    L family) rather than left at its raw linear-interpolation value: an unsnapped intermediate
    value would not be edge-aligned on this module's mesh, silently reintroducing the staircase
    error Task A's own mesh-alignment requirement exists to avoid."""
    import numpy as np

    mid = max_n // 2
    left = np.linspace(low_l_mm, center_l_mm, mid + 1)
    right = np.linspace(center_l_mm, high_l_mm, max_n - mid)
    raw_row = list(left[:-1]) + list(right)

    def snap(v):
        return min(family_mm, key=lambda f: abs(f - v))

    row = [snap(v) for v in raw_row]
    row[mid] = snap(center_l_mm)
    return [[(row[j], row[j]) for j in range(max_n)] for _ in range(max_n)]


def build_checkerboard_pattern_mm(center_l_mm, opposite_l_mm, max_n):
    """Phase 11's own maximally-demanding checkerboard(centre, opposite-extreme), reused for direct
    grounded-vs-ungrounded comparability (the SAME pattern definition Phase 11 itself used)."""
    return [
        [(center_l_mm, center_l_mm) if (i + j) % 2 == 0 else (opposite_l_mm, opposite_l_mm) for j in range(max_n)]
        for i in range(max_n)
    ]


def run_coupling_radius_series(full_grid_mm, sizes=SIZES, mesh_cells_per_subcell=MESH_CELLS_PER_SUBCELL):
    """As ``phase11.coupling_radius.run_coupling_radius_series``, grounded: every ``truncate_centered``
    (reused, unmodified) size in ``sizes`` is solved as its own periodic grounded structure and the
    TRUE centre cell's response read out."""
    max_n = max(sizes)
    points = []
    for n in sorted(sizes):
        grid_mm = truncate_centered(full_grid_mm, n)
        grid_m = [[(l * 1e-3, w * 1e-3) for (l, w) in row] for row in grid_mm]
        supercell, system, currents = solve_grounded_grid(n, grid_m, mesh_cells_per_subcell)
        response = _center_response(supercell, system, currents, n)
        points.append(CouplingRadiusPoint(n, response, system.basis.n_unknowns))
    return points


@dataclass(frozen=True)
class CouplingSeriesReport:
    name: str
    points: tuple
    increments: tuple
    convergence: object
    decay_observation: object
    decay_fit: object


def analyze_coupling_series(name, points):
    """Sec. C: the full Phase-11-style diagnostic bundle for one series, every function reused
    unmodified from ``phase11``."""
    increments = tuple(compute_shell_increments(points))
    convergence = classify_convergence(increments, center_name=name)
    pairs = [(inc.n_to, inc.delta_relative) for inc in increments]
    decay_obs = observe_decay(pairs)
    decay_fit = fit_exponential_if_monotonic(decay_obs)
    return CouplingSeriesReport(name, tuple(points), increments, convergence, decay_obs, decay_fit)
