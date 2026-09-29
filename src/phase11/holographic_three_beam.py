"""Sec. 24 (Phase 11S): the three-beam holographic test, using the project's approved beam
specification (theta=30 deg, phi=0/120/240 deg) and the SAME reference-wave/object-wave/phase-
difference-hologram machinery Phase 8B/8E already validated -- reused, not reimplemented. Geometry
per aperture cell is assigned by ``phase11.inverse``'s deterministic (non-ML) isolated-cell search,
never an inverted/learned mapping; the ASSEMBLED heterogeneous aperture is then solved ONCE with the
real, direct MoM solver (Sec. 24: "heterogeneous physical geometry -> direct MoM -> aperture
currents -> far field"), never inferred from a lookup table.

This is a disclosed, bounded (not per-cell-MoM-verified) search: ranking every interior cell by its
ISOLATED-CELL response is itself a real physics computation (zero fitting), and the single MoM
solve of the fully assembled aperture that follows is the actual verification Sec. 24 asks for --
individually re-verifying each of the 25 interior cells with its own separate heterogeneous MoM
solve was judged not worth the ~25x solve-time cost for this diagnostic test (a disclosed scoping
choice, the same kind Phase 9/10 made repeatedly).
"""

from dataclasses import dataclass

import numpy as np

from ..phase8b.object_wave import Beam, combine_beams
from ..phase8b.holographic_modulation import phase_difference_hologram
from ..phase8e.beam_reconstruction import evaluate_beam_reconstruction
from ..phase8e.reference_field_input import spherical_feed_reference_field
from ..phase10.environment_geometry import CENTER_CELLS_MM
from . import beam_validation as bv
from .inverse import _grid_values_mm
from .reference_solver import solve_isolated

THREE_BEAMS = (
    Beam(theta_deg=30.0, phi_deg=0.0),
    Beam(theta_deg=30.0, phi_deg=120.0),
    Beam(theta_deg=30.0, phi_deg=240.0),
)
DEMO_FEED_HEIGHT_M = 0.060


def _representative_isolated_magnitude(period_m, substrate, metal, frequency_hz, mesh_cells_per_subcell=8):
    """A physically-grounded target magnitude (never an arbitrary unit assumption): the mean
    |isolated-cell response| over this project's own 5 named center geometries."""
    mags = []
    for lw_mm in CENTER_CELLS_MM.values():
        resp, _ = solve_isolated(period_m, lw_mm[0] * 1e-3, lw_mm[1] * 1e-3, substrate, metal,
                                  frequency_hz, mesh_cells_per_subcell)
        mags.append(abs(resp))
    return float(np.mean(mags))


def build_holographic_target_grid(n, period_m, substrate, metal, frequency_hz, margin=1,
                                   mesh_cells_per_subcell=8):
    """Sec. 24: reference wave + object wave -> per-cell required complex response, over the
    interior (n - 2*margin) x (n - 2*margin) region of an n x n aperture."""
    xs_1d = np.array([-n * period_m / 2 + (j + 0.5) * period_m for j in range(n)])
    ys_1d = np.array([-n * period_m / 2 + (i + 0.5) * period_m for i in range(n)])
    x_grid, y_grid = np.meshgrid(xs_1d, ys_1d)

    ref_field = spherical_feed_reference_field(x_grid, y_grid, DEMO_FEED_HEIGHT_M, frequency_hz)
    ex_obj, ey_obj = combine_beams(x_grid, y_grid, THREE_BEAMS, frequency_hz)
    required_phase_deg = phase_difference_hologram(ref_field.Ex, ex_obj)

    magnitude = _representative_isolated_magnitude(period_m, substrate, metal, frequency_hz,
                                                     mesh_cells_per_subcell)
    target_complex = magnitude * np.exp(1j * np.radians(required_phase_deg))

    targets = {}
    for i in range(margin, n - margin):
        for j in range(margin, n - margin):
            targets[(i, j)] = complex(target_complex[i, j])
    return targets


def _isolated_response_table(period_m, substrate, metal, frequency_hz, mesh_cells_per_subcell,
                              candidate_lw_values_mm):
    """Every candidate (L, W)'s isolated-cell response is independent of which target cell it will
    be matched against, so it is computed ONCE here and reused for every cell (still zero fitting --
    it is exactly Baseline 1's own physical quantity from Phase 10, just not re-solved redundantly)."""
    table = {}
    for l_mm in candidate_lw_values_mm:
        for w_mm in candidate_lw_values_mm:
            resp, _ = solve_isolated(period_m, l_mm * 1e-3, w_mm * 1e-3, substrate, metal,
                                      frequency_hz, mesh_cells_per_subcell)
            table[(l_mm, w_mm)] = resp
    return table


def assign_geometry_by_isolated_search(targets, response_table):
    """Sec. 21's permitted deterministic method: for each cell, ranks every candidate purely by its
    (already-tabulated) ISOLATED-CELL response distance to that cell's target -- no fitting, no
    learned parameters."""
    assigned = {}
    for pos, target in targets.items():
        best_lw = min(response_table, key=lambda lw: abs(response_table[lw] - target))
        assigned[pos] = best_lw
    return assigned


@dataclass(frozen=True)
class ThreeBeamResult:
    grid_mm: list
    targets: dict
    assigned: dict
    achieved: dict
    metrics: object
    far_field: object


def run_three_beam_test(center_lw_mm, n, period_m, substrate, metal, frequency_hz, margin=1,
                         mesh_cells_per_subcell=8, candidate_lw_values_mm=None):
    """Sec. 24's full pipeline: reference+object wave -> per-cell target -> deterministic
    geometry assignment -> ONE direct heterogeneous MoM solve of the assembled aperture -> far
    field -> beam-reconstruction metrics against the real three-beam target."""
    targets = build_holographic_target_grid(n, period_m, substrate, metal, frequency_hz, margin,
                                             mesh_cells_per_subcell)
    values = candidate_lw_values_mm or _grid_values_mm()
    response_table = _isolated_response_table(period_m, substrate, metal, frequency_hz,
                                               mesh_cells_per_subcell, values)
    assigned = assign_geometry_by_isolated_search(targets, response_table)

    grid_mm = [[center_lw_mm] * n for _ in range(n)]
    for (i, j), lw in assigned.items():
        grid_mm[i][j] = lw

    sg, points = bv.solve_aperture(grid_mm, period_m, substrate, metal, frequency_hz, mesh_cells_per_subcell)
    achieved = {(p.i, p.j): p.response for p in points}
    far_field = bv.far_field_from_aperture(points, period_m, frequency_hz, n)
    metrics = evaluate_beam_reconstruction(far_field, bv.THETA_GRID_DEG, bv.PHI_GRID_DEG, THREE_BEAMS)

    return ThreeBeamResult(grid_mm, targets, assigned, achieved, metrics, far_field)
