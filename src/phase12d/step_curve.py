"""Task A: fill in the local-periodicity step-size curve with a finer L grid, and report the
design-relevant reflection-phase/amplitude error alongside Phase 9's own dipole-moment metric.

Mesh cost note (measured, not assumed -- see ``PHASE12D_LOCAL_PERIODICITY.md`` Sec. 2 for the full
timing table). The brief asks for a 0.05-0.1 mm L grid near resonance, aligned via the Phase 12B
helper. A 3x3 EMBEDDED supercell's cost is dominated by its dense Z matrix
(``n_unknowns**2 * 16`` bytes) and scales roughly as the 6th power of the per-subcell mesh count
once floquet order tracks it (Phase 12C's own finding, confirmed again here): measured 76.8 s at
24 cells/subcell, 196.4 s at 28, and a 3x3 at 40 cells/subcell (the resolution a 0.05 mm-aligned L
family near this design's ~6 mm resonance would actually need is far finer even than that) did not
finish in 14+ minutes and was killed -- the SAME infeasibility Phase 12C's own Task A already found
for n=7 supercells at ITS 20-cell/subcell mesh. **24 cells/subcell (aligned step = 2*period/24 =
0.667 mm) is used here as the largest embedded-solve resolution this project can afford for a
multi-case study** -- finer than Phase 12C's own 20-cell/0.8 mm family, but still far coarser than
the brief's 0.05-0.1 mm aspiration, which is genuinely computationally unreachable for an EMBEDDED
(supercell) solve at this project's validated floquet convention within a practical time budget.
The isolated-cell (single-patch) curve, by contrast, is cheap at any resolution used here.
"""

from dataclasses import dataclass

from ..discretization import feature_edges_aligned
from ..phase import wrapped_difference_deg
from ..phase12c import design_point as dp
from ..phase12c.grounded_supercell import (
    build_grounded_supercell,
    solve_grounded_supercell_center,
    solve_isolated_grounded_cell,
)
from ..phase12c.locality_study import (
    CurvePoint,
    build_1d_row_grid,
    build_2d_gradient_grid,
    classify_centers,
)
from .reflection_metric import reflection_phase_amplitude_error

EMBEDDED_MESH_CELLS_PER_SUBCELL = 24
_UNIT_MM = 2.0 * (dp.PERIOD_M * 1e3) / EMBEDDED_MESH_CELLS_PER_SUBCELL   # aligned step, mm
L_FAMILY_MM = tuple(round(_UNIT_MM * k, 10) for k in range(2, 12))       # ~1.33-7.33 mm

for _l in L_FAMILY_MM:
    assert feature_edges_aligned(dp.PERIOD_M, EMBEDDED_MESH_CELLS_PER_SUBCELL, _l * 1e-3, 0.0), _l


def floquet_for(total_cells):
    """The validated ``2 * total_cells`` convention (Phase 9/12C), reused: Task A's embedded
    resolution is chosen to keep THIS affordable, not the other way round."""
    return 2 * total_cells


def build_reference_curve(mesh_cells_per_subcell=EMBEDDED_MESH_CELLS_PER_SUBCELL, l_values_mm=L_FAMILY_MM):
    from ..basis_functions import MeshSettings, RooftopBasis
    from ..sparameters import compute_s_parameters

    floquet = floquet_for(mesh_cells_per_subcell)
    points = []
    for l_mm in l_values_mm:
        cell = dp.unit_cell_at(l_mm * 1e-3)
        mesh = MeshSettings(mesh_cells_per_subcell, mesh_cells_per_subcell, floquet)
        s = compute_s_parameters(cell, dp.DESIGN_FREQUENCY_HZ, mesh, dp.REFERENCE_PLANES,
                                  warn_grating=False, warn_surface_wave=False)
        n_unknowns = RooftopBasis(cell, mesh_cells_per_subcell, mesh_cells_per_subcell).n_unknowns
        points.append(CurvePoint(l_mm, s.magnitude("S11_xx"), s.phase_deg("S11_xx"), n_unknowns))
    return tuple(points)


@dataclass(frozen=True)
class StepCase:
    name: str
    kind: str              # "1d_row" | "2d_gradient"
    center_lw_mm: tuple
    grid_mm: tuple
    k_separation: int       # how many L_FAMILY_MM steps away the neighbour(s) are
    achieved_step_deg: float
    isolated_px: complex
    embedded_px: complex
    dipole_relative_complex_error: float   # Phase 9's own metric, for comparison
    reflection_phase_error_deg: float      # the design-relevant metric (Task A's own contribution)
    reflection_magnitude_error_db: float
    n_unknowns_isolated: int
    n_unknowns_embedded: int


def _kth_neighbor(curve, center, k, direction):
    idx = curve.index(center)
    j = idx + k if direction == "higher_l" else idx - k
    return curve[j] if 0 <= j < len(curve) else None


def _solve_case(name, kind, center, grid_mm, achieved_step_deg, k_separation,
                 mesh_cells_per_subcell=EMBEDDED_MESH_CELLS_PER_SUBCELL):
    from ..phase9.local_periodicity import compute_lpa_error

    n = len(grid_mm)
    floquet = floquet_for(n * mesh_cells_per_subcell)
    center_l_m = center.l_mm * 1e-3
    isolated = solve_isolated_grounded_cell(
        dp.PERIOD_M, dp.PERIOD_M, center_l_m, center_l_m, dp.SUBSTRATE, dp.METAL,
        dp.DESIGN_FREQUENCY_HZ, mesh_cells_per_subcell, floquet_for(mesh_cells_per_subcell),
    )
    grid_m = [[(l * 1e-3, w * 1e-3) for (l, w) in row] for row in grid_mm]
    supercell = build_grounded_supercell(n, dp.PERIOD_M, dp.PERIOD_M, dp.SUBSTRATE, dp.METAL, grid_m)
    embedded = solve_grounded_supercell_center(
        supercell, (center_l_m, center_l_m), dp.DESIGN_FREQUENCY_HZ, mesh_cells_per_subcell, floquet,
    )
    _, _, _, dipole_rel = compute_lpa_error(isolated.dipole, embedded.dipole)
    phase_err, mag_db_err, _, _ = reflection_phase_amplitude_error(
        isolated.dipole.px, embedded.dipole.px, dp.PERIOD_M, dp.SUBSTRATE, dp.DESIGN_FREQUENCY_HZ,
        dp.REFERENCE_PLANES,
    )
    return StepCase(
        name=name, kind=kind, center_lw_mm=(center.l_mm, center.l_mm),
        grid_mm=tuple(tuple(row) for row in grid_mm), k_separation=k_separation,
        achieved_step_deg=achieved_step_deg, isolated_px=isolated.dipole.px, embedded_px=embedded.dipole.px,
        dipole_relative_complex_error=dipole_rel, reflection_phase_error_deg=phase_err,
        reflection_magnitude_error_db=mag_db_err, n_unknowns_isolated=isolated.n_unknowns,
        n_unknowns_embedded=embedded.n_unknowns,
    )


def build_step_cases(curve, center_names_and_points, max_k=3,
                      mesh_cells_per_subcell=EMBEDDED_MESH_CELLS_PER_SUBCELL):
    """For each named center, both 1-D-row and 2-D-gradient cases at every available
    k=1..max_k family-step separation in whichever of the two (east/west, i.e. higher/lower L)
    directions the family supports -- an honest sweep over whatever achieved steps this family
    and mesh can actually produce (brief: "report the achieved steps honestly"), not a forced set
    of target values."""
    cases = []
    for center_name, center in center_names_and_points:
        for k in range(1, max_k + 1):
            for direction in ("higher_l", "lower_l"):
                neighbor = _kth_neighbor(curve, center, k, direction)
                if neighbor is None:
                    continue
                step_deg = wrapped_difference_deg(neighbor.s11_phase_deg, center.s11_phase_deg)
                row_grid = build_1d_row_grid(
                    (center.l_mm, center.l_mm),
                    (neighbor.l_mm, neighbor.l_mm) if direction == "lower_l" else (center.l_mm, center.l_mm),
                    (neighbor.l_mm, neighbor.l_mm) if direction == "higher_l" else (center.l_mm, center.l_mm),
                )
                cases.append((f"{center_name}_1d_row_k{k}_{direction}", "1d_row", center, row_grid,
                              abs(step_deg), k))
                grid2d, _clamped = build_2d_gradient_grid((center.l_mm, center.l_mm), curve, center, step_deg)
                cases.append((f"{center_name}_2d_gradient_k{k}_{direction}", "2d_gradient", center, grid2d,
                              abs(step_deg), k))
    return [
        _solve_case(name, kind, center, grid, step_deg, k, mesh_cells_per_subcell)
        for name, kind, center, grid, step_deg, k in cases
    ]
