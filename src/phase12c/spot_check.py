"""Phase 12C Task D: spot-check Phases 9-11 on a finer ALIGNED mesh -- read-only for earlier
results (nothing under ``results/phase9_validation/``, ``results/phase10_*/`` or
``results/phase11_heterogeneous_mom/`` is edited or re-run; every output here goes under
``results/phase12c_*``). Reuses Phase 9/10/11's own frozen geometry-construction functions
verbatim (``phase9.local_periodicity.build_3x3_grid``, ``phase9.supercell_geometry.build_supercell``,
``phase9.dipole_extraction.local_dipole_moment``, ``phase10.environment_geometry.CENTER_CELLS_MM``,
``phase11.patterns.checkerboard``/``abrupt_boundary``) -- only the mesh resolution and the solve
call change (a direct ``solver.assemble_system``/``solve`` call with an EXPLICIT, cost-bounded
floquet order, rather than Phase 9/11's own ``floquet_factor=2`` convention, which would make an
n>1 supercell at this mesh categorically infeasible; see module-level cost note).

Domain: Phase 9/10/11's own P = 4.8 mm, RT/Duroid 5880 (1.575 mm) stack, ungrounded (this task is
about mesh resolution, not the grounded/8 mm design point the rest of Phase 12C studies).

Mesh: ``discretization.smallest_aligned_cells`` for each specific (small) set of L values actually
re-checked here -- NOT Phase 9/10/11's original full 18-value, 0.2 mm-step family's own 48-cell
figure applied blindly to every case; a smaller specific subset can align more cheaply (measured:
Phase 10's C2/C3 pair aligns at 24 cells, cheaper than the full family's 48).

Cost note (measured during this phase's own development): assemble+solve time for a supercell grows
roughly as the 6th power of the per-side mesh-cell count once floquet_order is scaled with it. At
48 cells/subcell, an n=3 supercell (144 total mesh cells/side) is a single feasible spot-check
(tens of minutes); Phase 11's own original n=7 pattern study at 48 cells/subcell (336 total) is not
(days, by the same scaling) -- Task D's own Case D re-check is therefore run at a REDUCED n=3
(disclosed below, not silently substituted for the original n=7 finding), and its coupling-radius
spot-check stops at n=3 (n=5 at this mesh, ~93 min extrapolated from measured n<=3 timings, was
deemed a disproportionate cost for a single spot-check; n=7+ is not attempted).
"""

from dataclasses import dataclass

from .. import solver as solver_mod
from ..basis_functions import MeshSettings
from ..discretization import metal_mask_key, smallest_aligned_cells
from ..excitation import PlaneWaveExcitation
from ..geometry import Geometry, UnitCell
from ..materials import Metal, Substrate
from ..periodic_bc import RectangularLattice
from ..phase9 import local_periodicity as p9_lpa
from ..phase9.dipole_extraction import local_dipole_moment
from ..phase9.supercell_geometry import build_supercell
from ..phase10 import environment_geometry as p10_eg
from ..phase11 import patterns as p11_patterns

PERIOD_M = 4.8e-3  # Phase 9/10/11's own period -- unchanged, this task is about mesh resolution
SUBSTRATE = Substrate("Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3)
METAL = Metal()
DESIGN_FREQUENCY_HZ = 12.5e9

REDUCED_CASE_D_N = 3  # Sec. 'Cost note': reduced from Phase 11's original PATTERN_N=7
COUPLING_SIZES = (1, 3)  # Sec. 'Cost note': n=5+ not attempted at this mesh


def floquet_for(total_cells):
    """``total_cells // 2 + 8``: a rank-sufficient order with real margin, deliberately NOT the
    fully floquet-converged ``2 * total_cells`` convention ``locality_study.floquet_for`` validates
    for Task B's fine-grained (few-percent) deviation measurements -- at 48 cells/subcell that
    convention is measured to cost tens of minutes to hours per solve (module cost note). Task D's
    own re-checks all concern LARGE, qualitative effects (Phase 9's strong-contrast ~99.7% relative
    error; Phase 11's Case D "substantial" >=50% threshold) where a residual few-percent
    floquet-truncation artifact -- the size measured for a similarly modest order in
    ``locality_study``'s own convergence scan -- cannot flip the qualitative conclusion, so the
    cheaper order is the appropriate, disclosed choice here."""
    return total_cells // 2 + 8


def aligned_cells_for(l_values_mm, max_cells=200):
    features = [(l * 1e-3, 0.0) for l in l_values_mm]
    n = smallest_aligned_cells(PERIOD_M, features, max_cells=max_cells)
    if n is None:
        raise ValueError(f"no aligned mesh found for {l_values_mm} within {max_cells} cells")
    return n


def _isolated(length_m, width_m, mesh_cells_per_side, floquet_order):
    lattice = RectangularLattice(PERIOD_M, PERIOD_M)
    geometry = Geometry(kind="rectangle", length=length_m, width=width_m)
    unit_cell = UnitCell(lattice, geometry, SUBSTRATE, METAL)
    excitation = PlaneWaveExcitation(frequency=DESIGN_FREQUENCY_HZ, polarization="x")
    mesh = MeshSettings(mesh_cells_per_side, mesh_cells_per_side, floquet_order)
    system = solver_mod.assemble_system(unit_cell, excitation, mesh)
    currents = solver_mod.solve(system.z, system.v)
    dipole = local_dipole_moment(system.basis, currents, 0.0, 0.0, PERIOD_M / 2, PERIOD_M / 2)
    return dipole, system.basis.n_unknowns, unit_cell


def _supercell_center(n, lw_grid_m, mesh_cells_per_subcell, floquet_order):
    supercell = build_supercell(n, PERIOD_M, PERIOD_M, SUBSTRATE, METAL, lw_grid_m)
    mesh_cells = n * mesh_cells_per_subcell
    excitation = PlaneWaveExcitation(frequency=DESIGN_FREQUENCY_HZ, polarization="x")
    mesh = MeshSettings(mesh_cells, mesh_cells, floquet_order)
    system = solver_mod.assemble_system(supercell, excitation, mesh)
    currents = solver_mod.solve(system.z, system.v)
    mid = n // 2
    cx, cy = supercell.sub_cell_center(mid, mid)
    dipole = local_dipole_moment(system.basis, currents, cx, cy, PERIOD_M / 2, PERIOD_M / 2)
    return dipole, system.basis.n_unknowns


# ------------------------------------------------------------- Phase 10: C2 vs C3 re-verification


@dataclass(frozen=True)
class C2C3Recheck:
    mesh_cells_per_side: int
    mask_identical: bool
    c2_isolated_px: complex
    c3_isolated_px: complex
    relative_difference: float


def recheck_phase10_c2_c3():
    """Sec. D: re-verify Phase 10's C2 = (2.0, 2.8) mm vs C3 = (2.8, 2.8) mm (found mesh-identical
    by Phase 12B at the solver's own 8-cell/side mesh) on the smallest mesh that aligns BOTH their
    edges -- a real, isolated-cell solve at that finer mesh, not just the mask check alone (which
    Phase 12B's audit already performed; this adds the same real-solve confirmation the brief
    asks for)."""
    c2_lw, c3_lw = p10_eg.CENTER_CELLS_MM["C2"], p10_eg.CENTER_CELLS_MM["C3"]
    n = aligned_cells_for([c2_lw[0], c2_lw[1], c3_lw[0], c3_lw[1]])
    floquet = floquet_for(n)

    def _cell(lw_mm):
        return UnitCell(
            RectangularLattice(PERIOD_M, PERIOD_M),
            Geometry(kind="rectangle", length=lw_mm[0] * 1e-3, width=lw_mm[1] * 1e-3),
            SUBSTRATE,
            METAL,
        )

    mask_c2, mask_c3 = (
        metal_mask_key(_cell(c2_lw), MeshSettings(n, n, floquet)),
        metal_mask_key(_cell(c3_lw), MeshSettings(n, n, floquet)),
    )
    dip_c2, _, _ = _isolated(c2_lw[0] * 1e-3, c2_lw[1] * 1e-3, n, floquet)
    dip_c3, _, _ = _isolated(c3_lw[0] * 1e-3, c3_lw[1] * 1e-3, n, floquet)
    rel = abs(dip_c2.px - dip_c3.px) / abs(dip_c3.px) if abs(dip_c3.px) > 0 else float("inf")
    return C2C3Recheck(n, mask_c2 == mask_c3, dip_c2.px, dip_c3.px, rel)


# ------------------------------------------------------------------ Phase 9: strong-contrast 3x3


@dataclass(frozen=True)
class Phase9SpotCheck:
    mesh_cells_per_subcell: int
    n: int
    isolated_px: complex
    embedded_px: complex
    relative_complex_error: float
    n_unknowns_isolated: int
    n_unknowns_embedded: int
    earlier_relative_complex_error: float  # PHASE9_VALIDATION.md's own reported 3x3-strong figure


def recheck_phase9_strong_contrast(earlier_relative_complex_error=0.997, mesh_cells_per_subcell=None):
    """``mesh_cells_per_subcell``: override for testing only (skips the real aligned-mesh
    computation, which is expensive for this specific value set -- Sec. 'Cost note'); the
    production call (``run.py``) always leaves this ``None`` and gets the real aligned mesh."""
    # NOTE: phase9.local_periodicity.build_3x3_grid returns METRES despite its "_mm" grid naming
    # internally (it converts via `_mm_to_m` before returning) -- use it directly as the solve
    # grid, and get the mm-valued alignment inputs from the module's own explicit mm constants.
    grid_m = p9_lpa.build_3x3_grid("strong", p9_lpa.CENTER_LW_MM)
    center_l_mm, center_w_mm = p9_lpa.CENTER_LW_MM
    if mesh_cells_per_subcell is None:
        all_values_mm = sorted({v for pair in p9_lpa.STRONG_NEIGHBORS_MM for v in pair} | {center_l_mm, center_w_mm})
        n_mesh = aligned_cells_for(all_values_mm)
    else:
        n_mesh = mesh_cells_per_subcell
    n = 3
    floquet = floquet_for(n * n_mesh)
    isolated, n_iso, _ = _isolated(center_l_mm * 1e-3, center_w_mm * 1e-3, n_mesh, floquet_for(n_mesh))
    embedded, n_emb = _supercell_center(n, grid_m, n_mesh, floquet)
    rel = abs(embedded.px - isolated.px) / abs(isolated.px) if abs(isolated.px) > 0 else float("inf")
    return Phase9SpotCheck(n_mesh, n, isolated.px, embedded.px, rel, n_iso, n_emb, earlier_relative_complex_error)


# --------------------------------------------------------------- Phase 11: Case D, reduced size


@dataclass(frozen=True)
class CaseDRecheck:
    mesh_cells_per_subcell: int
    n: int
    checkerboard_px: complex
    block_px: complex
    relative_difference: float
    substantial: bool


def recheck_phase11_case_d(n=REDUCED_CASE_D_N, mesh_cells_per_subcell=None):
    """Sec. D: Phase 11's own Case D comparison (checkerboard vs a contiguous block, same two
    geometries, same 50/50 split), reduced to n (default 3, from the original PATTERN_N=7 --
    Sec. 'Cost note') so the aligned mesh stays feasible. ``mesh_cells_per_subcell``: override for
    testing only, as in ``recheck_phase9_strong_contrast``."""
    a_mm, b_mm = p11_patterns.A_DEFAULT_MM, p11_patterns.B_DEFAULT_MM
    n_mesh = mesh_cells_per_subcell if mesh_cells_per_subcell is not None else aligned_cells_for([a_mm[0], b_mm[0]])
    floquet = floquet_for(n * n_mesh)
    checker_mm = p11_patterns.checkerboard(n, a_mm, b_mm)
    block_mm = p11_patterns.abrupt_boundary(n, a_mm, b_mm)
    checker_m = [[(l * 1e-3, w * 1e-3) for l, w in row] for row in checker_mm]
    block_m = [[(l * 1e-3, w * 1e-3) for l, w in row] for row in block_mm]
    checker_dipole, _ = _supercell_center(n, checker_m, n_mesh, floquet)
    block_dipole, _ = _supercell_center(n, block_m, n_mesh, floquet)
    rel = (
        abs(block_dipole.px - checker_dipole.px) / abs(checker_dipole.px)
        if abs(checker_dipole.px) > 0
        else float("inf")
    )
    return CaseDRecheck(n_mesh, n, checker_dipole.px, block_dipole.px, rel, rel >= 0.5)


# ---------------------------------------------------------- Phase 11: coupling-radius spot-check


@dataclass(frozen=True)
class CouplingSpotCheckPoint:
    n: int
    response: complex
    n_unknowns: int


def recheck_phase11_coupling_radius(center_lw_mm=None, sizes=COUPLING_SIZES, mesh_cells_per_subcell=None):
    """Sec. D: one Phase 11 coupling-radius series re-solved at a finer aligned mesh, to the
    largest feasible size (Sec. 'Cost note': n=5+ not attempted). ``mesh_cells_per_subcell``:
    override for testing only, as in ``recheck_phase9_strong_contrast`` (the real "opposite"
    contrast value is always one of the domain's 0.1 mm-granularity extremes, so even a custom
    ``center_lw_mm`` still forces the expensive 48-cell alignment unless this is given)."""
    from ..phase10.environment_geometry import strong_neighbor_low_high
    from ..phase11.coupling_radius import truncate_centered

    center_lw_mm = center_lw_mm or p9_lpa.CENTER_LW_MM
    low, high = strong_neighbor_low_high(center_lw_mm)
    opposite = (low, low) if abs(center_lw_mm[0] - high) < abs(center_lw_mm[0] - low) else (high, high)
    max_n = max(sizes)
    full_grid_mm = p11_patterns.checkerboard(max_n, a_mm=center_lw_mm, b_mm=opposite)
    n_mesh = (
        mesh_cells_per_subcell
        if mesh_cells_per_subcell is not None
        else aligned_cells_for([center_lw_mm[0], opposite[0]])
    )
    points = []
    for n in sorted(sizes):
        grid_mm = truncate_centered(full_grid_mm, n)
        grid_m = [[(l * 1e-3, w * 1e-3) for l, w in row] for row in grid_mm]
        floquet = floquet_for(n * n_mesh)
        dipole, n_unknowns = _supercell_center(n, grid_m, n_mesh, floquet)
        points.append(CouplingSpotCheckPoint(n, dipole.px, n_unknowns))
    return n_mesh, tuple(points)
