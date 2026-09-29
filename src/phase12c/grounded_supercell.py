"""Phase 12C Task A: grounded heterogeneous supercells.

Wraps (never edits) Phase 9's frozen ``phase9.supercell_geometry.SupercellUnitCell`` with a
``backing`` field -- exactly the composition pattern Phase 12A itself established
(``getattr(unit_cell, "backing", "air")`` in ``solver.assemble_system``/
``sparameters.compute_s_parameters``): ``SupercellUnitCell`` has no ``backing`` attribute, so it
already defaults to ``"air"`` through that ``getattr``; this module adds a thin, duck-typed wrapper
that DOES have one, set to ``"pec"``, delegating every other attribute/method the solve pipeline
actually reads (grep-verified, same set Phase 9's own module docstring lists: ``.lattice``,
``.substrate``, ``.is_metal(x, y)``, plus ``.metal_area``/``.sub_cell_center``/``.n``/
``.sub_period_x``/``.sub_period_y``/``.patches`` for the dipole-extraction read-out) to the wrapped
object. ``supercell_geometry.py`` is not modified.
"""

from dataclasses import dataclass

from .. import solver as solver_mod
from ..basis_functions import MeshSettings
from ..excitation import PlaneWaveExcitation
from ..geometry import Geometry, UnitCell
from ..periodic_bc import RectangularLattice
from ..phase9.dipole_extraction import local_dipole_moment
from ..phase9.supercell_geometry import SupercellUnitCell, build_supercell


@dataclass(frozen=True)
class GroundedSupercellUnitCell:
    """Duck-typed grounded wrapper around a :class:`SupercellUnitCell`. ``backing="pec"`` is the
    only new information; every other attribute is the wrapped cell's own."""

    inner: SupercellUnitCell
    backing: str = "pec"

    @property
    def lattice(self):
        return self.inner.lattice

    @property
    def substrate(self):
        return self.inner.substrate

    @property
    def metal(self):
        return self.inner.metal

    @property
    def n(self):
        return self.inner.n

    @property
    def sub_period_x(self):
        return self.inner.sub_period_x

    @property
    def sub_period_y(self):
        return self.inner.sub_period_y

    @property
    def patches(self):
        return self.inner.patches

    @property
    def metal_area(self):
        return self.inner.metal_area

    def sub_cell_center(self, i, j):
        return self.inner.sub_cell_center(i, j)

    def is_metal(self, x, y):
        return self.inner.is_metal(x, y)


def build_grounded_supercell(n, sub_period_x, sub_period_y, substrate, metal, lw_grid):
    """As ``phase9.supercell_geometry.build_supercell``, grounded (``backing="pec"``)."""
    return GroundedSupercellUnitCell(build_supercell(n, sub_period_x, sub_period_y, substrate, metal, lw_grid))


@dataclass(frozen=True)
class SolvedDipole:
    dipole: object  # phase9.dipole_extraction.DipoleMoment
    n_unknowns: int


def solve_isolated_grounded_cell(
    period_x, period_y, length_m, width_m, substrate, metal, frequency_hz, mesh_cells_per_side, floquet_order
):
    """The grounded reference: Phase 9's own ``supercell_solver.solve_isolated_cell``, generalized
    with ``backing="pec"`` -- a real capability (a grounded ``UnitCell`` solved through the exact
    same, unmodified ``solver.assemble_system``/``solve``), not a new EM formulation."""
    lattice = RectangularLattice(period_x, period_y)
    geometry = Geometry(kind="rectangle", length=length_m, width=width_m)
    unit_cell = UnitCell(lattice, geometry, substrate, metal, backing="pec")
    excitation = PlaneWaveExcitation(frequency=frequency_hz, polarization="x")
    mesh = MeshSettings(mesh_cells_per_side, mesh_cells_per_side, floquet_order)
    system = solver_mod.assemble_system(unit_cell, excitation, mesh)
    currents = solver_mod.solve(system.z, system.v)
    dipole = local_dipole_moment(system.basis, currents, 0.0, 0.0, period_x / 2, period_y / 2)
    return SolvedDipole(dipole, system.basis.n_unknowns)


def solve_grounded_supercell_center(
    supercell: GroundedSupercellUnitCell, center_lw, frequency_hz, mesh_cells_per_subcell, floquet_order
):
    """As Phase 9's own ``supercell_solver.solve_supercell_center``, for a
    :class:`GroundedSupercellUnitCell`."""
    inner = supercell.inner
    i_center = j_center = inner.n // 2
    if inner.n % 2 == 0:
        raise ValueError("solve_grounded_supercell_center: n must be odd (a unique centre sub-cell)")
    center_geom = inner.patches[i_center][j_center]
    if (center_geom.length, center_geom.width) != tuple(center_lw):
        raise ValueError(
            f"supercell centre patch is {(center_geom.length, center_geom.width)}, expected {tuple(center_lw)}"
        )
    mesh_cells = inner.n * mesh_cells_per_subcell
    excitation = PlaneWaveExcitation(frequency=frequency_hz, polarization="x")
    mesh = MeshSettings(mesh_cells, mesh_cells, floquet_order)
    system = solver_mod.assemble_system(supercell, excitation, mesh)
    currents = solver_mod.solve(system.z, system.v)
    cx, cy = inner.sub_cell_center(i_center, j_center)
    dipole = local_dipole_moment(
        system.basis,
        currents,
        cx,
        cy,
        inner.sub_period_x / 2,
        inner.sub_period_y / 2,
    )
    return SolvedDipole(dipole, system.basis.n_unknowns)


# ------------------------------------------------ Task A: correctness proof, cost, convergence

import time


@dataclass(frozen=True)
class CorrectnessCheckResult:
    n: int
    mesh_cells_per_subcell: int
    isolated_px: complex
    embedded_px: complex
    relative_error: float
    n_unknowns_isolated: int
    n_unknowns_embedded: int
    runtime_s: float
    estimated_z_bytes: int  # dominant memory cost: n_unknowns^2 * 16 bytes (complex128 dense Z)


def homogeneous_correctness_check(length_m, n, mesh_cells_per_subcell, period_m, substrate, metal, frequency_hz):
    """Sec. A: a homogeneous n x n grounded supercell (every sub-cell identical) MUST reproduce the
    isolated grounded cell's own dipole moment exactly (Phase 9's own "9 identical cells repeated
    periodically == 1 cell repeated periodically" argument, unchanged by grounding), to the
    solver's own floating-point noise floor -- the load-bearing proof that the wrapper introduces
    no numerical error of its own. Also reports the real, measured runtime and the dominant
    (dense Z matrix) memory cost, as Task A's own brief requires."""
    t0 = time.perf_counter()
    # floquet_order = 2 * total_mesh_cells (Phase 9's own convention), applied to EACH solve's OWN
    # total mesh cells (the isolated cell's is ``mesh_cells_per_subcell``, the embedded supercell's
    # is ``n * mesh_cells_per_subcell``) -- measured (module docstring in locality_study.py) to be
    # the order that actually reaches machine precision here; a merely rank-rule-sufficient, smaller
    # order was tried first and measured to leave a genuine ~1% floquet-truncation residual, i.e.
    # this is not a cost-driven approximation but the validated value.
    floquet_isolated = 2 * mesh_cells_per_subcell
    floquet_embedded = 2 * mesh_cells_per_subcell * n
    isolated = solve_isolated_grounded_cell(
        period_m,
        period_m,
        length_m,
        length_m,
        substrate,
        metal,
        frequency_hz,
        mesh_cells_per_subcell,
        floquet_isolated,
    )
    grid = [[(length_m, length_m)] * n for _ in range(n)]
    supercell = build_grounded_supercell(n, period_m, period_m, substrate, metal, grid)
    embedded = solve_grounded_supercell_center(
        supercell, (length_m, length_m), frequency_hz, mesh_cells_per_subcell, floquet_embedded
    )
    runtime = time.perf_counter() - t0
    rel = abs(embedded.dipole.px - isolated.dipole.px) / abs(isolated.dipole.px)
    z_bytes = embedded.n_unknowns**2 * 16
    return CorrectnessCheckResult(
        n,
        mesh_cells_per_subcell,
        isolated.dipole.px,
        embedded.dipole.px,
        rel,
        isolated.n_unknowns,
        embedded.n_unknowns,
        runtime,
        z_bytes,
    )


@dataclass(frozen=True)
class MeshConvergenceResult:
    n: int
    base_mesh_cells_per_subcell: int
    finer_mesh_cells_per_subcell: int
    base_px: complex
    finer_px: complex
    relative_difference: float
    base_runtime_s: float
    finer_runtime_s: float


def supercell_mesh_convergence_check(
    length_m, n, base_mesh_cells_per_subcell, finer_mesh_cells_per_subcell, period_m, substrate, metal, frequency_hz
):
    """Sec. A: "verify on a 3x3 case that the result converges against a finer aligned mesh" --
    the SAME homogeneous n x n case at ``base_mesh_cells_per_subcell`` and (doubled, so every
    original mesh line stays a mesh line, Phase 12B's own alignment-preserving-doubling argument)
    ``finer_mesh_cells_per_subcell``."""
    grid = [[(length_m, length_m)] * n for _ in range(n)]

    def _run(mesh_cells_per_subcell):
        t0 = time.perf_counter()
        floquet = 2 * mesh_cells_per_subcell * n  # validated convention, see homogeneous_correctness_check
        supercell = build_grounded_supercell(n, period_m, period_m, substrate, metal, grid)
        embedded = solve_grounded_supercell_center(
            supercell, (length_m, length_m), frequency_hz, mesh_cells_per_subcell, floquet
        )
        return embedded.dipole.px, time.perf_counter() - t0

    base_px, base_t = _run(base_mesh_cells_per_subcell)
    finer_px, finer_t = _run(finer_mesh_cells_per_subcell)
    rel = abs(finer_px - base_px) / abs(base_px) if abs(base_px) > 0 else float("inf")
    return MeshConvergenceResult(
        n, base_mesh_cells_per_subcell, finer_mesh_cells_per_subcell, base_px, finer_px, rel, base_t, finer_t
    )
