"""A Phase-11-owned solve wrapper, identical to ``phase10.grid_solver.solve_grid`` except that the
excitation polarization is a parameter (Sec. 15's polarization-dependence study needs y-polarized
solves; Phase 10's own wrapper hardcodes "x"). Phase 10's file is not modified -- this mirrors
exactly how Phase 10 itself added a new wrapper (`grid_solver.py`) rather than modifying Phase 9's
`supercell_solver.py` when it needed generalized behavior the frozen file did not expose.

No new EM formulation, no new Green's function: `solver.assemble_system`/`solver.solve` and
`phase9.dipole_extraction.local_dipole_moment` are reused verbatim.
"""

from dataclasses import dataclass

from .. import solver as solver_mod
from ..basis_functions import MeshSettings
from ..excitation import PlaneWaveExcitation
from ..phase9.dipole_extraction import local_dipole_moment
from ..phase9.supercell_geometry import SupercellUnitCell, build_supercell

FLOQUET_FACTOR = 2  # matches phase9.run.FLOQUET_FACTOR / phase10.grid_solver.FLOQUET_FACTOR


@dataclass(frozen=True)
class SolvedGrid:
    supercell: SupercellUnitCell
    basis: object
    coeffs: object
    n_unknowns: int
    mesh_cells_per_subcell: int
    polarization: str

    def dipole_at(self, i, j):
        cx, cy = self.supercell.sub_cell_center(i, j)
        return local_dipole_moment(
            self.basis,
            self.coeffs,
            cx,
            cy,
            self.supercell.sub_period_x / 2,
            self.supercell.sub_period_y / 2,
        )

    def response_at(self, i, j):
        """The co-polarized dipole-moment component for this solve's excitation polarization
        (Sec. 15: the physically meaningful "center response" when comparing x- vs. y-polarized
        excitation is the component aligned with the incident field, not always px)."""
        d = self.dipole_at(i, j)
        return d.px if self.polarization == "x" else d.py


def solve_grid(
    n,
    period_m,
    substrate,
    metal,
    lw_grid_m,
    frequency_hz,
    mesh_cells_per_subcell=8,
    floquet_factor=FLOQUET_FACTOR,
    polarization="x",
):
    supercell = build_supercell(n, period_m, period_m, substrate, metal, lw_grid_m)
    mesh_cells = n * mesh_cells_per_subcell
    floquet_order = floquet_factor * mesh_cells
    excitation = PlaneWaveExcitation(frequency=frequency_hz, polarization=polarization)
    mesh = MeshSettings(mesh_cells, mesh_cells, floquet_order)
    system = solver_mod.assemble_system(supercell, excitation, mesh)
    currents = solver_mod.solve(system.z, system.v)
    return SolvedGrid(supercell, system.basis, currents, system.basis.n_unknowns, mesh_cells_per_subcell, polarization)


def solve_isolated(
    period_m,
    length_m,
    width_m,
    substrate,
    metal,
    frequency_hz,
    mesh_cells_per_side=8,
    floquet_factor=FLOQUET_FACTOR,
    polarization="x",
):
    """The n=1 (isolated periodic cell) reference, at a configurable polarization -- Phase 9's own
    ``solve_isolated_cell`` hardcodes x; reused conceptually, re-derived here only to add the
    polarization parameter (same physics, same solver calls)."""
    from ..geometry import Geometry, UnitCell
    from ..periodic_bc import RectangularLattice

    lattice = RectangularLattice(period_m, period_m)
    geometry = Geometry(kind="rectangle", length=length_m, width=width_m)
    unit_cell = UnitCell(lattice, geometry, substrate, metal)
    excitation = PlaneWaveExcitation(frequency=frequency_hz, polarization=polarization)
    mesh = MeshSettings(mesh_cells_per_side, mesh_cells_per_side, floquet_factor * mesh_cells_per_side)
    system = solver_mod.assemble_system(unit_cell, excitation, mesh)
    currents = solver_mod.solve(system.z, system.v)
    dipole = local_dipole_moment(system.basis, currents, 0.0, 0.0, period_m / 2, period_m / 2)
    response = dipole.px if polarization == "x" else dipole.py
    return response, system.basis.n_unknowns
