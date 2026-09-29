"""Sec. 20: solves the isolated single-cell reference and the heterogeneous supercell through the
*exact same*, unmodified `solver.assemble_system`/`solver.solve` MoM pipeline (Phase 1-3) and
extracts each one's centre-cell dipole moment (`dipole_extraction.py`). No new EM formulation, no
new Green's function, no coupled-cell solver: the only new code is the composite geometry
(`supercell_geometry.py`) and the dipole-moment read-out.
"""

from dataclasses import dataclass

from .. import solver as solver_mod
from ..basis_functions import MeshSettings
from ..excitation import PlaneWaveExcitation
from ..geometry import Geometry, UnitCell
from ..periodic_bc import RectangularLattice
from .dipole_extraction import local_dipole_moment
from .supercell_geometry import SupercellUnitCell


@dataclass(frozen=True)
class SolvedDipole:
    dipole: object  # DipoleMoment
    n_unknowns: int


def solve_isolated_cell(
    period_x, period_y, length_m, width_m, substrate, metal, frequency_hz, mesh_cells_per_side, floquet_order
):
    """The reference: this geometry alone, periodic, exactly the Phase 1-7 single-cell pipeline
    (`solver.assemble_system`, unmodified), at the same mesh resolution the supercell's own
    sub-cells use (so the comparison is apples-to-apples, Sec. 20-21)."""
    lattice = RectangularLattice(period_x, period_y)
    geometry = Geometry(kind="rectangle", length=length_m, width=width_m)
    unit_cell = UnitCell(lattice, geometry, substrate, metal)
    excitation = PlaneWaveExcitation(frequency=frequency_hz, polarization="x")
    mesh = MeshSettings(mesh_cells_per_side, mesh_cells_per_side, floquet_order)
    system = solver_mod.assemble_system(unit_cell, excitation, mesh)
    currents = solver_mod.solve(system.z, system.v)
    dipole = local_dipole_moment(system.basis, currents, 0.0, 0.0, period_x / 2, period_y / 2)
    return SolvedDipole(dipole, system.basis.n_unknowns)


def solve_supercell_center(
    supercell: SupercellUnitCell, center_lw, frequency_hz, mesh_cells_per_subcell, floquet_order
):
    """The heterogeneous case: solves the whole supercell (one MoM solve, `solver.assemble_system`,
    unmodified) and extracts the centre sub-cell's dipole moment. ``center_lw`` is not used for the
    solve (the supercell's own ``patches`` already define every sub-cell's geometry, including the
    centre); it is accepted only so the caller's intent is explicit and checked."""
    i_center = j_center = supercell.n // 2
    if supercell.n % 2 == 0:
        raise ValueError("solve_supercell_center: n must be odd (a unique centre sub-cell)")
    center_geom = supercell.patches[i_center][j_center]
    if (center_geom.length, center_geom.width) != tuple(center_lw):
        raise ValueError(
            f"supercell centre patch is {(center_geom.length, center_geom.width)}, expected {tuple(center_lw)}"
        )
    mesh_cells = supercell.n * mesh_cells_per_subcell
    excitation = PlaneWaveExcitation(frequency=frequency_hz, polarization="x")
    mesh = MeshSettings(mesh_cells, mesh_cells, floquet_order)
    system = solver_mod.assemble_system(supercell, excitation, mesh)
    currents = solver_mod.solve(system.z, system.v)
    cx, cy = supercell.sub_cell_center(i_center, j_center)
    dipole = local_dipole_moment(
        system.basis,
        currents,
        cx,
        cy,
        supercell.sub_period_x / 2,
        supercell.sub_period_y / 2,
    )
    return SolvedDipole(dipole, system.basis.n_unknowns)
