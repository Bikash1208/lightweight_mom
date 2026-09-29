"""A thin, Phase-10-owned wrapper around the same reuse chain
``phase9.supercell_geometry.build_supercell`` -> ``solver.assemble_system``/``solver.solve`` ->
``phase9.dipole_extraction.local_dipole_moment`` that ``phase9.supercell_solver`` already
established, generalized in one respect: ``phase9.supercell_solver.solve_supercell_center`` only
extracts the one center sub-cell of an odd-sized grid, but several Phase 10 studies (Sec. 30-32's
gradient/abrupt-discontinuity tests, Sec. 28-29's beam-level aperture) need many sub-cells' local
responses from the *same* solved system (one MoM solve, many read-outs -- solving once per pixel
would be wasteful and is not how the physical current distribution works: it is one self-consistent
solution over the whole supercell). No new EM formulation, no new Green's function, no coupled-cell
solver -- exactly Phase 9's own reuse discipline, extended to more than one read-out point.
"""

from dataclasses import dataclass

from .. import solver as solver_mod
from ..basis_functions import MeshSettings
from ..excitation import PlaneWaveExcitation
from ..phase9.dipole_extraction import local_dipole_moment
from ..phase9.supercell_geometry import SupercellUnitCell, build_supercell

FLOQUET_FACTOR = 2  # matches phase9.run.FLOQUET_FACTOR; kept identical so every solve in Phase 10
                     # uses the same, already-verified mesh/Floquet convergence ratio as Phase 9's
                     # homogeneous-exact-match proof.


@dataclass(frozen=True)
class SolvedGrid:
    supercell: SupercellUnitCell
    basis: object
    coeffs: object
    n_unknowns: int
    mesh_cells_per_subcell: int

    def dipole_at(self, i, j):
        """The local induced-current dipole moment of sub-cell (row i, column j)."""
        cx, cy = self.supercell.sub_cell_center(i, j)
        return local_dipole_moment(
            self.basis, self.coeffs, cx, cy,
            self.supercell.sub_period_x / 2, self.supercell.sub_period_y / 2,
        )


def solve_grid(n, period_m, substrate, metal, lw_grid_m, frequency_hz, mesh_cells_per_subcell=8,
               floquet_factor=FLOQUET_FACTOR):
    """Builds and solves one n x n heterogeneous supercell (``lw_grid_m[i][j] = (L_m, W_m)``),
    returning a ``SolvedGrid`` from which any sub-cell's dipole moment can be read out cheaply
    (no re-solve)."""
    supercell = build_supercell(n, period_m, period_m, substrate, metal, lw_grid_m)
    mesh_cells = n * mesh_cells_per_subcell
    floquet_order = floquet_factor * mesh_cells
    excitation = PlaneWaveExcitation(frequency=frequency_hz, polarization="x")
    mesh = MeshSettings(mesh_cells, mesh_cells, floquet_order)
    system = solver_mod.assemble_system(supercell, excitation, mesh)
    currents = solver_mod.solve(system.z, system.v)
    return SolvedGrid(supercell, system.basis, currents, system.basis.n_unknowns, mesh_cells_per_subcell)
