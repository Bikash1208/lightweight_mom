"""Build and solve the MoM system  Z I = V  (Phase 3: no S-parameters yet).

    geometry -> mesh / rooftop basis -> excitation vector V -> impedance matrix Z

``build_system`` runs that chain from a validated ``SimulationConfig``;
``solve`` performs the dense linear solve.
"""

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .basis_functions import RooftopBasis
from .excitation import stack_tangential_field
from .green_function import LayeredGreenFunction
from .mom_matrix import assemble_excitation_vector, assemble_impedance_matrix

if TYPE_CHECKING:       # config imports the sweep machinery, which imports this module
    from .config import SimulationConfig


@dataclass
class MomSystem:
    basis: RooftopBasis
    green: LayeredGreenFunction
    kx0: float
    ky0: float
    e_stack: np.ndarray      # tangential E of the unpatterned stack at z = 0, shape (2,)
    z: np.ndarray            # (n, n)
    v: np.ndarray            # (n,)


def assemble_system(unit_cell, excitation, mesh) -> MomSystem:
    """geometry -> mesh/basis -> excitation vector V -> impedance matrix Z."""
    basis = RooftopBasis(unit_cell, mesh.cells_x, mesh.cells_y)
    # getattr, not unit_cell.backing: duck-typed "UnitCell-like" objects that predate Phase 12A
    # (e.g. phase9.supercell_geometry.SupercellUnitCell) have no such attribute -- absent means
    # "air", exactly like every other pre-Phase-12A context/config (CONVENTIONS.md).
    green = LayeredGreenFunction(
        unit_cell.substrate, excitation.frequency, backing=getattr(unit_cell, "backing", "air")
    )
    kx0, ky0 = excitation.transverse_wavenumber
    e_stack = stack_tangential_field(green, excitation)
    v = assemble_excitation_vector(basis, e_stack, kx0, ky0)
    z = assemble_impedance_matrix(basis, green, kx0, ky0, mesh.floquet_order)
    return MomSystem(basis, green, kx0, ky0, e_stack, z, v)


def build_system(config: "SimulationConfig") -> MomSystem:
    """As ``assemble_system`` from a validated configuration. The reference planes in the
    configuration do not enter Z or V."""
    return assemble_system(config.unit_cell, config.excitation, config.mesh)


def solve(z, v):
    """Solve Z I = V (dense LU) for one right-hand side (n,) or several (n, k).

    Raises ``numpy.linalg.LinAlgError`` if Z is singular. Z does not depend on the excitation,
    so all polarizations/ports of a cell share one factorization.
    """
    z = np.asarray(z, dtype=complex)
    v = np.asarray(v, dtype=complex)
    if z.ndim != 2 or z.shape[0] != z.shape[1] or v.ndim not in (1, 2) or v.shape[0] != z.shape[0]:
        raise ValueError(f"incompatible shapes Z {z.shape}, V {v.shape}")
    current = np.linalg.solve(z, v)
    if not np.all(np.isfinite(current)):
        raise np.linalg.LinAlgError("solution contains non-finite values")
    return current
