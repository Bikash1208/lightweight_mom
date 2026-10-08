"""Phase 19 -- RCWA unit-cell solver for z-invariant dielectric (3D-printed) pillars.

New package; does not touch the periodic-MoM core.  See ``docs/theory/DIELECTRIC_CELL_RCWA.md``.
"""

from .cell import CellResponse, OctagonPillarCell, solve_cell
from .rcwa import PolygonLayer, UniformLayer, solve_stack

__all__ = [
    "CellResponse",
    "OctagonPillarCell",
    "PolygonLayer",
    "UniformLayer",
    "solve_cell",
    "solve_stack",
]
