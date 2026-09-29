"""Sec. 16-19: the heterogeneous N x N supercell as one larger-period `UnitCell`-compatible
object. Duck-typed against the exact interface the frozen solver internals actually use
(grep-verified: `basis_functions.RooftopBasis` needs only `.lattice` and `.is_metal(x, y)`;
`solver.assemble_system`/`sparameters.compute_s_parameters` additionally need `.substrate`) --
`geometry.py`/`solver.py`/`basis_functions.py`/`green_function.py`/`mom_matrix.py` are not
modified; this is a new, separate object solved through their existing, unmodified public
functions, exactly the way a real `geometry.UnitCell` is.
"""

from dataclasses import dataclass

import numpy as np

from ..geometry import Geometry
from ..materials import Metal, Substrate
from ..periodic_bc import RectangularLattice

_BOUNDS_RTOL = 1e-9


@dataclass(frozen=True)
class SupercellUnitCell:
    """An n x n grid of independent rectangular-patch sub-cells, each ``sub_period_x`` x
    ``sub_period_y`` in size, tiled into one big period (``n*sub_period_x``, ``n*sub_period_y``).
    ``patches[i][j]``: the `Geometry` (kind='rectangle', centred at (0, 0) -- placement within the
    supercell is handled by ``is_metal`` below, not by the patch's own center) for sub-cell row
    ``i`` (y-index), column ``j`` (x-index)."""

    lattice: RectangularLattice
    substrate: Substrate
    metal: Metal
    n: int
    sub_period_x: float
    sub_period_y: float
    patches: tuple      # tuple[tuple[Geometry, ...], ...], shape (n, n)

    def __post_init__(self):
        if self.n < 1:
            raise ValueError(f"n must be >= 1, got {self.n!r}")
        if len(self.patches) != self.n or any(len(row) != self.n for row in self.patches):
            raise ValueError(f"patches must be an n x n grid, n={self.n}")
        if abs(self.lattice.period_x - self.n * self.sub_period_x) > _BOUNDS_RTOL * self.lattice.period_x:
            raise ValueError("lattice.period_x must equal n * sub_period_x")
        if abs(self.lattice.period_y - self.n * self.sub_period_y) > _BOUNDS_RTOL * self.lattice.period_y:
            raise ValueError("lattice.period_y must equal n * sub_period_y")
        tol_x, tol_y = _BOUNDS_RTOL * self.sub_period_x, _BOUNDS_RTOL * self.sub_period_y
        for row in self.patches:
            for geom in row:
                if geom.center_x != 0.0 or geom.center_y != 0.0:
                    raise ValueError("each patch's Geometry must be centred at (0, 0); placement is by grid index")
                corners = geom.corners()
                if np.any(np.abs(corners[:, 0]) > self.sub_period_x / 2 + tol_x) or np.any(
                    np.abs(corners[:, 1]) > self.sub_period_y / 2 + tol_y
                ):
                    raise ValueError(f"a patch extends outside its {self.sub_period_x}x{self.sub_period_y} sub-cell")

    @property
    def metal_area(self):
        return sum(geom.feature_area for row in self.patches for geom in row)

    def sub_cell_center(self, i, j):
        """Physical (x, y) centre of sub-cell (row i = y-index, column j = x-index)."""
        px, py = self.lattice.period_x, self.lattice.period_y
        cx = -px / 2 + (j + 0.5) * self.sub_period_x
        cy = -py / 2 + (i + 0.5) * self.sub_period_y
        return cx, cy

    def is_metal(self, x, y):
        px, py = self.lattice.period_x, self.lattice.period_y
        xw = np.asarray(x, dtype=float)
        yw = np.asarray(y, dtype=float)
        xw = xw - px * np.floor(xw / px + 0.5)
        yw = yw - py * np.floor(yw / py + 0.5)
        j_idx = np.clip(np.floor((xw + px / 2) / self.sub_period_x).astype(int), 0, self.n - 1)
        i_idx = np.clip(np.floor((yw + py / 2) / self.sub_period_y).astype(int), 0, self.n - 1)

        result = np.zeros(xw.shape, dtype=bool)
        for i in range(self.n):
            for j in range(self.n):
                mask = (i_idx == i) & (j_idx == j)
                if not np.any(mask):
                    continue
                cx, cy = self.sub_cell_center(i, j)
                result[mask] = self.patches[i][j].in_feature(xw[mask] - cx, yw[mask] - cy)
        return result


def build_supercell(n, sub_period_x, sub_period_y, substrate, metal, lw_grid):
    """``lw_grid``: an n x n array-like of (L_m, W_m) pairs (row i = y-index, column j = x-index).
    Builds the `SupercellUnitCell` (period = n*sub_period_x, n*sub_period_y)."""
    lattice = RectangularLattice(n * sub_period_x, n * sub_period_y)
    patches = tuple(
        tuple(Geometry(kind="rectangle", length=lw_grid[i][j][0], width=lw_grid[i][j][1]) for j in range(n))
        for i in range(n)
    )
    return SupercellUnitCell(lattice, substrate, metal, n, sub_period_x, sub_period_y, patches)
