"""S_C, the isolated-periodic-cell reference response, is a physical quantity directly available
from the existing frozen solver (`phase9.supercell_solver.solve_isolated_cell`) for ANY (L, W) in
the validated continuous range -- it needs no fitting/training and generalizes automatically to a
center geometry never seen in the dataset (e.g. a continuous, non-grid-aligned L used by
`generalization.py`'s gradient test, or a fully held-out center like C5). Baseline 1, Baseline 2
and the surrogate (`surrogate.py`) all share ONE memoized cache of this reference rather than each
maintaining their own (inconsistent, training-data-limited) lookup table.
"""

from dataclasses import dataclass, field

from ..phase9.supercell_solver import solve_isolated_cell


@dataclass
class IsolatedReferenceCache:
    period_m: float
    substrate: object
    metal: object
    mesh_cells_per_subcell: int = 8
    _cache: dict = field(default_factory=dict)

    def get(self, center_lw_mm, frequency_hz):
        key = (round(center_lw_mm[0], 6), round(center_lw_mm[1], 6), round(frequency_hz, 3))
        if key not in self._cache:
            l_m, w_m = center_lw_mm[0] * 1e-3, center_lw_mm[1] * 1e-3
            solved = solve_isolated_cell(self.period_m, self.period_m, l_m, w_m, self.substrate,
                                          self.metal, frequency_hz, self.mesh_cells_per_subcell,
                                          2 * self.mesh_cells_per_subcell)
            self._cache[key] = solved.dipole.px
        return self._cache[key]
