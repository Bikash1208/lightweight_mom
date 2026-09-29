"""Sec. 21-22 (Phase 11P/Q): deterministic, non-ML geometry inversion. Candidates are ranked purely
by a physically-computed quantity -- the ISOLATED-CELL reference response at each candidate
geometry (Sec. 21's own permitted "interpolation between physically calculated points" / "bounded
local search", zero fitting, zero learned parameters) -- and the winner is always verified against
the REAL heterogeneous MoM solve in the actual fixed environment (Sec. 21: "every candidate must
ultimately be verified by the direct MoM solver"). No regression, no neural network, no random
forest: `EXHAUSTIVE DISCRETE SEARCH` over the validated 18-value sweep grid, full stop.
"""

from dataclasses import dataclass

from ..phase9.local_periodicity import SUBSTANTIAL_DEVIATION_RELATIVE
from ..phase10.environment_geometry import (DOMAIN_L_MAX_MM, DOMAIN_L_MIN_MM, DOMAIN_STEP_MM,
                                             grid_to_m)
from .reference_solver import solve_grid, solve_isolated

TARGET_UNREACHABLE = "TARGET_UNREACHABLE"
FOUND = "FOUND"


def _grid_values_mm():
    n = int(round((DOMAIN_L_MAX_MM - DOMAIN_L_MIN_MM) / DOMAIN_STEP_MM)) + 1
    return [round(DOMAIN_L_MIN_MM + i * DOMAIN_STEP_MM, 6) for i in range(n)]


@dataclass(frozen=True)
class InverseResult:
    target: complex
    best_candidate_lw_mm: object
    isolated_predicted_at_candidate: object
    mom_achieved_at_candidate: object
    mom_relative_error: float
    status: str


def invert_target_direct(target, environment_grid_mm, center_i, center_j, period_m, substrate, metal,
                          frequency_hz, mesh_cells_per_subcell=8, candidate_lw_values_mm=None,
                          substantial=SUBSTANTIAL_DEVIATION_RELATIVE):
    """``environment_grid_mm``: the FIXED surrounding pattern (the center cell's own current value
    is ignored/overwritten by each candidate). Ranks every candidate (L, W) on the sweep grid by
    its ISOLATED-CELL response distance to ``target`` (a pure, deterministic physics computation,
    no fitting), then MoM-verifies the single best candidate in the real, fixed heterogeneous
    environment."""
    values = candidate_lw_values_mm or _grid_values_mm()
    best_lw, best_err, best_pred = None, float("inf"), None
    for l_mm in values:
        for w_mm in values:
            pred, _ = solve_isolated(period_m, l_mm * 1e-3, w_mm * 1e-3, substrate, metal, frequency_hz,
                                      mesh_cells_per_subcell)
            err = abs(pred - target)
            if err < best_err:
                best_err, best_lw, best_pred = err, (l_mm, w_mm), pred

    if best_lw is None:
        return InverseResult(target, None, None, None, float("inf"), TARGET_UNREACHABLE)

    n = len(environment_grid_mm)
    grid_mm = [row[:] for row in environment_grid_mm]
    grid_mm[center_i][center_j] = best_lw
    grid_m = grid_to_m(grid_mm)
    sg = solve_grid(n, period_m, substrate, metal, grid_m, frequency_hz, mesh_cells_per_subcell)
    achieved = sg.response_at(center_i, center_j)

    relative = abs(achieved - target) / abs(target) if abs(target) > 0 else float("inf")
    status = TARGET_UNREACHABLE if relative >= substantial else FOUND
    return InverseResult(target, best_lw, best_pred, achieved, relative, status)
