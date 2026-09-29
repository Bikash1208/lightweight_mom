"""Sec. 38-40: first-stage geometry inversion. Given a FIXED environment (8 neighbor geometries)
and a target complex response, search the center cell's own geometry over the validated 18x18
sweep grid for the value the surrogate predicts comes closest to the target (Sec. 38's J = |S_target
- S_hat|^2), then VERIFY the winning candidate with the real heterogeneous MoM solver (Sec. 39: "the
surrogate must NOT be the sole judge of its own inverse predictions") -- exactly the same
`grid_solver.solve_grid` pipeline used everywhere else in Phase 10, not a new solver.
"""

from dataclasses import dataclass

from . import environment_geometry as eg
from .dataset import ALL_POSITIONS, EnvironmentSample
from .grid_solver import solve_grid

# Sec. 40: "no physical geometry can reproduce a target" is reported, not forced -- classified via
# the SAME self-referential "at least half the target's own magnitude" test this project already
# uses (`phase9.local_periodicity.SUBSTANTIAL_DEVIATION_RELATIVE`), reused rather than a new
# invented inversion-specific tolerance.
from ..phase9.local_periodicity import SUBSTANTIAL_DEVIATION_RELATIVE

TARGET_UNREACHABLE = "TARGET_UNREACHABLE"


@dataclass(frozen=True)
class InverseResult:
    target: complex
    neighbors_lw_mm: dict
    best_candidate_lw_mm: tuple
    surrogate_predicted_at_candidate: object
    mom_achieved_at_candidate: complex
    mom_relative_error: float
    status: str   # "REACHED" | TARGET_UNREACHABLE


def _grid_values_mm():
    import numpy as np
    return [float(v) for v in
            (round(eg.DOMAIN_L_MIN_MM + i * eg.DOMAIN_STEP_MM, 6)
             for i in range(int(round((eg.DOMAIN_L_MAX_MM - eg.DOMAIN_L_MIN_MM) / eg.DOMAIN_STEP_MM)) + 1))]


def invert_target(target, neighbors_lw_mm, frequency_hz, surrogate, period_m, substrate, metal,
                   mesh_cells_per_subcell=8, candidate_lw_values_mm=None):
    """Sec. 38-39: proposes a center geometry via surrogate search, then checks it against the
    real heterogeneous MoM solver in the same fixed environment."""
    values = candidate_lw_values_mm or _grid_values_mm()
    best_lw, best_err = None, float("inf")
    for l_mm in values:
        for w_mm in values:
            sample = EnvironmentSample("candidate", "n/a", (l_mm, w_mm), neighbors_lw_mm,
                                        "inverse_candidate", frequency_hz, 0j, "inverse")
            pred = surrogate.predict(sample)
            if not isinstance(pred, complex):
                continue
            err = abs(pred - target)
            if err < best_err:
                best_err, best_lw = err, (l_mm, w_mm)

    if best_lw is None:
        return InverseResult(target, neighbors_lw_mm, None, None, None, float("inf"), TARGET_UNREACHABLE)

    sample = EnvironmentSample("candidate", "n/a", best_lw, neighbors_lw_mm,
                                "inverse_candidate", frequency_hz, 0j, "inverse")
    predicted_at_candidate = surrogate.predict(sample)

    cardinal = {d: neighbors_lw_mm[d] for d in ("N", "E", "S", "W")}
    diagonal = {d: neighbors_lw_mm[d] for d in ("NE", "NW", "SE", "SW")}
    grid_mm = eg.build_9cell_environment(best_lw, cardinal_mm=cardinal, diagonal_mm=diagonal)
    grid_m = eg.grid_to_m(grid_mm)
    sg = solve_grid(3, period_m, substrate, metal, grid_m, frequency_hz, mesh_cells_per_subcell)
    achieved = sg.dipole_at(1, 1).px

    relative = abs(achieved - target) / abs(target) if abs(target) > 0 else float("inf")
    status = TARGET_UNREACHABLE if relative >= SUBSTANTIAL_DEVIATION_RELATIVE else "REACHED"
    return InverseResult(target, neighbors_lw_mm, best_lw, predicted_at_candidate, achieved, relative, status)
