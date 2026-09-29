"""Full-supercell MoM verification of a periodic-gradient lookup-table design: reuses
`phase12d.aperture_test.run_aperture_test` (already validated on known-answer cases, Phase 12D)
and `phase13.order_power_prediction.predict_order_powers` (Part 0's idealized, no-coupling
prediction from the SAME lookup table) unmodified -- this module only adds a cost estimate/refusal
gate and combines the two results into one report.

**Cost model**: fit from Phase 12D's own three MEASURED runtimes at this reduced-floquet
convention (N=6: 7.0s, N=8: 18.5s, N=12: 124.2s -- `results/phase12d_local_periodicity/task_b.json`,
not re-measured here). The fit is a plain power law through the two largest points (the steepest,
most conservative segment) extrapolated forward; N=12 is the largest size Phase 12D actually
solved, so anything larger is an ESTIMATE, clearly labelled as such, not a new measurement.
"""

from dataclasses import dataclass

import numpy as np

from ..phase12d.aperture_test import run_aperture_test
from ..phase12d.run import DEFAULT_APERTURE_L_FAMILY_MM as DEFAULT_L_FAMILY_MM
from ..phase12d.run import DEFAULT_APERTURE_MESH_CELLS_PER_SUBCELL as DEFAULT_MESH_CELLS_PER_SUBCELL
from ..phase13.order_power_prediction import predict_order_powers

# Phase 12D's own measured (N, runtime_s) points at the reduced-floquet convention.
_MEASURED_RUNTIMES = ((6, 6.997816), (8, 18.514908), (12, 124.244273))
MAX_MEASURED_N = 12         # largest N Phase 12D actually solved -- anything above is extrapolated
HARD_REFUSE_ABOVE_N = 20    # beyond this, extrapolation itself is not trusted; refuse outright


@dataclass(frozen=True)
class CostEstimate:
    n_cells: int
    estimated_runtime_s: float
    is_measured: bool          # True if n_cells is one of the 3 points Phase 12D actually measured
    is_extrapolated: bool
    refused: bool
    refusal_reason: "str | None"


def estimate_cost(n_cells):
    for measured_n, measured_t in _MEASURED_RUNTIMES:
        if n_cells == measured_n:
            return CostEstimate(n_cells, measured_t, True, False, False, None)

    if n_cells > HARD_REFUSE_ABOVE_N:
        return CostEstimate(n_cells, float("inf"), False, False, True,
                             f"n_cells={n_cells} > {HARD_REFUSE_ABOVE_N}: no measured data this large, "
                             f"and extrapolation this far is not trusted (Phase 12D's own largest solved "
                             f"size was N={MAX_MEASURED_N}).")

    # Power-law fit through the two largest measured points (steepest segment -> a conservative,
    # not optimistic, extrapolation): t = a * N^p.
    (n1, t1), (n2, t2) = _MEASURED_RUNTIMES[1], _MEASURED_RUNTIMES[2]
    p = np.log(t2 / t1) / np.log(n2 / n1)
    a = t1 / n1**p
    estimated_t = a * n_cells**p
    return CostEstimate(n_cells, float(estimated_t), False, True, False, None)


@dataclass(frozen=True)
class VerificationResult:
    n_cells: int
    aperture_result: object          # phase12d.aperture_test.ApertureTestResult
    predicted_order_powers: object   # phase13.order_power_prediction.PredictedOrderPowers
    predicted_vs_measured_intended_fraction_diff: float
    cost: CostEstimate


def run_verification(n_cells, mesh_cells_per_subcell=DEFAULT_MESH_CELLS_PER_SUBCELL,
                      l_values_mm=DEFAULT_L_FAMILY_MM, force=False):
    """Runs the full-supercell MoM check plus the idealized lookup-table prediction for the same
    N/mesh/L-family. Raises `ValueError` if the cost estimate says refuse, unless `force=True`
    (an explicit override the app UI gates behind a confirmation, never a default)."""
    cost = estimate_cost(n_cells)
    if cost.refused and not force:
        raise ValueError(cost.refusal_reason)

    aperture_result = run_aperture_test(n_cells, mesh_cells_per_subcell, l_values_mm)
    predicted = predict_order_powers(n_cells, mesh_cells_per_subcell, l_values_mm)
    intended_m = aperture_result.intended_order[0]
    predicted_fraction = predicted.fraction_at(intended_m)
    diff = abs(predicted_fraction - aperture_result.fraction_intended)

    return VerificationResult(
        n_cells=n_cells, aperture_result=aperture_result, predicted_order_powers=predicted,
        predicted_vs_measured_intended_fraction_diff=diff, cost=cost,
    )
