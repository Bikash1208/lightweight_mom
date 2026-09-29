"""Part 0 task 3 (PHASE13_WEBAPP_PROMPT.md): the missing comparison Phase 12D's Task B left out
-- lookup-table **predicted** per-Floquet-order power vs the full-supercell MoM ground truth
already measured and stored in ``results/phase12d_local_periodicity/task_b.json``, for N = 6/8/12.

**Why.** Phase 12D found that no tested N reaches 90% of reflected power in the intended
(blazed-grating) order (N=6: 0.529, N=8: 0.411, N=12: 0.475) and root-caused it to this design's
limited achievable phase span (186 deg of the 360 deg a full anomalous reflector needs). That
measured the DESIGN, not the lookup-table forward MODEL: a lookup table built from a curve that
cannot reach a full 2*pi cycle would produce this shortfall even if local periodicity were a
perfect physical assumption. This module supplies the missing control: an idealized, purely
analytic prediction of per-order power made from the SAME lookup-table entries used to build the
gradient (no coupling, no full-wave solve) -- if it is close to the measured full-supercell MoM
result, the lookup-table/local-periodicity forward model is vindicated (the design, not the model,
is the limit); if it is far off, the model itself is inadequate on top of the design limitation.

**Method.** The standard array-factor/grating argument: model the aperture's reflectivity as a
periodic (period = N cells), piecewise-constant function of x, r(x) = r_j for x in cell j, where
r_j = |S11_j| * exp(j * phase_j) is cell j's ISOLATED reflection coefficient (the same curve entry
`phase12d.aperture_test.build_lookup_curve` already looked up to build the gradient -- reused here
unmodified, so the grid is bit-identical to what was actually solved). Its Floquet expansion in the
supercell's own lattice (CONVENTIONS.md's forward-transform sign, `f_mn = (1/N) sum_j f_j
exp(+j*2*pi*m*j/N)`, matching `RooftopBasis.current_coefficients`'s own convention) gives the
predicted field amplitude of order m directly; squaring gives the predicted power fraction, exactly
matching `floquet_power.floquet_order_powers`'s own normalization at (0,0) for a uniform array
(N=1, or every r_j identical: all power predicted in (0,0), reproducing |S11|**2 exactly -- this
module's own trivial-case check, `tests/test_phase13_order_power_prediction.py`).

Does not edit `phase12d.aperture_test.ApertureTestResult` (whose `predicted_intended_phase_deg`
field is left `None` there, a Phase 12D acknowledged gap) -- this is new, additive analysis that
reuses its curve-building and grid-building functions unmodified and reads its already-measured
`task_b.json` output rather than re-running the expensive full-supercell MoM solve.
"""

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..phase12c.locality_study import build_periodic_phase_gradient_supercell
from ..phase12d.aperture_test import build_lookup_curve
from ..phase12d.run import DEFAULT_APERTURE_L_FAMILY_MM, DEFAULT_APERTURE_MESH_CELLS_PER_SUBCELL

TASK_B_JSON = Path(__file__).resolve().parent.parent.parent / "results" / "phase12d_local_periodicity" / "task_b.json"


@dataclass(frozen=True)
class PredictedOrderPowers:
    n_cells: int
    orders: tuple  # int, the supercell-lattice order index m (y is always 0 -- 1-D gradient)
    predicted_fraction: np.ndarray  # same length/order as `orders`
    row_l_mm: tuple  # the per-cell L sequence actually used (bit-identical to the MoM grid)

    def fraction_at(self, m):
        return float(self.predicted_fraction[self.orders.index(m)])


def predict_order_powers(
    n_cells,
    mesh_cells_per_subcell=DEFAULT_APERTURE_MESH_CELLS_PER_SUBCELL,
    l_values_mm=DEFAULT_APERTURE_L_FAMILY_MM,
    max_order=None,
):
    """The idealized (local-periodicity, no-coupling) prediction for the SAME N-cell gradient
    Phase 12D actually solved (reuses `build_lookup_curve` and
    `build_periodic_phase_gradient_supercell` unmodified, so `row_l_mm` reproduces Phase 12D's
    own grid bit-for-bit)."""
    curve = build_lookup_curve(mesh_cells_per_subcell, l_values_mm)
    grid_mm, _, _, _ = build_periodic_phase_gradient_supercell(curve, n_cells, mesh_cells_per_subcell)
    row_l_mm = tuple(l for (l, _) in grid_mm[0])

    lookup = {round(c.l_mm, 9): c for c in curve}
    r = np.array(
        [
            lookup[round(l_mm, 9)].s11_mag * np.exp(1j * np.radians(lookup[round(l_mm, 9)].s11_phase_deg))
            for l_mm in row_l_mm
        ]
    )
    n = len(r)
    max_order = n if max_order is None else max_order
    orders = tuple(range(-max_order, max_order + 1))
    j = np.arange(n)
    fractions = []
    for m in orders:
        c_m = np.sum(r * np.exp(1j * 2.0 * np.pi * m * j / n)) / n
        fractions.append(abs(c_m) ** 2)
    return PredictedOrderPowers(n_cells, orders, np.array(fractions), row_l_mm)


def load_measured_task_b(path=TASK_B_JSON):
    """The already-measured full-supercell MoM result Phase 12D wrote (not re-run here)."""
    return json.loads(Path(path).read_text())


@dataclass(frozen=True)
class ModelVsDesignVerdict:
    n_cells: int
    predicted_intended_fraction: float
    measured_intended_fraction: float
    absolute_difference: float
    predicted_specular_fraction: float
    measured_specular_fraction: float


def compare_predicted_vs_measured(measured_task_b=None):
    """For every N in `task_b.json`'s `requested_sizes` that was actually solved: predicted vs
    measured fraction in the intended order and in the specular order. Close agreement (this
    module's own finding, not assumed) says the ~50% shortfall Phase 12D measured is a design
    (phase-span) limit that even an idealized, coupling-free lookup-table prediction reproduces --
    i.e. the lookup-table forward MODEL is not the bottleneck for apertures like this one."""
    data = measured_task_b if measured_task_b is not None else load_measured_task_b()
    verdicts = []
    for entry in data["results"]:
        n_cells = entry["n_cells"]
        predicted = predict_order_powers(n_cells, entry["mesh_cells_per_subcell"])
        intended_m = entry["intended_order"][0]  # (m, 0); gradient is along x
        pred_intended = predicted.fraction_at(intended_m)
        pred_specular = predicted.fraction_at(0)
        verdicts.append(
            ModelVsDesignVerdict(
                n_cells=n_cells,
                predicted_intended_fraction=pred_intended,
                measured_intended_fraction=entry["fraction_intended"],
                absolute_difference=abs(pred_intended - entry["fraction_intended"]),
                predicted_specular_fraction=pred_specular,
                measured_specular_fraction=entry["fraction_specular"],
            )
        )
    return tuple(verdicts)
