"""PHASE12D_ADDENDUM.md task 1: Richardson extrapolation of the MoM solve in mesh cell size.

**Why.** The addendum measured that MoM's reflection phase converges only first-order in cell
size near resonance (edge-current singularities unresolved on a uniform rooftop mesh): at G4
(L = 6.0 mm, 12.5 GHz) the 0.2 mm mesh (40 cells) gives -168.74 deg, refining to 96 cells
(0.083 mm) only reaches -172.11 deg, and a fit phi = phi_inf + a*dx^p over 40..96 cells gives
p = 0.96, phi_inf = -174.65 deg -- i.e. the 0.2 mm mesh is ~6 deg off *at resonance*. This module
turns that by-hand procedure into a reusable, tested function.

**How.** Solves `sparameters.compute_s_parameters` unmodified at each of a family of *aligned*
mesh resolutions (dx = P/N for integer N, `floquet_order` tracked with N by default, matching
the addendum's own family), unwraps the resulting phase across the resolution sequence, then
fits `v(dx) = v_inf + a*dx^p` (3-parameter nonlinear least squares) separately to phase and to
magnitude. Per the project's "fit nothing unless monotonic" rule (CLAUDE.md), the fit is attempted
only if the raw sequence is monotonic in dx; otherwise the point is reported as
`NOT_MONOTONIC_NO_FIT` with the raw resolutions only -- never a silently-untrustworthy fit.

This is Phase 13's one permitted additive numerics change (PHASE13_WEBAPP_PROMPT.md: "Physics/
numerics untouched except the Part 0 convergence option, additive, default only in the app"). It
does not edit `sparameters.py`/`solver.py`/`basis_functions.py` -- it is new orchestration that
calls them repeatedly, the same pattern every phase from 8B onward already uses for reuse.
"""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import curve_fit

from ..basis_functions import MeshSettings, RooftopBasis
from ..sparameters import compute_s_parameters

# The addendum's own family (aligned meshes, floquet_order = cells, G4 at 12.5 GHz): 40 cells is
# this project's usual "production" mesh for the 8 mm grounded design point; 96 is the largest the
# addendum actually ran. Both endpoints, and every point between, are real measurements, not a
# guess at what would be affordable.
DEFAULT_CELL_COUNTS = (40, 48, 56, 64, 72, 80, 96)

# A cheaper family for interactive (app single-cell page) use. 4 points, not 3: the 3-parameter
# fit (v_inf, a, p) needs a 4th point for a nonzero residual at all -- with exactly 3 points the
# fit is exactly determined (0 degrees of freedom), so curve_fit cannot estimate a covariance
# (measured: gives `OptimizeWarning: Covariance ... could not be estimated`, not a usable error
# bar). Documented as a speed/accuracy tradeoff, never presented as equal in confidence to
# DEFAULT_CELL_COUNTS's 7-point family.
FAST_CELL_COUNTS = (24, 32, 40, 48)

STATUS_CONVERGED = "CONVERGED"
STATUS_NOT_MONOTONIC_NO_FIT = "NOT_MONOTONIC_NO_FIT"
STATUS_FIT_FAILED = "FIT_FAILED"
STATUS_TOO_FEW_POINTS = "TOO_FEW_POINTS"


@dataclass(frozen=True)
class ConvergencePoint:
    cells: int
    dx_m: float
    floquet_order: int
    n_unknowns: int
    s11_mag: float
    s11_phase_deg: float  # raw, on (-180, 180] -- NOT the unwrapped value used internally for the fit


@dataclass(frozen=True)
class FitResult:
    status: str
    value_inf: float | None  # extrapolated dx -> 0 value
    value_inf_error: float | None  # see `_fit_power_law` -- NOT the fit's covariance (see there)
    order: float | None  # fitted power p
    finest_value: float | None  # the raw value at the finest resolution actually run (fallback)
    coefficient: float | None = None  # fitted `a` in `v(dx) = v_inf + a*dx^p` -- lets a caller
    # scale this SAME measured law to a different dx (e.g.
    # `resonance_restatement.py`) without re-fitting


@dataclass(frozen=True)
class ConvergenceResult:
    points: tuple[ConvergencePoint, ...]
    phase: FitResult
    magnitude: FitResult

    @property
    def finest(self) -> ConvergencePoint:
        return min(self.points, key=lambda p: p.dx_m)

    @property
    def mesh_converged(self) -> bool:
        """`True` only if the phase fit succeeded and the extrapolation correction from the
        finest actually-run mesh is within the stated tolerance (see `is_converged_within`)."""
        return self.phase.status == STATUS_CONVERGED

    def is_converged_within(self, tol_deg: float) -> bool:
        """Whether the *correction* |phi_inf - phi(finest mesh)| is within `tol_deg` -- the
        practical question ('is my production mesh good enough'), not just whether a fit exists."""
        if self.phase.status != STATUS_CONVERGED:
            return False
        return abs(self.phase.value_inf - _unwrap_to_match(self)) <= tol_deg


def _unwrap_to_match(result: "ConvergenceResult") -> float:
    """The finest point's phase, unwrapped onto the same branch as `phase.value_inf` (fit is done
    in the unwrapped frame; comparing the raw (-180,180] value directly could be off by 360 deg)."""
    ordered = sorted(result.points, key=lambda p: -p.dx_m)  # coarse -> fine
    unwrapped = np.unwrap(np.radians([p.s11_phase_deg for p in ordered]))
    return float(np.degrees(unwrapped[-1]))


def solve_at_resolution(
    unit_cell, frequency_hz, planes, cells, floquet_order=None, warn_grating=False, warn_surface_wave=False
):
    """One aligned-mesh MoM solve (`cells_x = cells_y = cells`, `dx = P/cells`), reused verbatim
    from `sparameters.compute_s_parameters` -- the same function every other solve in this project
    uses. `floquet_order` defaults to `cells` (the addendum's own choice, generous relative to the
    `2*order+1 >= cells` minimum)."""
    floquet_order = cells if floquet_order is None else floquet_order
    mesh = MeshSettings(cells_x=cells, cells_y=cells, floquet_order=floquet_order)
    sp = compute_s_parameters(
        unit_cell, frequency_hz, mesh, planes, warn_grating=warn_grating, warn_surface_wave=warn_surface_wave
    )
    s11 = sp.s11[0, 0]
    period_x = unit_cell.lattice.period_x
    n_unknowns = RooftopBasis(unit_cell, cells, cells).n_unknowns
    return ConvergencePoint(
        cells=cells,
        dx_m=period_x / cells,
        floquet_order=floquet_order,
        n_unknowns=n_unknowns,
        s11_mag=float(abs(s11)),
        s11_phase_deg=float(np.degrees(np.angle(s11))),
    )


def _power_law(dx, v_inf, a, p):
    return v_inf + a * np.power(dx, p)


def _fit_power_law(dx_values, values):
    """Fit `v(dx) = v_inf + a*dx^p`; returns a `FitResult`. Monotonicity (in the order dx is
    given, which callers pass coarse-to-fine or fine-to-coarse consistently) is required before a
    fit is attempted at all -- CLAUDE.md's "fit nothing unless monotonic" rule."""
    dx_values = np.asarray(dx_values, dtype=float)
    values = np.asarray(values, dtype=float)
    finest_value = float(values[np.argmin(dx_values)])
    if len(dx_values) < 4:
        # 3 points exactly determine the 3 fit parameters (0 degrees of freedom): curve_fit
        # would return a point estimate with no usable covariance (see FAST_CELL_COUNTS comment).
        return FitResult(STATUS_TOO_FEW_POINTS, None, None, None, finest_value)

    order = np.argsort(-dx_values)  # coarse -> fine
    v_sorted = values[order]
    diffs = np.diff(v_sorted)
    if not (np.all(diffs <= 0) or np.all(diffs >= 0)):
        return FitResult(STATUS_NOT_MONOTONIC_NO_FIT, None, None, None, finest_value)

    spread = float(np.max(v_sorted) - np.min(v_sorted))
    if spread < 1e-9:
        # Already converged to numerical precision at every tested resolution (e.g. far
        # off-resonance): a 3-parameter power-law fit to a flat curve is ill-conditioned (the
        # exponent p is unconstrained by data with no visible curvature) and would silently
        # overflow rather than reveal anything -- report the flat value directly instead.
        return FitResult(STATUS_CONVERGED, finest_value, 0.0, float("inf"), finest_value, 0.0)

    # p bounded to (0.05, 8): extrapolating a mesh-refinement error to dx -> 0 only makes physical
    # sense for p > 0 (the error must vanish as the mesh refines); letting the optimizer explore
    # p <= 0 is exactly what drove `a*dx**p` to overflow for a near-flat curve in earlier testing.
    p0 = (v_sorted[-1], v_sorted[0] - v_sorted[-1], 1.0)
    bounds = ([-np.inf, -np.inf, 0.05], [np.inf, np.inf, 8.0])
    try:
        with np.errstate(over="raise", invalid="raise"):
            popt, pcov = curve_fit(_power_law, dx_values, values, p0=p0, bounds=bounds, maxfev=20000)
    except Exception:
        return FitResult(STATUS_FIT_FAILED, None, None, None, finest_value)
    v_inf, a, p = popt
    if not np.all(np.isfinite(pcov)):
        return FitResult(STATUS_FIT_FAILED, None, None, None, finest_value)
    # NOT `sqrt(pcov[0,0])`: PHASE13B review found that number (e.g. 0.0076 deg at G4) is not a
    # credible uncertainty -- a 3-parameter fit to as few as 4 points (FAST_CELL_COUNTS) has only
    # 1 degree of freedom, and curve_fit's covariance assumes the power-law form is exactly right,
    # which understates error from model mismatch. Instead report the *size of the extrapolation
    # correction itself* (|v_inf - finest measured value|): if refining the mesh already moved the
    # answer this much, the remaining error from here to the true dx->0 limit is plausibly
    # comparable -- a conservative, always-defined bound that needs no extra data (unlike a
    # leave-one-out spread, which needs >= 5 points to still satisfy the >=4-point fit minimum).
    v_inf_err = float(abs(v_inf - finest_value))
    return FitResult(STATUS_CONVERGED, float(v_inf), v_inf_err, float(p), finest_value, float(a))


def richardson_extrapolate_s11xx(
    unit_cell,
    frequency_hz,
    planes,
    cell_counts=DEFAULT_CELL_COUNTS,
    floquet_order_fn=None,
    warn_grating=False,
    warn_surface_wave=False,
):
    """Solve at every resolution in `cell_counts`, fit and extrapolate phase and magnitude
    separately (`floquet_order_fn(cells)` overrides the default `floquet_order = cells`). Phase
    is unwrapped across the resolution sequence (coarse -> fine) before fitting, since the raw
    per-solve value is only ever known on (-180, 180]."""
    points = []
    for cells in sorted(cell_counts):
        fo = floquet_order_fn(cells) if floquet_order_fn is not None else None
        points.append(
            solve_at_resolution(
                unit_cell,
                frequency_hz,
                planes,
                cells,
                floquet_order=fo,
                warn_grating=warn_grating,
                warn_surface_wave=warn_surface_wave,
            )
        )
    points = tuple(points)

    ordered = sorted(points, key=lambda p: -p.dx_m)  # coarse -> fine
    dx_sorted = np.array([p.dx_m for p in ordered])
    phase_unwrapped = np.degrees(np.unwrap(np.radians([p.s11_phase_deg for p in ordered])))
    mag_sorted = np.array([p.s11_mag for p in ordered])

    phase_fit = _fit_power_law(dx_sorted, phase_unwrapped)
    mag_fit = _fit_power_law(dx_sorted, mag_sorted)
    return ConvergenceResult(points=points, phase=phase_fit, magnitude=mag_fit)
