"""The standard reflectarray required-phase-map formula and cell selection.

**Required phase** (e.g. Huang & Encinar, *Reflectarray Antennas*, Sec. 1.3): a feed at
`(xf, yf, zf)` illuminates the array in the z=0 plane; cell `(x, y)` must impose a reflection
phase that cancels the feed's own spherical-wave phase delay AND adds the progressive phase of the
desired planar output beam `(theta0, phi0)`:

    phase_required(x, y) = wrap( k0 * [ R_feed(x, y) - (x*sin(theta0)*cos(phi0) + y*sin(theta0)*sin(phi0)) ] )
    R_feed(x, y) = sqrt((x-xf)^2 + (y-yf)^2 + zf^2)

(a constant additive phase reference is free -- any single global offset produces the same
far-field pattern; not fixed here since no project convention picks one).

Feed illumination taper: the standard `cos^q(theta_feed)` amplitude pattern, `theta_feed` = the
angle between the feed boresight (assumed -z, pointing at the array centre) and the ray to
`(x, y, 0)`.
"""

from dataclasses import dataclass

import numpy as np

from ..frequency import free_space_wavenumber
from ..labels import LOCAL_PERIODICITY_APPROXIMATION, phase_coverage_shortfall, steep_gradient
from ..phase import wrap_deg
from ..phase8d_rev.response_matching import nearest_response_match

DESIGN_STEP_LIMIT_DEG_PER_CELL = 5.8  # Phase 12D's own measured, off-resonance, design-tolerance limit

# PHASE13B A4 (live review): a library was reporting `reachable: true` for every cell just because
# SOME candidate existed, even at 25-50deg of phase error -- there is no measured physical limit
# for "how much phase-quantization error is too much" the way DESIGN_STEP_LIMIT_DEG_PER_CELL above
# is a measured cell-to-cell gradient tolerance; this is a documented design-quality policy default
# (roughly a coarse 12-level phase quantizer's worst-case bin, 360/12), not a fitted or measured
# number -- override per call if a design needs a stricter (or looser) bar.
PHASE_ERROR_UNREACHABLE_THRESHOLD_DEG = 30.0


def feed_distance_m(x_m, y_m, feed_xyz_m):
    xf, yf, zf = feed_xyz_m
    return np.sqrt((x_m - xf) ** 2 + (y_m - yf) ** 2 + zf**2)


def required_phase_deg(x_m, y_m, frequency_hz, feed_xyz_m, theta0_deg, phi0_deg):
    """`phase_required(x, y)` [deg], wrapped to (-180, 180]."""
    k0 = free_space_wavenumber(frequency_hz)
    r_feed = feed_distance_m(x_m, y_m, feed_xyz_m)
    theta0, phi0 = np.radians(theta0_deg), np.radians(phi0_deg)
    steering_term = x_m * np.sin(theta0) * np.cos(phi0) + y_m * np.sin(theta0) * np.sin(phi0)
    phase_rad = k0 * (r_feed - steering_term)
    return wrap_deg(np.degrees(phase_rad))


def feed_illumination_taper(x_m, y_m, feed_xyz_m, q=1.0):
    """`cos^q(theta_feed)` amplitude taper, `theta_feed` measured from the feed's boresight
    (assumed to point at the array's geometric centre, i.e. straight down -z from the feed)."""
    xf, yf, zf = feed_xyz_m
    r_feed = feed_distance_m(x_m, y_m, feed_xyz_m)
    cos_theta_feed = zf / r_feed  # boresight along -z: cos(theta_feed) = zf / |r|
    cos_theta_feed = np.clip(cos_theta_feed, 0.0, 1.0)
    return cos_theta_feed**q


@dataclass(frozen=True)
class CellAssignment:
    i: int
    j: int
    x_m: float
    y_m: float
    required_phase_deg: float
    chosen_l_mm: "float | None"
    chosen_w_mm: "float | None"
    achieved_phase_deg: "float | None"
    phase_error_deg: "float | None"
    library_id: "str | None"
    phase_coverage_ok: bool  # False if no candidate existed, or its phase_error_deg exceeds threshold


@dataclass(frozen=True)
class PhaseErrorStats:
    mean_deg: float
    rms_deg: float
    max_deg: float
    n_cells: int
    n_unreachable: int  # phase_coverage_ok is False


@dataclass(frozen=True)
class ReflectarrayDesign:
    x_m: np.ndarray
    y_m: np.ndarray
    required_phase_deg: np.ndarray
    assignments: "tuple[CellAssignment, ...]"
    max_phase_gradient_deg_per_cell: float
    phase_error_stats: "PhaseErrorStats | None"
    labels: tuple


def assign_cells_from_library(
    x_m, y_m, required_phase_deg_grid, library_candidates, threshold_deg=PHASE_ERROR_UNREACHABLE_THRESHOLD_DEG
):
    """`library_candidates`: an iterable of `(record_id, L_mm, W_mm, library_id, complex_S11)` --
    the same schema `phase8d_rev.response_matching.nearest_response_match` expects. For each grid
    point, picks the nearest-phase candidate (magnitude-aware: `nearest_response_match`'s own
    complex-plane distance already penalises a magnitude mismatch, not phase alone) and returns
    every assignment plus the achieved-vs-required error. `library_candidates` is materialized
    once and reused for every cell (the caller's library doesn't change per cell). A cell is
    `phase_coverage_ok` only if a candidate existed AND its phase error is within `threshold_deg`
    (PHASE13B A4 -- previously any candidate at all, however large its error, counted as
    "reachable")."""
    candidates = list(library_candidates)
    ny, nx = np.shape(x_m)
    assignments = []
    for i in range(ny):
        for j in range(nx):
            required = required_phase_deg_grid[i, j]
            required_complex = np.exp(1j * np.radians(required))  # unit magnitude: phase-only match
            matches = nearest_response_match(required_complex, candidates)
            best = matches[0] if matches else None
            assignments.append(
                CellAssignment(
                    i=i,
                    j=j,
                    x_m=float(x_m[i, j]),
                    y_m=float(y_m[i, j]),
                    required_phase_deg=float(required),
                    chosen_l_mm=best.L_mm if best else None,
                    chosen_w_mm=best.W_mm if best else None,
                    achieved_phase_deg=float(np.degrees(np.angle(best.complex_value))) if best else None,
                    phase_error_deg=best.phase_error_deg if best else None,
                    library_id=best.library_id if best else None,
                    phase_coverage_ok=best is not None and best.phase_error_deg <= threshold_deg,
                )
            )
    return tuple(assignments)


def phase_error_stats(assignments):
    """Summary stats over every assignment with a defined `phase_error_deg` (PHASE13B A4)."""
    errors = np.array([a.phase_error_deg for a in assignments if a.phase_error_deg is not None])
    if errors.size == 0:
        return None
    return PhaseErrorStats(
        mean_deg=float(np.mean(errors)),
        rms_deg=float(np.sqrt(np.mean(errors**2))),
        max_deg=float(np.max(errors)),
        n_cells=len(assignments),
        n_unreachable=sum(1 for a in assignments if not a.phase_coverage_ok),
    )


def max_phase_gradient_deg_per_cell(phase_deg_grid):
    """Largest |wrapped difference| between any two adjacent cells (row- or column-neighbours) --
    the same per-cell step-size quantity Phase 12D's local-periodicity study measured."""
    from ..phase import wrapped_difference_deg

    grid = np.asarray(phase_deg_grid, dtype=float)
    dx = wrapped_difference_deg(grid[:, 1:], grid[:, :-1])
    dy = wrapped_difference_deg(grid[1:, :], grid[:-1, :])
    return float(max(np.max(np.abs(dx)) if dx.size else 0.0, np.max(np.abs(dy)) if dy.size else 0.0))


def design_reflectarray(
    x_m,
    y_m,
    frequency_hz,
    feed_xyz_m,
    theta0_deg,
    phi0_deg,
    library_candidates,
    phase_error_threshold_deg=PHASE_ERROR_UNREACHABLE_THRESHOLD_DEG,
):
    required = required_phase_deg(x_m, y_m, frequency_hz, feed_xyz_m, theta0_deg, phi0_deg)
    assignments = assign_cells_from_library(x_m, y_m, required, library_candidates, phase_error_threshold_deg)
    achieved_grid = np.array([a.achieved_phase_deg for a in assignments]).reshape(np.shape(x_m))
    max_grad = max_phase_gradient_deg_per_cell(achieved_grid)
    stats = phase_error_stats(assignments)
    labels = [LOCAL_PERIODICITY_APPROXIMATION]
    if max_grad > DESIGN_STEP_LIMIT_DEG_PER_CELL:
        labels.append(steep_gradient(DESIGN_STEP_LIMIT_DEG_PER_CELL))
    if stats is not None and stats.n_unreachable > 0:
        worst = max(a.phase_error_deg for a in assignments if not a.phase_coverage_ok)
        labels.append(phase_coverage_shortfall(worst))
    return ReflectarrayDesign(
        x_m=x_m,
        y_m=y_m,
        required_phase_deg=required,
        assignments=assignments,
        max_phase_gradient_deg_per_cell=max_grad,
        phase_error_stats=stats,
        labels=tuple(labels),
    )
