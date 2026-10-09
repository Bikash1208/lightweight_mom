"""Turns a parsed `.fld` file (`fld_parser.FldGrid`) into a
`phase8e.reference_field_input.ReferenceFieldGrid` (reused unmodified), plus the validation and
sanity-check machinery PHASE13_WEBAPP_PROMPT.md's "Reference-field import" section asks for:
grid-regularity, local-wavenumber estimate vs the analytic TM0 baseline, and a cells-per-fringe
aliasing warning.

**2-D phase unwrapping**: `numpy.unwrap` applied along each axis in turn -- a standard, practical
approach, NOT a branch-cut-optimal algorithm. Adequate for a smooth, quasi-radial fringe pattern
like this project's own example file (empirically verified: reproduces the brief's own cited
~-225 rad/m plateau, `tests/test_reference_field_grid.py`), but can fail near strong local phase
singularities away from the origin -- documented, not hidden. The single feed-point NaN (at the
origin, by construction of a monopole/point-feed hologram) is filled by nearest-neighbor
interpolation *only* to stop it corrupting `numpy.unwrap`'s sequential accumulation along its row
and column; the filled value is never reported back to the caller (`FEED_REGION_MASKED`).
"""

from dataclasses import dataclass, replace

import numpy as np
from scipy.interpolate import griddata

from ..dispersion.analytic_tm0 import tm0_beta
from ..frequency import free_space_wavenumber
from ..phase8e.reference_field_input import ReferenceFieldGrid, phase_only_reference_field
from .fld_parser import FldGrid

RECOMMENDED_MIN_CELLS_PER_FRINGE = 6.0  # brief: "surface-wave holograms normally use ~ lambda/6-lambda/10"


def fld_to_reference_field_grid(phase_fld: FldGrid, magnitude_fld: "FldGrid | None" = None):
    """`phase_fld` must be a `CmplxPhase(...)` file. Without `magnitude_fld`, returns a phase-only
    grid (`REFERENCE_AMPLITUDE_ASSUMED_UNIFORM` is the caller's label to apply, not this
    function's -- it only reports what data it actually has). With a matching `magnitude_fld`
    (same grid), returns a full complex field."""
    if phase_fld.header.expression.kind != "phase_rad":
        raise ValueError(
            f"fld_to_reference_field_grid: phase_fld must be a CmplxPhase(...) export, "
            f"got kind={phase_fld.header.expression.kind!r}"
        )
    phase_deg = np.degrees(phase_fld.value)
    dx_m = phase_fld.x_m[0, 1] - phase_fld.x_m[0, 0] if phase_fld.x_m.shape[1] > 1 else None
    dy_m = phase_fld.y_m[1, 0] - phase_fld.y_m[0, 0] if phase_fld.y_m.shape[0] > 1 else None
    z_reference_m = phase_fld.header.min_m[2]

    if magnitude_fld is None:
        grid = phase_only_reference_field(phase_fld.x_m, phase_fld.y_m, phase_deg)
        return replace(grid, dx_m=dx_m, dy_m=dy_m, z_reference_m=z_reference_m, source=f"FILE:{phase_fld.path}")

    if magnitude_fld.value.shape != phase_fld.value.shape or not np.allclose(magnitude_fld.x_m, phase_fld.x_m):
        raise ValueError("fld_to_reference_field_grid: phase_fld and magnitude_fld are not on the same grid")
    ex = magnitude_fld.value * np.exp(1j * phase_fld.value)
    return ReferenceFieldGrid(
        x_m=phase_fld.x_m,
        y_m=phase_fld.y_m,
        Ex=ex,
        Ey=np.zeros_like(ex),
        amplitude_available=True,
        source=f"FILE:{phase_fld.path}+{magnitude_fld.path}",
        dx_m=dx_m,
        dy_m=dy_m,
        z_reference_m=z_reference_m,
        phase_deg=phase_deg,
    )


@dataclass(frozen=True)
class GridRegularityReport:
    regular: bool
    dx_m: "float | None"
    dy_m: "float | None"
    dx_relative_spread: float
    dy_relative_spread: float
    nx: int
    ny: int


def check_grid_regularity(fld_grid: FldGrid, rtol=1e-6):
    dx_steps = np.diff(fld_grid.x_m[0, :]) if fld_grid.x_m.shape[1] > 1 else np.array([0.0])
    dy_steps = np.diff(fld_grid.y_m[:, 0]) if fld_grid.y_m.shape[0] > 1 else np.array([0.0])
    dx_spread = float(np.std(dx_steps) / abs(np.mean(dx_steps))) if np.mean(dx_steps) else 0.0
    dy_spread = float(np.std(dy_steps) / abs(np.mean(dy_steps))) if np.mean(dy_steps) else 0.0
    regular = dx_spread <= rtol and dy_spread <= rtol
    return GridRegularityReport(
        regular=regular,
        dx_m=float(np.mean(dx_steps)),
        dy_m=float(np.mean(dy_steps)),
        dx_relative_spread=dx_spread,
        dy_relative_spread=dy_spread,
        nx=fld_grid.value.shape[1],
        ny=fld_grid.value.shape[0],
    )


def _fill_nan_for_unwrap(value):
    filled = value.copy()
    nan_mask = np.isnan(filled)
    if nan_mask.any():
        good_y, good_x = np.nonzero(~nan_mask)
        bad_y, bad_x = np.nonzero(nan_mask)
        filled[nan_mask] = griddata((good_y, good_x), filled[~nan_mask], (bad_y, bad_x), method="nearest")
    return filled, nan_mask


def unwrap_phase_2d(phase_rad):
    """Fills any NaN (see module docstring), then `numpy.unwrap` along axis 1 (x) then axis 0
    (y). Returns `(unwrapped, nan_mask)`; NaN positions are NOT re-masked into the unwrapped
    array (the fill is needed for every downstream gradient sample near that point too) -- callers
    that need to exclude the feed point use `nan_mask`."""
    filled, nan_mask = _fill_nan_for_unwrap(phase_rad)
    unwrapped = np.unwrap(filled, axis=1)
    unwrapped = np.unwrap(unwrapped, axis=0)
    return unwrapped, nan_mask


def local_wavenumber_field(fld_grid: FldGrid):
    """beta(x,y) = |grad(unwrapped phase)| [rad/m], a 2-D array the same shape as `fld_grid.value`.
    Requires a regular grid (`check_grid_regularity`) for `numpy.gradient`'s spacing arguments to
    be meaningful."""
    unwrapped, nan_mask = unwrap_phase_2d(fld_grid.value)
    dy_m = fld_grid.y_m[1, 0] - fld_grid.y_m[0, 0]
    dx_m = fld_grid.x_m[0, 1] - fld_grid.x_m[0, 0]
    dphi_dy, dphi_dx = np.gradient(unwrapped, dy_m, dx_m)
    beta = np.hypot(dphi_dx, dphi_dy)
    beta[nan_mask] = np.nan
    return beta


def radial_beta_profile(fld_grid: FldGrid, feed_xy_m=(0.0, 0.0), n_bins=30):
    """Median beta(r) in `n_bins` radial bins from the feed point -- the brief's own worked check
    ("constant from r=10 to 140mm"). Returns `(r_bin_centers_m, median_beta)`, NaN-safe."""
    beta = local_wavenumber_field(fld_grid)
    r = np.hypot(fld_grid.x_m - feed_xy_m[0], fld_grid.y_m - feed_xy_m[1])
    r_max = np.nanmax(r)
    edges = np.linspace(0.0, r_max, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    medians = np.full(n_bins, np.nan)
    for i in range(n_bins):
        mask = (r >= edges[i]) & (r < edges[i + 1]) & ~np.isnan(beta)
        if mask.any():
            medians[i] = float(np.median(beta[mask]))
    return centers, medians


@dataclass(frozen=True)
class TM0Comparison:
    beta_measured: float
    beta_tm0: float
    n_eff_measured: "float | None"
    n_eff_tm0: "float | None"
    relative_difference: float


def compare_beta_to_tm0(beta_measured, frequency_hz, eps_r, thickness_m):
    """Sanity check only (brief: "not a requirement") -- never gates anything, just reports how
    far the imported field's own measured slope is from the bare-slab TM0 baseline."""
    beta_tm0 = tm0_beta(frequency_hz, eps_r, thickness_m)
    k0 = free_space_wavenumber(frequency_hz)
    return TM0Comparison(
        beta_measured=beta_measured,
        beta_tm0=beta_tm0,
        n_eff_measured=beta_measured / k0,
        n_eff_tm0=beta_tm0 / k0,
        relative_difference=abs(beta_measured - beta_tm0) / beta_tm0,
    )


def cells_per_fringe(beta_rad_per_m, period_m):
    """lambda_sw / period, where lambda_sw = 2*pi/beta is the surface-wave fringe period. The
    brief's own worked example: beta ~= 225 rad/m, period = 8mm -> ~3.5 cells/fringe (undersampled
    against the recommended >= 6)."""
    if beta_rad_per_m <= 0.0:
        raise ValueError(f"beta_rad_per_m must be > 0, got {beta_rad_per_m!r}")
    lambda_sw = 2.0 * np.pi / beta_rad_per_m
    return lambda_sw / period_m


def _cell_centers(vmin, vmax, period_m):
    """`n` centres of `period_m`-spaced cells, symmetric about the midpoint, spanning no more than
    `[vmin, vmax]` -- `n = floor((vmax-vmin)/period_m)` (e.g. a 300mm aperture at an 8mm period
    gives 37 cells, matching the PHASE13B review's own worked example, not 38). The floor has a
    relative tolerance (1e-6 cells): 300.8mm / 6.4mm is 46.999... in floating point, and a
    whole-cell side must not silently lose a row and column (Phase 21 C1)."""
    span = vmax - vmin
    n = max(1, int(np.floor(span / period_m + 1e-6)))
    center = 0.5 * (vmin + vmax)
    offsets = (np.arange(n) - (n - 1) / 2.0) * period_m
    return center + offsets, n


def _complex_interpolate(x_1d, y_1d, values, x_new, y_new):
    """Bilinear interpolation of a complex-valued field, done on the real and imaginary parts
    SEPARATELY (never on magnitude/phase -- interpolating phase directly is exactly the wrapping
    hazard PHASE13_WEBAPP_PROMPT.md's "complex (not phase) interpolation" note warns about, e.g.
    averaging -179deg and +179deg naively gives ~0deg instead of +-180deg)."""
    from scipy.interpolate import RegularGridInterpolator

    real_interp = RegularGridInterpolator((y_1d, x_1d), values.real, bounds_error=False, fill_value=None)
    imag_interp = RegularGridInterpolator((y_1d, x_1d), values.imag, bounds_error=False, fill_value=None)
    pts = np.stack([y_new.ravel(), x_new.ravel()], axis=-1)
    real_new = real_interp(pts).reshape(x_new.shape)
    imag_new = imag_interp(pts).reshape(x_new.shape)
    return real_new + 1j * imag_new


def resample_to_cell_grid(reference_field, period_m):
    """PHASE13B A2: resample an imported reference field (typically a fine field-calculator grid,
    e.g. 2mm pitch) onto the cell-CENTRE grid of a `period_m` unit cell, covering the same
    aperture extent -- so hologram synthesis runs once per physical unit cell, not once per raw
    field sample. Complex interpolation throughout (see `_complex_interpolate`); works for both a
    phase-only field (interpolates the unit-magnitude phasor, keeps `amplitude_available=False`)
    and a complex one (interpolates `Ex` itself, so magnitude is resampled too)."""
    if period_m <= 0.0:
        raise ValueError(f"period_m must be > 0, got {period_m!r}")
    x_1d = reference_field.x_m[0, :]
    y_1d = reference_field.y_m[:, 0]
    x_centers, _ = _cell_centers(float(x_1d.min()), float(x_1d.max()), period_m)
    y_centers, _ = _cell_centers(float(y_1d.min()), float(y_1d.max()), period_m)
    x_new, y_new = np.meshgrid(x_centers, y_centers, indexing="xy")

    source = f"{reference_field.source}+RESAMPLED({period_m * 1e3:.3g}mm)"
    if reference_field.amplitude_available:
        ex_new = _complex_interpolate(x_1d, y_1d, reference_field.Ex, x_new, y_new)
        return replace(
            reference_field,
            x_m=x_new,
            y_m=y_new,
            Ex=ex_new,
            Ey=np.zeros_like(ex_new),
            dx_m=period_m,
            dy_m=period_m,
            phase_deg=np.degrees(np.angle(ex_new)),
            source=source,
        )

    phasor = np.exp(1j * np.radians(reference_field.phase()))
    phasor_new = _complex_interpolate(x_1d, y_1d, phasor, x_new, y_new)
    return replace(
        reference_field,
        x_m=x_new,
        y_m=y_new,
        dx_m=period_m,
        dy_m=period_m,
        phase_deg=np.degrees(np.angle(phasor_new)),
        source=source,
    )


def feed_position_from_nan_mask(x_m, y_m, nan_mask):
    """PHASE15 Part A: "feed position (x0, y0) is taken from the reference field's NaN/singularity
    location (default)" -- the centroid of the masked (NaN) cells, since a monopole/point-feed
    reference field's own NaN sits exactly at the feed (`fld_to_reference_field_grid`'s module
    docstring: "the single feed-point NaN, at the origin, by construction"). Returns `(0.0, 0.0)`
    (this project's existing origin-centered feed convention, e.g. `phase8b.reference_wave.
    monopole_reference_field`) if `nan_mask` has no True entries at all -- a synthetic or already-
    masked field carries no singularity to detect, not an error."""
    nan_mask = np.asarray(nan_mask, dtype=bool)
    if not nan_mask.any():
        return 0.0, 0.0
    return float(np.mean(x_m[nan_mask])), float(np.mean(y_m[nan_mask]))


def aliasing_warning(cells_per_fringe_value, minimum=RECOMMENDED_MIN_CELLS_PER_FRINGE):
    if cells_per_fringe_value >= minimum:
        return None
    return (
        f"only {cells_per_fringe_value:.1f} cells/fringe at this period -- surface-wave holograms "
        f"normally need >= {minimum:.0f} (~lambda_sw/6); consider a smaller unit-cell period."
    )
