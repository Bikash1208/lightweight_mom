"""Phase 15 Part C: multi-feed ("multipole") reference-field construction. Methodology source:
`docs/history/briefs/PHASE15_POLARIZATION_VORTEX_MULTIFEED.md` Part C (paper eq. 5 for superposition).

Each feed contributes a reference field either supplied directly (an uploaded `.fld`, already a
`phase8e.reference_field_input.ReferenceFieldGrid`) or built analytically about that feed's own
position: `phase8b.reference_wave.monopole_reference_field` already computes the standard
`rho^(-1/2) exp(-j*k_swg*rho)` cylindrical wave about the ORIGIN, unmodified -- `_monopole_field_at`
below reuses it exactly as-is, just shifting the coordinates it is evaluated at, rather than
duplicating or generalizing that function.

Two strategies, both producing a single `ReferenceFieldGrid` on the caller's `(x_m, y_m)` grid so
the existing single-reference-field `hologram.hmsa_synthesis.synthesize_formulation_a/b` need no
changes to consume either one:

- **Superposition** (paper eq. 5): `E_ref(x,y) = sum_k amplitude_k * E_ref,k(x,y)`, each feed's field
  resampled to the common grid first (complex interpolation, never on magnitude/phase separately --
  same reasoning as `fields.reference_field_grid._complex_interpolate`, whose exact technique this
  module's `_interpolate_complex_to_grid` duplicates rather than importing a sibling module's
  underscore-prefixed helper, matching how `pattern.vector_aperture_pattern` already duplicates
  `aperture_pattern._directivity_estimate` for the same reason). Known literature caveat (the brief:
  "report the aperture-efficiency penalty noted in the literature for superposition"): summing N
  independent cylindrical waves spreads the total feed power over a wider effective illumination
  than any one feed alone, typically costing aperture efficiency relative to segmentation -- this
  project has no independent measurement of that penalty's SIZE (no citation-backed number is
  invented here); Part B's vector pattern model is what actually measures the real difference for a
  given design, per the brief's own "let the vector model show the result".
- **Segmentation**: each feed owns the aperture region nearest to it (a Voronoi partition by
  straight-line distance, `feed_owner_index`); the combined field is PIECEWISE -- each cell uses only
  its owning feed's own field, no summation, so cells near a "seam" between two feeds' regions do not
  see the interference the superposition strategy has there.

Either way, Part A's polarization compensation needs to know which feed a given BEAM's current is
radial about: `feed_xy_m_for_beams` builds the `(x0_sequence, y0_sequence)` format
`hologram.polarization_compensation.combine_beams_with_polarization`'s `feed_xy_m` parameter already
accepts (Part A shipped that generality specifically so Part C would need no changes there).
"""

from dataclasses import dataclass

import numpy as np

from ..phase8b.reference_wave import monopole_reference_field
from ..phase8e.reference_field_input import ReferenceFieldGrid


@dataclass(frozen=True)
class FeedSpec:
    """One feed: its aperture position, a relative complex excitation `amplitude` (magnitude and
    phase, applied on top of whichever field this feed contributes), and EITHER an uploaded
    `reference_field` (already a `ReferenceFieldGrid`, any source grid -- resampled onto the target
    grid by this module) or `None` to build an analytic TM0 cylindrical wave about `(x_m, y_m)`
    (needs `k_swg` at combination time -- `superposed_reference_field`/
    `segmented_reference_field_and_regions` raise if a feed needs it and none was supplied)."""

    x_m: float
    y_m: float
    amplitude: complex = 1.0
    reference_field: "ReferenceFieldGrid | None" = None

    def __post_init__(self):
        if not np.isfinite(self.x_m) or not np.isfinite(self.y_m):
            raise ValueError(f"FeedSpec position must be finite, got ({self.x_m!r}, {self.y_m!r})")


def _monopole_field_at(x_m, y_m, feed: FeedSpec, k_swg, amplitude_at_1m=1.0):
    """`monopole_reference_field` (unmodified) evaluated about `feed`'s own position, by shifting
    coordinates -- that function's own `rho = sqrt(x**2+y**2)` already IS "distance from wherever
    (0,0) of the arrays passed in is", so `x_m - feed.x_m, y_m - feed.y_m` is the whole
    generalization needed."""
    return monopole_reference_field(x_m - feed.x_m, y_m - feed.y_m, k_swg, amplitude_at_1m)


def _interpolate_complex_to_grid(source_x_m, source_y_m, values, x_m, y_m):
    """Bilinear interpolation of a complex field onto an arbitrary target grid, real and imaginary
    parts separately (never magnitude/phase, which would mishandle the +-180deg wrap) -- the same
    technique `fields.reference_field_grid._complex_interpolate` uses, duplicated here rather than
    reaching into that module's private helper (module docstring)."""
    from scipy.interpolate import RegularGridInterpolator

    x_1d = np.asarray(source_x_m)[0, :]
    y_1d = np.asarray(source_y_m)[:, 0]
    real_interp = RegularGridInterpolator((y_1d, x_1d), values.real, bounds_error=False, fill_value=None)
    imag_interp = RegularGridInterpolator((y_1d, x_1d), values.imag, bounds_error=False, fill_value=None)
    pts = np.stack([np.asarray(y_m).ravel(), np.asarray(x_m).ravel()], axis=-1)
    real_new = real_interp(pts).reshape(np.shape(x_m))
    imag_new = imag_interp(pts).reshape(np.shape(x_m))
    return real_new + 1j * imag_new


def _feed_complex_field(feed: FeedSpec, x_m, y_m, k_swg):
    """`feed`'s own contribution, resampled/evaluated onto `(x_m, y_m)`, already scaled by its
    `amplitude` -- shared by both strategies."""
    if feed.reference_field is not None:
        rf = feed.reference_field
        source_values = rf.Ex if rf.amplitude_available else np.exp(1j * np.radians(rf.phase()))
        values = _interpolate_complex_to_grid(rf.x_m, rf.y_m, source_values, x_m, y_m)
    else:
        if k_swg is None:
            raise ValueError(
                f"feed at ({feed.x_m}, {feed.y_m}) has no reference_field and no k_swg was supplied "
                "for the analytic TM0 fallback"
            )
        values = _monopole_field_at(x_m, y_m, feed, k_swg)
    return feed.amplitude * values


def superposed_reference_field(feeds, x_m, y_m, k_swg=None):
    """Paper eq. 5: `E_ref = sum_k amplitude_k * E_ref,k`, every feed resampled/evaluated onto the
    common `(x_m, y_m)` grid first. Returns a `ReferenceFieldGrid` (amplitude available, always
    complex -- a sum of complex fields has no meaningful "phase-only" form)."""
    if not feeds:
        raise ValueError("superposed_reference_field: at least one feed is required")
    x_m, y_m = np.asarray(x_m, dtype=float), np.asarray(y_m, dtype=float)
    total = np.zeros(np.broadcast(x_m, y_m).shape, dtype=complex)
    for feed in feeds:
        total = total + _feed_complex_field(feed, x_m, y_m, k_swg)
    return ReferenceFieldGrid(
        x_m=x_m,
        y_m=y_m,
        Ex=total,
        Ey=np.zeros_like(total),
        amplitude_available=True,
        source=f"MULTI_FEED_SUPERPOSITION({len(feeds)} feeds)",
        dx_m=None,
        dy_m=None,
        z_reference_m=None,
        phase_deg=np.degrees(np.angle(total)),
    )


def feed_owner_index(feeds, x_m, y_m):
    """Voronoi partition of `(x_m, y_m)` by nearest feed (straight-line distance) -- which feed
    each cell belongs to under the segmentation strategy."""
    if not feeds:
        raise ValueError("feed_owner_index: at least one feed is required")
    x_m, y_m = np.asarray(x_m, dtype=float), np.asarray(y_m, dtype=float)
    distances = np.stack([np.hypot(x_m - f.x_m, y_m - f.y_m) for f in feeds], axis=-1)
    return np.argmin(distances, axis=-1)


def segmented_reference_field_and_regions(feeds, x_m, y_m, k_swg=None):
    """Each feed owns its own nearest-distance aperture region (`feed_owner_index`); the combined
    field is PIECEWISE (no summation across the region boundary), each cell holding only its owning
    feed's own field there. Returns `(ReferenceFieldGrid, owner_index_array)`."""
    if not feeds:
        raise ValueError("segmented_reference_field_and_regions: at least one feed is required")
    x_m, y_m = np.asarray(x_m, dtype=float), np.asarray(y_m, dtype=float)
    owner = feed_owner_index(feeds, x_m, y_m)
    combined = np.zeros(np.broadcast(x_m, y_m).shape, dtype=complex)
    for i, feed in enumerate(feeds):
        mask = owner == i
        if not mask.any():
            continue
        combined = np.where(mask, _feed_complex_field(feed, x_m, y_m, k_swg), combined)
    grid = ReferenceFieldGrid(
        x_m=x_m,
        y_m=y_m,
        Ex=combined,
        Ey=np.zeros_like(combined),
        amplitude_available=True,
        source=f"MULTI_FEED_SEGMENTATION({len(feeds)} feeds)",
        dx_m=None,
        dy_m=None,
        z_reference_m=None,
        phase_deg=np.degrees(np.angle(combined)),
    )
    return grid, owner


def feed_xy_m_for_beams(feeds, beam_feed_indices):
    """Builds the `(x0_sequence, y0_sequence)` format `hologram.polarization_compensation.
    combine_beams_with_polarization`'s `feed_xy_m` parameter accepts, from a per-beam list of
    indices into `feeds` (the brief's "each beam can be assigned to a feed", superposition
    strategy). `beam_feed_indices[i]` is the index into `feeds` that beam `i`'s current is radial
    about."""
    xs = [feeds[i].x_m for i in beam_feed_indices]
    ys = [feeds[i].y_m for i in beam_feed_indices]
    return xs, ys
