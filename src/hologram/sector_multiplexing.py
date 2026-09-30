"""Phase 15 Part A: sector (quadrant) multiplexing, the brief's eq. 17 option -- each object beam
owns its own sector of the aperture (rather than the default full-aperture multiplexing, where every
beam's hologram term is superposed everywhere, eq. 16). No project document picks a sector SHAPE for
an arbitrary beam set, so this module fixes the natural, simplest one and documents it: a nearest-
beam-azimuth partition of the aperture (about the feed), i.e. each cell belongs to whichever beam's
own (theta_b, phi_b) direction its own azimuth phi' is closest to -- this reduces to equal quadrants
for 4 beams spaced 90 degrees apart (the brief's own example), and to equal wedges for any other
even spacing, while still being well-defined for an arbitrary, uneven set of beam azimuths.

`hologram.polarization_compensation.combine_beams_with_polarization`'s `sector_owner_index`
parameter is exactly this module's `sector_index_map` output -- this module only computes the
partition; the actual per-sector hologram restriction happens there.
"""

import numpy as np


def sector_index_map(x_m, y_m, feed_xy_m, beams):
    """For each `(x_m, y_m)` cell, the index into `beams` of the beam whose own `phi_deg` is
    angularly closest (about `feed_xy_m`) -- a Voronoi partition of azimuth by beam direction.
    `beams`: any sequence exposing `.phi_deg` (a `phase8b.object_wave.Beam`, typically)."""
    if not beams:
        raise ValueError("sector_index_map: at least one beam is required")
    x_m, y_m = np.asarray(x_m, dtype=float), np.asarray(y_m, dtype=float)
    x0, y0 = feed_xy_m
    phi_prime_deg = np.degrees(np.arctan2(y_m - y0, x_m - x0))

    beam_phis_deg = np.array([b.phi_deg for b in beams])
    # angular distance on a circle: wrap the difference into (-180, 180] before comparing magnitude
    diffs = phi_prime_deg[..., None] - beam_phis_deg[None, ...]
    wrapped = (diffs + 180.0) % 360.0 - 180.0
    return np.argmin(np.abs(wrapped), axis=-1)


def sector_boundaries_deg(feed_xy_m, beams):
    """The angular midpoints (deg, sorted) between consecutive beam azimuths -- the sector edges
    `sector_index_map` actually applies, exposed separately for reporting/plotting (e.g. drawing the
    sector lines on the layout preview)."""
    if not beams:
        raise ValueError("sector_boundaries_deg: at least one beam is required")
    phis = sorted(b.phi_deg % 360.0 for b in beams)
    n = len(phis)
    boundaries = []
    for i in range(n):
        a, b = phis[i], phis[(i + 1) % n]
        span = (b - a) % 360.0
        boundaries.append((a + span / 2.0) % 360.0)
    return sorted(boundaries)
