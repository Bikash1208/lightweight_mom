"""Phase 15 Part C (gap-fill): multi-port sub-apertures. Methodology source: Mishra, Barad, Sarkar,
IEEE TVT 2026 (`docs/history/briefs/PHASE15_POLARIZATION_VORTEX_MULTIFEED.md` Part C) -- five coaxial
monopoles on one grounded board, each feeding its OWN rectangular sub-aperture (ports 1-4: 5*lambda
x 5*lambda single-beam; port 5: 10*lambda x 10*lambda quad-beam), rather than several feeds sharing
one hologram (`hologram.multi_feed`'s superposition/segmentation strategies, still available as the
brief's own "Part C'" alternative for users who want that instead).

A `PortRegion` is just an axis-aligned rectangle (in the SAME global aperture coordinates every
other Phase 15 module uses) plus the feed position inside it -- `multi_port_synthesis.py` runs the
existing single-aperture `hologram.hmsa_synthesis.synthesize_formulation_a/b` independently on each
port's own cropped sub-grid (Part A's own `feed_xy_m` already makes a single synthesis call radial
about an arbitrary feed, so a port needs no new synthesis machinery of its own -- only the cropping
and layout logic below).
"""

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class PortRegion:
    """One sub-aperture: an axis-aligned rectangle `[x_min_m, x_max_m] x [y_min_m, y_max_m]` in the
    board's global coordinates, and the feed position inside it (the brief's default: the
    sub-aperture's own centre; `custom_layout` allows overriding it, e.g. for an uploaded `.fld`
    whose own measured singularity is not exactly centred)."""

    port_id: int
    x_min_m: float
    x_max_m: float
    y_min_m: float
    y_max_m: float
    feed_x_m: float
    feed_y_m: float

    def __post_init__(self):
        for name in ("x_min_m", "x_max_m", "y_min_m", "y_max_m", "feed_x_m", "feed_y_m"):
            v = getattr(self, name)
            if not math.isfinite(v):
                raise ValueError(f"PortRegion.{name} must be finite, got {v!r}")
        if self.x_max_m <= self.x_min_m:
            raise ValueError(f"port {self.port_id}: x_max_m ({self.x_max_m}) must be > x_min_m ({self.x_min_m})")
        if self.y_max_m <= self.y_min_m:
            raise ValueError(f"port {self.port_id}: y_max_m ({self.y_max_m}) must be > y_min_m ({self.y_min_m})")

    @property
    def width_m(self):
        return self.x_max_m - self.x_min_m

    @property
    def height_m(self):
        return self.y_max_m - self.y_min_m

    @property
    def area_m2(self):
        return self.width_m * self.height_m

    def contains(self, x_m, y_m):
        return (self.x_min_m <= x_m <= self.x_max_m) and (self.y_min_m <= y_m <= self.y_max_m)


def layout_preset_2x2_plus_1_large(small_side_m, gap_m=0.0):
    """Mishra's own 5-port layout (Fig. 1 of the paper): four `small_side_m` x `small_side_m`
    active sub-apertures (ports 0-3, a 2x2 block) beside one `2*small_side_m` x `2*small_side_m`
    active sub-aperture (port 4). With `gap_m=0` the board tiles exactly. `gap_m` (Part C's guard
    gap) is added between neighbouring ports, outside every active side."""
    if small_side_m <= 0.0:
        raise ValueError(f"small_side_m must be > 0, got {small_side_m!r}")
    if gap_m < 0.0:
        raise ValueError(f"gap_m must be >= 0, got {gap_m!r}")
    s, g = small_side_m, gap_m
    raw = [
        (-2 * s - g, -s - g, g, s + g),  # port 0: top-left small
        (-s, 0.0, g, s + g),  # port 1: top-right small
        (-2 * s - g, -s - g, -s, 0.0),  # port 2: bottom-left small
        (-s, 0.0, -s, 0.0),  # port 3: bottom-right small
        (g, 2 * s + g, -s, s),  # port 4: the large quad-beam port
    ]
    regions = []
    for i, (x0, x1, y0, y1) in enumerate(raw):
        if i < 4:
            y0, y1 = y0 - g / 2.0, y1 - g / 2.0  # the rows' own gap is centred on the origin
        # Port 4 (the large port) already spans [-s, s], centred on the origin -- the same centre
        # the 2x2 block above ends up at once its own row gap is centred. Shifting it by -g/2 too
        # (PHASE16 follow-ups bug) pulled its feed off that centre by g/2.
        regions.append(PortRegion(i, x0, x1, y0, y1, (x0 + x1) / 2.0, (y0 + y1) / 2.0))
    return tuple(regions)


def layout_preset_1xn(n, cell_side_m, gap_m=0.0):
    """`n` equal `cell_side_m` x `cell_side_m` active sub-apertures in a single row, centred on the
    origin, with `gap_m` between neighbours (outside each side) -- the brief's "1xN strip" preset."""
    if not isinstance(n, int) or n < 1:
        raise ValueError(f"n must be a positive integer, got {n!r}")
    if cell_side_m <= 0.0:
        raise ValueError(f"cell_side_m must be > 0, got {cell_side_m!r}")
    if gap_m < 0.0:
        raise ValueError(f"gap_m must be >= 0, got {gap_m!r}")
    s = cell_side_m
    half_total = (n * s + (n - 1) * gap_m) / 2.0
    regions = []
    for i in range(n):
        x0 = -half_total + i * (s + gap_m)
        regions.append(PortRegion(i, x0, x0 + s, -s / 2.0, s / 2.0, x0 + s / 2.0, 0.0))
    return tuple(regions)


def custom_layout(rectangles_m):
    """`rectangles_m`: an iterable of `(x_min_m, x_max_m, y_min_m, y_max_m)` or
    `(x_min_m, x_max_m, y_min_m, y_max_m, feed_x_m, feed_y_m)` tuples (the brief's "draw
    sub-apertures" option) -- the feed defaults to each rectangle's own centre when omitted."""
    regions = []
    for i, r in enumerate(rectangles_m):
        if len(r) == 4:
            x0, x1, y0, y1 = r
            fx, fy = (x0 + x1) / 2.0, (y0 + y1) / 2.0
        elif len(r) == 6:
            x0, x1, y0, y1, fx, fy = r
        else:
            raise ValueError(f"custom_layout: each rectangle needs 4 or 6 values, got {len(r)}")
        regions.append(PortRegion(i, x0, x1, y0, y1, fx, fy))
    return tuple(regions)


def _rect_gap_m(a: PortRegion, b: PortRegion):
    """Euclidean gap [m] between two axis-aligned rectangles: 0 if they touch or overlap along
    both axes (the standard separating-axis distance for axis-aligned boxes)."""
    dx = max(a.x_min_m - b.x_max_m, b.x_min_m - a.x_max_m, 0.0)
    dy = max(a.y_min_m - b.y_max_m, b.y_min_m - a.y_max_m, 0.0)
    return math.hypot(dx, dy)


@dataclass(frozen=True)
class IsolationFlag:
    """A first-order, NOT full-wave, isolation estimate (the brief: "flag adjacent-edge distance
    in lambda_sw and recommend a gap/guard row; say plainly that true isolation needs HFSS") --
    `gap_lambda_sw` is the measured edge-to-edge gap divided by the surface-wave wavelength this
    design's own feed uses; `flagged` is True when that gap is below `min_lambda_sw` by more than
    `_FLAG_TOLERANCE` (a tight seam, where surface-wave coupling between the two ports' currents is
    least suppressed by distance alone -- not a claim of a specific achievable isolation in dB,
    which this project cannot compute without a full-wave solve). This is INFORMATION, not a fault:
    isolation comes mainly from the holograms themselves radiating the surface wave away rather
    than from edge distance alone, and should be confirmed with an HFSS S-parameter solve either
    way."""

    port_a: int
    port_b: int
    gap_m: float
    gap_lambda_sw: float
    flagged: bool
    guard_row_suggestion_m: "float | None"  # extra gap (in addition to what's already there) to
    # reach `min_lambda_sw`, only when `flagged`


_FLAG_TOLERANCE = 0.01  # lambda_sw, PHASE16 follow-ups: the app's own default gap is lambda_sw
# rounded to 0.1 mm for display, so a design left at that default could come back ~0.1-0.2%
# under 1.00 lambda_sw from rounding alone and get flagged as "too tight" at its own recommended
# minimum; this tolerance comfortably covers that rounding without hiding a genuinely tight gap.


def isolation_estimate(port_regions, k_swg, min_lambda_sw=1.0):
    """Pairwise `IsolationFlag` for every pair of ports -- `k_swg` [rad/m] the guided surface-wave
    propagation constant feeding this board (the same quantity `hologram.multi_feed` already needs
    for an analytic feed), `lambda_sw = 2*pi/k_swg`."""
    if k_swg is None or k_swg <= 0.0:
        raise ValueError(f"isolation_estimate needs a positive k_swg, got {k_swg!r}")
    lambda_sw = 2.0 * math.pi / k_swg
    flags = []
    for i in range(len(port_regions)):
        for j in range(i + 1, len(port_regions)):
            gap_m = _rect_gap_m(port_regions[i], port_regions[j])
            gap_lsw = gap_m / lambda_sw
            flagged = gap_lsw < min_lambda_sw - _FLAG_TOLERANCE
            suggestion = (min_lambda_sw - gap_lsw) * lambda_sw if flagged else None
            flags.append(
                IsolationFlag(port_regions[i].port_id, port_regions[j].port_id, gap_m, gap_lsw, flagged, suggestion)
            )
    return tuple(flags)
