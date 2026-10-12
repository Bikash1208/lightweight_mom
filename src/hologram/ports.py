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

import numpy as np


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
    # Phase 23: `None` = the axis-aligned rectangle above; otherwise a convex polygon, counter-
    # clockwise `((x, y), ...)` in the same global coordinates, lying inside the box above.
    vertices_m: "tuple | None" = None

    def __post_init__(self):
        for name in ("x_min_m", "x_max_m", "y_min_m", "y_max_m", "feed_x_m", "feed_y_m"):
            v = getattr(self, name)
            if not math.isfinite(v):
                raise ValueError(f"PortRegion.{name} must be finite, got {v!r}")
        if self.x_max_m <= self.x_min_m:
            raise ValueError(f"port {self.port_id}: x_max_m ({self.x_max_m}) must be > x_min_m ({self.x_min_m})")
        if self.y_max_m <= self.y_min_m:
            raise ValueError(f"port {self.port_id}: y_max_m ({self.y_max_m}) must be > y_min_m ({self.y_min_m})")
        if self.vertices_m is not None:
            verts = tuple((float(x), float(y)) for x, y in self.vertices_m)
            _check_convex_ccw(self.port_id, verts)
            slack = 1e-9
            if (
                min(v[0] for v in verts) < self.x_min_m - slack
                or max(v[0] for v in verts) > self.x_max_m + slack
                or min(v[1] for v in verts) < self.y_min_m - slack
                or max(v[1] for v in verts) > self.y_max_m + slack
            ):
                raise ValueError(f"port {self.port_id}: polygon vertices lie outside the region's bounding box")
            object.__setattr__(self, "vertices_m", verts)

    @property
    def width_m(self):
        return self.x_max_m - self.x_min_m

    @property
    def height_m(self):
        return self.y_max_m - self.y_min_m

    @property
    def is_polygon(self):
        return self.vertices_m is not None

    @property
    def outline_m(self):
        """The region's boundary as counter-clockwise `((x, y), ...)`: the polygon, or the rectangle's corners."""
        if self.vertices_m is not None:
            return self.vertices_m
        return (
            (self.x_min_m, self.y_min_m),
            (self.x_max_m, self.y_min_m),
            (self.x_max_m, self.y_max_m),
            (self.x_min_m, self.y_max_m),
        )

    @property
    def area_m2(self):
        if self.vertices_m is not None:
            return polygon_area_m2(self.vertices_m)
        return self.width_m * self.height_m

    def contains(self, x_m, y_m):
        if self.vertices_m is not None:
            return _point_in_convex_polygon(self.vertices_m, x_m, y_m)
        return (self.x_min_m <= x_m <= self.x_max_m) and (self.y_min_m <= y_m <= self.y_max_m)

    def contains_array(self, x_m, y_m):
        """Boolean array, same shape as `x_m`: which points lie inside (edges inclusive, as `contains`)."""
        x_m, y_m = np.asarray(x_m, dtype=float), np.asarray(y_m, dtype=float)
        if self.vertices_m is None:
            return (x_m >= self.x_min_m) & (x_m <= self.x_max_m) & (y_m >= self.y_min_m) & (y_m <= self.y_max_m)
        inside = np.ones(x_m.shape, dtype=bool)
        v = self.vertices_m
        for k in range(len(v)):
            (x0, y0), (x1, y1) = v[k], v[(k + 1) % len(v)]
            cross = (x1 - x0) * (y_m - y0) - (y1 - y0) * (x_m - x0)
            inside &= cross >= -_POLYGON_EDGE_TOLERANCE_M * math.hypot(x1 - x0, y1 - y0)
        return inside


_POLYGON_EDGE_TOLERANCE_M = 1e-9  # a point within 1 nm of an edge counts as on it (inclusive)


def polygon_area_m2(vertices_m):
    """Shoelace area of a simple polygon (positive for counter-clockwise vertices)."""
    v = vertices_m
    twice = 0.0
    for k in range(len(v)):
        (x0, y0), (x1, y1) = v[k], v[(k + 1) % len(v)]
        twice += x0 * y1 - x1 * y0
    return 0.5 * abs(twice)


def _check_convex_ccw(port_id, verts):
    if len(verts) < 3:
        raise ValueError(f"port {port_id}: a polygon needs at least 3 vertices, got {len(verts)}")
    for x, y in verts:
        if not (math.isfinite(x) and math.isfinite(y)):
            raise ValueError(f"port {port_id}: polygon vertices must be finite")
    n = len(verts)
    for k in range(n):
        (x0, y0), (x1, y1), (x2, y2) = verts[k], verts[(k + 1) % n], verts[(k + 2) % n]
        cross = (x1 - x0) * (y2 - y1) - (y1 - y0) * (x2 - x1)
        if cross < -1e-18:
            raise ValueError(f"port {port_id}: polygon must be convex with counter-clockwise vertices")
    if polygon_area_m2(verts) <= 0.0:
        raise ValueError(f"port {port_id}: polygon has zero area")


def _point_in_convex_polygon(verts, x_m, y_m):
    n = len(verts)
    for k in range(n):
        (x0, y0), (x1, y1) = verts[k], verts[(k + 1) % n]
        cross = (x1 - x0) * (y_m - y0) - (y1 - y0) * (x_m - x0)
        if cross < -_POLYGON_EDGE_TOLERANCE_M * math.hypot(x1 - x0, y1 - y0):
            return False
    return True


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


def layout_preset_grid(n_rows, n_cols, side_m, gap_m=0.0):
    """Phase 21: `n_rows` x `n_cols` equal `side_m` x `side_m` sub-apertures, centred on the origin,
    `gap_m` between neighbours (outside each side), one feed at each sub-aperture centre.

    Port order is a serpentine from the top-right: the top row runs right -> left, the next row
    left -> right, and so on. For 2x2 that is P1 (+,+), P2 (-,+), P3 (-,-), P4 (+,-) -- the
    counter-clockwise order of the hand-built four-port reference board."""
    for name, n in (("n_rows", n_rows), ("n_cols", n_cols)):
        if not isinstance(n, int) or isinstance(n, bool) or n < 1:
            raise ValueError(f"{name} must be a positive integer, got {n!r}")
    if side_m <= 0.0:
        raise ValueError(f"side_m must be > 0, got {side_m!r}")
    if gap_m < 0.0:
        raise ValueError(f"gap_m must be >= 0, got {gap_m!r}")
    pitch = side_m + gap_m
    half_w = (n_cols * side_m + (n_cols - 1) * gap_m) / 2.0
    half_h = (n_rows * side_m + (n_rows - 1) * gap_m) / 2.0
    regions = []
    for r in range(n_rows):  # r = 0 is the top row
        cols = range(n_cols - 1, -1, -1) if r % 2 == 0 else range(n_cols)
        for c in cols:
            x0 = -half_w + c * pitch
            y1 = half_h - r * pitch
            regions.append(PortRegion(len(regions), x0, x0 + side_m, y1 - side_m, y1, x0 + side_m / 2.0, y1 - side_m / 2.0))
    return tuple(regions)


def regular_hexagon_vertices(cx_m, cy_m, side_m, orientation="pointy"):
    """Counter-clockwise vertices of a regular hexagon of side (= circumradius) `side_m` centred on
    `(cx_m, cy_m)`. `"pointy"`: corners at top and bottom (vertical edges left/right); `"flat"`: corners
    left and right (horizontal edges top/bottom)."""
    if orientation not in ("pointy", "flat"):
        raise ValueError(f"orientation must be 'pointy' or 'flat', got {orientation!r}")
    start = 90.0 if orientation == "pointy" else 0.0
    return tuple(
        (cx_m + side_m * math.cos(math.radians(start + 60.0 * k)), cy_m + side_m * math.sin(math.radians(start + 60.0 * k)))
        for k in range(6)
    )


def _floor_to_lattice(v, period_m):
    return math.floor(v / period_m + 1e-9) * period_m


def _ceil_to_lattice(v, period_m):
    return math.ceil(v / period_m - 1e-9) * period_m


def layout_preset_hex4(side_m=None, flat_to_flat_m=None, gap_m=0.0, orientation="pointy", period_m=None):
    """Phase 23: four regular hexagonal sub-apertures in the compact cluster (two offset rows of two),
    centred on the origin, one monopole at each hexagon centre.

    `side_m` is the hexagon side = circumradius (flat-to-flat = sqrt(3) * side, corner-to-corner =
    2 * side); `flat_to_flat_m` is an alternative way to give the size (exactly one of the two).
    `gap_m` is the distance between facing edges (0 = the hexagons tile). `orientation`: `"pointy"`
    (corners up/down, two horizontal rows) or `"flat"` (the same cluster turned by 90 degrees: two
    columns). `period_m`: when given, each region's cropping box is snapped OUTWARD onto the lattice of
    multiples of `period_m` from the origin, so every port's cells sit on one shared lattice (no
    off-grid seams between neighbouring hexagons); the polygon itself is unchanged.

    Port order (pointy): P1 upper row, right; P2 upper row, left; P3 lower row, left; P4 lower row,
    right -- the counter-clockwise order of `layout_preset_grid`, so the default beams
    phi = 45/135/225/315 deg point away from the cluster centre (the ports' own directions from the
    centre are 60/150/240/330 deg). `"flat"` rotates the cluster by 90 deg and re-cycles the order so
    P1 is again the port nearest 45 deg (directions 60/150/240/330 deg)."""
    if (side_m is None) == (flat_to_flat_m is None):
        raise ValueError("give exactly one of side_m and flat_to_flat_m")
    if side_m is None:
        side_m = flat_to_flat_m / math.sqrt(3.0)
    if not side_m > 0.0:
        raise ValueError(f"hexagon size must be > 0, got side_m={side_m!r}")
    if gap_m < 0.0:
        raise ValueError(f"gap_m must be >= 0, got {gap_m!r}")
    if orientation not in ("pointy", "flat"):
        raise ValueError(f"orientation must be 'pointy' or 'flat', got {orientation!r}")
    if period_m is not None and not period_m > 0.0:
        raise ValueError(f"period_m must be > 0, got {period_m!r}")
    pitch = math.sqrt(3.0) * side_m + gap_m  # centre-to-centre of neighbours (edge-normal direction)
    row_y = pitch * math.sqrt(3.0) / 4.0
    centres = [  # pointy-top, port order P1..P4
        (0.25 * pitch, +row_y),
        (-0.75 * pitch, +row_y),
        (-0.25 * pitch, -row_y),
        (0.75 * pitch, -row_y),
    ]
    if orientation == "flat":
        rotated = [(-y, x) for x, y in centres]  # +90 deg: directions 150/240/330/60 deg
        centres = rotated[3:] + rotated[:3]  # P1 = the one at 60 deg
    regions = []
    for k, (cx, cy) in enumerate(centres):
        verts = regular_hexagon_vertices(cx, cy, side_m, orientation)
        x0, x1 = min(v[0] for v in verts), max(v[0] for v in verts)
        y0, y1 = min(v[1] for v in verts), max(v[1] for v in verts)
        if period_m is not None:
            x0, x1 = _floor_to_lattice(x0, period_m), _ceil_to_lattice(x1, period_m)
            y0, y1 = _floor_to_lattice(y0, period_m), _ceil_to_lattice(y1, period_m)
        regions.append(PortRegion(k, x0, x1, y0, y1, cx, cy, vertices_m=verts))
    return tuple(regions)


def snap_side_to_cells(side_m, period_m):
    """`(snapped_side_m, n_cells)`: `side_m` rounded to the nearest whole number of `period_m`
    cells (at least 1). Off-grid sides make the combined board resample cells by nearest neighbour."""
    if side_m <= 0.0 or period_m <= 0.0:
        raise ValueError(f"side_m and period_m must be > 0, got {side_m!r}, {period_m!r}")
    n = max(1, round(side_m / period_m))
    return n * period_m, n


def _clip_convex(subject, clip):
    """Sutherland-Hodgman: the polygon `subject` clipped to the convex counter-clockwise polygon `clip`."""
    out = list(subject)
    n = len(clip)
    for k in range(n):
        if not out:
            break
        (ax, ay), (bx, by) = clip[k], clip[(k + 1) % n]
        inp, out = out, []

        def side(p, ax=ax, ay=ay, bx=bx, by=by):
            return (bx - ax) * (p[1] - ay) - (by - ay) * (p[0] - ax)

        for m in range(len(inp)):
            cur, nxt = inp[m], inp[(m + 1) % len(inp)]
            s_cur, s_nxt = side(cur), side(nxt)
            if s_cur >= 0.0:
                out.append(cur)
            if (s_cur >= 0.0) != (s_nxt >= 0.0):
                t = s_cur / (s_cur - s_nxt)
                out.append((cur[0] + t * (nxt[0] - cur[0]), cur[1] + t * (nxt[1] - cur[1])))
    return out


def overlap_area_m2(a: PortRegion, b: PortRegion):
    """Area of the intersection of two regions (rectangles or convex polygons)."""
    if not a.is_polygon and not b.is_polygon:
        ox = max(0.0, min(a.x_max_m, b.x_max_m) - max(a.x_min_m, b.x_min_m))
        oy = max(0.0, min(a.y_max_m, b.y_max_m) - max(a.y_min_m, b.y_min_m))
        return ox * oy
    clipped = _clip_convex(a.outline_m, b.outline_m)
    return polygon_area_m2(clipped) if len(clipped) >= 3 else 0.0


def check_no_overlap(port_regions, tolerance_m=1e-9):
    """Raises ValueError naming the first pair of ports whose regions overlap with positive area
    (touching edges are fine -- a gap of 0 is a valid board). Rectangles: both overlaps along x and
    y must exceed `tolerance_m`. Polygons: the intersection area must exceed a `tolerance_m` x 1 mm
    sliver."""
    for i, a in enumerate(port_regions):
        for b in port_regions[i + 1 :]:
            if a.is_polygon or b.is_polygon:
                if overlap_area_m2(a, b) > tolerance_m * 1e-3:
                    raise ValueError(f"ports {a.port_id + 1} and {b.port_id + 1} overlap")
                continue
            ox = min(a.x_max_m, b.x_max_m) - max(a.x_min_m, b.x_min_m)
            oy = min(a.y_max_m, b.y_max_m) - max(a.y_min_m, b.y_min_m)
            if ox > tolerance_m and oy > tolerance_m:
                raise ValueError(f"ports {a.port_id + 1} and {b.port_id + 1} overlap")


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


def _point_segment_distance_m(px, py, ax, ay, bx, by):
    dx, dy = bx - ax, by - ay
    length2 = dx * dx + dy * dy
    t = 0.0 if length2 == 0.0 else max(0.0, min(1.0, ((px - ax) * dx + (py - ay) * dy) / length2))
    return math.hypot(px - (ax + t * dx), py - (ay + t * dy))


def region_gap_m(a: PortRegion, b: PortRegion):
    """Euclidean edge-to-edge gap [m] between two regions: the axis-aligned rectangle formula when both
    are rectangles (unchanged), otherwise the minimum vertex-to-edge distance of the two convex
    polygons (0 when they touch or overlap)."""
    if not a.is_polygon and not b.is_polygon:
        return _rect_gap_m(a, b)
    if overlap_area_m2(a, b) > 0.0:
        return 0.0
    best = math.inf
    for src, dst in ((a.outline_m, b.outline_m), (b.outline_m, a.outline_m)):
        for px, py in src:
            for k in range(len(dst)):
                (ax, ay), (bx, by) = dst[k], dst[(k + 1) % len(dst)]
                best = min(best, _point_segment_distance_m(px, py, ax, ay, bx, by))
    return 0.0 if best < 1e-12 else best


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
            gap_m = region_gap_m(port_regions[i], port_regions[j])
            gap_lsw = gap_m / lambda_sw
            flagged = gap_lsw < min_lambda_sw - _FLAG_TOLERANCE
            suggestion = (min_lambda_sw - gap_lsw) * lambda_sw if flagged else None
            flags.append(
                IsolationFlag(port_regions[i].port_id, port_regions[j].port_id, gap_m, gap_lsw, flagged, suggestion)
            )
    return tuple(flags)
