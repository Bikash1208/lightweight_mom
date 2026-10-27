"""Phase 28 item 3: candidate diode-switched CP rotating-axis cells (own code; additive, nothing existing is modified).

`ShortenableStarCell` -- a circular patch with a star of four through-slots (axes 0/45/90/135 deg, 8 arms). Every arm carries a conducting WALL across the slot
at radius `wall_radius` (metal over the whole slot width, a lumped diode in its middle, `diode_span` long): diode ON = the wall conducts, the arm is closed
beyond the wall (the arm is SHORTENED to `wall_radius`); diode OFF = the wall is open (capacitive), the arm keeps its full length. The geometry is the same in
every state; only the diode impedances change. State k (k = 0..3) = the arms of the slot PERPENDICULAR to axis k are shortened, so the long arms of axis k make
the cell's principal axis (a rotating-axis tensor sheet); `single_sided=True` puts a wall on one arm of each slot only (4 diodes instead of 8).

The wall is the "diode bridging the slot" of the Phase 20 R1 loads (`LumpedLoad`, rotated loads of R3): the load region lies inside the wall metal so that
rooftop functions cross it (`loads.check_load_geometry`). Diagonal walls use `transverse_conducts=True` (the R3 caveat for 45-degree strips).
Labels: MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET, VIA_OMITTED_APPROXIMATION.
"""

import math
from dataclasses import dataclass

import numpy as np

from ..geometry import UnitCell
from ..materials import Metal
from ..periodic_bc import RectangularLattice
from ..ris.cells import OFF, ON, RT5880
from ..ris.diodes import diode_load

ANGLES_DEG = (0.0, 45.0, 90.0, 135.0)


def _rot(dx, dy, ang):
    c, s = math.cos(ang), math.sin(ang)
    return dx * c + dy * s, -dx * s + dy * c  # (along, across)


@dataclass(frozen=True)
class StarPatchGeometry:
    """Duck-typed `geometry.Geometry`: circle minus the star slots plus the walls (mask based)."""

    diameter: float
    slot_length: float  # through-length of every slot [m] (arm length = half)
    slot_width: float
    walls: tuple  # ((arm_angle_rad, radius, wall_width), ...)

    kind = "patch_circle_star_walls"
    polarity = "metal"
    rotation = 0.0

    @property
    def feature_area(self):
        return math.pi * self.diameter**2 / 4.0

    def corners(self):
        r = self.diameter / 2.0
        return np.array([[r, r], [-r, r], [-r, -r], [r, -r]])

    def in_feature(self, x, y):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        inside = np.hypot(x, y) <= self.diameter / 2.0
        slot = np.zeros(x.shape, dtype=bool)
        for a in np.radians(ANGLES_DEG):
            u, v = _rot(x, y, a)
            slot |= (np.abs(u) <= self.slot_length / 2) & (np.abs(v) <= self.slot_width / 2)
        metal = inside & ~slot
        for ang, r, w in self.walls:
            u, v = _rot(x, y, ang)  # u along the arm
            metal |= (np.abs(u - r) <= w / 2) & (np.abs(v) <= self.slot_width / 2)
        return metal


@dataclass(frozen=True)
class ShortenableStarCell:
    period: float = 6.4e-3
    diameter: float = 6.0e-3
    slot_length: float = 5.2e-3
    slot_width: float = 1.6e-3
    wall_radius: float = 1.4e-3
    wall_width: float = 0.8e-3
    diode_span: float = 0.8e-3
    single_sided: bool = False

    kind = "star_shortenable"

    def arms_deg(self):
        arms = []
        for a in ANGLES_DEG:
            arms.append(a)
            if not self.single_sided:
                arms.append(a + 180.0)
        return arms

    @property
    def n_diodes(self):
        return len(self.arms_deg())

    def geometry(self):
        walls = tuple((math.radians(a), self.wall_radius, self.wall_width) for a in self.arms_deg())
        return StarPatchGeometry(self.diameter, self.slot_length, self.slot_width, walls)

    def states(self):
        """Four states: state k shortens the arms of the slot perpendicular to axis k (so axis k stays long)."""
        out = []
        for k in range(4):
            perp = (ANGLES_DEG[k] + 90.0) % 180.0
            out.append(tuple(ON if (a % 180.0) == perp else OFF for a in self.arms_deg()))
        return out

    def diode_sites(self):
        sites = []
        for a in self.arms_deg():
            ar = math.radians(a)
            x, y = self.wall_radius * math.cos(ar), self.wall_radius * math.sin(ar)
            sites.append((x, y, self.diode_span, self.wall_width, (a + 90.0) % 180.0))
        return sites

    def unit_cell(self, state, kit):
        if len(state) != self.n_diodes:
            raise ValueError(f"state must have {self.n_diodes} entries")
        loads = []
        for s, (x, y, length, width, direction) in zip(state, self.diode_sites(), strict=True):
            diagonal = abs(math.sin(math.radians(2 * direction))) > 1e-9
            loads.append(
                diode_load(kit.model(s), x, y, length, width, direction=direction, transverse_conducts=diagonal)
            )
        return UnitCell(
            RectangularLattice(self.period, self.period),
            self.geometry(),
            RT5880,
            Metal(),
            backing="pec",
            loads=tuple(loads),
        )


# ---------------------------------------------------------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class SectorHubGeometry:
    """Duck-typed `geometry.Geometry`: a central hub disc and eight 45-degree annular sectors (outer radius `radius`, separated from each other by a straight gap
    `sector_gap` and from the hub by the radial gap `radial_gap`), each sector joined to the hub by a neck of width `neck_width` that carries the diode."""

    radius: float
    hub_radius: float
    radial_gap: float
    sector_gap: float
    neck_width: float
    radius_diag: "float | None" = None  # outer radius of the four diagonal sectors (45, 135, ... deg); None = `radius`
    half_cell: "float | None" = None  # metal is clipped to |x|, |y| <= half_cell (keeps neighbouring cells apart)

    kind = "sector_hub"
    polarity = "metal"
    rotation = 0.0

    @property
    def feature_area(self):
        return math.pi * self.radius**2 * 0.85  # nominal; only used by generic reports

    def corners(self):
        r = self.radius
        return np.array([[r, r], [-r, r], [-r, -r], [r, -r]])

    def in_feature(self, x, y):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        r = np.hypot(x, y)
        th = np.arctan2(y, x)
        k = np.round(th / (math.pi / 4))
        off = th - k * math.pi / 4  # in [-pi/8, pi/8]
        r_in = self.hub_radius + self.radial_gap
        r_out = np.where((k % 2) == 0, self.radius, self.radius if self.radius_diag is None else self.radius_diag)
        sector = (r >= r_in) & (r <= r_out) & (r * np.sin(math.pi / 8 - np.abs(off)) >= self.sector_gap / 2)
        hub = r <= self.hub_radius
        u = r * np.cos(off)
        v = r * np.sin(off)
        neck = (u >= self.hub_radius - 1e-9) & (u <= r_in + 1e-9) & (np.abs(v) <= self.neck_width / 2)
        out = sector | hub | neck
        if self.half_cell is not None:
            out &= (np.abs(x) <= self.half_cell) & (np.abs(y) <= self.half_cell)
        return out


@dataclass(frozen=True)
class SectorHubCell:
    """Hub-and-sector cell: eight diode necks (one per sector); state k (axis k * 45 deg, k = 0..3) = the two necks on that axis ON (hub + two opposite sectors
    form a bow-tie along the axis), the other six OFF. `single_sided=True` switches only the neck at the axis angle ON (four necks are used, the opposite
    sectors are not connected)."""

    period: float = 6.4e-3
    radius: float = 3.0e-3
    hub_radius: float = 0.5e-3
    radial_gap: float = 0.6e-3
    sector_gap: float = 0.3e-3
    neck_width: float = 0.6e-3
    single_sided: bool = False
    radius_diag: "float | None" = None
    mode: str = "bowtie"  # "bowtie": the two necks on the axis ON; "notch": every neck ON except the two on the perpendicular axis (higher mean loading, more forward diodes)
    edge_margin: float = 0.2e-3  # metal stays this far from the cell boundary (neighbour gap = 2 x margin)

    kind = "sector_hub"
    n_diodes = 8

    def geometry(self):
        return SectorHubGeometry(
            self.radius,
            self.hub_radius,
            self.radial_gap,
            self.sector_gap,
            self.neck_width,
            self.radius_diag,
            self.period / 2 - self.edge_margin,
        )

    def arms_deg(self):
        return [45.0 * k for k in range(8)]

    def states(self):
        out = []
        for k in range(4):
            if self.mode == "notch":
                cut = {(45.0 * k + 90.0) % 180.0}
                out.append(tuple(OFF if (a % 180.0) in cut else ON for a in self.arms_deg()))
                continue
            if self.mode != "bowtie":
                raise ValueError(f"unknown mode {self.mode!r}")
            on = {45.0 * k} if self.single_sided else {45.0 * k, 45.0 * k + 180.0}
            out.append(tuple(ON if a in on else OFF for a in self.arms_deg()))
        return out

    def diode_sites(self):
        r_mid = self.hub_radius + self.radial_gap / 2
        sites = []
        for a in self.arms_deg():
            ar = math.radians(a)
            sites.append((r_mid * math.cos(ar), r_mid * math.sin(ar), self.radial_gap, self.neck_width, a % 180.0))
        return sites

    def unit_cell(self, state, kit):
        if len(state) != 8:
            raise ValueError("state must have 8 entries")
        loads = []
        for s, (x, y, length, width, direction) in zip(state, self.diode_sites(), strict=True):
            diagonal = abs(math.sin(math.radians(2 * direction))) > 1e-9
            loads.append(
                diode_load(kit.model(s), x, y, length, width, direction=direction, transverse_conducts=diagonal)
            )
        return UnitCell(
            RectangularLattice(self.period, self.period),
            self.geometry(),
            RT5880,
            Metal(),
            backing="pec",
            loads=tuple(loads),
        )


# ---------------------------------------------------------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class ManhattanHubGeometry:
    """Duck-typed `geometry.Geometry`, EVERY edge axis-aligned (so the rooftop mesh represents it without staircasing when the dimensions are multiples of the mesh
    step): a square hub (half-size `hub_half`), four arm bars along +-x and +-y (half-width `arm_half_width`, from the hub gap to `outer`), and four corner squares
    (from `arm_half_width + separation` to `outer` in the quadrant, spaced `neck_gap` from the hub). The arms are joined to the hub by necks along their axis; a corner
    square is joined to the hub by a neck whose direction rotates by 90 degrees from corner to corner (C4 symmetry), through the overlap of the corner square's inner
    edge with the hub's edge. Dimensions in metres."""

    hub_half: float
    neck_gap: float
    neck_width: float
    arm_half_width: float
    separation: float
    outer: float

    kind = "manhattan_hub"
    polarity = "metal"
    rotation = 0.0

    @property
    def feature_area(self):
        return 4 * self.outer**2

    def corners(self):
        r = self.outer
        return np.array([[r, r], [-r, r], [-r, -r], [r, -r]])

    # -- the pieces of the +x / ++ family; the other three are rotations by 90 degrees about the origin
    def _base_rects(self):
        h, g = self.hub_half, self.neck_gap
        wa, s, X = self.arm_half_width, self.separation, self.outer
        y0 = wa + s
        yn = 0.5 * (y0 + h)  # centre of the corner neck (inside the overlap [y0, h])
        wn = min(self.neck_width, h - y0)
        return {
            "hub": (-h, h, -h, h),
            "arm": (h + g, X, -wa, wa),
            "arm_neck": (h, h + g, -self.neck_width / 2, self.neck_width / 2),
            "corner": (h + g, X, y0, X),
            "corner_neck": (h, h + g, yn - wn / 2, yn + wn / 2),
        }

    def neck_rect(self, kind, k):
        """Neck `kind` ('arm_neck' or 'corner_neck') rotated by `k` x 90 degrees: `(x0, x1, y0, y1)`."""
        x0, x1, y0, y1 = self._base_rects()[kind]
        for _ in range(k % 4):
            x0, x1, y0, y1 = -y1, -y0, x0, x1  # (x, y) -> (-y, x)
        return x0, x1, y0, y1

    def in_feature(self, x, y):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        out = np.zeros(x.shape, dtype=bool)
        for name, (x0, x1, y0, y1) in self._base_rects().items():
            rots = [0] if name == "hub" else [0, 1, 2, 3]
            for k in rots:
                a0, a1, b0, b1 = (x0, x1, y0, y1)
                for _ in range(k):
                    a0, a1, b0, b1 = -b1, -b0, a0, a1
                out |= (x >= a0 - 1e-12) & (x <= a1 + 1e-12) & (y >= b0 - 1e-12) & (y <= b1 + 1e-12)
        return out


@dataclass(frozen=True)
class ManhattanHubCell:
    """Manhattan hub cell, eight diode necks. States: 0 = the two arm necks on +-x ON (axis 0 deg), 1 = the corner necks of the ++ and -- corners ON (axis 45 deg), 2 = arm
    necks on +-y (axis 90 deg), 3 = corner necks of the -+ and +- corners ON (axis 135 deg). Rotating a state by 90 degrees gives the state two steps on EXACTLY; the 45 degree
    steps are equal only to the extent the corner pieces are tuned to the arms (the balance the search scores)."""

    period: float = 6.4e-3
    hub_half: float = 1.2e-3
    neck_gap: float = 0.4e-3
    neck_width: float = 0.4e-3
    arm_half_width: float = 0.4e-3
    separation: float = 0.4e-3
    outer: float = 2.8e-3
    mode: str = "bowtie"  # "bowtie": the two necks of the state's axis ON; "notch": every neck ON except the two of the perpendicular state (more forward diodes, higher mean loading)

    kind = "manhattan_hub"
    n_diodes = 8

    def geometry(self):
        return ManhattanHubGeometry(
            self.hub_half, self.neck_gap, self.neck_width, self.arm_half_width, self.separation, self.outer
        )

    def _sites(self):
        """Eight diode sites `(kind, k, x, y, length, width, direction)`; arms k = 0..3 (0, 90, 180, 270 deg), corners k = 0..3 (++, -+, --, +-)."""
        geo = self.geometry()
        sites = []
        for kind in ("arm_neck", "corner_neck"):
            for k in range(4):
                x0, x1, y0, y1 = geo.neck_rect(kind, k)
                horizontal = (
                    k % 2 == 0
                )  # the base necks are horizontal; each 90-degree rotation swaps the direction (no float comparison of equal extents)
                length = (x1 - x0) if horizontal else (y1 - y0)
                width = (y1 - y0) if horizontal else (x1 - x0)
                sites.append((kind, k, 0.5 * (x0 + x1), 0.5 * (y0 + y1), length, width, "x" if horizontal else "y"))
        return sites

    def states(self):
        out = []
        for st in range(4):
            if st == 0:
                on = {("arm_neck", 0), ("arm_neck", 2)}
            elif st == 2:
                on = {("arm_neck", 1), ("arm_neck", 3)}
            elif st == 1:
                on = {("corner_neck", 0), ("corner_neck", 2)}
            else:
                on = {("corner_neck", 1), ("corner_neck", 3)}
            if self.mode == "notch":
                perp = {
                    0: {("arm_neck", 1), ("arm_neck", 3)},
                    2: {("arm_neck", 0), ("arm_neck", 2)},
                    1: {("corner_neck", 1), ("corner_neck", 3)},
                    3: {("corner_neck", 0), ("corner_neck", 2)},
                }[st]
                on = {(s[0], s[1]) for s in self._sites()} - perp
            elif self.mode != "bowtie":
                raise ValueError(f"unknown mode {self.mode!r}")
            out.append(tuple(ON if (s[0], s[1]) in on else OFF for s in self._sites()))
        return out

    def diode_sites(self):
        return [(x, y, length, width, d) for (_, _, x, y, length, width, d) in self._sites()]

    def unit_cell(self, state, kit):
        if len(state) != 8:
            raise ValueError("state must have 8 entries")
        loads = [
            diode_load(kit.model(s), x, y, length, width, direction=d, transverse_conducts=False)
            for s, (x, y, length, width, d) in zip(state, self.diode_sites(), strict=True)
        ]
        return UnitCell(
            RectangularLattice(self.period, self.period),
            self.geometry(),
            RT5880,
            Metal(),
            backing="pec",
            loads=tuple(loads),
        )


# ---------------------------------------------------------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class RationalSectorGeometry:
    """Duck-typed `geometry.Geometry`: the hub-and-sector idea with every edge either axis-aligned or of RATIONAL slope (0, 1/2, 1, 2, infinity), so that the rooftop
    mesh represents it by a self-similar staircase (no rasterisation jitter between meshes when the dimensions are multiples of the mesh step; the radial wedges of
    `SectorHubGeometry`, at 22.5 degrees, are not). With u = |x|, v = |y| (all pieces are symmetric under x -> -x, y -> -y and x <-> y):

    hub:        max(u, v) <= hub_half
    arm (+-x):  r_in <= u <= arm_outer and v <= u / 2 - gap_t           (+-y: x <-> y)
    diagonal:   r_in <= max(u, v) <= diag_outer and v >= u / 2 + gap_t and u >= v / 2 + gap_t        (one in each quadrant)
    necks:      arm necks along the axis (width `neck_width`); the diagonal piece of the ++ quadrant is joined to the hub's right edge by a horizontal neck, the other three
                by its 90-degree rotations (C4 symmetry, a pinwheel).
    `r_in = hub_half + neck_gap`. The whole pattern is clipped to |x|, |y| <= clip."""

    hub_half: float
    neck_gap: float
    neck_width: float
    gap_t: float
    arm_outer: float
    diag_outer: float
    clip: float
    neck_width_diag: "float | None" = None  # neck width of the diagonal pieces; None = `neck_width`

    kind = "rational_sector_hub"
    polarity = "metal"
    rotation = 0.0

    @property
    def r_in(self):
        return self.hub_half + self.neck_gap

    @property
    def feature_area(self):
        return 4 * self.clip**2

    def corners(self):
        return np.array(
            [[self.clip, self.clip], [-self.clip, self.clip], [-self.clip, -self.clip], [self.clip, -self.clip]]
        )

    def diag_neck_base(self):
        """Horizontal neck of the ++ diagonal piece: `(x0, x1, y0, y1)` inside the allowed interval of the piece at u = r_in and the hub's right edge."""
        r_in, h, t, w = self.r_in, self.hub_half, self.gap_t, self.neck_width_diag or self.neck_width
        lo = r_in / 2 + t
        hi = min(h, 2 * (r_in - t))
        yc = 0.5 * (lo + hi)
        wn = min(w, hi - lo)
        return (h, r_in, yc - wn / 2, yc + wn / 2)

    def neck_rect(self, kind, k):
        x0, x1, y0, y1 = (
            (self.hub_half, self.r_in, -self.neck_width / 2, self.neck_width / 2)
            if kind == "arm_neck"
            else self.diag_neck_base()
        )
        for _ in range(k % 4):
            x0, x1, y0, y1 = -y1, -y0, x0, x1
        return x0, x1, y0, y1

    def in_feature(self, x, y):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        u, v = np.abs(x), np.abs(y)
        r_in, t = self.r_in, self.gap_t
        m = np.maximum(u, v)
        hub = m <= self.hub_half
        arm_x = (u >= r_in) & (u <= self.arm_outer) & (v <= u / 2 - t)
        arm_y = (v >= r_in) & (v <= self.arm_outer) & (u <= v / 2 - t)
        diag = (m >= r_in) & (m <= self.diag_outer) & (v >= u / 2 + t) & (u >= v / 2 + t)
        out = hub | arm_x | arm_y | diag
        for kind in ("arm_neck", "corner_neck"):
            for k in range(4):
                x0, x1, y0, y1 = self.neck_rect("arm_neck" if kind == "arm_neck" else "diag", k)
                out |= (x >= x0 - 1e-12) & (x <= x1 + 1e-12) & (y >= y0 - 1e-12) & (y <= y1 + 1e-12)
        return out & (np.abs(x) <= self.clip + 1e-12) & (np.abs(y) <= self.clip + 1e-12)


@dataclass(frozen=True)
class RationalSectorCell:
    """Cell on `RationalSectorGeometry` with eight diode necks; states as `ManhattanHubCell` (0: arm necks +-x, 1: diagonal necks of the ++ and -- pieces, 2: arm necks +-y, 3:
    the other diagonal pair; `mode` 'bowtie' = those two necks ON, 'notch' = all ON except the perpendicular two)."""

    period: float = 8.0e-3
    hub_half: float = 1.2e-3
    neck_gap: float = 0.4e-3
    neck_width: float = 0.4e-3
    gap_t: float = 0.2e-3
    arm_outer: float = 3.6e-3
    diag_outer: float = 3.6e-3
    clip: float = 3.8e-3
    mode: str = "bowtie"
    neck_width_diag: "float | None" = None

    kind = "rational_sector_hub"
    n_diodes = 8

    def geometry(self):
        return RationalSectorGeometry(
            self.hub_half,
            self.neck_gap,
            self.neck_width,
            self.gap_t,
            self.arm_outer,
            self.diag_outer,
            self.clip,
            self.neck_width_diag,
        )

    def _sites(self):
        geo = self.geometry()
        sites = []
        for kind, name in (("arm_neck", "arm_neck"), ("corner_neck", "diag")):
            for k in range(4):
                x0, x1, y0, y1 = geo.neck_rect(name, k)
                horizontal = k % 2 == 0
                length = (x1 - x0) if horizontal else (y1 - y0)
                width = (y1 - y0) if horizontal else (x1 - x0)
                sites.append((kind, k, 0.5 * (x0 + x1), 0.5 * (y0 + y1), length, width, "x" if horizontal else "y"))
        return sites

    def states(self):
        out = []
        pairs = {
            0: {("arm_neck", 0), ("arm_neck", 2)},
            1: {("corner_neck", 0), ("corner_neck", 2)},
            2: {("arm_neck", 1), ("arm_neck", 3)},
            3: {("corner_neck", 1), ("corner_neck", 3)},
        }
        perp = {0: pairs[2], 1: pairs[3], 2: pairs[0], 3: pairs[1]}
        for st in range(4):
            if self.mode == "notch":
                on = {(s[0], s[1]) for s in self._sites()} - perp[st]
            elif self.mode == "bowtie":
                on = pairs[st]
            else:
                raise ValueError(f"unknown mode {self.mode!r}")
            out.append(tuple(ON if (s[0], s[1]) in on else OFF for s in self._sites()))
        return out

    def diode_sites(self):
        return [(x, y, length, width, d) for (_, _, x, y, length, width, d) in self._sites()]

    def unit_cell(self, state, kit):
        if len(state) != 8:
            raise ValueError("state must have 8 entries")
        loads = [
            diode_load(kit.model(s), x, y, length, width, direction=d, transverse_conducts=False)
            for s, (x, y, length, width, d) in zip(state, self.diode_sites(), strict=True)
        ]
        return UnitCell(
            RectangularLattice(self.period, self.period),
            self.geometry(),
            RT5880,
            Metal(),
            backing="pec",
            loads=tuple(loads),
        )
