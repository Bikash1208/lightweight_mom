"""Phase 30 item 2: elliptical / super-elliptical patch cell family ("grain of rice"), optionally with a centred rectangular slot
or a cross slot -- additive geometry that duck-types `geometry.Geometry` (as `slotted_cells` does), so `UnitCell`,
`RooftopBasis` and the MoM need no change.

`EllipticalPatch(major, minor, rotation, slot_length, slot_width, slot_angle, cross)`: full axes `major >= minor`
(diameters), the major axis at `rotation` from +x; the slot (if `slot_length`, `slot_width` > 0) is centred, its
long side at `slot_angle` RELATIVE to the major axis (the whole cell, slot included, then turns with `rotation`);
`cross=True` adds the orthogonal second slot. `exponent` n is the super-ellipse exponent
(|u/a|^n + |v/b|^n <= 1: n = 2 ellipse, n = 4 squircle, large n -> rectangle; default 2 = the plain ellipse, unchanged). With `major == minor`, no slot and `rotation == 0` the mask is
identical to `slotted_cells.SlottedCircularPatch` without a slot (tested).
"""

import math
from dataclasses import dataclass

import numpy as np

from .._validate import finite_real
from .slotted_cells import _in_slot


@dataclass(frozen=True)
class EllipticalPatch:
    major: float
    minor: float
    rotation: float = 0.0  # [rad], major axis from +x
    slot_length: float = 0.0
    slot_width: float = 0.0
    slot_angle: float = 0.0  # [rad], slot long side from the major axis
    cross: bool = False
    exponent: float = 2.0
    center_x: float = 0.0
    center_y: float = 0.0

    kind = "patch_ellipse_slot"

    def __post_init__(self):
        for name in ("major", "minor", "rotation", "slot_length", "slot_width", "slot_angle", "exponent", "center_x", "center_y"):
            object.__setattr__(self, name, finite_real(name, getattr(self, name)))
        if self.exponent < 1.0:
            raise ValueError("exponent must be >= 1")
        if self.minor <= 0.0 or self.major < self.minor:
            raise ValueError(f"need major >= minor > 0, got major={self.major}, minor={self.minor}")
        if self.slot_length < 0.0 or self.slot_width < 0.0:
            raise ValueError("slot_length and slot_width must be >= 0")

    @property
    def polarity(self):
        return "metal"

    @property
    def _slotted(self):
        return self.slot_length > 0.0 and self.slot_width > 0.0

    @property
    def feature_area(self):
        """Analytic area of the ideal outline (slot overlap with the ellipse is taken as full: a slot inside the
        outline, as every sensible cell has)."""
        n = self.exponent
        area = self.major * self.minor * math.gamma(1.0 + 1.0 / n) ** 2 / math.gamma(1.0 + 2.0 / n)  # n = 2: pi a b / 4
        if self._slotted:
            slot = self.slot_length * self.slot_width
            area -= (2.0 * slot - self.slot_width**2) if self.cross else slot
        return area

    def corners(self):
        r = self.major / 2.0
        return np.array([[r, r], [-r, r], [-r, -r], [r, -r]]) + np.array([self.center_x, self.center_y])

    def in_feature(self, x, y):
        dx = np.asarray(x, dtype=float) - self.center_x
        dy = np.asarray(y, dtype=float) - self.center_y
        c, s = math.cos(self.rotation), math.sin(self.rotation)
        u = dx * c + dy * s  # along the major axis
        v = -dx * s + dy * c
        if self.exponent == 2.0:
            inside = (u / (self.major / 2.0)) ** 2 + (v / (self.minor / 2.0)) ** 2 <= 1.0
        else:
            inside = np.abs(u / (self.major / 2.0)) ** self.exponent + np.abs(v / (self.minor / 2.0)) ** self.exponent <= 1.0
        if self._slotted:
            inside &= ~_in_slot(u, v, self.slot_angle, self.slot_length, self.slot_width)
            if self.cross:
                inside &= ~_in_slot(u, v, self.slot_angle + math.pi / 2.0, self.slot_length, self.slot_width)
        return inside


@dataclass(frozen=True)
class MeanderSlotPatch:
    """Circular metal disc with one continuous square-wave cut.

    ``tooth_count`` is the number of turns, giving ``tooth_count + 1``
    horizontal slot runs. ``tooth_width`` is the copper finger between runs.
    An optional removed tooth is replaced by an electrically isolated island;
    this is geometry only and has no diode or switching model.
    All dimensions are metres and ``orientation`` is radians.
    """

    diameter: float
    slot_span: float
    tooth_count: int
    tooth_width: float
    slot_width: float
    orientation: float = 0.0
    removed_tooth: int | None = None
    island_length: float = 0.0
    island_clearance: float = 0.2e-3
    center_x: float = 0.0
    center_y: float = 0.0

    kind = "patch_circle_meander_slot"

    def __post_init__(self):
        for name in ("diameter", "slot_span", "tooth_width", "slot_width", "orientation",
                     "island_length", "island_clearance", "center_x", "center_y"):
            object.__setattr__(self, name, finite_real(name, getattr(self, name)))
        if isinstance(self.tooth_count, bool) or not isinstance(self.tooth_count, int) or self.tooth_count < 1:
            raise ValueError("tooth_count must be a positive integer")
        if self.slot_width < 0.2e-3 or self.tooth_width < 0.2e-3:
            raise ValueError("slot width and copper finger gap must each be >= 0.2 mm")
        if self.slot_span <= 2 * self.slot_width:
            raise ValueError("slot_span must exceed twice the slot width")
        if self.diameter <= 0:
            raise ValueError("diameter must be positive")
        half_y = (self.tooth_count * (self.tooth_width + self.slot_width) + self.slot_width) / 2
        if math.hypot(self.slot_span / 2, half_y) + 0.2e-3 > self.diameter / 2:
            raise ValueError("meander cut must clear the disc rim by at least 0.2 mm")
        if self.removed_tooth is not None:
            if isinstance(self.removed_tooth, bool) or not isinstance(self.removed_tooth, int) or not (0 <= self.removed_tooth < self.tooth_count):
                raise ValueError("removed_tooth must index a copper finger between slot runs")
            if self.island_clearance < 0.2e-3 or self.tooth_width <= 2 * self.island_clearance:
                raise ValueError("island needs >= 0.2 mm clearance on both sides")
            if not (0 < self.island_length < self.slot_span - 2 * self.slot_width - 2 * self.island_clearance):
                raise ValueError("island_length must fit within the removed tooth")
        elif self.island_length != 0:
            raise ValueError("island_length requires removed_tooth")

    @property
    def polarity(self):
        return "metal"

    @property
    def rotation(self):
        return self.orientation

    @property
    def feature_area(self):
        pitch = self.tooth_width + self.slot_width
        n = self.tooth_count + 1
        cut = n * self.slot_span * self.slot_width + (n - 1) * self.slot_width * (pitch - self.slot_width)
        if self.removed_tooth is not None:
            opening = self.island_length + 2 * self.island_clearance
            cut += opening * self.tooth_width - self.island_length * (self.tooth_width - 2 * self.island_clearance)
        return math.pi * (self.diameter / 2) ** 2 - cut

    def corners(self):
        r = self.diameter / 2
        return np.array([[r, r], [-r, r], [-r, -r], [r, -r]]) + np.array([self.center_x, self.center_y])

    def in_feature(self, x, y):
        dx = np.asarray(x, dtype=float) - self.center_x
        dy = np.asarray(y, dtype=float) - self.center_y
        c, s = math.cos(self.orientation), math.sin(self.orientation)
        u, v = dx * c + dy * s, -dx * s + dy * c
        metal = np.hypot(u, v) <= self.diameter / 2
        pitch = self.tooth_width + self.slot_width
        first = -self.tooth_count * pitch / 2
        slot = np.zeros(np.broadcast_shapes(np.shape(u), np.shape(v)), dtype=bool)
        for i in range(self.tooth_count + 1):
            slot |= (np.abs(u) <= self.slot_span / 2) & (np.abs(v - (first + i * pitch)) <= self.slot_width / 2)
        for i in range(self.tooth_count):
            end = self.slot_span / 2 - self.slot_width / 2
            ux = end if i % 2 == 0 else -end
            mid = first + (i + 0.5) * pitch
            slot |= (np.abs(u - ux) <= self.slot_width / 2) & (np.abs(v - mid) <= pitch / 2)
        metal &= ~slot
        if self.removed_tooth is not None:
            mid = first + (self.removed_tooth + 0.5) * pitch
            opening = self.island_length + 2 * self.island_clearance
            metal &= ~((np.abs(u) <= opening / 2) & (np.abs(v - mid) <= self.tooth_width / 2))
            island = (np.abs(u) <= self.island_length / 2) & (np.abs(v - mid) <= (self.tooth_width - 2 * self.island_clearance) / 2)
            metal |= island
        return metal


@dataclass(frozen=True)
class VerticalLegMeanderSlotPatch:
    """Continuous square-wave slot: vertical legs pitched along local x.

    ``leg_depth`` is the full vertical run. Adjacent legs join alternately at
    their top and bottom ends. The mask removes the whole slot strip, leaving
    no metal inside it. Local x rotates by ``orientation`` from global +x.
    """

    diameter: float
    leg_count: int
    leg_pitch: float
    leg_depth: float
    slot_width: float
    orientation: float = 0.0
    center_x: float = 0.0
    center_y: float = 0.0

    kind = "patch_circle_vertical_leg_meander_slot"

    def __post_init__(self):
        for name in ("diameter", "leg_pitch", "leg_depth", "slot_width", "orientation", "center_x", "center_y"):
            object.__setattr__(self, name, finite_real(name, getattr(self, name)))
        if isinstance(self.leg_count, bool) or not isinstance(self.leg_count, int) or self.leg_count < 2:
            raise ValueError("leg_count must be an integer >= 2")
        if self.slot_width < 0.2e-3 or self.leg_pitch - self.slot_width < 0.2e-3:
            raise ValueError("slot width and copper finger gap must each be >= 0.2 mm")
        if self.leg_depth <= self.slot_width or self.diameter <= 0:
            raise ValueError("leg depth and disc diameter must be positive; depth must exceed slot width")
        half_x = (self.leg_count - 1) * self.leg_pitch / 2 + self.slot_width / 2
        half_y = self.leg_depth / 2 + self.slot_width / 2
        if math.hypot(half_x, half_y) + 0.2e-3 > self.diameter / 2:
            raise ValueError("slot must clear disc rim by at least 0.2 mm")

    @property
    def polarity(self):
        return "metal"

    @property
    def rotation(self):
        return self.orientation

    @property
    def feature_area(self):
        n = self.leg_count
        w = self.slot_width
        return math.pi * self.diameter**2 / 4 - n * self.leg_depth * w - (n - 1) * w * (self.leg_pitch - w)

    def corners(self):
        r = self.diameter / 2
        return np.array([[r, r], [-r, r], [-r, -r], [r, -r]]) + np.array([self.center_x, self.center_y])

    def in_feature(self, x, y):
        dx = np.asarray(x, dtype=float) - self.center_x
        dy = np.asarray(y, dtype=float) - self.center_y
        c, s = math.cos(self.orientation), math.sin(self.orientation)
        u, v = dx * c + dy * s, -dx * s + dy * c
        slot = np.zeros(np.broadcast_shapes(np.shape(u), np.shape(v)), dtype=bool)
        first = -(self.leg_count - 1) * self.leg_pitch / 2
        for i in range(self.leg_count):
            slot |= (np.abs(u - (first + i * self.leg_pitch)) <= self.slot_width / 2) & (np.abs(v) <= self.leg_depth / 2)
        for i in range(self.leg_count - 1):
            mid = first + (i + 0.5) * self.leg_pitch
            end = (1 if i % 2 == 0 else -1) * (self.leg_depth / 2 - self.slot_width / 2)
            slot |= (np.abs(u - mid) <= self.leg_pitch / 2) & (np.abs(v - end) <= self.slot_width / 2)
        return (np.hypot(u, v) <= self.diameter / 2) & ~slot
