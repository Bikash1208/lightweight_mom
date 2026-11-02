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
