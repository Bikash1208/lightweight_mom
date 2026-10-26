"""Phase 26 B2: fixed-outline slotted patches -- geometries that duck-type `geometry.Geometry` (as
`AnnularSlotPatch` does), so `UnitCell`, `RooftopBasis` and the MoM need no change.

* `SlottedCircularPatch`: a circular patch (diameter sets the mean loading) with a centred rectangular slot
  (length / width set the anisotropy, `slot_angle` the axis). The outline is rotation-invariant: every slot
  angle exists at every diameter, and the metal mask rotates with the slot only.
* `SlottedSquarePatch`: a square patch with the same slot (comparison family; its outline is NOT rotation
  invariant -- only the slot turns).
"""

import math
from dataclasses import dataclass

import numpy as np

from .._validate import finite_real


def _in_slot(dx, dy, angle, length, width):
    c, s = math.cos(angle), math.sin(angle)
    u = dx * c + dy * s
    v = -dx * s + dy * c
    return (np.abs(u) <= length / 2) & (np.abs(v) <= width / 2)


@dataclass(frozen=True)
class SlottedCircularPatch:
    diameter: float
    slot_length: float
    slot_width: float
    slot_angle: float = 0.0  # [rad], axis of the slot's long side from +x
    center_x: float = 0.0
    center_y: float = 0.0

    kind = "patch_circle_rect_slot"

    def __post_init__(self):
        for name in ("diameter", "slot_length", "slot_width", "slot_angle", "center_x", "center_y"):
            object.__setattr__(self, name, finite_real(name, getattr(self, name)))
        if self.diameter <= 0.0:
            raise ValueError("diameter must be > 0")
        if self.slot_length < 0.0 or self.slot_width < 0.0:
            raise ValueError("slot_length and slot_width must be >= 0")

    @property
    def polarity(self):
        return "metal"

    @property
    def rotation(self):
        return self.slot_angle

    @property
    def feature_area(self):
        # circle minus the part of the slot inside it (slot assumed to fit when length <= diameter)
        area = math.pi * self.diameter**2 / 4.0
        if self.slot_length > 0 and self.slot_width > 0:
            area -= self.slot_length * self.slot_width
        return area

    def corners(self):
        r = self.diameter / 2.0
        return np.array([[r, r], [-r, r], [-r, -r], [r, -r]]) + np.array([self.center_x, self.center_y])

    def in_feature(self, x, y):
        dx = np.asarray(x, dtype=float) - self.center_x
        dy = np.asarray(y, dtype=float) - self.center_y
        inside = np.hypot(dx, dy) <= self.diameter / 2.0
        if self.slot_length > 0.0 and self.slot_width > 0.0:
            inside &= ~_in_slot(dx, dy, self.slot_angle, self.slot_length, self.slot_width)
        return inside


@dataclass(frozen=True)
class SlottedSquarePatch:
    side: float
    slot_length: float
    slot_width: float
    slot_angle: float = 0.0
    center_x: float = 0.0
    center_y: float = 0.0

    kind = "patch_square_rect_slot"

    def __post_init__(self):
        for name in ("side", "slot_length", "slot_width", "slot_angle", "center_x", "center_y"):
            object.__setattr__(self, name, finite_real(name, getattr(self, name)))
        if self.side <= 0.0:
            raise ValueError("side must be > 0")

    @property
    def polarity(self):
        return "metal"

    @property
    def rotation(self):
        return 0.0

    @property
    def feature_area(self):
        area = self.side**2
        if self.slot_length > 0 and self.slot_width > 0:
            area -= self.slot_length * self.slot_width
        return area

    def corners(self):
        r = self.side / 2.0
        return np.array([[r, r], [-r, r], [-r, -r], [r, -r]]) + np.array([self.center_x, self.center_y])

    def in_feature(self, x, y):
        dx = np.asarray(x, dtype=float) - self.center_x
        dy = np.asarray(y, dtype=float) - self.center_y
        inside = (np.abs(dx) <= self.side / 2.0) & (np.abs(dy) <= self.side / 2.0)
        if self.slot_length > 0.0 and self.slot_width > 0.0:
            inside &= ~_in_slot(dx, dy, self.slot_angle, self.slot_length, self.slot_width)
        return inside


@dataclass(frozen=True)
class CrossSlottedCircularPatch:
    """Circular patch with TWO orthogonal centred rectangular slots (a cross): the 4-fold-symmetric, isotropic-
    tensor way of adding the slot loading (Phase 26 V2 launch pad). `slot_angle` rotates the cross (no effect on
    the tensor, kept for the exporter)."""

    diameter: float
    slot_length: float
    slot_width: float
    slot_angle: float = 0.0
    center_x: float = 0.0
    center_y: float = 0.0

    kind = "patch_circle_cross_slot"

    def __post_init__(self):
        for name in ("diameter", "slot_length", "slot_width", "slot_angle", "center_x", "center_y"):
            object.__setattr__(self, name, finite_real(name, getattr(self, name)))
        if self.diameter <= 0.0:
            raise ValueError("diameter must be > 0")
        if self.slot_length < 0.0 or self.slot_width < 0.0:
            raise ValueError("slot_length and slot_width must be >= 0")

    @property
    def polarity(self):
        return "metal"

    @property
    def rotation(self):
        return 0.0

    @property
    def feature_area(self):
        area = math.pi * self.diameter**2 / 4.0
        if self.slot_length > 0 and self.slot_width > 0:
            area -= 2.0 * self.slot_length * self.slot_width - self.slot_width**2
        return area

    def corners(self):
        r = self.diameter / 2.0
        return np.array([[r, r], [-r, r], [-r, -r], [r, -r]]) + np.array([self.center_x, self.center_y])

    def in_feature(self, x, y):
        dx = np.asarray(x, dtype=float) - self.center_x
        dy = np.asarray(y, dtype=float) - self.center_y
        inside = np.hypot(dx, dy) <= self.diameter / 2.0
        if self.slot_length > 0.0 and self.slot_width > 0.0:
            inside &= ~_in_slot(dx, dy, self.slot_angle, self.slot_length, self.slot_width)
            inside &= ~_in_slot(dx, dy, self.slot_angle + math.pi / 2.0, self.slot_length, self.slot_width)
        return inside


@dataclass(frozen=True)
class CellRow:
    """A row of `cells` (each a slotted patch centred in its own `pitch`-wide cell) as ONE supercell pattern:
    cell i sits at `x = -N pitch / 2 + (i + 1/2) pitch`. Duck-types `Geometry` like the cells it holds, so a
    `UnitCell(RectangularLattice(N pitch, pitch), CellRow(...))` is the supercell of a ramp."""

    cells: tuple
    pitch: float

    kind = "cell_row"
    polarity = "metal"
    rotation = 0.0

    def _centres(self):
        n = len(self.cells)
        return [-n * self.pitch / 2.0 + (i + 0.5) * self.pitch for i in range(n)]

    @property
    def feature_area(self):
        return sum(c.feature_area for c in self.cells)

    def corners(self):
        n = len(self.cells)
        w = n * self.pitch / 2.0
        h = max(max(abs(c.corners()[:, 1]).max() for c in self.cells), 0.0)
        return np.array([[w, h], [-w, h], [-w, -h], [w, -h]]) * 0.999

    def in_feature(self, x, y):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        out = np.zeros(np.broadcast(x, y).shape, dtype=bool)
        for cx, cell in zip(self._centres(), self.cells, strict=True):
            out |= cell.in_feature(x - cx, y)
        return out
