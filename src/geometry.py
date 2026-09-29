"""Parametric unit-cell geometry.

Cell coordinates: origin at the cell centre, x in [-Px/2, Px/2], y in
[-Py/2, Py/2]. The metallic pattern lies in the z = 0 plane.

All supported features are rectangles of ``length`` (along the feature's local
u axis) and ``width`` (along local v), centred at (center_x, center_y) and
rotated counter-clockwise by ``rotation`` [rad] about their centre:

    u =  (x-cx)*cos(rot) + (y-cy)*sin(rot)
    v = -(x-cx)*sin(rot) + (y-cy)*cos(rot)
    inside feature  <=>  |u| <= length/2  and  |v| <= width/2

=========  ========  ==========================================================
kind       polarity  meaning
=========  ========  ==========================================================
rectangle  metal     metal rectangle, no shape restriction
patch      metal     metal rectangle (same maths as ``rectangle``; kept as a
                     separate name for design intent)
dipole     metal     metal strip, requires length >= width
slot       aperture  rectangular hole in an otherwise solid metal cell,
                     requires length >= width
=========  ========  ==========================================================

Features must lie inside the cell (touching the cell edge is allowed; a feature
touching opposite edges forms a continuous conductor/aperture across cells).
Wrap-around features are not supported in Phase 1.

This module only defines the continuous geometry. Meshing into MoM basis
functions is a later phase.
"""

import math
from dataclasses import dataclass

import numpy as np

from ._validate import finite_real
from .materials import Metal, Substrate
from .periodic_bc import RectangularLattice

POLARITY = {
    "rectangle": "metal",
    "patch": "metal",
    "dipole": "metal",
    "slot": "aperture",
}
_ELONGATED = {"dipole", "slot"}

# Relative tolerance (of the larger period) for "inside the cell" checks, so a
# feature exactly filling the cell is not rejected by rounding in cos/sin.
_BOUNDS_RTOL = 1e-12


@dataclass(frozen=True)
class Geometry:
    kind: str
    length: float
    width: float
    center_x: float = 0.0
    center_y: float = 0.0
    rotation: float = 0.0  # [rad]

    def __post_init__(self):
        if self.kind not in POLARITY:
            raise ValueError(f"unknown geometry type {self.kind!r}; supported: {sorted(POLARITY)}")
        for name in ("length", "width", "center_x", "center_y", "rotation"):
            object.__setattr__(self, name, finite_real(name, getattr(self, name)))
        if self.length <= 0.0 or self.width <= 0.0:
            raise ValueError(
                f"length and width must be > 0, got length={self.length}, width={self.width}"
            )
        if self.kind in _ELONGATED and self.length < self.width:
            raise ValueError(
                f"{self.kind} requires length >= width, got length={self.length}, width={self.width}"
            )

    @property
    def polarity(self):
        return POLARITY[self.kind]

    @property
    def feature_area(self):
        return self.length * self.width

    def corners(self):
        """(4, 2) array of feature corner coordinates, counter-clockwise."""
        hl, hw = self.length / 2, self.width / 2
        local = np.array([[hl, hw], [-hl, hw], [-hl, -hw], [hl, -hw]])
        c, s = math.cos(self.rotation), math.sin(self.rotation)
        rot = np.array([[c, -s], [s, c]])
        return local @ rot.T + np.array([self.center_x, self.center_y])

    def in_feature(self, x, y):
        """True where (x, y) is inside the rectangle (closed set). No periodic wrap."""
        dx = np.asarray(x, dtype=float) - self.center_x
        dy = np.asarray(y, dtype=float) - self.center_y
        c, s = math.cos(self.rotation), math.sin(self.rotation)
        u = dx * c + dy * s
        v = -dx * s + dy * c
        return (np.abs(u) <= self.length / 2) & (np.abs(v) <= self.width / 2)


_BACKINGS = ("air", "pec")


@dataclass(frozen=True)
class UnitCell:
    """Lattice + pattern + materials: everything that defines one periodic cell.

    ``backing``: the stack's bottom termination at z = -h (Phase 12A). ``"air"`` (the default,
    unchanged since Phase 1) is the ungrounded stack air | pattern | substrate | air. ``"pec"``
    adds a PEC ground plane at z = -h (a reflectarray/grounded-metasurface stack); the bottom
    port then does not exist (GREEN_FUNCTION.md, S_PARAMETERS.md). This field only selects which
    closed form ``LayeredGreenFunction`` uses; it does not change ``Geometry`` or metal-pattern
    handling at z = 0.
    """

    lattice: RectangularLattice
    geometry: Geometry
    substrate: Substrate
    metal: Metal
    backing: str = "air"

    def __post_init__(self):
        if self.backing not in _BACKINGS:
            raise ValueError(f"backing must be one of {_BACKINGS}, got {self.backing!r}")
        px, py = self.lattice.period_x, self.lattice.period_y
        tol = _BOUNDS_RTOL * max(px, py)
        xy = self.geometry.corners()
        if (
            np.any(np.abs(xy[:, 0]) > px / 2 + tol)
            or np.any(np.abs(xy[:, 1]) > py / 2 + tol)
        ):
            raise ValueError(
                f"{self.geometry.kind} feature extends outside the unit cell "
                f"(Px={px}, Py={py}); corners={xy.tolist()}"
            )
        if self.metal_area <= _BOUNDS_RTOL * self.lattice.cell_area:
            raise ValueError("pattern leaves no metal in the unit cell")

    @property
    def metal_area(self):
        a = self.geometry.feature_area
        if self.geometry.polarity == "metal":
            return a
        return self.lattice.cell_area - a

    @property
    def metal_fill_fraction(self):
        return self.metal_area / self.lattice.cell_area

    def is_metal(self, x, y):
        """Boolean array: True where the (periodically wrapped) point is metal.

        Points are wrapped into [-P/2, P/2) so the function is exactly
        periodic with the lattice.
        """
        px, py = self.lattice.period_x, self.lattice.period_y
        xw = np.asarray(x, dtype=float)
        yw = np.asarray(y, dtype=float)
        xw = xw - px * np.floor(xw / px + 0.5)
        yw = yw - py * np.floor(yw / py + 0.5)
        inside = self.geometry.in_feature(xw, yw)
        return inside if self.geometry.polarity == "metal" else ~inside
