"""Composite patch geometry (Phase 20 R1 item 3): a union of axis-aligned metal rectangles minus a union of axis-aligned
aperture rectangles, so a patch with gaps (and arms of different lengths) can be described.

Duck-types `geometry.Geometry` (`kind`, `polarity`, `feature_area`, `corners()`, `in_feature()`), so `UnitCell`, `RooftopBasis`,
the MoM solve and `discretization` need no changes -- the same extension point `AnnularSlotPatch` uses. `in_feature(x, y)` for a
single metal rectangle performs exactly the operations of `Geometry.in_feature` at rotation 0, so a composite of one rectangle gives
the SAME metal mask and bit-identical S-parameters as the `rectangle` geometry (tested).

A rectangle is `(cx, cy, length, width)`: centre, extent along x, extent along y [m]. Rectangles may overlap (union). Metal =
union(metal) minus union(apertures); apertures are cut out of the metal only (they do not add metal).
"""

from dataclasses import dataclass

import numpy as np

from ._validate import finite_real


@dataclass(frozen=True)
class AxisRect:
    cx: float
    cy: float
    length: float  # extent along x [m]
    width: float  # extent along y [m]

    def __post_init__(self):
        for name in ("cx", "cy", "length", "width"):
            object.__setattr__(self, name, finite_real(name, getattr(self, name)))
        if self.length <= 0.0 or self.width <= 0.0:
            raise ValueError(f"rectangle length and width must be > 0, got {self.length}, {self.width}")

    @property
    def x_edges(self):
        return self.cx - self.length / 2.0, self.cx + self.length / 2.0

    @property
    def y_edges(self):
        return self.cy - self.width / 2.0, self.cy + self.width / 2.0

    def contains(self, x, y):
        dx = np.asarray(x, dtype=float) - self.cx
        dy = np.asarray(y, dtype=float) - self.cy
        return (np.abs(dx) <= self.length / 2) & (np.abs(dy) <= self.width / 2)


@dataclass(frozen=True)
class CompositeGeometry:
    metal: tuple  # AxisRect, ...
    apertures: tuple = ()

    kind = "composite"
    polarity = "metal"
    rotation = 0.0  # axis-aligned only

    def __post_init__(self):
        object.__setattr__(self, "metal", tuple(self.metal))
        object.__setattr__(self, "apertures", tuple(self.apertures))
        if not self.metal:
            raise ValueError("a composite geometry needs at least one metal rectangle")
        if not all(isinstance(r, AxisRect) for r in self.metal + self.apertures):
            raise ValueError("metal and apertures must be AxisRect instances")

    # ---------------------------------------------------------------- geometry.Geometry's interface
    @property
    def length(self):  # bounding box (used only by generic tools)
        x0, x1, _, _ = self.bbox()
        return x1 - x0

    @property
    def width(self):
        _, _, y0, y1 = self.bbox()
        return y1 - y0

    def bbox(self):
        x0 = min(r.x_edges[0] for r in self.metal)
        x1 = max(r.x_edges[1] for r in self.metal)
        y0 = min(r.y_edges[0] for r in self.metal)
        y1 = max(r.y_edges[1] for r in self.metal)
        return x0, x1, y0, y1

    def corners(self):
        x0, x1, y0, y1 = self.bbox()
        return np.array([[x1, y1], [x0, y1], [x0, y0], [x1, y0]])

    def in_feature(self, x, y):
        m = self.metal[0].contains(x, y)
        for r in self.metal[1:]:
            m = m | r.contains(x, y)
        for r in self.apertures:
            m = m & ~r.contains(x, y)
        return m

    @property
    def feature_area(self):
        """Exact area of union(metal) minus union(apertures) (coordinate compression)."""
        rects = self.metal + self.apertures
        xs = np.unique([e for r in rects for e in r.x_edges])
        ys = np.unique([e for r in rects for e in r.y_edges])
        cx, cy = 0.5 * (xs[:-1] + xs[1:]), 0.5 * (ys[:-1] + ys[1:])
        gx, gy = np.meshgrid(cx, cy, indexing="ij")
        inside = self.in_feature(gx, gy)
        return float(np.sum(inside * np.outer(np.diff(xs), np.diff(ys))))

    # ---------------------------------------------------------------- mesh alignment
    def edge_positions(self):
        """`(x_edges, y_edges)` of every rectangle (metal and aperture), the lines the mesh should contain."""
        rects = self.metal + self.apertures
        return [e for r in rects for e in r.x_edges], [e for r in rects for e in r.y_edges]

    def misaligned_edges(self, lattice, mesh, rtol=1e-6):
        """Edges (axis, position) that do not lie on a mesh line of the `mesh.cells_x x mesh.cells_y` mesh."""
        from .discretization import is_position_mesh_aligned

        xs, ys = self.edge_positions()
        bad = [("x", p) for p in xs if not is_position_mesh_aligned(p, lattice.period_x, mesh.cells_x, rtol)]
        bad += [("y", p) for p in ys if not is_position_mesh_aligned(p, lattice.period_y, mesh.cells_y, rtol)]
        return bad

    def edges_aligned(self, lattice, mesh, rtol=1e-6):
        return not self.misaligned_edges(lattice, mesh, rtol)
