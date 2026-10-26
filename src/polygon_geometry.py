"""Polygon-composite patch geometry (Phase 20 R3 item 2): a union of simple polygons / axis-aligned rectangles minus a union of apertures,
for non-Manhattan RIS cells (an octagonal patch with 45-degree delay lines).

Duck-types `geometry.Geometry` exactly like `composite_geometry.CompositeGeometry` (`kind`, `polarity`, `feature_area`, `corners()`,
`in_feature()`), so `UnitCell`, `RooftopBasis` and the MoM solve need no changes. Edges that are not axis-aligned are staircased by the
rooftop mesh (as the rotated rectangles of `geometry.Geometry` always were): `edges_aligned` is True only when every edge is axis-aligned
and on a mesh line, and `non_manhattan_edges()` lists the others so a caller can report the staircasing honestly.

`Poly.contains` is the even-odd rule with the boundary inclusive (a point within `BOUNDARY_TOL` of an edge is inside), vectorised over arrays.
`feature_area` is evaluated on a fine grid (no exact polygon-clipping dependency): relative error about 1e-3 for typical cells.
"""

from dataclasses import dataclass

import numpy as np

from ._validate import finite_real
from .composite_geometry import AxisRect

BOUNDARY_TOL = 1e-12  # [m]


@dataclass(frozen=True)
class Poly:
    vertices: tuple  # ((x, y), ...) simple polygon, >= 3 vertices, either orientation

    def __post_init__(self):
        v = tuple((finite_real("x", a), finite_real("y", b)) for a, b in self.vertices)
        if len(v) < 3:
            raise ValueError("a polygon needs at least 3 vertices")
        object.__setattr__(self, "vertices", v)
        if abs(self.signed_area) <= 0.0:
            raise ValueError("polygon has zero area")

    @property
    def array(self):
        return np.asarray(self.vertices, dtype=float)

    @property
    def signed_area(self):
        a = self.array
        x, y = a[:, 0], a[:, 1]
        return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))

    @property
    def area(self):
        return abs(self.signed_area)

    def bbox(self):
        a = self.array
        return float(a[:, 0].min()), float(a[:, 0].max()), float(a[:, 1].min()), float(a[:, 1].max())

    def edges(self):
        a = self.array
        return list(zip(a, np.roll(a, -1, axis=0), strict=True))

    def is_manhattan(self, tol=1e-9):
        return all(abs(p[0] - q[0]) < tol or abs(p[1] - q[1]) < tol for p, q in self.edges())

    def contains(self, x, y):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        shape = np.broadcast(x, y).shape
        px, py = np.broadcast_to(x, shape).ravel(), np.broadcast_to(y, shape).ravel()
        inside = np.zeros(px.shape, dtype=bool)
        on_edge = np.zeros(px.shape, dtype=bool)
        for (x0, y0), (x1, y1) in self.edges():
            dx, dy = x1 - x0, y1 - y0
            cond = (y0 > py) != (y1 > py)
            with np.errstate(divide="ignore", invalid="ignore"):
                xc = x0 + (py - y0) * dx / dy
            inside ^= cond & (px < xc)
            seg2 = dx * dx + dy * dy
            t = np.clip(((px - x0) * dx + (py - y0) * dy) / seg2, 0.0, 1.0)
            on_edge |= np.hypot(px - (x0 + t * dx), py - (y0 + t * dy)) <= BOUNDARY_TOL
        return (inside | on_edge).reshape(shape)


def regular_polygon(cx, cy, circumradius, n=8, rotation=0.0):
    """Regular `n`-gon with the given circumradius; `rotation` [rad] turns the first vertex from +x."""
    ang = rotation + 2.0 * np.pi * np.arange(n) / n
    return Poly(tuple((cx + circumradius * np.cos(a), cy + circumradius * np.sin(a)) for a in ang))


def oriented_rect(cx, cy, length, width, angle):
    """Rectangle of `length` along the direction `angle` [rad] (from +x) and `width` across it, centred at (cx, cy)."""
    if length <= 0.0 or width <= 0.0:
        raise ValueError("length and width must be > 0")
    u = np.array([np.cos(angle), np.sin(angle)])
    v = np.array([-np.sin(angle), np.cos(angle)])
    c = np.array([cx, cy])
    pts = [c + sl * length / 2 * u + sw * width / 2 * v for sl, sw in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
    return Poly(tuple((float(p[0]), float(p[1])) for p in pts))


def _contains(shape, x, y):
    return shape.contains(x, y)  # Poly and AxisRect both provide it


@dataclass(frozen=True)
class PolyComposite:
    metal: tuple  # Poly | AxisRect, ...
    apertures: tuple = ()

    kind = "poly_composite"
    polarity = "metal"
    rotation = 0.0

    def __post_init__(self):
        object.__setattr__(self, "metal", tuple(self.metal))
        object.__setattr__(self, "apertures", tuple(self.apertures))
        if not self.metal:
            raise ValueError("a polygon composite needs at least one metal shape")
        if not all(isinstance(r, Poly | AxisRect) for r in self.metal + self.apertures):
            raise ValueError("metal and apertures must be Poly or AxisRect instances")

    # ---------------------------------------------------------------- geometry.Geometry's interface
    def bbox(self):
        bs = [r.bbox() if isinstance(r, Poly) else (*r.x_edges, *r.y_edges) for r in self.metal]
        return min(b[0] for b in bs), max(b[1] for b in bs), min(b[2] for b in bs), max(b[3] for b in bs)

    @property
    def length(self):
        x0, x1, _, _ = self.bbox()
        return x1 - x0

    @property
    def width(self):
        _, _, y0, y1 = self.bbox()
        return y1 - y0

    def corners(self):
        x0, x1, y0, y1 = self.bbox()
        return np.array([[x1, y1], [x0, y1], [x0, y0], [x1, y0]])

    def in_feature(self, x, y):
        m = _contains(self.metal[0], x, y)
        for r in self.metal[1:]:
            m = m | _contains(r, x, y)
        for r in self.apertures:
            m = m & ~_contains(r, x, y)
        return m

    @property
    def feature_area(self):
        """Area of union(metal) minus union(apertures) on a 1500 x 1500 grid over the bounding box (relative error ~1e-3)."""
        x0, x1, y0, y1 = self.bbox()
        n = 1500
        xs = x0 + (np.arange(n) + 0.5) * (x1 - x0) / n
        ys = y0 + (np.arange(n) + 0.5) * (y1 - y0) / n
        gx, gy = np.meshgrid(xs, ys, indexing="ij")
        return float(np.count_nonzero(self.in_feature(gx, gy)) * (x1 - x0) * (y1 - y0) / (n * n))

    # ---------------------------------------------------------------- mesh alignment
    def non_manhattan_edges(self):
        return [
            (i, p) for i, p in enumerate(self.metal + self.apertures) if isinstance(p, Poly) and not p.is_manhattan()
        ]

    def _axis_lines(self):
        xs, ys = [], []
        for r in self.metal + self.apertures:
            if isinstance(r, AxisRect):
                xs += list(r.x_edges)
                ys += list(r.y_edges)
            else:
                for p, q in r.edges():
                    if abs(p[0] - q[0]) < 1e-9:
                        xs.append(float(p[0]))
                    elif abs(p[1] - q[1]) < 1e-9:
                        ys.append(float(p[1]))
        return xs, ys

    def misaligned_edges(self, lattice, mesh, rtol=1e-6):
        from .discretization import is_position_mesh_aligned

        xs, ys = self._axis_lines()
        bad = [("x", p) for p in xs if not is_position_mesh_aligned(p, lattice.period_x, mesh.cells_x, rtol)]
        bad += [("y", p) for p in ys if not is_position_mesh_aligned(p, lattice.period_y, mesh.cells_y, rtol)]
        return bad

    def edges_aligned(self, lattice, mesh, rtol=1e-6):
        """True only if there is no non-axis-aligned edge and every axis-aligned edge is on a mesh line."""
        return not self.non_manhattan_edges() and not self.misaligned_edges(lattice, mesh, rtol)
