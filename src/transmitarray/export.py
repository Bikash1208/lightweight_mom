"""STEP / STL / layout export of the finished lens as ONE fused solid named "Lens" (needs the optional cadquery extra).

Coordinates are millimetres, z = 0 is the exit (bottom) face, z = total thickness is the feed-side face; the lens is
centred on the axis, octagon flats on the x/y axes.  Octagon vertices are rounded to ``ROUND_MM`` (no 4.9999999 values),
and the volume check uses the *same rounded vertices*.

Geometry: base plate (outline prism, 0..plate) + per-cell upper layer (plate..total): pillar cells add an octagonal
prism, hole cells a full cell square; every hole cell's octagon is then cut straight through plate and lens
(through-holes); finally everything is intersected with the octagonal outline prism.
"""

from __future__ import annotations

import math

import numpy as np

ROUND_MM = 1e-4
LENS_NAME = "Lens"


def _require_cadquery():
    try:
        import cadquery as cq
    except ImportError as e:  # pragma: no cover
        raise RuntimeError("cadquery is required (pip install -r requirements-cad.txt)") from e
    return cq


def octagon_pts(cx, cy, radius_mm, rotation_deg=22.5):
    ang = np.deg2rad(rotation_deg + 45.0 * np.arange(8))
    return [(round(cx + radius_mm * math.cos(a), 4), round(cy + radius_mm * math.sin(a), 4)) for a in ang]


def outline_pts(d_flats_mm):
    return octagon_pts(0.0, 0.0, d_flats_mm / 2 / math.cos(math.radians(22.5)))


def _prism(cq, pts, z0, z1):
    return cq.Workplane("XY").workplane(offset=z0).polyline(pts).close().extrude(z1 - z0)


def build_lens(cells, d_flats_mm, pitch_mm, plate_mm, total_mm):
    """cells: iterable of dicts with x_mm, y_mm, kind, radius_mm.  Returns a cadquery Workplane (one solid)."""
    cq = _require_cadquery()
    out = outline_pts(d_flats_mm)
    plate = _prism(cq, out, 0.0, plate_mm).val()
    adds, cuts = [], []
    for c in cells:
        x, y, r = c["x_mm"], c["y_mm"], c["radius_mm"]
        if c["kind"] == "pillar":
            if r > 0:
                adds.append(_prism(cq, octagon_pts(x, y, r), plate_mm, total_mm).val())
        else:
            h = pitch_mm / 2
            sq = [(x - h, y - h), (x + h, y - h), (x + h, y + h), (x - h, y + h)]
            adds.append(_prism(cq, sq, plate_mm, total_mm).val())
            if r > 0:
                cuts.append(_prism(cq, octagon_pts(x, y, r), -1.0, total_mm + 1.0).val())
    body = plate.fuse(*adds) if adds else plate
    if cuts:
        body = body.cut(*cuts)
    body = body.intersect(_prism(cq, out, 0.0, total_mm).val())
    return cq.Workplane("XY").add(body.clean())


def save_step_named(workplane, path, name=LENS_NAME):
    cq = _require_cadquery()
    assy = cq.Assembly(workplane, name=name)
    assy.save(str(path), exportType="STEP")


def save_stl(workplane, path, tolerance=0.005, angular=0.1):
    cq = _require_cadquery()
    cq.exporters.export(workplane, str(path), tolerance=tolerance, angularTolerance=angular)


def raster_volume_mm3(cells, d_flats_mm, pitch_mm, plate_mm, total_mm, step_mm=0.04):
    """Independent volume estimate from a pixel raster of the layout (no CAD kernel involved)."""
    a = d_flats_mm / 2
    n = int(round(d_flats_mm / step_mm))
    xs = -a + (np.arange(n) + 0.5) * step_mm
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    inside = (np.abs(X) + np.abs(Y)) <= a * math.sqrt(2)
    upper = np.zeros_like(inside)
    holes = np.zeros_like(inside)
    h = pitch_mm / 2
    for c in cells:
        x, y, r = c["x_mm"], c["y_mm"], c["radius_mm"]
        ix0 = max(int((x - h + a) / step_mm) - 1, 0)
        ix1 = min(int((x + h + a) / step_mm) + 2, n)
        iy0 = max(int((y - h + a) / step_mm) - 1, 0)
        iy1 = min(int((y + h + a) / step_mm) + 2, n)
        sl = (slice(ix0, ix1), slice(iy0, iy1))
        px, py = X[sl], Y[sl]
        poly = np.ones(px.shape, bool)
        if r > 0:
            pts = octagon_pts(x, y, r)
            for i in range(8):
                p0, p1 = pts[i], pts[(i + 1) % 8]
                poly &= (p1[0] - p0[0]) * (py - p0[1]) - (p1[1] - p0[1]) * (px - p0[0]) >= 0
        else:
            poly[:] = False
        if c["kind"] == "pillar":
            upper[sl] |= poly
        else:
            upper[sl] |= (np.abs(px - x) <= h) & (np.abs(py - y) <= h)
            holes[sl] |= poly
    cell_area = step_mm**2
    plate_area = np.sum(inside & ~holes) * cell_area
    upper_area = np.sum(inside & upper & ~holes) * cell_area
    return float(plate_area * plate_mm + upper_area * (total_mm - plate_mm))


def verify_step(path, cells, d_flats_mm, pitch_mm, plate_mm, total_mm, vol_tol=0.004):
    """Re-import the STEP from disk and check one valid solid, bbox, volume vs raster, and pillar-top-face count."""
    cq = _require_cadquery()
    from OCP.BRepCheck import BRepCheck_Analyzer

    shape = cq.importers.importStep(str(path)).val()
    solids = shape.Solids()
    bb = shape.BoundingBox()
    vol = shape.Volume()
    ref = raster_volume_mm3(cells, d_flats_mm, pitch_mm, plate_mm, total_mm)
    text = pathlib_read(path)
    out = outline_pts(d_flats_mm)
    # pillar tops: one planar +z single-wire face at z = total whose centroid lies within 1.5 mm of each pillar centre
    # (a clipped pillar keeps its top face; clipped hole-cell fragments are >= 2.5 mm from every pillar centre)
    pillars = [c for c in cells if c["kind"] == "pillar" and c["radius_mm"] > 0]
    tops = [
        f.Center()
        for f in shape.Faces()
        if abs(f.Center().z - total_mm) < 1e-6
        and f.geomType() == "PLANE"
        and len(f.Wires()) == 1
        and f.Area() <= pitch_mm**2 * 1.0001
    ]
    matched = sum(1 for c in pillars if sum(1 for t in tops if math.hypot(t.x - c["x_mm"], t.y - c["y_mm"]) < 1.5) == 1)
    return dict(
        n_solids=len(solids),
        valid=all(BRepCheck_Analyzer(s.wrapped).IsValid() for s in solids),
        bbox_mm=[bb.xmax - bb.xmin, bb.ymax - bb.ymin, bb.zmax - bb.zmin],
        bbox_expected_mm=[d_flats_mm, d_flats_mm, total_mm],
        volume_mm3=vol,
        raster_volume_mm3=ref,
        volume_rel_diff=abs(vol - ref) / ref,
        volume_ok=abs(vol - ref) / ref < vol_tol,
        name_in_step=LENS_NAME in text,
        pillar_cells=len(pillars),
        pillar_top_faces_found=matched,
        outline_vertices=out,
    )


def pathlib_read(path):
    import pathlib

    return pathlib.Path(path).read_text(errors="ignore")
