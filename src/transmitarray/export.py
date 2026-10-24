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

    return pathlib.Path(path).read_text(errors="ignore", encoding="utf-8")


# ----------------------------------------------------------------------------------------------------------------------
# Exact coordinates (Phase 19 Part D-bis).  OpenCASCADE prints STEP reals with 12 significant digits and its boolean kernel
# leaves ~1e-16..1e-10 noise on values that should be exact (0.085 -> 8.500000000001E-02, 0 -> 1.110223024625E-16), so the
# written file is post-processed: every coordinate within SNAP_TOL_MM of a value the builder deliberately placed is set to
# exactly that double, and every real is re-serialised with Python's shortest round-trip repr.
# ----------------------------------------------------------------------------------------------------------------------
SNAP_TOL_MM = 1e-7
UNCERTAINTY_MM = 1e-6


def intended_coordinates(cells, d_flats_mm, pitch_mm, plate_mm, total_mm):
    """Values the builder placed on purpose, taken from the layout table / design inputs (never recomputed by a kernel).

    Returns dict(xy=sorted array, z=sorted array, dir=sorted array).  ``xy`` holds: the cell-square corners and centres
    (centre +- pitch/2), the rounded outline vertices, every octagon vertex coordinate (rounded to ROUND_MM), and the
    intersections of the 45-degree outline flats with the cell-grid lines (+-(c -+ g), c = |x|+|y| of an outline flat) --
    those are exact decimals because the outline vertices are.  What is NOT in the set (irrational, cannot be exact in
    decimal): intersections of pillar / pocket octagon edges with the outline flats, and the edge directions of the rounded
    octagons.
    """
    h = pitch_mm / 2
    xy = {0.0}
    out = outline_pts(d_flats_mm)
    for px, py in out:
        xy.update((px, py))
    c = round(max(abs(px) + abs(py) for px, py in out), 4)
    grid = set()
    for cell in cells:
        for v in (cell["x_mm"], cell["y_mm"]):
            grid.update((v, round(v - h, 4), round(v + h, 4)))
        if cell["radius_mm"] > 0:
            for px, py in octagon_pts(cell["x_mm"], cell["y_mm"], cell["radius_mm"]):
                xy.update((px, py))
    xy.update(grid)
    for g in list(grid):
        for s in (+1, -1):
            xy.update((round(c - s * g, 4), round(-(c - s * g), 4)))
    xy.update(-v for v in list(xy))
    z = {0.0, float(plate_mm), float(total_mm), round(total_mm - plate_mm, 4)}
    z.update(-v for v in list(z))
    root_half = math.sqrt(0.5)
    return dict(
        xy=np.array(sorted(xy)),
        z=np.array(sorted(z)),
        dir=np.array(sorted({0.0, 1.0, -1.0, root_half, -root_half})),
    )


def _snap_value(v, table, tol=SNAP_TOL_MM):
    i = int(np.searchsorted(table, v))
    best = None
    for j in (i - 1, i):
        if 0 <= j < len(table) and abs(v - table[j]) <= tol and (best is None or abs(v - table[j]) < abs(v - best)):
            best = float(table[j])
    return best


def step_real(v):
    """STEP REAL for double ``v`` using Python's shortest round-trip repr: 4.0 -> '4.', 0.085 -> '0.085', 1e-16 -> '1.E-16'."""
    if v == 0:
        return "0."
    r = repr(float(v))
    if "e" in r or "E" in r:
        mant, exp = r.lower().split("e")
        if "." not in mant:
            mant += "."
        return f"{mant}E{int(exp)}"
    if "." in r:
        r = r.rstrip("0")
    return r if "." in r else r + "."


_POINT = None


def _regexes():
    import re

    global _POINT
    if _POINT is None:
        _POINT = dict(
            point=re.compile(r"CARTESIAN_POINT\('([^']*)',\s*\(([^()]*)\)\s*\)"),
            direction=re.compile(r"DIRECTION\('([^']*)',\s*\(([^()]*)\)\s*\)"),
            vector=re.compile(r"VECTOR\('([^']*)',\s*(#\d+)\s*,\s*([^()]*?)\s*\)"),
            uncertainty=re.compile(r"(UNCERTAINTY_MEASURE_WITH_UNIT\(LENGTH_MEASURE\()([^()]*)(\))"),
            description=re.compile(r"FILE_DESCRIPTION\(\(.*?\),\s*'[^']*'\);", re.S),
        )
    return _POINT


def snap_step_text(text, intended, uncertainty_mm=UNCERTAINTY_MM, tol=SNAP_TOL_MM, description=None):
    """Snap + re-serialise every CARTESIAN_POINT / DIRECTION / VECTOR magnitude.  Returns (new_text, stats)."""
    rx = _regexes()
    union = np.array(sorted(set(intended["xy"]) | set(intended["z"])))
    stats = dict(n_tokens=0, n_snapped=0, max_snap_mm=0.0, vector_magnitudes={})

    def fix(tok, table):
        v = float(tok)
        stats["n_tokens"] += 1
        s = _snap_value(v, table, tol)
        if s is not None:
            if s != v:
                stats["n_snapped"] += 1
                stats["max_snap_mm"] = max(stats["max_snap_mm"], abs(s - v))
            v = s
        return step_real(v)

    def point(m):
        toks = [t.strip() for t in m.group(2).split(",")]
        # 3-D points use the per-axis sets; 2-D parameter-space points on a face / edge only the union of design values
        tabs = (intended["xy"], intended["xy"], intended["z"]) if len(toks) == 3 else (union,) * len(toks)
        return f"CARTESIAN_POINT('{m.group(1)}',({','.join(fix(t, tb) for t, tb in zip(toks, tabs, strict=True))}))"

    def direction(m):
        toks = [t.strip() for t in m.group(2).split(",")]
        return f"DIRECTION('{m.group(1)}',({','.join(fix(t, intended['dir']) for t in toks)}))"

    def vector(m):
        mag = fix(m.group(3), np.array([1.0]))
        stats["vector_magnitudes"][mag] = stats["vector_magnitudes"].get(mag, 0) + 1
        return f"VECTOR('{m.group(1)}',{m.group(2)},{mag})"

    text = rx["point"].sub(point, text)
    text = rx["direction"].sub(direction, text)
    text = rx["vector"].sub(vector, text)
    text, n = rx["uncertainty"].subn(lambda m: m.group(1) + step_real(uncertainty_mm) + m.group(3), text)
    stats["uncertainty_entities_set"] = n
    if description:
        text = rx["description"].sub(f"FILE_DESCRIPTION(('{description}'),'2;1');", text, count=1)
    return text, stats


def inexact_literals(text, intended, tol=SNAP_TOL_MM):
    """Coordinates in ``text`` that sit within ``tol`` of a design value but are not exactly that double.

    This is the regression predicate: an empty list means every near-design-value literal parses to exactly the
    intended double (doubles compared, not rounded strings).  A single-ulp error makes it non-empty.
    """
    rx = _regexes()
    union = np.array(sorted(set(intended["xy"]) | set(intended["z"])))
    bad = []

    def check(tok, table, kind):
        v = float(tok)
        s = _snap_value(v, table, tol)
        if s is not None and s != v:
            bad.append((kind, tok.strip(), s))

    for m in rx["point"].finditer(text):
        toks = m.group(2).split(",")
        tabs = (intended["xy"], intended["xy"], intended["z"]) if len(toks) == 3 else (union,) * len(toks)
        for t, tb in zip(toks, tabs, strict=True):
            check(t, tb, "CARTESIAN_POINT")
    for m in rx["direction"].finditer(text):
        for t in m.group(2).split(","):
            check(t, intended["dir"], "DIRECTION")
    for m in rx["vector"].finditer(text):
        check(m.group(3), np.array([1.0]), "VECTOR")
    return bad


def _entity_points(text):
    import re

    pts = {}
    for m in re.finditer(r"#(\d+)\s*=\s*CARTESIAN_POINT\('[^']*',\s*\(([^()]*)\)\s*\)", text):
        c = tuple(float(t) for t in m.group(2).split(","))
        if len(c) == 3:
            pts[int(m.group(1))] = c
    return pts


def verify_snap(raw_text, snapped_text, distinct_tol=SNAP_TOL_MM):
    """Snap safety: (1) no two previously distinct vertices became coincident, (2) no edge collapsed.

    'Previously distinct' = raw positions differing by more than ``distinct_tol``.  Vertices that merge from a raw
    separation <= distinct_tol are kernel noise and are only counted.
    """
    import re
    from collections import defaultdict

    raw, snp = _entity_points(raw_text), _entity_points(snapped_text)
    verts = {int(a): int(b) for a, b in re.findall(r"#(\d+)\s*=\s*VERTEX_POINT\('[^']*',#(\d+)\)", snapped_text)}
    groups = defaultdict(list)
    for pid in verts.values():
        groups[snp[pid]].append(raw[pid])
    violations, noise_merges, max_dist = [], 0, 0.0
    for s, raws in groups.items():
        uniq = set(raws)
        if len(uniq) > 1:
            noise_merges += 1
            lst = list(uniq)
            d = max(math.dist(lst[i], lst[j]) for i in range(len(lst)) for j in range(i + 1, len(lst)))
            max_dist = max(max_dist, d)
            if d > distinct_tol:
                violations.append((s, d))
    collapsed = 0
    for a, b in re.findall(r"EDGE_CURVE\('[^']*',#(\d+),#(\d+),", snapped_text):
        pa, pb = verts[int(a)], verts[int(b)]
        if snp[pa] == snp[pb] and raw[pa] != raw[pb] and math.dist(raw[pa], raw[pb]) > 0:
            collapsed += 1
    # smallest separation between distinct snapped vertex positions (sets the usable header tolerance)
    arr = np.array(sorted(groups))
    min_sep, n_close = float("inf"), 0
    if len(arr) > 1:
        from scipy.spatial import cKDTree

        tree = cKDTree(arr)
        d, _ = tree.query(arr, k=2)
        min_sep = float(d[:, 1].min())
        n_close = len(tree.query_pairs(0.6))  # vertex pairs closer than the 0.6 mm minimum feature
    return dict(
        n_vertices=len(verts),
        n_distinct_snapped_vertices=len(groups),
        merged_groups_from_noise=noise_merges,
        max_raw_separation_in_merged_group_mm=max_dist,
        distinct_merge_violations=violations,
        collapsed_edges=collapsed,
        min_distinct_vertex_separation_mm=min_sep,
        n_vertex_pairs_closer_than_0p6mm=n_close,
    )


def self_intersection_free(shape):
    """OpenCASCADE argument analyzer (self-interference, small edges, vertex/edge merge candidates, tangency): True iff clean.

    Sanity-checked in the tests against a deliberately self-overlapping compound, which it must flag.
    """
    from OCP.BOPAlgo import BOPAlgo_ArgumentAnalyzer

    an = BOPAlgo_ArgumentAnalyzer()
    an.SetShape1(shape.wrapped)
    for mode in (
        "SelfInterMode",
        "SmallEdgeMode",
        "MergeVertexMode",
        "MergeEdgeMode",
        "TangentMode",
        "CurveOnSurfaceMode",
    ):
        setattr(an, mode, True)
    an.Perform()
    return not an.HasFaulty()


def export_exact_step(workplane, path, cells, d_flats_mm, pitch_mm, plate_mm, total_mm, check_self_intersection=True):
    """Write ``path`` with exact design coordinates; returns the verification report (also asserts nothing silently).

    Steps: kernel STEP (raw) -> snap + shortest-round-trip rewrite + explicit uncertainty -> re-import both ->
    vertex-coincidence check, validity, optional self-intersection, volume before/after.  No tolerance is loosened to pass:
    the report carries the raw numbers and a boolean per assertion.
    """
    import pathlib
    import tempfile

    cq = _require_cadquery()
    from OCP.BRepCheck import BRepCheck_Analyzer

    path = pathlib.Path(path)
    intended = intended_coordinates(cells, d_flats_mm, pitch_mm, plate_mm, total_mm)
    with tempfile.TemporaryDirectory() as td:
        raw_path = pathlib.Path(td) / "raw.step"
        save_step_named(workplane, raw_path)
        raw_text = raw_path.read_text(encoding="utf-8")
        desc = (
            f"Nylon dual-beam lens, millimetres; coordinates snapped to design values within {SNAP_TOL_MM:g} mm; "
            f"uncertainty {UNCERTAINTY_MM:g} mm"
        )
        snapped_text, stats = snap_step_text(raw_text, intended, description=desc)
        path.write_text(snapped_text, encoding="utf-8", newline="\n")
        raw_shape = cq.importers.importStep(str(raw_path)).val()
        new_shape = cq.importers.importStep(str(path)).val()
    vol_raw, vol_new = raw_shape.Volume(), new_shape.Volume()
    vsnap = verify_snap(raw_text, snapped_text)
    report = dict(
        snap=stats,
        vertex_check=vsnap,
        no_distinct_vertices_merged=not vsnap["distinct_merge_violations"] and vsnap["collapsed_edges"] == 0,
        n_solids=len(new_shape.Solids()),
        valid=all(BRepCheck_Analyzer(s.wrapped).IsValid() for s in new_shape.Solids()),
        volume_before_mm3=vol_raw,
        volume_after_mm3=vol_new,
        volume_rel_diff=abs(vol_new - vol_raw) / vol_raw,
        volume_within_1e9=abs(vol_new - vol_raw) / vol_raw < 1e-9,
        inexact_literals_after=len(inexact_literals(snapped_text, intended)),
        inexact_literals_before=len(inexact_literals(raw_text, intended)),
        uncertainty_mm=UNCERTAINTY_MM,
        snap_tol_mm=SNAP_TOL_MM,
    )
    if check_self_intersection:
        report["self_intersection_free"] = self_intersection_free(new_shape)
    return report, intended
