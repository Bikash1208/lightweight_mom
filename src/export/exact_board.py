"""Phase 21 D1: an EXACT multi-port board model and its STEP/CSV writers.

A multi-port design's cells are placed at each port's OWN exact positions (`x = x_mid + (j - (nx-1)/2)
* period`, Decimal arithmetic -- no shared lattice, no float residue), every dimension is rounded to
the 0.001 mm grid, and the STEP is built from exact vertex coordinates (no transforms). After
writing, `verify_step_exact` parses every CARTESIAN_POINT and fails the export if a coordinate is
off the 0.001 mm grid or prints with float residue (`1.9999...`, `2.0000001`).

Rectangular patches only: the exact path has no annular-slot solid (use `export.array_export` for
those). Methodology source: `scripts/adhoc/square4/build_sq_export.py` (the hand-built reference)."""

import re
from dataclasses import dataclass
from decimal import ROUND_HALF_EVEN, Decimal
from pathlib import Path

from ..errors import cad_unavailable, config_invalid

GRID_MM = Decimal("0.001")


def dec(value):
    """Exact Decimal of a float via its shortest repr (5.85 -> Decimal('5.85'), not 5.8499999...)."""
    return Decimal(repr(float(value))) if not isinstance(value, Decimal) else value


def quantize(value):
    return dec(value).quantize(GRID_MM, rounding=ROUND_HALF_EVEN)


def fmt(value):
    """A Decimal as a plain, trimmed string ('5.85', '-14', '0')."""
    text = format(value.normalize(), "f")
    return "0" if text in ("-0", "") else text


@dataclass(frozen=True)
class ExactPatch:
    port: int  # 1-based
    i: int
    j: int
    x0: Decimal  # lower-left corner, mm
    y0: Decimal
    x1: Decimal
    y1: Decimal


@dataclass(frozen=True)
class ExactBoard:
    patches: "tuple[ExactPatch, ...]"
    feeds_mm: "tuple[tuple[Decimal, Decimal], ...]"  # one per port, port order
    port_rects_mm: "tuple[tuple[Decimal, Decimal, Decimal, Decimal], ...]"  # x0, x1, y0, y1 per port
    x_min: Decimal
    x_max: Decimal
    y_min: Decimal
    y_max: Decimal
    thickness_mm: Decimal
    period_mm: Decimal
    max_rounding_mm: float  # largest change made by snapping to the 0.001 mm grid
    # Phase 23: per-port polygon outline `((x, y), ...)` in mm (counter-clockwise), `None` for a
    # rectangular port; `()` when no port is a polygon. `trim_board`: substrate/ground follow the
    # outlines (union) instead of the bounding rectangle.
    port_outlines_mm: tuple = ()
    trim_board: bool = False

    @property
    def width_mm(self):
        return self.x_max - self.x_min

    @property
    def height_mm(self):
        return self.y_max - self.y_min


def exact_board_from_ports_result(result, period_mm, thickness_mm):
    """Builds an `ExactBoard` from a `/synthesize_ports` job result (`ports[*].patches_mm`)."""
    ports = result.get("ports") or []
    if not ports or any("patches_mm" not in p for p in ports):
        raise config_invalid(message="the source job has no per-port cell layout ('patches_mm') -- exact export needs a multi-port design")
    period = dec(period_mm)
    patches, feeds, rects = [], [], []
    worst = Decimal(0)

    def snap(v):
        nonlocal worst
        q = quantize(v)
        worst = max(worst, abs(q - dec(v)))
        return q

    for p in ports:
        layout = p["patches_mm"]
        x_mid, y_mid = snap(layout["x_mid_mm"]), snap(layout["y_mid_mm"])
        nx, ny = layout["nx"], layout["ny"]
        rects.append(
            (snap(p["x_min_mm"]), snap(p["x_max_mm"]), snap(p["y_min_mm"]), snap(p["y_max_mm"]))
        )
        feeds.append((snap(p["feed_x_mm"]), snap(p["feed_y_mm"])))
        for c in layout["cells"]:
            geometry = c["geometry"]
            if isinstance(geometry, dict):
                raise config_invalid(
                    message="exact STEP export supports rectangular patches only; this design has annular-slot cells "
                    "(use the standard STEP export)"
                )
            w, h = dec(geometry[0]), dec(geometry[1])
            cx = x_mid + (Decimal(c["j"]) - Decimal(nx - 1) / 2) * period
            cy = y_mid + (Decimal(c["i"]) - Decimal(ny - 1) / 2) * period
            x0, y0 = snap(cx - w / 2), snap(cy - h / 2)
            patches.append(ExactPatch(p["port_id"] + 1, c["i"], c["j"], x0, y0, snap(cx + w / 2), snap(cy + h / 2)))
    outlines = ()
    trim = False
    if any("outline_mm" in p for p in ports):
        # a shared hexagon corner is snapped from numerically equal inputs, so neighbours still coincide
        outlines = tuple(
            tuple((snap(x), snap(y)) for x, y in p["outline_mm"]) if p.get("outline_mm") else None for p in ports
        )
        trim = bool(result.get("trim_board", False))
    x_min, x_max = min(r[0] for r in rects), max(r[1] for r in rects)
    y_min, y_max = min(r[2] for r in rects), max(r[3] for r in rects)
    if trim:
        pts = [pt for o, r in zip(outlines, rects, strict=True) for pt in (o if o else ((r[0], r[2]), (r[1], r[3])))]
        x_min, x_max = min(q[0] for q in pts), max(q[0] for q in pts)
        y_min, y_max = min(q[1] for q in pts), max(q[1] for q in pts)
    return ExactBoard(
        patches=tuple(patches),
        feeds_mm=tuple(feeds),
        port_rects_mm=tuple(rects),
        x_min=x_min,
        x_max=x_max,
        y_min=y_min,
        y_max=y_max,
        thickness_mm=quantize(thickness_mm),
        period_mm=period,
        max_rounding_mm=float(worst),
        port_outlines_mm=outlines,
        trim_board=trim,
    )


def write_layout_csv(board, path):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("port,row,col,x_min_mm,y_min_mm,x_max_mm,y_max_mm\n")
        for p in board.patches:
            f.write(",".join([str(p.port), str(p.i), str(p.j), fmt(p.x0), fmt(p.y0), fmt(p.x1), fmt(p.y1)]) + "\n")
    return Path(path)


DEFAULT_R_PIN_MM = "0.635"  # monopole pin radius -- the same defaults as the AEDT script's design variables
DEFAULT_MONO_H_MM = "7.1"  # pin height from the ground (z = 0)
FEED_MARKER_HALF_MM = Decimal("0.3")  # feed markers are 0.6 x 0.6 mm square sheets centred on the feed

_POINT_LINE_RE = re.compile(r"#(\d+)\s*=\s*CARTESIAN_POINT\('[^']*',\(([^)]*)\)\)")
_VERTEX_RE = re.compile(r"VERTEX_POINT\('[^']*',#(\d+)\)")


def _step_points(text):
    """`{entity id: (floats...)}` of every CARTESIAN_POINT, and the set of ids that are 3-D vertex points."""
    points = {int(m.group(1)): tuple(float(t) for t in m.group(2).split(",")) for m in _POINT_LINE_RE.finditer(text)}
    vertices = {int(m.group(1)) for m in _VERTEX_RE.finditer(text)}
    return points, vertices


def _bad_coordinates(coords):
    raw = [repr(c) for c in coords]
    off_grid = [c for c in coords if abs(c * 1000 - round(c * 1000)) > 1e-9]
    ugly = [t for t in raw if re.search(r"(9999|0000\d)", t)]
    return off_grid, ugly


def verify_step_exact(path, strict_points=True, pins=()):
    """Raises ValueError if a STEP coordinate is off the 0.001 mm grid or prints with float residue.

    By default (`strict_points=True`, a board of planar faces) every 3-D CARTESIAN_POINT is checked.
    Cylinders (monopole pins) add circle geometry whose 2-D pcurve points and surface-placement
    parameters (angles such as 6.28318530718) can never be on the grid, so with `strict_points=False`
    only the 3-D VERTEX/edge-end points are checked. `pins`: `(x_mm, y_mm, radius_mm, height_mm)` per
    pin -- each must have vertices at exactly z = 0 and z = height, at exactly `radius` from its axis.
    Returns the number of coordinates checked."""
    text = Path(path).read_text()
    points, vertices = _step_points(text)
    if strict_points:
        checked = [p for p in points.values() if len(p) == 3]
    else:
        checked = [points[i] for i in vertices if i in points and len(points[i]) == 3]
    off_grid, ugly = [], []
    for p in checked:
        o, u = _bad_coordinates(p)
        off_grid += o
        ugly += u
    if off_grid or ugly:
        raise ValueError(
            f"STEP export is not exact: {len(off_grid)} coordinate(s) off the 0.001 mm grid "
            f"(e.g. {off_grid[:3]}), {len(ugly)} printed with float residue (e.g. {ugly[:3]})"
        )
    vertex_xyz = [points[i] for i in vertices if i in points and len(points[i]) == 3]
    for x, y, r, h in pins:
        near = [v for v in vertex_xyz if abs(((v[0] - x) ** 2 + (v[1] - y) ** 2) ** 0.5 - r) < 1e-9]
        if not any(v[2] == 0.0 for v in near) or not any(v[2] == float(h) for v in near):
            raise ValueError(f"monopole pin at ({x}, {y}) is not exact: needs vertices at z=0 and z={h}, radius {r}")
    return len(checked)


def write_exact_step(board, path, include_pins=False, r_pin_mm=DEFAULT_R_PIN_MM, mono_h_mm=DEFAULT_MONO_H_MM):
    """Substrate box, Ground (Face.makePlane at the origin when the board is centred there, else an
    exact polygon face), Patches (one compound of exact rectangular faces at z = thickness) and
    Feed1..FeedN: one 0.6 x 0.6 mm square sheet per port at z = 0, corners at feed +-0.3 mm (exact
    vertices -- a circle would leave pcurve residues). `include_pins=True` adds PEC pins
    Monopole1..MonopoleN: cylinders of radius `r_pin_mm` and height `mono_h_mm` from z = 0 at each feed."""
    try:
        import cadquery as cq
    except ImportError as exc:  # pragma: no cover - cadquery is an optional dependency
        raise cad_unavailable(exc) from exc

    def vec(x, y, z):
        return cq.Vector(float(x), float(y), float(z))

    def rect(x0, y0, x1, y1, z):
        wire = cq.Wire.makePolygon([vec(x0, y0, z), vec(x1, y0, z), vec(x1, y1, z), vec(x0, y1, z)], close=True)
        return cq.Face.makeFromWires(wire)

    def polygon_face(points, z):
        # The plane origin is pinned to the first vertex: `Face.makeFromWires` would put it at the
        # polygon's centroid, a CARTESIAN_POINT that is not on the 0.001 mm grid.
        from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace
        from OCP.gp import gp_Dir, gp_Pln, gp_Pnt

        wire = cq.Wire.makePolygon([vec(x, y, z) for x, y in points], close=True)
        plane = gp_Pln(gp_Pnt(float(points[0][0]), float(points[0][1]), float(z)), gp_Dir(0, 0, 1))
        return cq.Face(BRepBuilderAPI_MakeFace(plane, wire.wrapped).Face())

    h = board.thickness_mm
    if board.trim_board and any(board.port_outlines_mm):
        # substrate = the port outlines extruded to the thickness, ground = their faces. Kept as
        # separate prisms in one compound (a boolean union would create new, off-grid vertices where
        # snapped edges meet) -- they touch edge to edge and carry the exact vertex coordinates.
        outlines = [o if o else ((r[0], r[2]), (r[1], r[2]), (r[1], r[3]), (r[0], r[3]))
                    for o, r in zip(board.port_outlines_mm, board.port_rects_mm, strict=True)]
        prisms = [cq.Solid.extrudeLinear(polygon_face(o, 0), cq.Vector(0, 0, float(h))) for o in outlines]
        substrate = cq.Compound.makeCompound(prisms)
        ground = cq.Compound.makeCompound([polygon_face(o, 0) for o in outlines])
    else:
        substrate = cq.Solid.makeBox(float(board.width_mm), float(board.height_mm), float(h), vec(board.x_min, board.y_min, 0))
        if board.x_min == -board.x_max and board.y_min == -board.y_max:
            ground = cq.Face.makePlane(float(board.width_mm), float(board.height_mm), cq.Vector(0, 0, 0), cq.Vector(0, 0, 1))
        else:
            ground = rect(board.x_min, board.y_min, board.x_max, board.y_max, 0)
    patches = cq.Compound.makeCompound([rect(p.x0, p.y0, p.x1, p.y1, h) for p in board.patches])
    assembly = cq.Assembly(name="HMSA_board")
    assembly.add(substrate, name="Substrate")
    assembly.add(ground, name="Ground")
    assembly.add(patches, name="Patches")
    for k, (fx, fy) in enumerate(board.feeds_mm, start=1):
        m = FEED_MARKER_HALF_MM
        assembly.add(rect(quantize(fx - m), quantize(fy - m), quantize(fx + m), quantize(fy + m), 0), name=f"Feed{k}")
    pins = []
    if include_pins:
        r_pin, mono_h = quantize(r_pin_mm), quantize(mono_h_mm)
        for k, (fx, fy) in enumerate(board.feeds_mm, start=1):
            pin = cq.Solid.makeCylinder(float(r_pin), float(mono_h), vec(fx, fy, 0), cq.Vector(0, 0, 1))
            assembly.add(pin, name=f"Monopole{k}")
            pins.append((float(fx), float(fy), float(r_pin), float(mono_h)))
    assembly.save(str(path))
    verify_step_exact(path, strict_points=not include_pins, pins=pins)
    return Path(path)


def write_feeds_csv(board, path):
    """`port,x_mm,y_mm` -- one row per port, in port order."""
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("port,x_mm,y_mm\n")
        for k, (fx, fy) in enumerate(board.feeds_mm, start=1):
            f.write(f"{k},{fmt(fx)},{fmt(fy)}\n")
    return Path(path)


def write_layout_json(board, path):
    """Patches (`port,row,col` and the rectangle corners, mm) and the feed table, as JSON."""
    import json

    data = {
        "units": "mm",
        "board": {"x_min": fmt(board.x_min), "x_max": fmt(board.x_max), "y_min": fmt(board.y_min), "y_max": fmt(board.y_max)},
        "feeds": [{"port": k, "x_mm": fmt(fx), "y_mm": fmt(fy)} for k, (fx, fy) in enumerate(board.feeds_mm, start=1)],
        "patches": [
            {"port": p.port, "row": p.i, "col": p.j, "x_min_mm": fmt(p.x0), "y_min_mm": fmt(p.y0),
             "x_max_mm": fmt(p.x1), "y_max_mm": fmt(p.y1)}
            for p in board.patches
        ],
    }
    if any(board.port_outlines_mm):
        data["port_outlines"] = [
            None if o is None else [[fmt(x), fmt(y)] for x, y in o] for o in board.port_outlines_mm
        ]
        data["trim_board"] = board.trim_board
    Path(path).write_text(json.dumps(data, indent=1), encoding="utf-8", newline="\n")
    return Path(path)


def write_layout_dxf(board, path):
    """DXF R12 (millimetres): patch outlines on layer PATCHES, the board outline on BOARD, and one marker
    (CIRCLE r = 0.3 mm plus a "P<k>" label) per feed on layer FEEDS. A board with polygon (hexagon)
    ports adds layer PORTS: every sub-aperture's outline (Phase 23)."""
    parts = ["0", "SECTION", "2", "HEADER", "9", "$ACADVER", "1", "AC1009", "0", "ENDSEC", "0", "SECTION", "2", "ENTITIES"]

    def line(layer, x1, y1, x2, y2):
        return ["0", "LINE", "8", layer, "10", fmt(x1), "20", fmt(y1), "30", "0.0", "11", fmt(x2), "21", fmt(y2), "31", "0.0"]

    def rectangle(layer, x0, y0, x1, y1):
        return (line(layer, x0, y0, x1, y0) + line(layer, x1, y0, x1, y1) + line(layer, x1, y1, x0, y1)
                + line(layer, x0, y1, x0, y0))

    def polygon(layer, pts):
        out = []
        for k in range(len(pts)):
            (xa, ya), (xb, yb) = pts[k], pts[(k + 1) % len(pts)]
            out += line(layer, xa, ya, xb, yb)
        return out

    has_outlines = any(board.port_outlines_mm)
    if has_outlines and board.trim_board:
        for o, r in zip(board.port_outlines_mm, board.port_rects_mm, strict=True):
            parts += polygon("BOARD", o) if o else rectangle("BOARD", r[0], r[2], r[1], r[3])
    else:
        parts += rectangle("BOARD", board.x_min, board.y_min, board.x_max, board.y_max)
    if has_outlines:  # Phase 23: every sub-aperture's outline on its own layer
        for o, r in zip(board.port_outlines_mm, board.port_rects_mm, strict=True):
            parts += polygon("PORTS", o) if o else rectangle("PORTS", r[0], r[2], r[1], r[3])
    for p in board.patches:
        parts += rectangle("PATCHES", p.x0, p.y0, p.x1, p.y1)
    for k, (fx, fy) in enumerate(board.feeds_mm, start=1):
        parts += ["0", "CIRCLE", "8", "FEEDS", "10", fmt(fx), "20", fmt(fy), "30", "0.0", "40", fmt(FEED_MARKER_HALF_MM)]
        parts += ["0", "TEXT", "8", "FEEDS", "10", fmt(fx + 1), "20", fmt(fy + 1), "30", "0.0", "40", "3.0", "1", f"P{k}"]
    parts += ["0", "ENDSEC", "0", "EOF"]
    Path(path).write_text("\n".join(parts) + "\n", encoding="utf-8", newline="\n")
    return Path(path)
