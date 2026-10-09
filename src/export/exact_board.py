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
    return ExactBoard(
        patches=tuple(patches),
        feeds_mm=tuple(feeds),
        port_rects_mm=tuple(rects),
        x_min=min(r[0] for r in rects),
        x_max=max(r[1] for r in rects),
        y_min=min(r[2] for r in rects),
        y_max=max(r[3] for r in rects),
        thickness_mm=quantize(thickness_mm),
        period_mm=period,
        max_rounding_mm=float(worst),
    )


def write_layout_csv(board, path):
    with open(path, "w") as f:
        f.write("port,row,col,x_min_mm,y_min_mm,x_max_mm,y_max_mm\n")
        for p in board.patches:
            f.write(",".join([str(p.port), str(p.i), str(p.j), fmt(p.x0), fmt(p.y0), fmt(p.x1), fmt(p.y1)]) + "\n")
    return Path(path)


_POINT_RE = re.compile(r"CARTESIAN_POINT\([^,]*,\(([^)]*)\)\)")


def verify_step_exact(path):
    """Parses every CARTESIAN_POINT of the STEP file at `path`; raises ValueError if any coordinate is
    off the 0.001 mm grid or prints with float residue. Returns the number of points checked."""
    text = Path(path).read_text()
    values = []
    for m in _POINT_RE.finditer(text):
        values += [t.strip() for t in m.group(1).split(",")]
    off_grid = [t for t in values if abs(float(t) * 1000 - round(float(t) * 1000)) > 1e-9]
    ugly = [t for t in values if re.search(r"(9999|0000\d)", t)]
    if off_grid or ugly:
        raise ValueError(
            f"STEP export is not exact: {len(off_grid)} coordinate(s) off the 0.001 mm grid "
            f"(e.g. {off_grid[:3]}), {len(ugly)} printed with float residue (e.g. {ugly[:3]})"
        )
    return len(values)


def write_exact_step(board, path):
    """Substrate box, Ground (Face.makePlane at the origin when the board is centred there, else an
    exact polygon face) and Patches (one compound of exact rectangular faces at z = thickness)."""
    try:
        import cadquery as cq
    except ImportError as exc:  # pragma: no cover - cadquery is an optional dependency
        raise cad_unavailable(exc) from exc

    def vec(x, y, z):
        return cq.Vector(float(x), float(y), float(z))

    def rect(x0, y0, x1, y1, z):
        wire = cq.Wire.makePolygon([vec(x0, y0, z), vec(x1, y0, z), vec(x1, y1, z), vec(x0, y1, z)], close=True)
        return cq.Face.makeFromWires(wire)

    h = board.thickness_mm
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
    assembly.save(str(path))
    verify_step_exact(path)
    return Path(path)
