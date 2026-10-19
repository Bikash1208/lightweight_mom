"""STEP export of the physical unit-cell geometry used by the MoM solver (Phase 6).

Builds solids **directly from the same `UnitCell`/`Geometry`/`Substrate` objects the solver
consumes** (`geometry.py`, Phase 1) — no parameters are re-entered or re-derived here, so there is
one source of truth for the geometry. This module does not read, assemble or solve `Z I = V`; it
has no effect on any EM result and is not imported by `solver.py` or `sparameters.py`.

Needs the optional `cadquery` dependency (`pip install -r requirements-cad.txt`); it is not part
of the base `requirements.txt` because it pulls in a full OpenCASCADE build and is only needed for
this export, never for solving. Importing this module without it raises a clear `ImportError`.

Coordinate mapping (SI internally, STEP in millimetres — see `units.py`; this module is an I/O
boundary and is the only place besides `units.py` that converts units)
    x, y as in `geometry.py`: origin at the cell centre, lattice-aligned axes.
    z: the solver's stack is air (z > 0) | substrate (−h < z < 0) | zero-thickness PEC pattern at
    z = 0 | air (z < −h) (`GREEN_FUNCTION.md`). The exported **substrate** solid therefore occupies
    z in [−h, 0]. The solver's metal is an idealized zero-thickness surface current; CAD solids
    need a finite thickness to be solid models, so the exported **metal** solid is a thin slab
    occupying z in [0, metal_export_thickness], flush with the substrate's top face and extending
    into the (modelled-as-semi-infinite) air region above it. **`metal_export_thickness` is a CAD
    convenience with no effect on any MoM/Green's-function/S-parameter equation** — it is not
    read anywhere else in the codebase; the default (`DEFAULT_METAL_EXPORT_THICKNESS_M`, 35 µm,
    a standard 1 oz copper thickness) is a plausible physical thickness only so the solid renders
    sensibly, not a claim that the solver models finite-thickness metal. See VALIDATION_EXPORT.md
    Sec. 4 for the exact HFSS reproduction steps this implies (assign a lossless/PerfE boundary;
    do not assign a lossy copper volume unless you are deliberately testing thickness sensitivity,
    which this project has not done).

Feature construction (`geometry.py::Geometry`, unchanged formulas): a `length` (local x) by
`width` (local y) rectangle centred at the origin is rotated **counter-clockwise about the origin**
by `rotation` [rad] (`cq.Workplane.rotate` about the +Z axis matches numpy's/`Geometry.corners()`'s
sign convention — verified in `tests/test_cad_export.py`), then translated to
`(center_x, center_y)`; this reproduces `Geometry.corners()` exactly (rotate-about-centroid then
translate is the same operation whether the translation is applied to the corners or, as here,
folded into the solid's placement, because rotation is applied to the centred shape first).

* `polarity == "metal"` (`rectangle`, `patch`, `dipole`): the feature itself is the metal solid.
* `polarity == "aperture"` (`slot`): a full-period metal plate with the feature subtracted
  (boolean cut), matching `UnitCell.is_metal`'s complement rule exactly.

`export_unit_cell_step` returns the ground-truth solid properties (from the *pre-export* cadquery
shapes) used to build `geometry.json`; `import_and_summarize` independently re-derives the same
quantities from the *STEP file on disk* by re-importing it and analysing the resulting solids'
bounding boxes/volumes — the two must agree for the geometry-consistency tests (Sec. 6 of the
Phase 6 brief) to be meaningful, since they exercise the actual STEP bytes, not just this module's
inputs.
"""

import math
from dataclasses import dataclass

from . import units
from ._validate import finite_real
from .geometry import UnitCell

try:
    import cadquery as cq
except ImportError as _e:  # pragma: no cover - exercised only without the optional dependency
    cq = None
    _IMPORT_ERROR = _e
else:
    _IMPORT_ERROR = None

# 35 um: a standard 1 oz/ft^2 copper thickness. CAD solid convenience only; see module docstring.
DEFAULT_METAL_EXPORT_THICKNESS_M = 35e-6

_PART_SUBSTRATE = "Substrate"
_PART_METAL = "Metal"
_PART_GROUND = "Ground"


def _require_cadquery():
    if cq is None:
        raise ImportError(
            "cadquery is required for STEP export (src/cad_export.py); install it with "
            "'pip install -r requirements-cad.txt'. It is not part of the base solver "
            "dependencies."
        ) from _IMPORT_ERROR


def _check_thickness(metal_export_thickness_m):
    t = finite_real("metal_export_thickness_m", metal_export_thickness_m)
    if t <= 0.0:
        raise ValueError(f"metal_export_thickness_m must be > 0, got {t!r}")
    return t


def _feature_solid(geometry, thickness_mm, z0_mm):
    """The Geometry's rectangle as a solid slab [z0, z0+thickness], its own construction (Sec.
    'Feature construction' above), all dimensions already in millimetres."""
    length_mm, width_mm = units.to_mm(geometry.length), units.to_mm(geometry.width)
    cx_mm, cy_mm = units.to_mm(geometry.center_x), units.to_mm(geometry.center_y)
    return (
        cq.Workplane("XY")
        .workplane(offset=z0_mm)
        .rect(length_mm, width_mm)
        .extrude(thickness_mm)
        .rotate((0, 0, z0_mm), (0, 0, z0_mm + 1), math.degrees(geometry.rotation))
        .translate((cx_mm, cy_mm, 0.0))
    )


def _axis_rect_solid(rect, t_mm, z0_mm):
    return (
        cq.Workplane("XY")
        .workplane(offset=z0_mm)
        .center(units.to_mm(rect.cx), units.to_mm(rect.cy))
        .rect(units.to_mm(rect.length), units.to_mm(rect.width))
        .extrude(t_mm)
    )


def _composite_metal_solid(geo, t_mm):
    """Phase 20 R1: `CompositeGeometry` as one Metal part: union of the metal rectangles minus the apertures. Disjoint pieces
    (e.g. the two arms of a gapped patch) become separate solids inside the one part."""
    metal = _axis_rect_solid(geo.metal[0], t_mm, 0.0)
    for r in geo.metal[1:]:
        metal = metal.union(_axis_rect_solid(r, t_mm, 0.0))
    for r in geo.apertures:
        metal = metal.cut(_axis_rect_solid(r, t_mm, 0.0))
    return metal


def _poly_solid(poly, t_mm, z0_mm):
    pts = [(units.to_mm(x), units.to_mm(y)) for x, y in poly.vertices]
    return cq.Workplane("XY").workplane(offset=z0_mm).polyline(pts).close().extrude(t_mm)


def _poly_composite_metal_solid(geo, t_mm):
    """Phase 20 R3: `PolyComposite` (polygons and axis-aligned rectangles, minus apertures) as one Metal part."""
    from .polygon_geometry import Poly

    def solid(shape):
        return _poly_solid(shape, t_mm, 0.0) if isinstance(shape, Poly) else _axis_rect_solid(shape, t_mm, 0.0)

    metal = solid(geo.metal[0])
    for r in geo.metal[1:]:
        metal = metal.union(solid(r))
    for r in geo.apertures:
        metal = metal.cut(solid(r))
    return metal


def build_load_sheets(unit_cell: UnitCell, metal_export_thickness_m=DEFAULT_METAL_EXPORT_THICKNESS_M):
    """Phase 20 R1: one zero-thickness sheet (a single face) per `unit_cell.loads` entry, in the plane z = t/2 (mid-height of the
    exported metal slab, so the sheet bridges the gap between the side faces of the adjacent metal bodies), covering the load
    rectangle. Named `Load_<k>` (k from 1) in the STEP assembly, so each can be assigned an RLC boundary in HFSS (the current
    direction is `loads[k-1].direction`). Returns a list of cadquery Workplanes holding one Face each."""
    _require_cadquery()
    t_mm = units.to_mm(_check_thickness(metal_export_thickness_m))
    sheets = []
    for ld in unit_cell.loads:
        z = t_mm / 2.0
        if ld.is_rotated:  # Phase 20 R3: a load on a non-axis-aligned conduction line is a rotated rectangle
            corners = [cq.Vector(units.to_mm(x), units.to_mm(y), z) for x, y in ld.corners()]
        else:
            x0, x1, y0, y1 = (units.to_mm(v) for v in ld.region())
            corners = [cq.Vector(x0, y0, z), cq.Vector(x1, y0, z), cq.Vector(x1, y1, z), cq.Vector(x0, y1, z)]
        wire = cq.Wire.makePolygon(corners, close=True)
        sheets.append(cq.Workplane("XY").add(cq.Face.makeFromWires(wire)))
    return sheets


def build_unit_cell_solids(unit_cell: UnitCell, metal_export_thickness_m=DEFAULT_METAL_EXPORT_THICKNESS_M):
    """(substrate_solid, metal_solid): cadquery `Workplane` solids in millimetres, built only from
    `unit_cell` (no independent parameter path — see module docstring). Unchanged by Phase 12A
    (the return type stays a 2-tuple for every existing caller); the optional ground-plane solid
    (``unit_cell.backing == "pec"``) is built separately by ``build_ground_solid``."""
    _require_cadquery()
    t_mm = units.to_mm(_check_thickness(metal_export_thickness_m))
    lat, geo, sub = unit_cell.lattice, unit_cell.geometry, unit_cell.substrate
    px_mm, py_mm, h_mm = units.to_mm(lat.period_x), units.to_mm(lat.period_y), units.to_mm(sub.thickness)

    substrate = cq.Workplane("XY").workplane(offset=-h_mm).rect(px_mm, py_mm).extrude(h_mm)

    if geo.kind == "composite":  # Phase 20 R1: union of metal rectangles minus aperture rectangles (one Metal part)
        metal = _composite_metal_solid(geo, t_mm)
    elif geo.kind == "poly_composite":  # Phase 20 R3: polygons (octagon, 45-degree lines) and rectangles
        metal = _poly_composite_metal_solid(geo, t_mm)
    elif geo.polarity == "metal":
        metal = _feature_solid(geo, t_mm, 0.0)
    else:  # "aperture" (slot): full-period plate minus the feature
        plate = cq.Workplane("XY").rect(px_mm, py_mm).extrude(t_mm)
        metal = plate.cut(_feature_solid(geo, t_mm, 0.0))
    for ld in unit_cell.loads:  # Phase 20 R1: the MoM treats a load region as conductor; in CAD it is the gap the Load_<k> sheet bridges
        if ld.is_rotated:
            pts = [(units.to_mm(x), units.to_mm(y)) for x, y in ld.corners()]
            cutter = cq.Workplane("XY").polyline(pts).close().extrude(t_mm)
        else:
            x0, x1, y0, y1 = (units.to_mm(v) for v in ld.region())
            cutter = cq.Workplane("XY").center(0.5 * (x0 + x1), 0.5 * (y0 + y1)).rect(x1 - x0, y1 - y0).extrude(t_mm)
        metal = metal.cut(cutter)
    return substrate, metal


def build_ground_solid(unit_cell: UnitCell, ground_thickness_m=DEFAULT_METAL_EXPORT_THICKNESS_M):
    """The PEC ground-plane solid at the bottom of the substrate (Phase 12A): a full-period slab
    flush with, and extending downward from, the substrate's bottom face at z = -h, i.e.
    z in [-h - ground_thickness_m, -h] — the mirror of the top metal pattern's CAD-convenience
    slab (module docstring): the solver's ground is an idealized zero-thickness PEC sheet at
    z = -h, and ``ground_thickness_m`` (default: the same nominal copper thickness as
    ``metal_export_thickness_m``, since a solder-mounted/plated ground plane is physically
    similar metal) exists only so the exported solid has nonzero volume; it enters no
    Green's-function/S-parameter equation.

    Returns ``None`` if ``unit_cell.backing != "pec"`` (no ground plane to export).
    """
    if unit_cell.backing != "pec":
        return None
    _require_cadquery()
    t_mm = units.to_mm(_check_thickness(ground_thickness_m))
    lat, sub = unit_cell.lattice, unit_cell.substrate
    px_mm, py_mm, h_mm = units.to_mm(lat.period_x), units.to_mm(lat.period_y), units.to_mm(sub.thickness)
    return cq.Workplane("XY").workplane(offset=-h_mm - t_mm).rect(px_mm, py_mm).extrude(t_mm)


@dataclass(frozen=True)
class SolidProperties:
    """Ground-truth properties of one exported solid, millimetres/mm^3, read directly from the
    cadquery shape (used both as the recorded ground truth and, after a STEP round trip, as the
    independently recovered values for the consistency tests)."""

    n_solids: int
    bbox_mm: tuple  # (xmin, xmax, ymin, ymax, zmin, zmax)
    volume_mm3: float
    centroid_xy_mm: tuple  # bounding-box centre (x, y); equals the rotation centre for a solid
    # rotated about its own centroid, independent of the rotation angle

    @classmethod
    def of(cls, solid):
        n = len(solid.val().Solids())
        bb = solid.val().BoundingBox()
        box = (bb.xmin, bb.xmax, bb.ymin, bb.ymax, bb.zmin, bb.zmax)
        return cls(n, box, solid.val().Volume(), ((bb.xmin + bb.xmax) / 2, (bb.ymin + bb.ymax) / 2))


@dataclass(frozen=True)
class StepExportResult:
    path: str
    metal_export_thickness_m: float
    substrate: SolidProperties
    metal: SolidProperties
    ground: SolidProperties | None = None  # Phase 12A: only for backing="pec"
    load_bboxes_mm: tuple = ()  # Phase 20 R1: (xmin, xmax, ymin, ymax, z) of each Load_<k> sheet


def export_unit_cell_step(unit_cell: UnitCell, path, metal_export_thickness_m=DEFAULT_METAL_EXPORT_THICKNESS_M):
    """Write one STEP file (`Substrate` + `Metal`, plus `Ground` when ``unit_cell.backing=="pec"``
    (Phase 12A), as named parts of one assembly) and return the ground-truth `StepExportResult`
    used to fill `geometry.json`."""
    substrate, metal = build_unit_cell_solids(unit_cell, metal_export_thickness_m)
    ground = build_ground_solid(unit_cell, metal_export_thickness_m)
    asm = cq.Assembly()
    asm.add(substrate, name=_PART_SUBSTRATE, color=cq.Color(0.25, 0.55, 0.25, 1.0))
    asm.add(metal, name=_PART_METAL, color=cq.Color(0.80, 0.55, 0.10, 1.0))
    if ground is not None:
        asm.add(ground, name=_PART_GROUND, color=cq.Color(0.60, 0.60, 0.65, 1.0))
    sheets = build_load_sheets(unit_cell, metal_export_thickness_m) if unit_cell.loads else []
    for k, sheet in enumerate(sheets):
        asm.add(sheet, name=f"Load_{k + 1}", color=cq.Color(0.85, 0.10, 0.10, 1.0))
    asm.export(str(path), exportType="STEP")
    return StepExportResult(
        str(path),
        _check_thickness(metal_export_thickness_m),
        SolidProperties.of(substrate),
        SolidProperties.of(metal),
        SolidProperties.of(ground) if ground is not None else None,
        tuple((b.xmin, b.xmax, b.ymin, b.ymax, b.zmin) for b in (sh.val().BoundingBox() for sh in sheets)),
    )


def build_supercell_solids(supercell, metal_export_thickness_m=DEFAULT_METAL_EXPORT_THICKNESS_M):
    """Phase 12C: STEP export for a heterogeneous supercell (``phase9.supercell_geometry.
    SupercellUnitCell``, or the grounded ``phase12c.grounded_supercell.GroundedSupercellUnitCell``
    wrapper -- duck-typed: only ``.n``, ``.sub_period_x``/``.sub_period_y``, ``.patches``,
    ``.sub_cell_center``, ``.lattice``, ``.substrate`` are read, so neither module is modified).

    Returns ``(substrate_solid, metal_solids)``: one substrate slab sized to the WHOLE supercell
    period, and one metal solid PER SUB-CELL PATCH (a tuple, row-major (i, j) order), each built by
    the same ``_feature_solid`` construction a single ``UnitCell``'s own patch uses, then
    translated to that sub-cell's own physical centre (``supercell.sub_cell_center`` -- the
    patches themselves are centred at (0, 0) by ``supercell_geometry``'s own invariant, so this
    translation, not the patch's own ``center_x``/``center_y``, is what places each one)."""
    _require_cadquery()
    t_mm = units.to_mm(_check_thickness(metal_export_thickness_m))
    lat, sub = supercell.lattice, supercell.substrate
    px_mm, py_mm, h_mm = units.to_mm(lat.period_x), units.to_mm(lat.period_y), units.to_mm(sub.thickness)
    substrate = cq.Workplane("XY").workplane(offset=-h_mm).rect(px_mm, py_mm).extrude(h_mm)

    metals = []
    for i in range(supercell.n):
        for j in range(supercell.n):
            geom = supercell.patches[i][j]
            if geom.polarity != "metal":
                raise ValueError(
                    f"build_supercell_solids: sub-cell ({i}, {j}) has polarity {geom.polarity!r}; "
                    f"only 'metal' (rectangle) sub-cell patches are supported"
                )
            cx_mm, cy_mm = (
                units.to_mm(supercell.sub_cell_center(i, j)[0]),
                units.to_mm(supercell.sub_cell_center(i, j)[1]),
            )
            metals.append(_feature_solid(geom, t_mm, 0.0).translate((cx_mm, cy_mm, 0.0)))
    return substrate, tuple(metals)


def build_supercell_ground_solid(supercell, ground_thickness_m=DEFAULT_METAL_EXPORT_THICKNESS_M):
    """As ``build_ground_solid``, sized to the supercell's WHOLE period (``supercell.lattice``).
    Returns ``None`` if ``supercell.backing != "pec"`` (absence, via ``getattr``, means "air" --
    the same duck-typed default every other Phase 12A consumer uses for a pre-Phase-12A,
    ungrounded ``SupercellUnitCell``)."""
    if getattr(supercell, "backing", "air") != "pec":
        return None
    _require_cadquery()
    t_mm = units.to_mm(_check_thickness(ground_thickness_m))
    lat, sub = supercell.lattice, supercell.substrate
    px_mm, py_mm, h_mm = units.to_mm(lat.period_x), units.to_mm(lat.period_y), units.to_mm(sub.thickness)
    return cq.Workplane("XY").workplane(offset=-h_mm - t_mm).rect(px_mm, py_mm).extrude(t_mm)


@dataclass(frozen=True)
class SupercellStepExportResult:
    path: str
    metal_export_thickness_m: float
    substrate: SolidProperties
    metals: tuple  # one SolidProperties per sub-cell patch, row-major (i, j)
    ground: SolidProperties | None = None


def export_supercell_step(supercell, path, metal_export_thickness_m=DEFAULT_METAL_EXPORT_THICKNESS_M):
    """Write one STEP file for a heterogeneous supercell (``Substrate`` + one ``Metal_<k>`` part
    per sub-cell + optional ``Ground``, as named parts of one assembly)."""
    substrate, metals = build_supercell_solids(supercell, metal_export_thickness_m)
    ground = build_supercell_ground_solid(supercell, metal_export_thickness_m)
    asm = cq.Assembly()
    asm.add(substrate, name=_PART_SUBSTRATE, color=cq.Color(0.25, 0.55, 0.25, 1.0))
    for idx, metal in enumerate(metals):
        asm.add(metal, name=f"{_PART_METAL}_{idx}", color=cq.Color(0.80, 0.55, 0.10, 1.0))
    if ground is not None:
        asm.add(ground, name=_PART_GROUND, color=cq.Color(0.60, 0.60, 0.65, 1.0))
    asm.export(str(path), exportType="STEP")
    return SupercellStepExportResult(
        str(path),
        _check_thickness(metal_export_thickness_m),
        SolidProperties.of(substrate),
        tuple(SolidProperties.of(m) for m in metals),
        SolidProperties.of(ground) if ground is not None else None,
    )


def solid_validity(solid):
    """True iff OpenCASCADE's own topology checker accepts every solid in ``solid`` (the
    'every exported STEP opens as a valid solid model' acceptance criterion, checked with the
    kernel's own validator rather than merely 'did not raise on import')."""
    from OCP.BRepCheck import BRepCheck_Analyzer

    return all(BRepCheck_Analyzer(s.wrapped).IsValid() for s in solid.val().Solids())


def import_and_summarize(path):
    """Re-import ``path`` from disk (independent of anything held in memory) and split its solids
    into substrate, metal and (Phase 12A, if present) ground, by their geometric z-range — the
    same rule the export uses to place them, evaluated on the recovered shape, not on any name or
    import order.

    Two solids (the pre-Phase-12A, ungrounded case, unchanged): substrate (bounding box entirely
    at z <= 0) and metal (entirely at z >= 0). Three solids (a grounded export): the substrate is
    still the one solid whose bounding box has zmax <= 0 AND zmin > the third solid's zmax (i.e.
    it does not extend to the most negative z); metal is still the one with zmin >= 0; ground is
    the remaining solid, entirely below the substrate's bottom face. Raises if the file does not
    contain exactly one solid of each expected kind (2 or 3 total) or any solid fails the
    OpenCASCADE validity check.
    """
    _require_cadquery()
    imported = cq.importers.importStep(str(path))
    if not solid_validity(imported):
        raise ValueError(f"{path}: STEP file contains an invalid solid (BRepCheck_Analyzer)")
    solids = imported.val().Solids()

    def props(s):
        bb = s.BoundingBox()
        box = (bb.xmin, bb.xmax, bb.ymin, bb.ymax, bb.zmin, bb.zmax)
        return SolidProperties(1, box, s.Volume(), ((bb.xmin + bb.xmax) / 2, (bb.ymin + bb.ymax) / 2))

    metal = [s for s in solids if s.BoundingBox().zmin >= -1e-9]
    below_z0 = [s for s in solids if s.BoundingBox().zmax <= 1e-9]
    if len(solids) == 3 and len(metal) == 1 and len(below_z0) == 2:
        # The two z<=0 solids are substrate and ground; ground is the one entirely below the
        # other's bottom face (its own zmax <= the other's zmin, with a tolerance for touching).
        a, b = below_z0
        if a.BoundingBox().zmax <= b.BoundingBox().zmin + 1e-9:
            ground, substrate = a, b
        elif b.BoundingBox().zmax <= a.BoundingBox().zmin + 1e-9:
            ground, substrate = b, a
        else:
            raise ValueError(f"{path}: the two below-z=0 solids do not stack (substrate, ground)")
        return StepExportResult(
            str(path),
            float("nan"),
            props(substrate),
            props(metal[0]),
            props(ground),
        )
    if len(solids) == 2 and len(metal) == 1 and len(below_z0) == 1:
        return StepExportResult(str(path), float("nan"), props(below_z0[0]), props(metal[0]))
    raise ValueError(
        f"{path}: expected exactly one substrate (z<=0) and one metal (z>=0) solid, optionally "
        f"plus one ground solid entirely below the substrate (Phase 12A), got {len(solids)} "
        f"solid(s) ({len(below_z0)} at/below z=0, {len(metal)} at/above)"
    )


def expected_feature_bbox_mm(geometry):
    """Closed-form (xmin, xmax, ymin, ymax) of a `length` x `width` rectangle rotated about its
    own centroid by `rotation` and placed at (center_x, center_y): for a rectangle the axis-aligned
    bounding box after rotation by theta has half-extents
    (L|cos theta| + W|sin theta|)/2 and (L|sin theta| + W|cos theta|)/2 — used by the tests as an
    independent (non-cadquery) check of the STEP-recovered metal bounding box."""
    length_mm, width_mm = units.to_mm(geometry.length), units.to_mm(geometry.width)
    cx_mm, cy_mm = units.to_mm(geometry.center_x), units.to_mm(geometry.center_y)
    c, s = abs(math.cos(geometry.rotation)), abs(math.sin(geometry.rotation))
    hx, hy = (length_mm * c + width_mm * s) / 2.0, (length_mm * s + width_mm * c) / 2.0
    return (cx_mm - hx, cx_mm + hx, cy_mm - hy, cy_mm + hy)


def import_named_parts(path):
    """Phase 20 R1: re-import a STEP assembly from disk with its part names; returns
    `{name: {"n_solids", "n_faces", "bbox_mm": (xmin, xmax, ymin, ymax, zmin, zmax), "volume_mm3", "area_mm2"}}`
    (e.g. `Substrate`, `Metal`, `Ground`, `Load_1` ...). Independent of anything held in memory."""
    _require_cadquery()
    asm = cq.Assembly.importStep(str(path))
    out = {}
    for child in asm.children:
        shape = child.obj.val() if hasattr(child.obj, "val") else child.obj
        bb = shape.BoundingBox()
        out[child.name] = {
            "n_solids": len(shape.Solids()),
            "n_faces": len(shape.Faces()),
            "bbox_mm": (bb.xmin, bb.xmax, bb.ymin, bb.ymax, bb.zmin, bb.zmax),
            "volume_mm3": float(shape.Volume()),
            "area_mm2": float(shape.Area()),
        }
    return out
