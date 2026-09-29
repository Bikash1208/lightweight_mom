"""Phase 7: one STEP file per successfully-solved geometry (req. 16-17).

Reuses `cad_export.py` (Phase 6, unmodified) end to end: the `UnitCell` exported here is rebuilt
by `database.Database.unit_cell` from the *same* solved database (lattice, substrate, geometry),
not re-entered by hand, so "the geometry must represent the same physical unit cell used by the
solver" (req. 16) holds by construction, not by a separate cross-check. The STEP export itself
does not alter the EM model (`cad_export.py`'s own docstring: the exported metal solid's thickness
is a CAD-only convenience with zero effect on any solved quantity).

The unit cell it produces already contains what req. 16 asks for: `Substrate` (dielectric layer),
`Metal` (the PEC patch), both solids referenced to one coordinate system (cell centre at the
origin, `geometry.py`), which is exactly `periodic_cell_boundary` (the substrate/metal solids'
own footprint *is* the unit cell boundary -- there is no separate boundary object to export,
since a periodic unit cell's boundary is implicit in the solids' extent). Nothing about
"periodic lattice pairing" itself is representable in a STEP file (STEP has no periodic-boundary
-condition concept); this limitation is stated here and in `LIBRARY_REPORT.md`, not hidden.
"""

from dataclasses import dataclass
from pathlib import Path

from .cad_export import DEFAULT_METAL_EXPORT_THICKNESS_M, export_unit_cell_step, expected_feature_bbox_mm
from . import units


def step_filename(l_mm, w_mm):
    """Deterministic name, req. 16: `UC_L1.0_W1.0.step` .. `UC_L4.4_W4.4.step`. `.1f` matches
    the locked sweep's one-decimal-mm values exactly (`library_schema.L_VALUES_MM`)."""
    return f"UC_L{l_mm:.1f}_W{w_mm:.1f}.step"


@dataclass(frozen=True)
class GeometryBoundsCheck:
    """Req. 17's explicit numeric verification, independent of (not merely trusting) the fact
    that `UnitCell.__post_init__` already raises if a feature extends outside its cell -- every
    inequality req. 17 lists is computed and checked here directly."""

    period_x_mm: float
    period_y_mm: float
    cell_x_min_mm: float
    cell_x_max_mm: float
    cell_y_min_mm: float
    cell_y_max_mm: float
    patch_x_min_mm: float
    patch_x_max_mm: float
    patch_y_min_mm: float
    patch_y_max_mm: float
    l_in_range: bool
    w_in_range: bool
    within_cell: bool

    @property
    def ok(self):
        return self.l_in_range and self.w_in_range and self.within_cell


def check_geometry_bounds(unit_cell):
    geo, lat = unit_cell.geometry, unit_cell.lattice
    px_mm, py_mm = units.to_mm(lat.period_x), units.to_mm(lat.period_y)
    l_mm, w_mm = units.to_mm(geo.length), units.to_mm(geo.width)
    cell_x_min, cell_x_max = -px_mm / 2.0, px_mm / 2.0
    cell_y_min, cell_y_max = -py_mm / 2.0, py_mm / 2.0
    patch_x_min, patch_x_max, patch_y_min, patch_y_max = expected_feature_bbox_mm(geo)

    return GeometryBoundsCheck(
        period_x_mm=px_mm, period_y_mm=py_mm,
        cell_x_min_mm=cell_x_min, cell_x_max_mm=cell_x_max,
        cell_y_min_mm=cell_y_min, cell_y_max_mm=cell_y_max,
        patch_x_min_mm=patch_x_min, patch_x_max_mm=patch_x_max,
        patch_y_min_mm=patch_y_min, patch_y_max_mm=patch_y_max,
        l_in_range=bool(0.0 < l_mm < px_mm),
        w_in_range=bool(0.0 < w_mm < py_mm),
        within_cell=bool(
            patch_x_min >= cell_x_min - 1e-9 and patch_x_max <= cell_x_max + 1e-9
            and patch_y_min >= cell_y_min - 1e-9 and patch_y_max <= cell_y_max + 1e-9
        ),
    )


@dataclass(frozen=True)
class CellGeometryExport:
    cell_id: str
    l_mm: float
    w_mm: float
    step_path: Path
    bounds: GeometryBoundsCheck
    step_result: object          # cad_export.StepExportResult


def export_library_geometries(database, config, out_dir, *, progress=None):
    """One STEP file per `status == 'ok'` cell in ``database`` into ``out_dir``
    (`UC_L{L}_W{W}.step`); returns `{cell_id: CellGeometryExport}`. Every export is bounds-checked
    (req. 17); a failing check raises rather than writing a misleading geometry (req. 16's "report
    the limitation rather than creating a misleading geometry" -- here the limitation would be
    this module's own bug, not the exporter's, since the solver already refused any such
    geometry at solve time, so this is a belt-and-braces invariant, not an expected path)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    results = {}
    ok_cells = [c for c in database.cells if c.status == "ok"]
    for i, cell in enumerate(ok_cells):
        l_mm = round(units.to_mm(cell.parameters["length_m"]), 6)
        w_mm = round(units.to_mm(cell.parameters["width_m"]), 6)
        unit_cell = database.unit_cell(cell.index)
        bounds = check_geometry_bounds(unit_cell)
        if not bounds.ok:
            raise ValueError(
                f"cell {cell.cell_id} (L={l_mm} mm, W={w_mm} mm) fails its own geometry bounds "
                f"check: {bounds}; the solver should have already rejected this geometry"
            )
        path = out_dir / step_filename(l_mm, w_mm)
        step_result = export_unit_cell_step(unit_cell, path, DEFAULT_METAL_EXPORT_THICKNESS_M)
        results[cell.cell_id] = CellGeometryExport(cell.cell_id, l_mm, w_mm, path, bounds, step_result)
        if progress is not None:
            progress({"done": i + 1, "total": len(ok_cells), "cell_id": cell.cell_id})
    return results
