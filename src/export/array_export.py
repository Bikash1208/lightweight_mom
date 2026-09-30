"""Full-array STEP export, layout CSV/JSON, and an HFSS-ready README.

**STEP export** reuses `phase12c.grounded_supercell.build_grounded_supercell` and
`cad_export.export_supercell_step` unmodified. Both require a SQUARE N x N grid (a limitation of
`phase9.supercell_geometry.build_supercell`, which this phase does not touch, per "physics/
numerics untouched"): a non-square design must be padded or cropped by the caller first, and this
module says so explicitly rather than silently truncating.
"""

import csv
import json
from dataclasses import dataclass
from pathlib import Path

from ..errors import cad_unavailable
from ..phase12c.grounded_supercell import build_grounded_supercell

_HFSS_README_TEMPLATE = """\
# HFSS setup notes for this export (lessons from Phase 12, PHASE12D_ADDENDUM.md)

- Floquet port >= 18 modes, air-box top at 5 mm above the pattern plane, no de-embed (matches this
  project's own MoM reference plane).
- PEC ground and PEC pattern reproduce the MoM reference to ~0.3 deg mean phase off-resonance
  (G1); a 2-mode-port / copper-ground setup was off by up to 4.7 deg -- do not use it.
- Patch thickness: this export is zero-thickness in the solver's own model
  (`ZERO_THICKNESS_METAL_MODEL`); a 35 micron drawn thickness in the CAD file changes the resonant
  phase by ~3.7 deg near resonance (measured, G4) -- expect a corresponding HFSS/MoM offset if you
  keep the drawn thickness.
- Mesh convergence at resonance is first-order in cell size on this project's own MoM (unresolved
  edge-current singularity); do not assume a coarse-mesh MoM comparison point is converged near
  resonance (see `docs/phases/PHASE13_WEBAPP.md`).
- 10-15 GHz, 401 points reproduces this project's own validated frequency sweep.

## Checking the surface-wave dispersion table (Phase 14) against HFSS

If you want a spot check of the eigenmode method's auto-built dispersion table (`docs/phases/
PHASE14_EIGENMODE.md`) against a real HFSS eigenmode solve, for one representative patch size:

- Master/slave periodic boundaries on the unit cell's four side walls, with a linear phase shift
  (not a Floquet port -- this is a source-free eigenmode solve, not an excited one) swept to find
  the frequency where a bound mode exists at the design frequency, or equivalently sweep the phase
  shift at the fixed design frequency and locate the eigenfrequency crossing it.
- Eigenmode solver (not the driven-modal/terminal solver the S-parameter exports above use).
- Air region above the patch tall enough that the bound mode's evanescent tail (decay length
  ~1/alpha, alpha = sqrt(beta^2 - k0^2)) has decayed well before the outer boundary -- a PML or an
  absorbing boundary a few decay lengths up is typically needed; a plain radiation boundary too
  close will perturb the found beta.
- Record beta (or the effective index n_eff = beta/k0) at the same patch size and frequency this
  project's own table used, and compare via the CSV upload path (Step 2 Advanced, "Upload an HFSS
  eigenmode table") or `dispersion.hfss_table.read_dispersion_csv` directly.
"""


@dataclass(frozen=True)
class ArrayExportResult:
    step_path: "str | None"
    layout_csv_path: str
    layout_json_path: str
    layout_dxf_path: str
    readme_path: str
    n_cells_per_side: int  # PHASE13B A8: was `n_cells` but held only N (the side length) -- an
    # N x N design's summary reported "N cells" instead of the real N*N, misleadingly small.
    total_cells: int
    file_size_bytes: "int | None"
    export_time_s: "float | None"


def _validate_square(lw_grid_m):
    n = len(lw_grid_m)
    if any(len(row) != n for row in lw_grid_m):
        raise ValueError(
            "export_array requires a square N x N grid (build_grounded_supercell's own "
            "limitation); pad or crop the design to a square aperture first"
        )
    return n


def write_layout_files(lw_grid_m, output_dir, basename="layout"):
    """Layout CSV (`i, j, x_m, y_m, length_m, width_m`) and the same data as JSON -- always
    written, independent of whether STEP export (needs the optional `cad` extra) succeeds."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    n = len(lw_grid_m)
    rows = []
    for i, row in enumerate(lw_grid_m):
        for j, (length_m, width_m) in enumerate(row):
            rows.append({"i": i, "j": j, "length_m": length_m, "width_m": width_m})

    csv_path = output_dir / f"{basename}.csv"
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["i", "j", "length_m", "width_m"])
        writer.writeheader()
        writer.writerows(rows)

    json_path = output_dir / f"{basename}.json"
    json_path.write_text(json.dumps({"n_cells_per_side": n, "total_cells": n * n, "cells": rows}, indent=2))
    return csv_path, json_path


def write_hfss_readme(output_dir, basename="HFSS_README"):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"{basename}.md"
    path.write_text(_HFSS_README_TEMPLATE)
    return path


def _rectangle_lines_mm(cx_mm, cy_mm, length_mm, width_mm):
    hl, hw = length_mm / 2.0, width_mm / 2.0
    corners = [
        (cx_mm - hl, cy_mm - hw),
        (cx_mm + hl, cy_mm - hw),
        (cx_mm + hl, cy_mm + hw),
        (cx_mm - hl, cy_mm + hw),
    ]
    return [(corners[k], corners[(k + 1) % 4]) for k in range(4)]


def write_layout_dxf(lw_grid_m, period_x_m, period_y_m, output_dir, basename="layout"):
    """PHASE13D item 4: a lightweight DXF of the patch layer (millimetres, layer ``PATCHES``) -- a
    PCB fab wants this or a Gerber, not a 100 MB STEP (PHASE13D_POLISH_AND_PHYSICS.md). Hand-written
    DXF R12 (``LINE`` entities only -- valid in every DXF reader since R12, unlike ``LWPOLYLINE``
    which needs R14+): no new dependency (`ezdxf` was never added -- see this module's own history).
    Same cell-centre convention as `phase9.supercell_geometry.SupercellUnitCell.sub_cell_center`
    (row i = y-index, column j = x-index), so the drawing lines up with the STEP/CSV export of the
    same design.
    """
    n = _validate_square(lw_grid_m)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    px_mm, py_mm = period_x_m * 1e3, period_y_m * 1e3
    sub_px_mm, sub_py_mm = px_mm / n, py_mm / n

    lines = []
    for i, row in enumerate(lw_grid_m):
        for j, (length_m, width_m) in enumerate(row):
            cx_mm = -px_mm / 2.0 + (j + 0.5) * sub_px_mm
            cy_mm = -py_mm / 2.0 + (i + 0.5) * sub_py_mm
            lines.extend(_rectangle_lines_mm(cx_mm, cy_mm, length_m * 1e3, width_m * 1e3))

    parts = ["0", "SECTION", "2", "HEADER", "9", "$ACADVER", "1", "AC1009", "0", "ENDSEC"]
    parts += ["0", "SECTION", "2", "ENTITIES"]
    for (x1, y1), (x2, y2) in lines:
        parts += [
            "0", "LINE", "8", "PATCHES",
            "10", f"{x1:.6f}", "20", f"{y1:.6f}", "30", "0.0",
            "11", f"{x2:.6f}", "21", f"{y2:.6f}", "31", "0.0",
        ]
    parts += ["0", "ENDSEC", "0", "EOF"]

    path = output_dir / f"{basename}.dxf"
    path.write_text("\n".join(parts) + "\n")
    return path


def write_hfss_step(lw_grid_m, period_x_m, period_y_m, thickness_m, path, decimals=3, size_step_mm=None):
    """HFSS-friendly full-array STEP: every coordinate/size rounded in millimetres (no
    3.9999999999 values from m->mm float conversion), patches and ground as zero-thickness sheets
    (the solver's own model; no 35 um solids that force tiny mesh elements), all patches in ONE
    body ("Patches") so a single Perfect E assignment covers them. Optional `size_step_mm`
    snaps patch sizes to a fabrication grid (fewer distinct geometries)."""
    import cadquery as cq

    def r(v):
        return round(float(v), decimals) + 0.0

    n_y = len(lw_grid_m)
    n_x = len(lw_grid_m[0])
    px, py = r(period_x_m * 1e3), r(period_y_m * 1e3)
    h = r(thickness_m * 1e3)
    lx, ly = r(n_x * px), r(n_y * py)
    faces = []
    for i, row in enumerate(lw_grid_m):
        for j, (l_m, w_m) in enumerate(row):
            l_mm, w_mm = l_m * 1e3, w_m * 1e3
            if size_step_mm:
                l_mm = round(l_mm / size_step_mm) * size_step_mm
                w_mm = round(w_mm / size_step_mm) * size_step_mm
            l_mm, w_mm = r(l_mm), r(w_mm)
            if l_mm <= 0 or w_mm <= 0:
                continue
            cx = r(-lx / 2 + (j + 0.5) * px)
            cy = r(-ly / 2 + (i + 0.5) * py)
            x0, x1 = r(cx - l_mm / 2), r(cx + l_mm / 2)
            y0, y1 = r(cy - w_mm / 2), r(cy + w_mm / 2)
            wire = cq.Wire.makePolygon(
                [cq.Vector(x0, y0, 0), cq.Vector(x1, y0, 0), cq.Vector(x1, y1, 0), cq.Vector(x0, y1, 0)], close=True
            )
            faces.append(cq.Face.makeFromWires(wire))
    patches = cq.Compound.makeCompound(faces)
    substrate = cq.Solid.makeBox(lx, ly, h, pnt=cq.Vector(r(-lx / 2), r(-ly / 2), -h))
    gw = cq.Wire.makePolygon(
        [cq.Vector(r(-lx / 2), r(-ly / 2), -h), cq.Vector(r(lx / 2), r(-ly / 2), -h),
         cq.Vector(r(lx / 2), r(ly / 2), -h), cq.Vector(r(-lx / 2), r(ly / 2), -h)], close=True
    )
    ground = cq.Face.makeFromWires(gw)
    asm = cq.Assembly()
    asm.add(substrate, name="Substrate", color=cq.Color(0.85, 0.85, 0.86, 1.0))
    asm.add(patches, name="Patches", color=cq.Color(0.91, 0.54, 0.23, 1.0))
    asm.add(ground, name="Ground", color=cq.Color(0.60, 0.60, 0.65, 1.0))
    asm.export(str(path), exportType="STEP")
    return {"n_patches": len(faces), "board_mm": (lx, ly, h)}


def export_array(
    lw_grid_m, period_x_m, period_y_m, substrate, metal, output_dir, metal_export_thickness_m=None, include_step=True
):
    """Writes layout CSV/JSON, an HFSS-ready README, and (if `include_step` and the `cad` extra is
    installed) a full-array STEP file. Raises `errors.LwmomError` (`CAD_UNAVAILABLE`) only if
    `include_step=True` and cadquery is missing -- layout/README are always written first, so a
    partial result is never silently lost."""
    import time

    n = _validate_square(lw_grid_m)
    output_dir = Path(output_dir)
    csv_path, json_path = write_layout_files(lw_grid_m, output_dir)
    dxf_path = write_layout_dxf(lw_grid_m, period_x_m, period_y_m, output_dir)
    readme_path = write_hfss_readme(output_dir)

    step_path, file_size, export_time = None, None, None
    if include_step:
        from .. import cad_export

        step_path = output_dir / "array.step"
        t0 = time.perf_counter()
        try:
            # HFSS-clean geometry: mm values rounded to 1 um (no 3.99999999 float artefacts),
            # zero-thickness patch/ground sheets (the solver's own model), one "Patches" body.
            write_hfss_step(lw_grid_m, period_x_m, period_y_m, substrate.thickness, step_path)
        except ImportError as exc:
            raise cad_unavailable(exc, message="STEP export needs the optional 'cad' extra (cadquery)") from exc
        export_time = time.perf_counter() - t0
        file_size = step_path.stat().st_size
        step_path = str(step_path)

    return ArrayExportResult(
        step_path=step_path,
        layout_csv_path=str(csv_path),
        layout_json_path=str(json_path),
        layout_dxf_path=str(dxf_path),
        readme_path=str(readme_path),
        n_cells_per_side=n,
        total_cells=n * n,
        file_size_bytes=file_size,
        export_time_s=export_time,
    )
