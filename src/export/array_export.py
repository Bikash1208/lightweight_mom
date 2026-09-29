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
"""


@dataclass(frozen=True)
class ArrayExportResult:
    step_path: "str | None"
    layout_csv_path: str
    layout_json_path: str
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
    readme_path = write_hfss_readme(output_dir)

    step_path, file_size, export_time = None, None, None
    if include_step:
        from .. import cad_export

        supercell = build_grounded_supercell(n, period_x_m, period_y_m, substrate, metal, lw_grid_m)
        step_path = output_dir / "array.step"
        t0 = time.perf_counter()
        try:
            kwargs = {} if metal_export_thickness_m is None else {"metal_export_thickness_m": metal_export_thickness_m}
            cad_export.export_supercell_step(supercell, step_path, **kwargs)
        except ImportError as exc:
            raise cad_unavailable(exc, message="STEP export needs the optional 'cad' extra (cadquery)") from exc
        export_time = time.perf_counter() - t0
        file_size = step_path.stat().st_size
        step_path = str(step_path)

    return ArrayExportResult(
        step_path=step_path,
        layout_csv_path=str(csv_path),
        layout_json_path=str(json_path),
        readme_path=str(readme_path),
        n_cells_per_side=n,
        total_cells=n * n,
        file_size_bytes=file_size,
        export_time_s=export_time,
    )
