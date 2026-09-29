"""Phase 13 Part 3, page 2 ("Library"): browse a registered unit-cell library (the read-only Phase
7 benchmark, or one this app built) and run a new L(/W) x frequency sweep for an arbitrary base
config. Reuses `library_sweep.py` (itself `runner.run_sweep`, Phase 5, unmodified) plus
`library_classifier`/`phase_extraction` (Phase 7, unmodified) for the resonance/coverage/
phase-usability analysis -- no physics lives in this router.
"""

from fastapi import APIRouter, HTTPException, Request
from pydantic import BaseModel

from ...config import ConfigError, config_from_dict
from ...errors import config_invalid
from ...library_sweep import (
    SweepSpec,
    classify_cell_backing_aware,
    extract_design_frequency,
    run_custom_sweep,
    write_app_library_json,
)
from ...phase_extraction import phase_coverage, phase_vs_geometry_table
from ...unit_cell_library import UnitCellLibrary
from .. import library_registry
from ..storage import job_output_dir

router = APIRouter(prefix="/api/library", tags=["library"])


@router.get("")
def list_libraries():
    return library_registry.list_libraries()


@router.get("/{library_id}")
def get_library(library_id: str, frequency_ghz: "float | None" = None):
    try:
        path = library_registry.resolve_path(library_id)
    except KeyError as exc:
        raise HTTPException(404, str(exc)) from exc
    lib = UnitCellLibrary.from_json(path)
    cells = lib.all_cells()
    body = {
        "library_id": library_id,
        "n_cells": len(cells),
        "cells": [
            {
                "cell_id": c["cell_id"],
                "L_mm": c["L_mm"],
                "W_mm": c["W_mm"],
                "geometry_status": c["geometry_status"],
                "solver_status": c["solver_status"],
            }
            for c in cells
        ],
    }
    if frequency_ghz is not None:
        table = []
        for c in cells:
            if c["solver_status"] != "OK":
                continue
            try:
                entry = lib.response(c["L_mm"], c["W_mm"], frequency_ghz, family="S11")
            except ValueError:
                continue  # this cell has no record at the requested frequency
            table.append(
                {
                    "L_mm": c["L_mm"],
                    "W_mm": c["W_mm"],
                    "phase_deg": entry["phase_deg"],
                    "magnitude_db": entry["mag_db"],
                }
            )
        body["frequency_ghz"] = frequency_ghz
        body["phase_vs_geometry"] = table
        body["coverage"] = phase_coverage(table) if table else None
    return body


class LibrarySweepRequest(BaseModel):
    library_id: str
    base_config: dict  # config.yaml shape; geometry.length/width are placeholders, overwritten per cell
    l_values_mm: "list[float]"
    w_values_mm: "list[float]"
    f_min_hz: float
    f_max_hz: float
    f_step_hz: float
    f_design_hz: float


def _run_sweep_job(job, library_id, base_config_dict, spec_kwargs):
    try:
        sim_config = config_from_dict(base_config_dict)
    except ConfigError as exc:
        raise config_invalid(exc, message=str(exc)) from exc
    spec = SweepSpec(**spec_kwargs)

    def _progress(event):
        job.progress = event["done"] / max(event["total"], 1)
        job.progress_message = f"{event['done']}/{event['total']} solved"

    out_dir = job_output_dir(job.id)
    result = run_custom_sweep(sim_config, spec, out_dir=out_dir, progress=_progress)
    path = write_app_library_json(result.database, out_dir / "unit_cell_library.json")

    unit_cell = sim_config.unit_cell
    library_registry.register(
        library_id,
        path,
        meta={
            "period_mm": unit_cell.lattice.period_x * 1e3,
            "backing": unit_cell.backing,
            "substrate": unit_cell.substrate.name,
            "n_geometries": len(result.database.cells),
            "f_min_hz": spec.f_min_hz,
            "f_max_hz": spec.f_max_hz,
        },
    )

    design_rows = extract_design_frequency(result.database, spec)
    table = phase_vs_geometry_table(design_rows)
    coverage = phase_coverage(table) if table else None
    by_id = {c.cell_id: c for c in result.database.cells}
    classifications = [
        classify_cell_backing_aware(by_id[row.cell_id], result.database, row, unit_cell) for row in design_rows
    ]

    return {
        "library_id": library_id,
        "download_url": f"/api/files/{job.id}/{path.name}",
        "n_geometries": len(result.database.cells),
        "n_ok": sum(1 for c in result.database.cells if c.status == "ok"),
        "table": table,
        "coverage": coverage,
        "classifications": classifications,
    }


@router.post("/sweep")
def sweep(body: LibrarySweepRequest, request: Request):
    spec_kwargs = dict(
        f_min_hz=body.f_min_hz,
        f_max_hz=body.f_max_hz,
        f_step_hz=body.f_step_hz,
        f_design_hz=body.f_design_hz,
        l_values_mm=tuple(body.l_values_mm),
        w_values_mm=tuple(body.w_values_mm),
    )
    job = request.app.state.jobs.submit("library_sweep", _run_sweep_job, body.library_id, body.base_config, spec_kwargs)
    return job.to_dict()
