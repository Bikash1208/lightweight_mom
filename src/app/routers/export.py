"""Phase 13 Part 3, page 6 ("Export"): full-array STEP + layout CSV/JSON + HFSS-ready README,
chaining off a completed reflectarray (or HMSA) design job's cell assignments. Thin glue over
`export.array_export.export_array` (itself `phase12c.grounded_supercell.build_grounded_supercell`
+ `cad_export.export_supercell_step`, unmodified) -- grounded-array-only, the same limitation
`export_array`'s own module docstring states (square N x N grid).
"""

from fastapi import APIRouter, HTTPException, Request
from pydantic import BaseModel

from ...errors import config_invalid
from ...export.array_export import export_array
from ...materials import Metal, Substrate
from ..jobs import JobStatus
from ..storage import job_output_dir

router = APIRouter(prefix="/api/export", tags=["export"])


class ExportRequest(BaseModel):
    source_job_id: str
    period_mm: float
    eps_r: float
    tan_delta: float
    thickness_mm: float
    include_step: bool = True


def _lw_grid_from_assignments(assignments, nx, ny):
    grid = [[None] * nx for _ in range(ny)]
    for a in assignments:
        if a.chosen_l_mm is None or a.chosen_w_mm is None:
            raise config_invalid(
                message=f"cell ({a.i},{a.j}) has no library match -- cannot export an incomplete array"
            )
        grid[a.i][a.j] = (a.chosen_l_mm * 1e-3, a.chosen_w_mm * 1e-3)
    return grid


def _run_export_job(job, source_result, params):
    assignments = source_result["assignments"]
    ny, nx = source_result["x_m"].shape
    lw_grid_m = _lw_grid_from_assignments(assignments, nx, ny)
    substrate = Substrate(
        name="user", eps_r=params["eps_r"], tan_delta=params["tan_delta"], thickness=params["thickness_mm"] * 1e-3
    )
    metal = Metal()
    period_m = params["period_mm"] * 1e-3
    out_dir = job_output_dir(job.id)
    result = export_array(lw_grid_m, period_m, period_m, substrate, metal, out_dir, include_step=params["include_step"])

    def _download(path):
        if path is None:
            return None
        from pathlib import Path

        return f"/api/files/{job.id}/{Path(path).name}"

    return {
        "n_cells": result.n_cells,
        "file_size_bytes": result.file_size_bytes,
        "export_time_s": result.export_time_s,
        "step_download_url": _download(result.step_path),
        "layout_csv_download_url": _download(result.layout_csv_path),
        "layout_json_download_url": _download(result.layout_json_path),
        "readme_download_url": _download(result.readme_path),
    }


@router.post("/array")
def export_reflectarray(body: ExportRequest, request: Request):
    source = request.app.state.jobs.get(body.source_job_id)
    if source is None:
        raise HTTPException(404, f"job {body.source_job_id!r} not found")
    if source.status != JobStatus.DONE:
        raise HTTPException(409, f"source job {body.source_job_id!r} is {source.status.value}, not done")
    if "assignments" not in (source.result or {}):
        raise HTTPException(422, f"source job {body.source_job_id!r} did not produce cell assignments")

    job = request.app.state.jobs.submit("array_export", _run_export_job, source.result, body.model_dump())
    return job.to_dict()
