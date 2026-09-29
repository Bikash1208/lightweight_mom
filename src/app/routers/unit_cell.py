"""Phase 13 Part 3, page 1 ("Unit cell"): thin HTTP glue over `diagnostics.instrumented_compute_s_parameters`
(a normal fixed-mesh solve) and `phase13.convergence.richardson_extrapolate_s11xx` (the app's
mesh-convergence-aware default, PHASE13_WEBAPP_PROMPT.md Part 0) -- no physics lives here. The
request body's ``config`` is the exact same mapping `config.config_from_dict` already validates
(the same shape as `configs/config.yaml`), so there is one schema and one validator, not two.
"""

from fastapi import APIRouter, Request
from pydantic import BaseModel

from ...basis_functions import MeshSettings
from ...cad_export import DEFAULT_METAL_EXPORT_THICKNESS_M, export_unit_cell_step
from ...config import ConfigError, config_from_dict
from ...diagnostics import instrumented_compute_s_parameters
from ...errors import cad_unavailable, config_invalid
from ...labels import ZERO_THICKNESS_METAL_MODEL
from ...phase13.convergence import FAST_CELL_COUNTS, richardson_extrapolate_s11xx
from ..jobs import config_hash
from ..storage import job_output_dir

router = APIRouter(prefix="/api/unit_cell", tags=["unit_cell"])

MESH_CONVERGED = "MESH_CONVERGED"
MESH_NOT_CONVERGED = "MESH_NOT_CONVERGED"


class UnitCellSolveRequest(BaseModel):
    config: dict
    frequencies_hz: "list[float] | None" = None
    mesh_mode: str = "fast"  # "fast" (use config.mesh as given) | "richardson" (Part 0 extrapolation)


def _build_sim_config(config_dict):
    try:
        return config_from_dict(config_dict)
    except ConfigError as exc:
        raise config_invalid(exc, message=str(exc)) from exc


def _solve_one_point(sim_config, frequency_hz, mesh_mode, run_id):
    labels = [ZERO_THICKNESS_METAL_MODEL]
    warnings = []
    convergence = None

    if mesh_mode == "richardson":
        convergence = richardson_extrapolate_s11xx(
            sim_config.unit_cell, frequency_hz, sim_config.reference_planes, cell_counts=FAST_CELL_COUNTS
        )
        labels.append(MESH_CONVERGED if convergence.mesh_converged else MESH_NOT_CONVERGED)
        finest = max(FAST_CELL_COUNTS)
        mesh = MeshSettings(cells_x=finest, cells_y=finest, floquet_order=finest)
    else:
        labels.append(MESH_NOT_CONVERGED)  # fixed mesh, never extrapolated -- honest by default
        mesh = sim_config.mesh

    result, diagnostics = instrumented_compute_s_parameters(
        sim_config.unit_cell, frequency_hz, mesh, sim_config.reference_planes, run_id=run_id
    )
    if diagnostics.grating_orders_propagating:
        warnings.append("a Floquet order other than (0,0) propagates at this frequency (grating lobe)")
    if diagnostics.surface_wave_margin is not None and diagnostics.surface_wave_margin.get("any_near_band"):
        warnings.append("frequency is close to a surface-wave onset for this stack")

    return {
        "frequency_hz": frequency_hz,
        "record": result.to_record(),
        "diagnostics": diagnostics,
        "convergence": convergence,
        "labels": labels,
        "warnings": warnings,
    }


def _run_solve_job(job, config_dict, frequencies_hz, mesh_mode):
    sim_config = _build_sim_config(config_dict)
    freqs = frequencies_hz or [sim_config.excitation.frequency]
    points = []
    for i, f in enumerate(freqs):
        if job.cancelled:
            break
        job.progress = i / len(freqs)
        job.progress_message = f"solving {i + 1}/{len(freqs)} ({f / 1e9:.3f} GHz)"
        points.append(_solve_one_point(sim_config, f, mesh_mode, run_id=job.id))
    job.progress = 1.0
    return {"mode": mesh_mode, "grounded": sim_config.unit_cell.backing == "pec", "points": points}


@router.post("/solve")
def solve(body: UnitCellSolveRequest, request: Request):
    key = config_hash({"kind": "unit_cell_solve", **body.model_dump()})
    job = request.app.state.jobs.submit(
        "unit_cell_solve",
        _run_solve_job,
        body.config,
        body.frequencies_hz,
        body.mesh_mode,
        cache_key=key,
    )
    return job.to_dict()


def _run_export_step_job(job, config_dict, metal_export_thickness_m):
    sim_config = _build_sim_config(config_dict)
    try:
        out_dir = job_output_dir(job.id)
        path = out_dir / "unit_cell.step"
        result = export_unit_cell_step(sim_config.unit_cell, path, metal_export_thickness_m)
    except ImportError as exc:
        raise cad_unavailable(exc) from exc
    return {
        "download_url": f"/api/files/{job.id}/{path.name}",
        "properties": result,
    }


class UnitCellExportStepRequest(BaseModel):
    config: dict
    metal_export_thickness_m: float = DEFAULT_METAL_EXPORT_THICKNESS_M


@router.post("/export_step")
def export_step(body: UnitCellExportStepRequest, request: Request):
    job = request.app.state.jobs.submit(
        "unit_cell_export_step", _run_export_step_job, body.config, body.metal_export_thickness_m
    )
    return job.to_dict()
