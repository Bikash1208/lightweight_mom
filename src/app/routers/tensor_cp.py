"""Phase 27 item 3: "Tensor CP (monopole-fed)" design mode.

`POST /api/tensor_cp/design` submits a background job (`hologram.tensor_cp_service.run_design`); poll `/api/jobs/{id}`.
`GET /api/tensor_cp/defaults` the default request per board; `GET /api/tensor_cp/families` the cell-family registry
(domains of validity, validation labels, known biases). A request outside the family's domain finishes with the named
state `NOT_AVAILABLE_FROM_CURRENT_DATA` and the reasons. The 1-bit HMSA flow (`/api/hmsa`) is untouched. Model-only."""

from typing import Literal

from fastapi import APIRouter, HTTPException, Request
from pydantic import BaseModel

from ...export.tensor_cp_export import export_app_result
from ...hologram import cell_family_registry as registry
from ...hologram import tensor_cp_service as service
from ..jobs import config_hash
from ..storage import job_output_dir

router = APIRouter(prefix="/api/tensor_cp", tags=["tensor_cp"])


class SubstrateBody(BaseModel):
    eps_r: float = 2.2
    tan_d: float = 0.0009
    h_mm: float = 1.575


class TensorCPRequest(BaseModel):
    board: Literal["single", "hex4_static", "hex4_combined"] = "single"
    beams_deg: "list[list[float]] | None" = None  # [[theta, phi], ...]: one per port (hex4_static) else one
    handedness: Literal["rhcp", "lhcp"] = "rhcp"
    pin_height_mm: float = 6.0  # pin top above the patterned surface
    frequency_ghz: float = 10.54
    period_mm: float = 6.4
    substrate: SubstrateBody = SubstrateBody()
    quick: bool = False  # coarser pattern grid (about 4x faster); not for reporting


def _run_job(job, body):
    out = job_output_dir(job.id)

    def progress(frac, msg):
        job.progress, job.progress_message = float(frac), msg

    result = service.run_design(body, out_dir=out, progress=progress, cancelled=lambda: job.cancelled)
    job.progress, job.progress_message = 1.0, "done"
    return result


class TensorCPExportRequest(BaseModel):
    job_id: str  # a finished /design job
    drive_port: int = 0  # four-hexagon boards: 0 = HFSS default sources, 1..4 = one port driven, the others 50 ohm
    with_step: bool = False  # single hexagon only: also write board.step (CadQuery; tens of MB)


def _export_job(job, source_result, body):
    job.progress_message = "writing the board files"
    names = export_app_result(
        source_result, job_output_dir(body["job_id"]), drive_port=body["drive_port"], with_step=body["with_step"]
    )
    job.progress = 1.0
    return {
        "job_id": body["job_id"],
        "files": names,
        "download": {n: f"/api/files/{body['job_id']}/{n}" for n in names},
    }


@router.post("/export")
def export(body: TensorCPExportRequest, request: Request):
    if body.drive_port not in (0, 1, 2, 3, 4):
        raise HTTPException(400, "drive_port must be 0..4")
    source = request.app.state.jobs.get(body.job_id)
    if source is None or source.kind != "tensor_cp_design":
        raise HTTPException(404, f"no finished tensor_cp_design job {body.job_id!r}")
    if source.status.value != "done" or (source.result or {}).get("state") != "OK":
        raise HTTPException(409, "the design job is not finished with state OK")
    return request.app.state.jobs.submit("tensor_cp_export", _export_job, source.result, body.model_dump()).to_dict()


@router.get("/defaults")
def defaults(board: str = "single"):
    if board not in service.BOARDS:
        raise HTTPException(400, f"board must be one of {service.BOARDS}")
    return service.default_request(board)


@router.get("/families")
def families():
    return {name: f.to_dict() for name, f in registry.list_families().items()}


@router.post("/design")
def design(body: TensorCPRequest, request: Request):
    payload = body.model_dump()
    if payload["beams_deg"] is None:
        payload.pop("beams_deg")
    try:
        norm = service.normalise_request(payload)
    except ValueError as exc:
        raise HTTPException(400, str(exc)) from exc
    job = request.app.state.jobs.submit(
        "tensor_cp_design", _run_job, norm, cache_key=config_hash({"kind": "tensor_cp_design", **norm})
    )
    return job.to_dict()
