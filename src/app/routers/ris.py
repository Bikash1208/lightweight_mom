"""Phase 20 R5-5: "RIS designer" (plane-wave reflecting RIS, CP-2 family). `POST /api/ris/design` submits a short background job (`ris.service.run_design`); poll `/api/jobs/{id}`.
`GET /api/ris/defaults`, `GET /api/ris/coupling` (array-level results merged from the server jobs), `GET /api/ris/families`. Requests outside the registry domain end in the named state
`NOT_AVAILABLE_FROM_CURRENT_DATA`. Model-only; no heavy compute in the app (supercell checks are server jobs)."""

from typing import Literal

from fastapi import APIRouter, HTTPException, Request
from pydantic import BaseModel

from ...hologram import cell_family_registry as registry
from ...ris import service

router = APIRouter(prefix="/api/ris", tags=["ris"])


class RisRequest(BaseModel):
    frequency_ghz: float = 10.5
    rows: int = 20
    cols: int = 20
    theta_s_deg: float = 30.0
    phi_s_deg: float = 0.0
    theta_i_deg: float = 0.0
    phi_i_deg: float = 0.0
    handedness: Literal["rhcp", "lhcp"] = "rhcp"
    variant: Literal["plain", "dither", "tilt", "refined", "auto"] = "plain"
    seed: int = 0
    sensitivity: bool = True
    library: Literal["v1", "v2"] = "v2"


def _run(job, body):
    def progress(frac, msg):
        job.progress, job.progress_message = float(frac), msg

    return service.run_design(body, progress=progress)


@router.get("/defaults")
def defaults():
    return service.default_request()


@router.get("/coupling")
def coupling():
    return service.coupling_status()


@router.get("/families")
def families():
    return {service.FAMILY: registry.get(service.FAMILY).to_dict()}


@router.post("/design")
def design(body: RisRequest, request: Request):
    try:
        norm = service.normalise_request(body.model_dump())
    except ValueError as exc:
        raise HTTPException(400, str(exc)) from exc
    return request.app.state.jobs.submit("ris_design", _run, norm).to_dict()
