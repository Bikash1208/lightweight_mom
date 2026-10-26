"""Phase 27 item 3: "Tensor CP (monopole-fed)" design mode.

`POST /api/tensor_cp/design` submits a background job (`hologram.tensor_cp_service.run_design`); poll `/api/jobs/{id}`.
`GET /api/tensor_cp/defaults` the default request per board; `GET /api/tensor_cp/families` the cell-family registry
(domains of validity, validation labels, known biases). A request outside the family's domain finishes with the named
state `NOT_AVAILABLE_FROM_CURRENT_DATA` and the reasons. The 1-bit HMSA flow (`/api/hmsa`) is untouched. Model-only."""

import json
import uuid
from types import SimpleNamespace
from typing import Literal

from fastapi import APIRouter, File, Form, HTTPException, Request, UploadFile
from pydantic import BaseModel

from ...export.tensor_cp_export import export_app_result
from ...hologram import cell_family_registry as registry
from ...hologram import tensor_cp_compare as compare_mod
from ...hologram import tensor_cp_design as design_mod
from ...hologram import tensor_cp_service as service
from ..jobs import config_hash
from ..jsonutil import to_jsonable
from ..storage import job_output_dir, upload_dir

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


async def _save(updir, upload):
    if upload is None or not upload.filename:
        return None
    path = updir / upload.filename
    path.write_bytes(await upload.read())
    return str(path)


def _pred_for_job(request, job_id):
    """Predictions in the Phase 26 `predictions.json` format for a finished single-hexagon design job."""
    job = request.app.state.jobs.get(job_id)
    res = None if job is None else job.result
    if not res or res.get("state") != "OK" or res["request"]["board"] != "single":
        raise HTTPException(404, f"no finished single-hexagon tensor_cp_design job {job_id!r}")
    req = res["request"]
    hand = +1 if req["handedness"] == "rhcp" else -1
    beam = tuple(req["beams_deg"][0])
    sel, _ = design_mod.single_selected(
        design_mod.SupercellTable.from_files(design_mod.TABLE_GLOB), beam=beam, hand=hand, quick=req.get("quick", False)
    )
    return design_mod.single_hex_predictions(sel, beam=beam, hand=hand)


@router.post("/compare")
async def compare(
    request: Request,
    board: Literal["single", "hex4"] = Form("single"),
    pin: str = Form("scripted"),
    design_job_id: "str | None" = Form(None),
    predictions_json: "UploadFile | None" = File(None),
    gain_csv: "UploadFile | None" = File(None),
    params_txt: "UploadFile | None" = File(None),
    s_param: "UploadFile | None" = File(None),
    s_port: int = Form(1),
    ar_csv: "UploadFile | None" = File(None),
    ez_csv: "UploadFile | None" = File(None),
    gain_ports: str = Form(""),
    gain_files: "list[UploadFile] | None" = File(None),
    params_ports: str = Form(""),
    params_files: "list[UploadFile] | None" = File(None),
):
    """Upload the HFSS exports of a tensor-CP board and get the Phase 26 comparison (named states, power budget). Synchronous:
    the comparison reads files, it solves nothing. `board=single`: gain_csv / params_txt / s_param / ar_csv / ez_csv.
    `board=hex4`: s_param (.s4p), and per driven port a gain CSV (`gain_files` + `gain_ports`, e.g. "3,2") and optional
    params files (`params_files` + `params_ports`)."""
    updir = upload_dir(str(uuid.uuid4()))
    pred_path = await _save(updir, predictions_json)
    if pred_path is None:
        if design_job_id and board == "single":
            pred_path = str(updir / "predictions.json")
            (updir / "predictions.json").write_text(json.dumps(_pred_for_job(request, design_job_id), default=float), encoding="utf-8", newline="\n")
        else:
            pred_path = str(compare_mod.PRED if board == "single" else compare_mod.PRED4)
    try:
        if board == "single":
            args = SimpleNamespace(
                gain=await _save(updir, gain_csv),
                s11=await _save(updir, s_param),
                port=s_port,
                params=await _save(updir, params_txt),
                ar=await _save(updir, ar_csv),
                ez=await _save(updir, ez_csv),
                pin=pin,
                pred=pred_path,
            )
            if not any([args.gain, args.params, args.s11, args.ar, args.ez]):
                raise HTTPException(400, "upload at least one HFSS export")
            rep = compare_mod.compare(args)
            text = compare_mod.markdown(rep) if args.gain else compare_mod.markdown_params_only(rep)
        else:
            gains = [p.strip() for p in gain_ports.split(",") if p.strip()]
            params = [p.strip() for p in params_ports.split(",") if p.strip()]
            if len(gains) != len(gain_files or []) or len(params) != len(params_files or []):
                raise HTTPException(400, "gain_ports / params_ports must list one port per uploaded file, in order")
            g_items = [f"{p}={await _save(updir, f)}" for p, f in zip(gains, gain_files or [], strict=True)]
            p_items = [f"{p}={await _save(updir, f)}" for p, f in zip(params, params_files or [], strict=True)]
            args = SimpleNamespace(s4p=await _save(updir, s_param), gain=g_items, params=p_items, pred=pred_path)
            if not (args.s4p or g_items):
                raise HTTPException(400, "upload a .s4p and/or at least one port's gain CSV")
            rep = compare_mod.compare_fourhex(args)
            text = compare_mod.markdown_fourhex(rep)
    except HTTPException:
        raise
    except (ValueError, KeyError, OSError) as exc:
        raise HTTPException(400, f"could not read the uploads: {exc}") from exc
    return to_jsonable({"report": rep, "markdown": text, "predictions_used": pred_path})
