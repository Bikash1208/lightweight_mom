"""Phase 13 Part 3, page 7 ("HFSS compare"): upload a Touchstone (.sNp), pick the matching case,
overlay against MoM and feed the on-disk comparison ledger. Thin glue over
`hfss_compare.service` (itself `phase9.hfss_reference.read_touchstone` +
`phase12d.hfss_ingest.compare_single_cell`, unmodified). Synchronous (no job queue): comparing
against an already-solved MoM reference CSV is fast, not a solve.
"""

import uuid

from fastapi import APIRouter, File, Form, HTTPException, UploadFile

from ...hfss_compare.service import (
    compare_uploaded_file,
    load_ledger,
    record_comparison,
    validate_touchstone_upload,
)
from ...phase9.hfss_reference import TouchstoneError
from ...phase12d.hfss_ingest import ALL_SINGLE_CELL_IDS, CASE_WITH_HFSS_DATA
from ..jsonutil import to_jsonable
from ..storage import upload_dir

router = APIRouter(prefix="/api/hfss_compare", tags=["hfss_compare"])


@router.get("/cases")
def cases():
    return {"case_ids": ALL_SINGLE_CELL_IDS, "case_with_real_hfss_data_on_disk": CASE_WITH_HFSS_DATA}


@router.get("/ledger")
def ledger():
    return load_ledger()


@router.post("/upload")
async def upload(
    touchstone_file: UploadFile = File(...),
    case_id: str = Form(...),
    top_port_name: str = Form("FloquetPort1"),
    bottom_port_name: str = Form("FloquetPort2"),
):
    updir = upload_dir(str(uuid.uuid4()))
    path = updir / touchstone_file.filename
    path.write_bytes(await touchstone_file.read())

    try:
        validate_touchstone_upload(path)
    except TouchstoneError as exc:
        raise HTTPException(400, str(exc)) from exc

    try:
        summary, comparison = compare_uploaded_file(
            case_id, path, top_port_name=top_port_name, bottom_port_name=bottom_port_name
        )
    except (ValueError, TouchstoneError) as exc:
        raise HTTPException(400, str(exc)) from exc

    ledger_path = record_comparison(summary)

    n = len(comparison.frequency_ghz)
    stride = max(1, n // 40)  # a manageable overlay table, not all 401 points
    overlay = [
        {
            "frequency_ghz": float(comparison.frequency_ghz[i]),
            "mom_mag": float(comparison.mom_mag[i]),
            "mom_phase_deg": float(comparison.mom_phase_deg[i]),
            "hfss_mode1_mag": float(comparison.hfss_mode1_mag[i]),
            "hfss_mode1_phase_deg": float(comparison.hfss_mode1_phase_deg[i]),
            "dphase_mode1_deg": float(comparison.dphase_mode1_deg[i]),
        }
        for i in range(0, n, stride)
    ]

    return to_jsonable({"summary": summary, "overlay": overlay, "ledger_path": str(ledger_path)})
