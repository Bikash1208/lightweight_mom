"""Phase 13 Part 3: downloads a file a job produced (STEP export, layout CSV/JSON, an app-built
library, a diagnostic bundle) -- everything under `storage.OUTPUT_ROOT`, one directory per job."""

from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse

from ..storage import output_root

router = APIRouter(prefix="/api/files", tags=["files"])


@router.get("/{job_id}/{filename}")
def download_file(job_id: str, filename: str):
    root = output_root().resolve()
    path = (root / job_id / filename).resolve()
    if root not in path.parents or not path.is_file():
        raise HTTPException(404, "file not found")
    return FileResponse(path, filename=filename)
