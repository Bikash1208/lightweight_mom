"""Phase 13 Part 3, page 4 ("Pattern prediction"): far field from an aperture, reusing
`pattern.aperture_pattern.compute_pattern` (itself built on `phase8e.beam_reconstruction.far_field_pattern`,
unmodified). Chains off a prior reflectarray (or, later, HMSA) design job's in-memory result
(`app/jobs.py`'s `Job.result` -- raw numpy, not the JSON-flattened HTTP response) rather than
re-uploading a potentially large aperture-field array through the browser.
"""

from fastapi import APIRouter, HTTPException, Request
from pydantic import BaseModel

from ...labels import LOCAL_PERIODICITY_APPROXIMATION
from ...pattern.aperture_pattern import compute_pattern
from ..jobs import JobStatus

router = APIRouter(prefix="/api/pattern", tags=["pattern"])


class PatternRequest(BaseModel):
    source_job_id: str
    frequency_hz: float
    theta_max_deg: float = 90.0
    n_theta: int = 91
    n_phi: int = 73


def _run_pattern_job(job, source_result, frequency_hz, theta_max_deg, n_theta, n_phi):
    x_m = source_result["x_m"]
    y_m = source_result["y_m"]
    aperture_field = source_result["aperture_field_real"] + 1j * source_result["aperture_field_imag"]
    period_m = source_result["period_mm"] * 1e-3
    result = compute_pattern(x_m, y_m, aperture_field, frequency_hz, period_m, period_m, theta_max_deg, n_theta, n_phi)
    return {"result": result, "labels": (LOCAL_PERIODICITY_APPROXIMATION,)}


@router.post("/predict")
def predict(body: PatternRequest, request: Request):
    source = request.app.state.jobs.get(body.source_job_id)
    if source is None:
        raise HTTPException(404, f"job {body.source_job_id!r} not found")
    if source.status != JobStatus.DONE:
        raise HTTPException(409, f"source job {body.source_job_id!r} is {source.status.value}, not done")
    if "aperture_field_real" not in (source.result or {}):
        raise HTTPException(422, f"source job {body.source_job_id!r} did not produce an aperture field")

    job = request.app.state.jobs.submit(
        "pattern_predict",
        _run_pattern_job,
        source.result,
        body.frequency_hz,
        body.theta_max_deg,
        body.n_theta,
        body.n_phi,
    )
    return job.to_dict()
