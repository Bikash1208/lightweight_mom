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
from ..plots import pattern_cut_png, pattern_uv_map_png
from ..storage import job_output_dir

router = APIRouter(prefix="/api/pattern", tags=["pattern"])


class BeamDirection(BaseModel):
    theta_deg: float
    phi_deg: float


class PatternRequest(BaseModel):
    source_job_id: str
    frequency_hz: float
    theta_max_deg: float = 90.0
    n_theta: int = 91
    n_phi: int = 73
    beams: "list[BeamDirection]" = []  # 13C Step 3: markers on the u-v map + one cut per beam


def _run_pattern_job(job, source_result, frequency_hz, theta_max_deg, n_theta, n_phi, beams):
    x_m = source_result["x_m"]
    y_m = source_result["y_m"]
    aperture_field = source_result["aperture_field_real"] + 1j * source_result["aperture_field_imag"]
    period_m = source_result["period_mm"] * 1e-3
    result = compute_pattern(x_m, y_m, aperture_field, frequency_hz, period_m, period_m, theta_max_deg, n_theta, n_phi)

    beams_theta_phi = [(b["theta_deg"], b["phi_deg"]) for b in beams]
    uv_png_url, cut_png_urls = None, []
    try:
        out_dir = job_output_dir(job.id)
        uv_path = pattern_uv_map_png(
            result.theta_deg, result.phi_deg, result.power_db, out_dir / "pattern_uv.png", beams_theta_phi
        )
        uv_png_url = f"/api/files/{job.id}/{uv_path.name}"
        cut_phis = [b[1] for b in beams_theta_phi] or [0.0]
        for i, phi_cut in enumerate(cut_phis):
            cut_path = pattern_cut_png(
                result.theta_deg, result.phi_deg, result.power_db, phi_cut, out_dir / f"pattern_cut_{i}.png"
            )
            cut_png_urls.append(f"/api/files/{job.id}/{cut_path.name}")
    except Exception:
        pass  # best-effort plots only; never fail a real pattern result over a plotting error

    return {
        "result": result,
        "labels": (LOCAL_PERIODICITY_APPROXIMATION,),
        "pattern_uv_png_url": uv_png_url,
        "pattern_cut_png_urls": cut_png_urls,
    }


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
        body.model_dump()["beams"],
    )
    return job.to_dict()
