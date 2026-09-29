"""Phase 13 Part 3: the one router every page shares -- job status polling and cancellation. Every
other router only ever creates jobs (via ``request.app.state.jobs.submit``); the frontend polls
this router's GET endpoint until the job is ``done``/``failed``/``cancelled``.
"""

from fastapi import APIRouter, HTTPException, Request

router = APIRouter(prefix="/api/jobs", tags=["jobs"])


@router.get("/{job_id}")
def get_job(job_id: str, request: Request):
    job = request.app.state.jobs.get(job_id)
    if job is None:
        raise HTTPException(404, f"job {job_id!r} not found")
    return job.to_dict(include_result=True)


@router.post("/{job_id}/cancel")
def cancel_job(job_id: str, request: Request):
    job = request.app.state.jobs.cancel(job_id)
    if job is None:
        raise HTTPException(404, f"job {job_id!r} not found")
    return job.to_dict()
