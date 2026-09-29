"""Phase 13 Part 3, page 5 ("Verify"): the full-supercell MoM check against the lookup-table
prediction, with the cost-estimate/refusal gate `verify.supercell_check` already provides
(reused unmodified) -- refuses sizes known to be infeasible with the reason, shown to the user
before they can force it.
"""

from fastapi import APIRouter, HTTPException, Request
from pydantic import BaseModel

from ...verify.supercell_check import estimate_cost, run_verification
from ..jsonutil import to_jsonable

router = APIRouter(prefix="/api/verify", tags=["verify"])


@router.get("/cost_estimate")
def cost_estimate(n_cells: int):
    return to_jsonable(estimate_cost(n_cells))


class VerifyRunRequest(BaseModel):
    n_cells: int
    mesh_cells_per_subcell: "int | None" = None
    force: bool = False


def _run_verification_job(job, n_cells, mesh_cells_per_subcell, force):
    kwargs = {} if mesh_cells_per_subcell is None else {"mesh_cells_per_subcell": mesh_cells_per_subcell}
    return run_verification(n_cells, force=force, **kwargs)


@router.post("/run")
def run(body: VerifyRunRequest, request: Request):
    cost = estimate_cost(body.n_cells)
    if cost.refused and not body.force:
        raise HTTPException(409, cost.refusal_reason)
    job = request.app.state.jobs.submit(
        "verify_run", _run_verification_job, body.n_cells, body.mesh_cells_per_subcell, body.force
    )
    return job.to_dict()
