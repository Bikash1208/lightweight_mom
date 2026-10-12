"""Phase 23 C: "Reduce sidelobes" for a coded 1-bit multi-port design.

`POST /api/sidelobes/optimize` takes a finished `/synthesize_ports` job and returns a NEW job result of
the same shape (so the existing board export, pattern and package endpoints take it as their source),
with the optimised cell states, the chosen port phases and a before / after table measured on a
held-out grid the optimiser never saw (`hologram.sll_board`). Model-only."""

from typing import Literal

from fastapi import APIRouter, HTTPException, Request
from pydantic import BaseModel

from ...errors import config_invalid
from ...hologram.sll_board import (
    HeldOutGrid,
    apply_levels,
    build_problem,
    cells_changed,
    evaluate,
    optimise,
)
from ...hologram.sll_optimizer import OptimiserConfig
from ..jobs import JobStatus
from ..storage import job_output_dir

router = APIRouter(prefix="/api/sidelobes", tags=["sidelobes"])


class SidelobeBeam(BaseModel):
    theta_deg: float
    phi_deg: float


class SidelobeRequest(BaseModel):
    source_job_id: str
    frequency_hz: float
    period_mm: float
    eps_r: float
    tan_delta: float
    thickness_mm: float
    port_beams: "list[list[SidelobeBeam]]"  # the requested beams of each port, port order
    polarization: str = "rhcp"
    mode: Literal["single", "all", "both"] = "single"
    max_gain_loss_db: float = 0.5
    max_flips: int = 600  # iteration budget (accepted single-cell changes)
    optimise_port_phases: bool = True  # all/both modes
    seed: int = 0
    frequency_check_percent: float = 2.0


def _round(x):
    return None if x is None else round(float(x), 2)


def _table(ev):
    """Per scenario, per beam: the numbers the user reads, rounded."""
    out = {}
    for name, sc in ev.items():
        out[name] = {
            "peak_directivity_dbi": _round(sc["peak_directivity_dbi"]),
            "level_spread_db": _round(sc["level_spread_db"]),
            "beams": [
                {
                    "port": b["port_id"] + 1,
                    "requested": [b["requested_theta_deg"], b["requested_phi_deg"]],
                    "achieved": [_round(b["theta_deg"]), _round(b["phi_deg"])],
                    "pointing_error_deg": _round(b["pointing_error_deg"]),
                    "directivity_dbi": _round(b["directivity_dbi"]),
                    "sll_global_db": _round(b["sll_global_db"]),
                    "sll_xcut_db": _round(b.get("sll_xcut_db")),
                    "sll_ycut_db": _round(b.get("sll_ycut_db")),
                    "xpd_db": _round(b["xpd_db"]),
                }
                for b in sc["beams"]
            ],
        }
    return out


def _run_optimize_job(job, source_result, body):
    n_ports = len(source_result["ports"])
    if len(body["port_beams"]) != n_ports:
        raise config_invalid(message=f"{len(body['port_beams'])} port beam list(s) for a {n_ports}-port design")
    if source_result.get("formulation") != "coded" or source_result.get("n_bits") != 1:
        raise config_invalid(message="Reduce sidelobes needs a coded 1-bit multi-port design")
    job.progress_message = "Building the cell-wise pattern model..."
    port_beams = [[(b["theta_deg"], b["phi_deg"]) for b in pb] for pb in body["port_beams"]]
    prob = build_problem(
        source_result,
        body["frequency_hz"],
        body["eps_r"],
        body["tan_delta"],
        body["thickness_mm"],
        body["period_mm"],
        port_beams,
        polarization=body["polarization"],
    )
    mode = body["mode"]
    if mode in ("all", "both") and n_ports < 2:
        raise config_invalid(message="the all-ports mode needs 2 or more ports")
    grid = HeldOutGrid()
    job.progress_message = "Measuring the starting design on the held-out grid..."
    before = evaluate(prob, grid=grid)
    cfg = OptimiserConfig(
        max_gain_loss_db=body["max_gain_loss_db"],
        max_flips=body["max_flips"],
        optimise_alpha=bool(body["optimise_port_phases"]) and mode in ("all", "both"),
        seed=body["seed"],
    )
    job.progress_message = f"Optimising cell states ({mode}, up to {cfg.max_flips} changes)..."
    levels, alpha, runs = optimise(prob, mode, cfg)
    job.progress_message = "Measuring the optimised design on the held-out grid..."
    after = evaluate(prob, levels, alpha, grid=grid)
    freq_checks = {}
    pct = body["frequency_check_percent"] / 100.0
    if pct > 0:
        for sign in (-1, 1):
            f = body["frequency_hz"] * (1 + sign * pct)
            freq_checks[f"{sign * body['frequency_check_percent']:+g}%"] = _table(
                evaluate(prob, levels, alpha, frequency_hz=f, grid=grid)
            )
    new = apply_levels(source_result, prob, levels, alpha, evaluation=after)
    new["sidelobe_optimisation"].update(
        {
            "mode": mode,
            "max_gain_loss_db": cfg.max_gain_loss_db,
            "iteration_budget": cfg.max_flips,
            "held_out_grid": {"theta_step_deg": grid.theta_step_deg, "phi_step_deg": grid.phi_step_deg},
            "before": _table(before),
            "after": _table(after),
            "frequency_check": freq_checks,
            "frequency_check_note": "same cell states re-evaluated at the shifted frequency (reference wave and "
            "far field moved, the two cells' own phases kept at the design frequency); not re-optimised",
            "cells_changed": cells_changed(prob, levels),
            "label": "MODEL_ONLY_SIDELOBE_OPTIMISATION",
        }
    )
    del runs
    try:  # a fresh layout picture of the optimised board (best effort)
        from ...export.exact_board import exact_board_from_ports_result
        from .export import _write_board_png

        board = exact_board_from_ports_result(new, body["period_mm"], body["thickness_mm"])
        path = _write_board_png(board, job_output_dir(job.id) / "patch_layout.png")
        if path:
            new["patch_layout_png_url"] = f"/api/files/{job.id}/patch_layout.png"
    except Exception:
        pass
    return new


@router.post("/optimize")
def optimize(body: SidelobeRequest, request: Request):
    source = request.app.state.jobs.get(body.source_job_id)
    if source is None:
        raise HTTPException(404, f"job {body.source_job_id!r} not found")
    if source.status != JobStatus.DONE:
        raise HTTPException(409, f"source job {body.source_job_id!r} is {source.status.value}, not done")
    result = source.result or {}
    if not result.get("ports") or "patches_mm" not in result["ports"][0]:
        raise HTTPException(422, f"source job {body.source_job_id!r} is not a multi-port design")
    job = request.app.state.jobs.submit("sidelobe_optimize", _run_optimize_job, result, body.model_dump())
    return job.to_dict()
