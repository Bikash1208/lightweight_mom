"""Phase 13 Part 3, page 4 ("Pattern prediction"): far field from an aperture, reusing
`pattern.aperture_pattern.compute_pattern` (itself built on `phase8e.beam_reconstruction.far_field_pattern`,
unmodified). Chains off a prior reflectarray (or, later, HMSA) design job's in-memory result
(`app/jobs.py`'s `Job.result` -- raw numpy, not the JSON-flattened HTTP response) rather than
re-uploading a potentially large aperture-field array through the browser.

PHASE15 Part D: a second, VECTOR pattern model (`pattern.vector_aperture_pattern`, Part B) is now
offered alongside the original scalar one. `PatternRequest.pattern_model` defaults to `"scalar"` --
the pre-Phase-15 behaviour, so any caller that omits the field (including every pre-Phase-15 test)
is completely unaffected; the app's own JS explicitly requests `"vector"` by default (a UI
convenience default, not a backend one -- the same split Part A's `Beam.desired_polarization`
already uses, for the same bit-identical-reproduction reason).
"""

import numpy as np
from fastapi import APIRouter, HTTPException, Request
from pydantic import BaseModel

from ...errors import config_invalid
from ...green_function import LayeredGreenFunction
from ...labels import LOCAL_PERIODICITY_APPROXIMATION, SCALAR_PATTERN_MODEL_LEGACY, VECTOR_PATTERN_MODEL
from ...materials import Substrate
from ...pattern.aperture_pattern import compute_pattern
from ...pattern.vector_aperture_pattern import (
    aperture_current_from_scalar,
    compute_vector_pattern,
    copol_component,
    measure_vortex_winding,
)
from ..jobs import JobStatus
from ..plots import pattern_cut_png, pattern_uv_map_png, vortex_phase_ring_png
from ..storage import job_output_dir

router = APIRouter(prefix="/api/pattern", tags=["pattern"])


class BeamDirection(BaseModel):
    theta_deg: float
    phi_deg: float
    # PHASE15 Part D: only meaningful for pattern_model="vector" -- which polarization/OAM charge
    # this beam was DESIGNED for (Part A), so the per-beam breakdown can report achieved vs
    # requested. Defaults match Beam's own (natural/0), so an old scalar-only request is unaffected.
    desired_polarization: str = "natural"
    oam_charge: int = 0


class PatternRequest(BaseModel):
    source_job_id: str
    frequency_hz: float
    theta_max_deg: float = 90.0
    n_theta: int = 91
    n_phi: int = 73
    beams: "list[BeamDirection]" = []  # 13C Step 3: markers on the u-v map + one cut per beam
    pattern_model: str = "scalar"  # "scalar" (default, pre-Phase-15 behaviour) | "vector"
    eps_r: "float | None" = None  # required for pattern_model="vector"
    tan_delta: "float | None" = None
    thickness_mm: "float | None" = None
    leakage_alpha_per_m: float = 0.0  # PHASE15 Part B: off (0) by default
    # PHASE15 Part C gap-fill: a multi-port synthesis job's result has no top-level x_m/aperture_field
    # of its own (`"ports"` is a LIST of per-port results instead) -- `port_index` selects which
    # port's own current this pattern call predicts from. `None` (default) is the ordinary
    # single-aperture path, completely unaffected.
    port_index: "int | None" = None


def _port_source_result(source_result, port_index):
    if port_index is None:
        return source_result
    ports = source_result.get("ports")
    if not ports:
        raise config_invalid(message="source job has no per-port results (not a multi-port synthesis job)")
    if not (0 <= port_index < len(ports)):
        raise config_invalid(message=f"port_index {port_index} out of range (0..{len(ports) - 1})")
    port = ports[port_index]
    return {
        "x_m": port["x_m"], "y_m": port["y_m"],
        "aperture_field_real": port["aperture_field_real"], "aperture_field_imag": port["aperture_field_imag"],
        "period_mm": source_result["period_mm"],
        "feeds_mm": [(port["feed_x_mm"], port["feed_y_mm"])],
    }


def _run_scalar_pattern_job(job, source_result, frequency_hz, theta_max_deg, n_theta, n_phi, beams):
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
        "labels": (LOCAL_PERIODICITY_APPROXIMATION, SCALAR_PATTERN_MODEL_LEGACY),
        "pattern_uv_png_url": uv_png_url,
        "pattern_cut_png_urls": cut_png_urls,
        "per_beam": [],
    }


def _nearest_grid_index(theta_grid_deg, phi_grid_deg, theta_deg, phi_deg):
    i = int(np.argmin(np.abs(theta_grid_deg[:, 0] - theta_deg)))
    j = int(np.argmin(np.abs((phi_grid_deg[0, :] - phi_deg + 180.0) % 360.0 - 180.0)))
    return i, j


_CROSS_LABEL = {"rhcp": "lhcp", "lhcp": "rhcp", "hp": "vp", "vp": "hp", "slant_p45": "slant_m45", "slant_m45": "slant_p45"}


def _per_beam_breakdown(job, x_m, y_m, jx, jy, frequency_hz, green, dx_m, dy_m, beams, pat, out_dir):
    """PHASE15 Part D: "Step 3 gains a per-beam result table (direction achieved, polarization
    achieved, XPD or axial ratio, 'Pencil check' or 'Vortex l=...')" -- one entry per REQUESTED
    beam, measured (not assumed) against the actual computed vector pattern."""
    peak_power = float(np.max(np.abs(pat.e_theta) ** 2 + np.abs(pat.e_phi) ** 2))
    entries = []
    for idx, b in enumerate(beams):
        pol = b["desired_polarization"]
        measure_label = pol if pol in _CROSS_LABEL else "rhcp"  # "natural" has no single co-pol
        i, j = _nearest_grid_index(pat.theta_deg, pat.phi_deg, b["theta_deg"], b["phi_deg"])
        e_theta_b, e_phi_b = pat.e_theta[i, j], pat.e_phi[i, j]

        xpd_db = None
        if pol in _CROSS_LABEL:
            co = copol_component(np.array([e_theta_b]), np.array([e_phi_b]), pol)[0]
            cross = copol_component(np.array([e_theta_b]), np.array([e_phi_b]), _CROSS_LABEL[pol])[0]
            if abs(cross) > 0.0:
                xpd_db = float(10.0 * np.log10(np.abs(co) ** 2 / np.abs(cross) ** 2))

        winding = measure_vortex_winding(
            x_m, y_m, jx, jy, frequency_hz, green, dx_m, dy_m,
            b["theta_deg"], b["phi_deg"], measure_label, pattern_peak_power=peak_power,
        )
        beam_type_requested = "Pencil" if b["oam_charge"] == 0 else f"Vortex l={b['oam_charge']}"
        looks_like_vortex = abs(winding.winding_number) > 0.5
        note = (
            "Vortex expected — Natural polarization keeps the monopole's radial currents."
            if pol == "natural" and looks_like_vortex
            else ("Vortex measured at this direction." if looks_like_vortex else "Pencil beam measured at this direction.")
        )

        phase_png_url = None
        try:
            path = vortex_phase_ring_png(
                winding.ring_alpha_deg, winding.ring_phase_deg, winding.winding_number,
                out_dir / f"vortex_phase_{idx}.png",
            )
            phase_png_url = f"/api/files/{job.id}/{path.name}"
        except Exception:
            pass

        entries.append(
            {
                "theta_deg": b["theta_deg"],
                "phi_deg": b["phi_deg"],
                "desired_polarization": pol,
                "beam_type_requested": beam_type_requested,
                "xpd_db": xpd_db,
                "null_depth_db": winding.null_depth_db,
                "winding_number": winding.winding_number,
                "note": note,
                "phase_ring_png_url": phase_png_url,
            }
        )
    return entries


def _run_vector_pattern_job(
    job, source_result, frequency_hz, theta_max_deg, n_theta, n_phi, beams, eps_r, tan_delta, thickness_mm, leakage_alpha_per_m
):
    if eps_r is None or tan_delta is None or thickness_mm is None:
        raise config_invalid(message="pattern_model='vector' needs eps_r, tan_delta and thickness_mm")
    x_m = source_result["x_m"]
    y_m = source_result["y_m"]
    aperture_field = source_result["aperture_field_real"] + 1j * source_result["aperture_field_imag"]
    period_m = source_result["period_mm"] * 1e-3
    feeds_mm = source_result.get("feeds_mm") or [(0.0, 0.0)]
    # PHASE15 Part D scoping (documented limit, PHASE15_POLARIZATION_VORTEX.md): a multi-feed
    # design's own per-region current direction is exact at SYNTHESIS time (Part A assigns each
    # beam its own feed); pattern prediction here vectorizes the combined aperture field about the
    # FIRST feed only -- a real per-region rho_hat model for the FAR-FIELD prediction itself is not
    # implemented. Correct (and the only feed) for every single-feed design.
    feed_xy_m = (feeds_mm[0][0] * 1e-3, feeds_mm[0][1] * 1e-3)

    substrate = Substrate(name="user", eps_r=eps_r, tan_delta=tan_delta, thickness=thickness_mm * 1e-3)
    green = LayeredGreenFunction(substrate, frequency_hz, backing="pec")
    jx, jy = aperture_current_from_scalar(x_m, y_m, aperture_field, feed_xy_m, leakage_alpha_per_m)
    pat = compute_vector_pattern(x_m, y_m, jx, jy, frequency_hz, green, period_m, period_m, theta_max_deg, n_theta, n_phi)

    beams_theta_phi = [(b["theta_deg"], b["phi_deg"]) for b in beams]
    uv_png_url, cut_png_urls = None, []
    out_dir = job_output_dir(job.id)
    try:
        uv_path = pattern_uv_map_png(
            pat.theta_deg, pat.phi_deg, pat.power_db, out_dir / "pattern_uv.png", beams_theta_phi
        )
        uv_png_url = f"/api/files/{job.id}/{uv_path.name}"
        cut_phis = [b[1] for b in beams_theta_phi] or [0.0]
        for i, phi_cut in enumerate(cut_phis):
            cut_path = pattern_cut_png(pat.theta_deg, pat.phi_deg, pat.power_db, phi_cut, out_dir / f"pattern_cut_{i}.png")
            cut_png_urls.append(f"/api/files/{job.id}/{cut_path.name}")
    except Exception:
        pass  # best-effort plots only

    # co/cross-pol u-v maps (the brief's Part B): only meaningful for ONE label at a time, so this
    # uses the FIRST beam's own (compensated) polarization if it has one; skipped (documented, not
    # silently wrong) for an all-"natural" beam set, which has no single co-pol to plot.
    copol_uv_url, cross_uv_url = None, None
    first_pol = beams[0]["desired_polarization"] if beams else "natural"
    if first_pol in _CROSS_LABEL:
        try:
            co = copol_component(pat.e_theta, pat.e_phi, first_pol)
            cross = copol_component(pat.e_theta, pat.e_phi, _CROSS_LABEL[first_pol])
            total_peak = np.max(np.abs(pat.e_theta) ** 2 + np.abs(pat.e_phi) ** 2)
            co_db = 10.0 * np.log10(np.maximum(np.abs(co) ** 2, 1e-30) / total_peak)
            cross_db = 10.0 * np.log10(np.maximum(np.abs(cross) ** 2, 1e-30) / total_peak)
            co_path = pattern_uv_map_png(
                pat.theta_deg, pat.phi_deg, co_db, out_dir / "pattern_uv_copol.png", beams_theta_phi,
                title=f"Co-pol ({first_pol})",
            )
            cross_path = pattern_uv_map_png(
                pat.theta_deg, pat.phi_deg, cross_db, out_dir / "pattern_uv_crosspol.png", beams_theta_phi,
                title=f"Cross-pol ({_CROSS_LABEL[first_pol]})",
            )
            copol_uv_url = f"/api/files/{job.id}/{co_path.name}"
            cross_uv_url = f"/api/files/{job.id}/{cross_path.name}"
        except Exception:
            pass

    per_beam = _per_beam_breakdown(job, x_m, y_m, jx, jy, frequency_hz, green, period_m, period_m, beams, pat, out_dir)

    return {
        "result": {
            "main_beam_theta_deg": pat.main_beam_theta_deg,
            "main_beam_phi_deg": pat.main_beam_phi_deg,
            "sll_db": pat.sll_db,
            "directivity_estimate_dbi": pat.directivity_estimate_dbi,
        },
        "labels": (LOCAL_PERIODICITY_APPROXIMATION, VECTOR_PATTERN_MODEL),
        "pattern_uv_png_url": uv_png_url,
        "pattern_cut_png_urls": cut_png_urls,
        "pattern_uv_copol_png_url": copol_uv_url,
        "pattern_uv_crosspol_png_url": cross_uv_url,
        "per_beam": per_beam,
    }


@router.post("/predict")
def predict(body: PatternRequest, request: Request):
    source = request.app.state.jobs.get(body.source_job_id)
    if source is None:
        raise HTTPException(404, f"job {body.source_job_id!r} not found")
    if source.status != JobStatus.DONE:
        raise HTTPException(409, f"source job {body.source_job_id!r} is {source.status.value}, not done")
    source_result = source.result or {}
    if "aperture_field_real" not in source_result and "ports" not in source_result:
        raise HTTPException(422, f"source job {body.source_job_id!r} did not produce an aperture field")
    try:
        result = _port_source_result(source_result, body.port_index)
    except Exception as exc:
        raise HTTPException(422, str(exc)) from exc

    beams = body.model_dump()["beams"]
    if body.pattern_model == "scalar":
        job = request.app.state.jobs.submit(
            "pattern_predict", _run_scalar_pattern_job,
            result, body.frequency_hz, body.theta_max_deg, body.n_theta, body.n_phi, beams,
        )
    elif body.pattern_model == "vector":
        job = request.app.state.jobs.submit(
            "pattern_predict", _run_vector_pattern_job,
            result, body.frequency_hz, body.theta_max_deg, body.n_theta, body.n_phi, beams,
            body.eps_r, body.tan_delta, body.thickness_mm, body.leakage_alpha_per_m,
        )
    else:
        raise HTTPException(400, f"pattern_model must be 'scalar' or 'vector', got {body.pattern_model!r}")
    return job.to_dict()
