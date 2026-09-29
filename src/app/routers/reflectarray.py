"""Phase 13 Part 3, page 3 ("Reflectarray designer"): required-phase-map -> cell selection from a
registered library -> phase-gradient/coverage labels. Thin glue over
`reflectarray.designer.design_reflectarray` (reused unmodified) and the same library query API the
Library page uses. The design job's result also carries a synthesized aperture field (feed taper x
achieved phase) so the Pattern page can chain off it without re-uploading the array through the
browser (`app/jobs.py`'s `Job.result` stays raw numpy in memory; only the HTTP response is
flattened by `to_jsonable`).
"""

import numpy as np
from fastapi import APIRouter, HTTPException, Request
from pydantic import BaseModel

from ...errors import library_frequency_unavailable
from ...reflectarray.designer import design_reflectarray, feed_illumination_taper
from ...unit_cell_library import UnitCellLibrary
from .. import library_registry
from ..library_query import available_frequencies_hz, library_candidates_at_frequency

router = APIRouter(prefix="/api/reflectarray", tags=["reflectarray"])


class ReflectarrayDesignRequest(BaseModel):
    library_id: str
    frequency_hz: float
    nx: int
    ny: int
    period_mm: float
    feed_xyz_mm: "tuple[float, float, float]"
    theta0_deg: float
    phi0_deg: float
    taper_q: float = 1.0


def _run_design_job(job, library_path, library_id, params):
    lib = UnitCellLibrary.from_json(library_path)
    candidates = library_candidates_at_frequency(lib, library_id, params["frequency_hz"])
    if not candidates:
        available_ghz = [f / 1e9 for f in available_frequencies_hz(lib)]
        raise library_frequency_unavailable(
            f"library {library_id!r} has no record at {params['frequency_hz'] / 1e9:g} GHz; "
            f"available: {available_ghz} GHz"
        )

    period_m = params["period_mm"] * 1e-3
    x_vals = (np.arange(params["nx"]) - (params["nx"] - 1) / 2) * period_m
    y_vals = (np.arange(params["ny"]) - (params["ny"] - 1) / 2) * period_m
    x_m, y_m = np.meshgrid(x_vals, y_vals)
    feed_xyz_m = tuple(v * 1e-3 for v in params["feed_xyz_mm"])

    design = design_reflectarray(
        x_m, y_m, params["frequency_hz"], feed_xyz_m, params["theta0_deg"], params["phi0_deg"], candidates
    )
    taper = feed_illumination_taper(x_m, y_m, feed_xyz_m, q=params["taper_q"])
    achieved = np.array(
        [a.achieved_phase_deg if a.achieved_phase_deg is not None else 0.0 for a in design.assignments]
    ).reshape(x_m.shape)
    aperture_field = taper * np.exp(1j * np.radians(achieved))

    return {
        "library_id": library_id,
        "x_m": x_m,
        "y_m": y_m,
        "period_mm": params["period_mm"],
        "required_phase_deg": design.required_phase_deg,
        "achieved_phase_deg": achieved,
        "assignments": design.assignments,
        "max_phase_gradient_deg_per_cell": design.max_phase_gradient_deg_per_cell,
        "phase_error_stats": design.phase_error_stats,
        "labels": design.labels,
        "aperture_field_real": aperture_field.real,
        "aperture_field_imag": aperture_field.imag,
    }


@router.post("/design")
def design(body: ReflectarrayDesignRequest, request: Request):
    try:
        library_path = library_registry.resolve_path(body.library_id)
    except KeyError as exc:
        raise HTTPException(404, str(exc)) from exc
    job = request.app.state.jobs.submit(
        "reflectarray_design", _run_design_job, library_path, body.library_id, body.model_dump()
    )
    return job.to_dict()
