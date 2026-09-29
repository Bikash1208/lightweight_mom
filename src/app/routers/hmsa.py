"""Phase 13 Part 3, page 8 ("HMSA designer"): import an HFSS field-calculator reference field,
combine it with one or more object beams, and synthesize a hologram via either of Phase 8B's two
formulations -- thin glue over `fields.fld_parser`/`fields.reference_field_grid`,
`dispersion.hfss_table`, and `hologram.hmsa_synthesis` (all reused unmodified).

Two synchronous upload endpoints (parsing is fast -- no job queue needed) plus one job-queued
synthesis endpoint (the object-wave combination + per-cell geometry lookup over the whole grid).
Formulation B's result also carries a synthesized aperture field (achieved phase from the library
match), so it can chain into the Pattern page exactly as a reflectarray design job does;
Formulation A's achieved quantity is a surface reactance, not a phase, so that chaining is not
offered for it here (documented limitation, not silently wrong).
"""

import uuid

import numpy as np
from fastapi import APIRouter, File, Form, HTTPException, Request, UploadFile
from pydantic import BaseModel

from ...dispersion.hfss_table import read_dispersion_csv
from ...errors import LwmomError, config_invalid
from ...fields.fld_parser import read_fld
from ...fields.reference_field_grid import (
    check_grid_regularity,
    compare_beta_to_tm0,
    fld_to_reference_field_grid,
    radial_beta_profile,
)
from ...hologram.hmsa_synthesis import synthesize_formulation_a, synthesize_formulation_b
from ...labels import FEED_REGION_MASKED
from ...phase8b.object_wave import Beam
from ...unit_cell_library import UnitCellLibrary
from .. import library_registry
from ..jsonutil import to_jsonable
from ..library_query import library_candidates_at_frequency
from ..storage import upload_dir

router = APIRouter(prefix="/api/hmsa", tags=["hmsa"])


@router.post("/upload_field")
async def upload_field(
    phase_file: UploadFile = File(...),
    magnitude_file: "UploadFile | None" = File(None),
    eps_r: "float | None" = Form(None),
    thickness_mm: "float | None" = Form(None),
    frequency_hz: "float | None" = Form(None),
):
    upload_id = str(uuid.uuid4())
    updir = upload_dir(upload_id)
    phase_path = updir / "phase.fld"
    phase_path.write_bytes(await phase_file.read())

    try:
        phase_fld = read_fld(phase_path)
    except (ValueError, LwmomError) as exc:
        raise HTTPException(400, str(exc)) from exc

    magnitude_path = None
    if magnitude_file is not None:
        magnitude_path = updir / "magnitude.fld"
        magnitude_path.write_bytes(await magnitude_file.read())

    regularity = check_grid_regularity(phase_fld)
    r_bins, beta_median = radial_beta_profile(phase_fld)

    tm0 = None
    if eps_r is not None and thickness_mm is not None and frequency_hz is not None:
        beta_typical = float(np.nanmedian(beta_median))
        if np.isfinite(beta_typical):
            tm0 = compare_beta_to_tm0(beta_typical, frequency_hz, eps_r, thickness_mm * 1e-3)

    return to_jsonable(
        {
            "upload_id": upload_id,
            "has_magnitude": magnitude_path is not None,
            "nx": phase_fld.value.shape[1],
            "ny": phase_fld.value.shape[0],
            "ordering": phase_fld.ordering,
            "unit": phase_fld.header.unit,
            "coordinate_system": phase_fld.header.coordinate_system,
            "nan_count": int(phase_fld.nan_mask.sum()),
            "regularity": regularity,
            "radial_beta_r_m": r_bins,
            "radial_beta_median_rad_per_m": beta_median,
            "tm0_comparison": tm0,
        }
    )


@router.post("/upload_dispersion")
async def upload_dispersion(csv_file: UploadFile = File(...), frequency_hz: float = Form(...)):
    upload_id = str(uuid.uuid4())
    updir = upload_dir(upload_id)
    csv_path = updir / "dispersion.csv"
    csv_path.write_bytes(await csv_file.read())
    table = read_dispersion_csv(csv_path, frequency_hz)
    return {
        "dispersion_id": upload_id,
        "n_samples": len(table.samples),
        "geometry_range": table.geometry_range,
        "has_beta": table.samples[0].beta_rad_per_m is not None,
        "has_x_s": table.samples[0].x_s_ohm is not None,
    }


class ObjectBeamSpec(BaseModel):
    theta_deg: float
    phi_deg: float
    weight_real: float = 1.0
    weight_imag: float = 0.0
    initial_phase_deg: float = 0.0
    polarization: str = "x"


class HmsaSynthesizeRequest(BaseModel):
    upload_id: str
    formulation: str  # "A" | "B"
    frequency_hz: float
    object_beams: "list[ObjectBeamSpec]"
    period_mm: float
    library_id: "str | None" = None  # formulation B
    dispersion_id: "str | None" = None  # formulation A
    x_avg_ohm: "float | None" = None  # formulation A
    modulation_index: "float | None" = None  # formulation A


def _load_reference_field(upload_id):
    updir = upload_dir(upload_id)
    phase_path = updir / "phase.fld"
    if not phase_path.exists():
        raise config_invalid(message=f"upload {upload_id!r} not found (or its files were cleaned up)")
    phase_fld = read_fld(phase_path)
    magnitude_path = updir / "magnitude.fld"
    magnitude_fld = read_fld(magnitude_path) if magnitude_path.exists() else None
    return fld_to_reference_field_grid(phase_fld, magnitude_fld)


def _mask_feed_region(reference_field):
    """The parsed field's phase can carry a NaN at the feed point (`fld_parser.py`'s own
    docstring: "the single feed-point NaN"); Phase 8B's hologram formulas reject non-finite input
    outright, so this fills the phase-only path by nearest-neighbor (the same technique
    `reference_field_grid._fill_nan_for_unwrap` uses for the wavenumber sanity check) purely so
    synthesis can run, and returns the masked grid indices so the caller can report
    `FEED_REGION_MASKED` honestly rather than silently interpolate and say nothing. Only the
    phase-only path is handled: a complex (amplitude-available) field with NaN in `Ex` would need
    `Ex` itself filled too, a path formulation A never reaches in practice here since the
    project's only real reference-field export is phase-only."""
    from dataclasses import replace

    phase = reference_field.phase()
    nan_mask = np.isnan(phase)
    if not nan_mask.any():
        return reference_field, []
    if reference_field.amplitude_available:
        return reference_field, []  # not handled for the complex-field path -- see docstring

    from scipy.interpolate import griddata

    good_y, good_x = np.nonzero(~nan_mask)
    bad_y, bad_x = np.nonzero(nan_mask)
    filled = phase.copy()
    filled[nan_mask] = griddata((good_y, good_x), phase[~nan_mask], (bad_y, bad_x), method="nearest")
    masked_indices = [{"i": int(i), "j": int(j)} for i, j in zip(bad_y.tolist(), bad_x.tolist(), strict=True)]
    return replace(reference_field, phase_deg=filled), masked_indices


def _run_synthesize_job(job, body: dict):
    reference_field = _load_reference_field(body["upload_id"])
    reference_field, masked_cells = _mask_feed_region(reference_field)
    feed_labels = (FEED_REGION_MASKED,) if masked_cells else ()
    beams = [
        Beam(
            theta_deg=b["theta_deg"],
            phi_deg=b["phi_deg"],
            weight=complex(b["weight_real"], b["weight_imag"]),
            initial_phase_deg=b["initial_phase_deg"],
            polarization=b["polarization"],
        )
        for b in body["object_beams"]
    ]
    period_m = body["period_mm"] * 1e-3

    if body["formulation"] == "A":
        if body["dispersion_id"] is None or body["x_avg_ohm"] is None or body["modulation_index"] is None:
            raise config_invalid(message="formulation A needs dispersion_id, x_avg_ohm and modulation_index")
        dispersion_path = upload_dir(body["dispersion_id"]) / "dispersion.csv"
        dispersion_table = read_dispersion_csv(dispersion_path, body["frequency_hz"])
        design = synthesize_formulation_a(
            reference_field,
            beams,
            body["frequency_hz"],
            body["x_avg_ohm"],
            body["modulation_index"],
            dispersion_table,
            period_m,
        )
        return {
            "formulation": design.formulation,
            "x_m": design.x_m,
            "y_m": design.y_m,
            "control_map": design.control_map,
            "cells": design.cells,
            "labels": design.labels + feed_labels,
            "unreachable_fraction": design.unreachable_fraction,
            "masked_cells": masked_cells,
        }

    if body["formulation"] == "B":
        if body["library_id"] is None:
            raise config_invalid(message="formulation B needs library_id")
        library_path = library_registry.resolve_path(body["library_id"])
        lib = UnitCellLibrary.from_json(library_path)
        candidates = library_candidates_at_frequency(lib, body["library_id"], body["frequency_hz"])
        design = synthesize_formulation_b(reference_field, beams, body["frequency_hz"], candidates, period_m)

        achieved = np.array([c.achieved_value if c.achieved_value is not None else 0.0 for c in design.cells]).reshape(
            design.x_m.shape
        )
        aperture_field = np.exp(1j * np.radians(achieved))
        return {
            "formulation": design.formulation,
            "x_m": design.x_m,
            "y_m": design.y_m,
            "period_mm": body["period_mm"],
            "control_map": design.control_map,
            "cells": design.cells,
            "labels": design.labels + feed_labels,
            "unreachable_fraction": design.unreachable_fraction,
            "masked_cells": masked_cells,
            "aperture_field_real": aperture_field.real,
            "aperture_field_imag": aperture_field.imag,
        }

    raise config_invalid(message=f"formulation must be 'A' or 'B', got {body['formulation']!r}")


@router.post("/synthesize")
def synthesize(body: HmsaSynthesizeRequest, request: Request):
    if body.formulation == "B" and body.library_id is not None:
        try:
            library_registry.resolve_path(body.library_id)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc
    job = request.app.state.jobs.submit("hmsa_synthesize", _run_synthesize_job, body.model_dump())
    return job.to_dict()
