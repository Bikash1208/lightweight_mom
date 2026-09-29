"""Phase 13 Part 3, page 8 ("HMSA designer"): import an HFSS field-calculator reference field,
combine it with one or more object beams, and synthesize a hologram via either of Phase 8B's two
formulations -- thin glue over `fields.fld_parser`/`fields.reference_field_grid`,
`dispersion.hfss_table`, and `hologram.hmsa_synthesis` (all reused unmodified).

Two synchronous upload endpoints (parsing is fast -- no job queue needed) plus one job-queued
synthesis endpoint (the object-wave combination + per-cell geometry lookup over the whole grid).
Formulation B's result also carries a synthesized aperture field (achieved phase from the library
match) AND a reflectarray-shaped `assignments` list (PHASE13B B2), so it can chain into the
Pattern AND Export pages exactly as a reflectarray design job does; Formulation A's achieved
quantity is a surface reactance, not a phase or an L/W pair, so neither chaining is offered for it
here (documented limitation, not silently wrong).
"""

import base64
import json
import uuid

import numpy as np
from fastapi import APIRouter, File, Form, HTTPException, Request, UploadFile
from pydantic import BaseModel

from ...dispersion.hfss_table import read_dispersion_csv
from ...errors import LwmomError, config_invalid, library_frequency_unavailable, phase_unreachable
from ...fields.fld_parser import read_fld
from ...fields.reference_field_grid import (
    check_grid_regularity,
    compare_beta_to_tm0,
    fld_to_reference_field_grid,
    radial_beta_profile,
)
from ...hologram.hmsa_synthesis import synthesize_formulation_a, synthesize_formulation_b
from ...labels import FEED_REGION_MASKED, reference_field_frequency_mismatch
from ...library_sweep import build_square_library_cells, write_square_library_json
from ...phase8b.object_wave import Beam
from ...reflectarray.designer import CellAssignment
from ...unit_cell_library import UnitCellLibrary
from .. import library_registry
from ..jobs import config_hash
from ..jsonutil import to_jsonable
from ..library_query import available_frequencies_hz, library_candidates_at_frequency
from ..plots import layout_map_png, reference_phase_map_png
from ..storage import job_output_dir, upload_dir

# PHASE13C: the wizard never sends the user to a separate Library page -- these are the fixed
# square-patch geometry family and mesh the auto-built board library always uses. A reasonable
# interactive-speed mesh, NOT alignment-searched (`discretization.smallest_aligned_cells` could
# need an unboundedly large mesh for an arbitrary user-chosen period/L combination -- deferred,
# documented in docs/phases/PHASE13_WEBAPP.md).
AUTO_LIBRARY_L_FRACTIONS_OF_PERIOD = (0.15, 0.30, 0.45, 0.60, 0.75, 0.90)
AUTO_LIBRARY_MESH_CELLS = 32

# PHASE13B A3: relative tolerance before an upload-vs-synthesis frequency difference is labeled a
# mismatch -- generous enough to absorb GHz/Hz round-trip and user-typed rounding, not a physics
# constant.
_FREQUENCY_MISMATCH_RTOL = 1e-4

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

    if frequency_hz is not None:
        # Persisted so `/synthesize` can carry this frequency forward and loudly label a mismatch
        # against whatever frequency synthesis is actually requested at (PHASE13B A3) -- the .fld
        # file itself carries no frequency (the app asks for it separately, docs/briefs/
        # PHASE13_WEBAPP_PROMPT.md), so this is the only place it's ever recorded.
        (updir / "meta.json").write_text(json.dumps({"frequency_hz": frequency_hz}))

    regularity = check_grid_regularity(phase_fld)
    r_bins, beta_median = radial_beta_profile(phase_fld)

    tm0 = None
    if eps_r is not None and thickness_mm is not None and frequency_hz is not None:
        beta_typical = float(np.nanmedian(beta_median))
        if np.isfinite(beta_typical):
            tm0 = compare_beta_to_tm0(beta_typical, frequency_hz, eps_r, thickness_mm * 1e-3)

    # 13C Step 1: "show immediately: a PNG of the reference phase map" -- inlined as base64 (no
    # job/download-route machinery needed for a one-off preview of a synchronous upload).
    phase_map_png_base64 = None
    try:
        png_path = reference_phase_map_png(
            phase_fld.x_m, phase_fld.y_m, np.degrees(phase_fld.value), updir / "phase_map.png"
        )
        phase_map_png_base64 = base64.b64encode(png_path.read_bytes()).decode("ascii")
    except Exception:
        pass  # best-effort preview only; never blocks the upload itself

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
            "phase_map_png_base64": phase_map_png_base64,
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
    rescale_to_beta_bar: bool = True  # formulation A (PHASE13B B3); on by default per the brief


def _load_reference_field(upload_id):
    updir = upload_dir(upload_id)
    phase_path = updir / "phase.fld"
    if not phase_path.exists():
        raise config_invalid(message=f"upload {upload_id!r} not found (or its files were cleaned up)")
    phase_fld = read_fld(phase_path)
    magnitude_path = updir / "magnitude.fld"
    magnitude_fld = read_fld(magnitude_path) if magnitude_path.exists() else None
    return fld_to_reference_field_grid(phase_fld, magnitude_fld)


def _uploaded_frequency_hz(upload_id):
    """The frequency recorded (if any) when this reference field was uploaded/analyzed -- see
    `upload_field`'s `meta.json`."""
    meta_path = upload_dir(upload_id) / "meta.json"
    if not meta_path.exists():
        return None
    return json.loads(meta_path.read_text()).get("frequency_hz")


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


def _layout_preview_url(job, x_m, y_m, value, label, title):
    """13C Step 3's "hologram/layout preview" PNG. Best-effort: a plotting failure must never
    fail a real synthesis result the user is waiting on."""
    try:
        out_dir = job_output_dir(job.id)
        path = layout_map_png(x_m, y_m, value, out_dir / "layout.png", title=title, label=label)
        return f"/api/files/{job.id}/{path.name}"
    except Exception:
        return None


def _run_synthesize_job(job, body: dict):
    reference_field = _load_reference_field(body["upload_id"])
    reference_field, masked_cells = _mask_feed_region(reference_field)
    extra_labels = (FEED_REGION_MASKED,) if masked_cells else ()
    upload_frequency_hz = _uploaded_frequency_hz(body["upload_id"])
    if upload_frequency_hz is not None and not np.isclose(
        upload_frequency_hz, body["frequency_hz"], rtol=_FREQUENCY_MISMATCH_RTOL
    ):
        # PHASE13B A3(b): the review found synthesis silently accepted a frequency that didn't
        # match the uploaded field's own (10.54 GHz field, 12.5 GHz synthesis -> 100% unreachable,
        # no explanation) -- carry the mismatch forward as a loud label rather than a silent no-op.
        extra_labels = (*extra_labels, reference_field_frequency_mismatch(upload_frequency_hz, body["frequency_hz"]))
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
            rescale_to_beta_bar=body["rescale_to_beta_bar"],
        )
        layout_png_url = _layout_preview_url(
            job, design.x_m, design.y_m, design.control_map, "reactance (ohm)", "Impedance modulation"
        )
        return {
            "formulation": design.formulation,
            "x_m": design.x_m,
            "y_m": design.y_m,
            "control_map": design.control_map,
            "cells": design.cells,
            "labels": design.labels + extra_labels,
            "unreachable_fraction": design.unreachable_fraction,
            "phase_error_stats": design.phase_error_stats,
            "masked_cells": masked_cells,
            "layout_png_url": layout_png_url,
        }

    if body["formulation"] == "B":
        if body["library_id"] is None:
            raise config_invalid(message="formulation B needs library_id")
        library_path = library_registry.resolve_path(body["library_id"])
        lib = UnitCellLibrary.from_json(library_path)
        candidates = library_candidates_at_frequency(lib, body["library_id"], body["frequency_hz"])
        if not candidates:
            # PHASE13B A3(a): report exactly what IS available instead of silently synthesizing
            # against zero candidates (every cell then has no achieved phase -- see A5 below).
            available_ghz = [f / 1e9 for f in available_frequencies_hz(lib)]
            raise library_frequency_unavailable(
                f"library {body['library_id']!r} has no record at {body['frequency_hz'] / 1e9:g} GHz; "
                f"available: {available_ghz} GHz"
            )
        design = synthesize_formulation_b(reference_field, beams, body["frequency_hz"], candidates, period_m)

        if any(c.achieved_value is None for c in design.cells):
            # PHASE13B A5: never hand Pattern prediction a placeholder (achieved_value defaulting
            # to 0deg -> aperture 1+0j) for a cell that has no real achieved phase -- refuse instead.
            n_missing = sum(1 for c in design.cells if c.achieved_value is None)
            raise phase_unreachable(
                f"{n_missing} of {len(design.cells)} cells have no achieved phase from library "
                f"{body['library_id']!r} at {body['frequency_hz'] / 1e9:g} GHz -- refusing to build "
                "a placeholder aperture field for Pattern prediction"
            )
        achieved = np.array([c.achieved_value for c in design.cells]).reshape(design.x_m.shape)
        aperture_field = design.aperture_field  # the RADIATED field, not the modulation phasor
        # alone -- see hologram.hmsa_synthesis.synthesize_formulation_b's docstring
        # PHASE13B B2: the same shape `export.py`'s `_lw_grid_from_assignments` reads off a
        # reflectarray design job (`.i/.j/.chosen_l_mm/.chosen_w_mm`) -- lets an HMSA design chain
        # straight into Export exactly as a reflectarray design already does.
        assignments = tuple(
            CellAssignment(
                i=c.i,
                j=c.j,
                x_m=c.x_m,
                y_m=c.y_m,
                required_phase_deg=c.control_value,
                chosen_l_mm=c.geometry_param,
                chosen_w_mm=c.geometry_param_w,
                achieved_phase_deg=c.achieved_value,
                phase_error_deg=c.phase_error_deg,
                library_id=body["library_id"],
                phase_coverage_ok=c.reachable,
            )
            for c in design.cells
        )
        layout_png_url = _layout_preview_url(
            job, design.x_m, design.y_m, achieved, "achieved phase (deg)", "Layout (achieved phase)"
        )
        return {
            "formulation": design.formulation,
            "x_m": design.x_m,
            "y_m": design.y_m,
            "period_mm": body["period_mm"],
            "control_map": design.control_map,
            "cells": design.cells,
            "assignments": assignments,
            "labels": design.labels + extra_labels,
            "unreachable_fraction": design.unreachable_fraction,
            "phase_error_stats": design.phase_error_stats,
            "masked_cells": masked_cells,
            "aperture_field_real": aperture_field.real,
            "aperture_field_imag": aperture_field.imag,
            "layout_png_url": layout_png_url,
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


class EnsureLibraryRequest(BaseModel):
    period_mm: float
    eps_r: float
    tan_delta: float
    thickness_mm: float
    frequency_hz: float


def _auto_library_id(body_dict):
    payload = {
        "kind": "hmsa_auto_library",
        "l_fractions": AUTO_LIBRARY_L_FRACTIONS_OF_PERIOD,
        "mesh_cells": AUTO_LIBRARY_MESH_CELLS,
        **body_dict,
    }
    return "auto_" + config_hash(payload)[:20]


def _run_ensure_library_job(job, body_dict, library_id):
    l_values_mm = tuple(round(body_dict["period_mm"] * f, 6) for f in AUTO_LIBRARY_L_FRACTIONS_OF_PERIOD)
    base_config = {
        "period_x": body_dict["period_mm"] * 1e-3,
        "period_y": body_dict["period_mm"] * 1e-3,
        "substrate": {
            "name": "user",
            "eps_r": body_dict["eps_r"],
            "tan_delta": body_dict["tan_delta"],
            "thickness": body_dict["thickness_mm"] * 1e-3,
            "backing": "pec",
        },
        "geometry": {"type": "rectangle", "length": l_values_mm[0], "width": l_values_mm[0]},  # placeholder
        "excitation": {"frequency": body_dict["frequency_hz"]},  # placeholder
        "mesh": {
            "cells_x": AUTO_LIBRARY_MESH_CELLS,
            "cells_y": AUTO_LIBRARY_MESH_CELLS,
            "floquet_order": AUTO_LIBRARY_MESH_CELLS,
        },
        "reference_planes": {"top_air_thickness": 5.0e-3},
    }

    def _progress(event):
        job.progress = event["done"] / max(event["total"], 1)
        job.progress_message = f"{event['done']}/{event['total']} cell sizes solved"

    cells = build_square_library_cells(base_config, l_values_mm, body_dict["frequency_hz"], progress=_progress)
    out_dir = job_output_dir(job.id)
    path = write_square_library_json(cells, out_dir / "unit_cell_library.json")
    library_registry.register(
        library_id,
        path,
        meta={
            "period_mm": body_dict["period_mm"],
            "backing": "pec",
            "substrate": "user",
            "n_geometries": len(cells),
            "auto_built": True,
            "frequency_hz": body_dict["frequency_hz"],
        },
    )
    return {"library_id": library_id, "n_geometries": len(cells), "l_values_mm": l_values_mm, "already_built": False}


@router.post("/ensure_library")
def ensure_library(body: EnsureLibraryRequest, request: Request):
    """PHASE13C: the wizard never sends the user to a separate Library page -- it builds (or
    reuses -- across requests via the job cache, and across app restarts via the persisted
    `library_registry`) a small square-patch library for the chosen board/frequency automatically.
    Always returns a job (even the "already built" fast path), so the frontend has one uniform
    poll-until-done code path."""
    library_id = _auto_library_id(body.model_dump())
    try:
        library_registry.resolve_path(library_id)
        job = request.app.state.jobs.submit(
            "hmsa_ensure_library", lambda job: {"library_id": library_id, "already_built": True}
        )
        return job.to_dict()
    except KeyError:
        pass
    key = config_hash({"kind": "hmsa_ensure_library", **body.model_dump()})
    job = request.app.state.jobs.submit(
        "hmsa_ensure_library", _run_ensure_library_job, body.model_dump(), library_id, cache_key=key
    )
    return job.to_dict()
