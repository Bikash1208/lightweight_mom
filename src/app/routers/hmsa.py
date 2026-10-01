"""Phase 13 Part 3, page 8 ("HMSA designer"): import an HFSS field-calculator reference field,
combine it with one or more object beams, and synthesize a hologram via either of Phase 8B's two
formulations -- thin glue over `fields.fld_parser`/`fields.reference_field_grid`,
`dispersion.hfss_table`, and `hologram.hmsa_synthesis` (all reused unmodified).

Two synchronous upload endpoints (parsing is fast -- no job queue needed) plus one job-queued
synthesis endpoint (the object-wave combination + per-cell geometry lookup over the whole grid).
Formulation B's result carries a synthesized aperture field (achieved phase from the library
match) AND a reflectarray-shaped `assignments` list (PHASE13B B2), so it can chain into the
Pattern AND Export pages exactly as a reflectarray design job does. Formulation A's result (PHASE14
Part C) also carries a synthesized aperture field (a documented leaky-wave-perturbation
approximation, `hologram.hmsa_synthesis.synthesize_formulation_a`'s own docstring) so it can chain
into Pattern too -- but its achieved quantity is a surface reactance, not an L/W pair, so Export
still isn't offered for it (a real geometry-inversion gap, not silently wrong).
"""

import base64
import json
import uuid
from typing import Literal

import numpy as np
from fastapi import APIRouter, File, Form, HTTPException, Request, UploadFile
from pydantic import BaseModel

from ...basis_functions import MeshSettings
from ...dispersion.analytic_tm0 import tm0_beta
from ...dispersion.hfss_table import read_dispersion_csv
from ...dispersion.mom_scaffold import square_patch_unit_cell
from ...dispersion.sheet_model import build_dispersion_table_from_sheet_model
from ...errors import LwmomError, config_invalid, library_frequency_unavailable, phase_unreachable
from ...feasibility_coverage import library_table_at_frequency, phase_density
from ...fields.fld_parser import read_fld
from ...fields.reference_field_grid import (
    check_grid_regularity,
    compare_beta_to_tm0,
    feed_position_from_nan_mask,
    fld_to_reference_field_grid,
    radial_beta_profile,
    resample_to_cell_grid,
)
from ...green_function import LayeredGreenFunction
from ...hologram.coded_synthesis import (
    continuous_reference_aperture_field,
    synthesize_coded_n_bit,
    synthesize_coded_one_bit,
)
from ...hologram.hmsa_synthesis import synthesize_formulation_a, synthesize_formulation_b
from ...hologram.multi_feed import (
    FeedSpec,
    feed_xy_m_for_beams,
    segmented_reference_field_and_regions,
    superposed_reference_field,
)
from ...hologram.multi_port_synthesis import combined_board_lw_grid_m, synthesize_ports
from ...hologram.ports import custom_layout, layout_preset_1xn, layout_preset_2x2_plus_1_large
from ...hologram.two_state_cells import (
    best_separated_pair,
    candidates_from_annular_slot_library,
    candidates_from_square_library,
    separation_warning,
    two_state_from_values,
)
from ...labels import (
    FEED_REGION_MASKED,
    MULTI_PORT_SHARED_METHOD,
    PORT_ISOLATION_NOT_FULL_WAVE,
    feed_position_auto_detected,
    isolation_edge_distance,
    quantization_loss_db,
    reference_field_frequency_mismatch,
    two_state_phase_separation,
)
from ...library_sweep import (
    build_annular_slot_library_cells,
    build_square_library_cells,
    write_annular_slot_library_json,
    write_square_library_json,
)
from ...materials import Metal, Substrate
from ...pattern.vector_aperture_pattern import aperture_current_from_scalar, compute_vector_pattern
from ...phase8b.object_wave import Beam
from ...reflectarray.designer import CellAssignment
from ...unit_cell_library import UnitCellLibrary
from .. import library_registry
from ..jobs import config_hash
from ..jsonutil import to_jsonable
from ..library_query import available_frequencies_hz, library_candidates_at_frequency
from ..plots import layout_map_png, patch_layout_png, reference_phase_map_png
from ..storage import job_output_dir, upload_dir

# PHASE13C/13D: the wizard never sends the user to a separate Library page -- these are the fixed
# square-patch geometry family and mesh the auto-built board library always uses. Mesh-aligned BY
# CONSTRUCTION (not a search, so no risk of `smallest_aligned_cells` needing an unboundedly large
# mesh): a centred feature of length L on an N-cell mesh over `period` has both edges on mesh
# lines iff `L / period == k / (N / 2)` for integer k (`discretization.feature_edges_aligned`'s
# condition, solved for L) -- i.e. every fraction is a multiple of `2 / N`.
# PHASE13D_POLISH_AND_PHYSICS.md's live walk-through found the old 6-point family (step 0.15, NOT
# a multiple of 2/32 -- unaligned, so subject to staircasing error on top of being sparse) too
# sparse: 4006/5476 cells unreachable at a 4 mm cell size. Both tiers below are aligned and denser
# (13 points, step 0.0625); "accurate" also uses a finer mesh (less staircasing error per solve --
# NOT Richardson extrapolation, which the auto-library builder does not currently support; an
# honest speed/resolution tradeoff, not claimed to be more than that).
AUTO_LIBRARY_MESH_CELLS_BY_ACCURACY = {"fast": 32, "accurate": 48}


def _l_fractions_of_period(mesh_cells):
    return tuple(k / (mesh_cells / 2) for k in range(2, 15))

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
    # PHASE15 Part A/D: desired_polarization/oam_charge map directly onto phase8b.object_wave.Beam's
    # OWN same-named fields (default "natural"/0 -- today's uncompensated behaviour, unchanged for
    # any caller that omits these, e.g. every pre-Phase-15 request body). "polarization" above is
    # the pre-existing, unrelated field (Sec. 18's Ex/Ey split) and is untouched.
    desired_polarization: str = "natural"
    oam_charge: int = 0
    feed_index: int = 0  # PHASE15 Part C: which entry of HmsaSynthesizeRequest.feeds this beam's
    # current is radial about (ignored for a single-feed design, body.feeds is None)


class FeedSpecRequest(BaseModel):
    """PHASE15 Part C: one feed of a multi-feed design -- `upload_id` (an `.fld` already analyzed
    via `/api/hmsa/upload_field`) for a real measured feed, or omitted for an analytic TM0
    cylindrical reference about `(x_mm, y_mm)` (needs the board's own eps_r/thickness, already
    known from the request's substrate-independent synthesize call -- PHASE15 Part C's
    `HmsaSynthesizeRequest` carries them only when `feeds` is given, see its own docstring)."""

    x_mm: float
    y_mm: float
    amplitude_real: float = 1.0
    amplitude_imag: float = 0.0
    upload_id: "str | None" = None


class HmsaSynthesizeRequest(BaseModel):
    upload_id: "str | None" = None  # None only when `feeds` is given (PHASE15 Part C)
    formulation: str  # "A" | "B"
    frequency_hz: float
    object_beams: "list[ObjectBeamSpec]"
    period_mm: float
    library_id: "str | None" = None  # formulation B
    dispersion_id: "str | None" = None  # formulation A
    x_avg_ohm: "float | None" = None  # formulation A
    modulation_index: "float | None" = None  # formulation A
    rescale_to_beta_bar: bool = True  # formulation A (PHASE13B B3); on by default per the brief
    # PHASE15 Part A
    cp_phase_offset_deg: float = 0.0
    multiplexing: str = "full_aperture"  # "full_aperture" | "sectors"
    feed_x_mm: "float | None" = None  # explicit feed-position override; None -> auto-detect from
    feed_y_mm: "float | None" = None  # the uploaded field's own NaN/singularity (Part A's default)
    # PHASE15 Part C: multi-feed. None (default) is the ordinary single-`upload_id` path, unchanged.
    feeds: "list[FeedSpecRequest] | None" = None
    multi_feed_strategy: str = "superposition"  # "superposition" | "segmentation"
    multi_feed_eps_r: "float | None" = None  # needed only for an analytic-TM0 feed (no upload_id)
    multi_feed_thickness_mm: "float | None" = None


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


def _resolve_feed_xy_m(body, raw_reference_field, raw_nan_mask):
    """PHASE15 Part A: "Feed position (x0, y0) is taken from the reference field's NaN/singularity
    location (default) or entered." Explicit `feed_x_mm`/`feed_y_mm` wins; otherwise
    `feed_position_from_nan_mask` on the RAW (pre-fill) field, matching where the actual monopole
    sits. Returns `(x0_m, y0_m, label_or_None)`."""
    if body["feed_x_mm"] is not None and body["feed_y_mm"] is not None:
        return body["feed_x_mm"] * 1e-3, body["feed_y_mm"] * 1e-3, None
    x0_m, y0_m = feed_position_from_nan_mask(raw_reference_field.x_m, raw_reference_field.y_m, raw_nan_mask)
    label = feed_position_auto_detected(x0_m * 1e3, y0_m * 1e3) if raw_nan_mask.any() else None
    return x0_m, y0_m, label


def _resolve_multi_feed(body, frequency_hz):
    """PHASE15 Part C: builds the combined `ReferenceFieldGrid` (superposed or segmented) plus the
    per-beam `(x0_sequence, y0_sequence)` `feed_xy_m` `synthesize_formulation_a/b` needs, from
    `HmsaSynthesizeRequest.feeds`. Returns `(reference_field, feed_xy_m, feeds_mm, extra_labels)`.
    `x_m`/`y_m` for the combined field are taken from the FIRST feed that carries an uploaded
    field (real grid extent); if no feed has an upload, a grid must come from elsewhere -- this
    project's only real reference-field export is phase-only per-feed, so at least one upload is
    required in that case (raised as `config_invalid`, not a silent synthetic grid guess)."""
    feed_specs = []
    k_swg = None
    uploaded_grid = None
    for f in body["feeds"]:
        rf = None
        if f["upload_id"] is not None:
            rf = _load_reference_field(f["upload_id"])
            rf, _masked = _mask_feed_region(rf)
            if uploaded_grid is None:
                uploaded_grid = (rf.x_m, rf.y_m)
        feed_specs.append(
            FeedSpec(x_m=f["x_mm"] * 1e-3, y_m=f["y_mm"] * 1e-3, amplitude=complex(f["amplitude_real"], f["amplitude_imag"]), reference_field=rf)
        )
    if any(fs.reference_field is None for fs in feed_specs):
        if body["multi_feed_eps_r"] is None or body["multi_feed_thickness_mm"] is None:
            raise config_invalid(
                message="a feed with no upload_id needs multi_feed_eps_r/multi_feed_thickness_mm "
                "for its analytic TM0 reference"
            )
        k_swg = tm0_beta(frequency_hz, body["multi_feed_eps_r"], body["multi_feed_thickness_mm"] * 1e-3)
    if uploaded_grid is None:
        raise config_invalid(message="multi-feed synthesis needs at least one feed with an uploaded reference field")
    x_m, y_m = uploaded_grid

    if body["multi_feed_strategy"] == "superposition":
        combined = superposed_reference_field(feed_specs, x_m, y_m, k_swg)
    elif body["multi_feed_strategy"] == "segmentation":
        combined, _owner = segmented_reference_field_and_regions(feed_specs, x_m, y_m, k_swg)
    else:
        raise config_invalid(message=f"multi_feed_strategy must be 'superposition' or 'segmentation', got {body['multi_feed_strategy']!r}")

    feed_indices = [b["feed_index"] for b in body["object_beams"]]
    feed_xy_m = feed_xy_m_for_beams(feed_specs, feed_indices)
    feeds_mm = [(fs.x_m * 1e3, fs.y_m * 1e3) for fs in feed_specs]
    extra_labels = (f"MULTI_FEED_{body['multi_feed_strategy'].upper()}({len(feed_specs)} feeds)",)
    return combined, feed_xy_m, feeds_mm, extra_labels


def _layout_preview_url(job, x_m, y_m, value, label, title):
    """13C Step 3's "hologram/layout preview" PNG. Best-effort: a plotting failure must never
    fail a real synthesis result the user is waiting on."""
    try:
        out_dir = job_output_dir(job.id)
        path = layout_map_png(x_m, y_m, value, out_dir / "layout.png", title=title, label=label)
        return f"/api/files/{job.id}/{path.name}"
    except Exception:
        return None


def _patch_layout_url(job, assignments, x_m, period_mm, title_suffix=""):
    """True-geometry copper/grey board drawing (PHASE13D item 5). Best-effort like the preview."""
    try:
        ny, nx = np.asarray(x_m).shape
        cells = [(a.x_m * 1e3, a.y_m * 1e3, a.chosen_l_mm, a.chosen_w_mm) for a in assignments]
        out_dir = job_output_dir(job.id)
        title = f"Patch layout — {nx} × {ny} cells, {period_mm:g} mm pitch{title_suffix}"
        path = patch_layout_png(cells, period_mm, nx, ny, out_dir / "patch_layout.png", title=title)
        return f"/api/files/{job.id}/{path.name}"
    except Exception:
        return None


def _geometry_param_to_mm(gp, period_mm):
    """Formulation A's dispersion-table geometry parameter is the square patch side: metres for
    the auto-built table (`ensure_dispersion`), millimetres for a typical HFSS upload. Anything
    smaller than 0.1 is taken as metres."""
    if gp is None:
        return None
    mm = gp * 1e3 if gp < 0.1 else gp
    return mm if 0.0 < mm <= period_mm else None


def _run_synthesize_job(job, body: dict):
    extra_labels = ()
    masked_cells = []
    feeds_mm = None

    if body["feeds"] is not None:
        reference_field, feed_xy_m, feeds_mm, multi_feed_labels = _resolve_multi_feed(body, body["frequency_hz"])
        extra_labels = (*extra_labels, *multi_feed_labels)
    else:
        if body["upload_id"] is None:
            raise config_invalid(message="synthesize needs either upload_id or feeds")
        raw_reference_field = _load_reference_field(body["upload_id"])
        raw_nan_mask = np.isnan(raw_reference_field.phase())
        reference_field, masked_cells = _mask_feed_region(raw_reference_field)
        if masked_cells:
            extra_labels = (*extra_labels, FEED_REGION_MASKED)
        upload_frequency_hz = _uploaded_frequency_hz(body["upload_id"])
        if upload_frequency_hz is not None and not np.isclose(
            upload_frequency_hz, body["frequency_hz"], rtol=_FREQUENCY_MISMATCH_RTOL
        ):
            # PHASE13B A3(b): the review found synthesis silently accepted a frequency that didn't
            # match the uploaded field's own (10.54 GHz field, 12.5 GHz synthesis -> 100% unreachable,
            # no explanation) -- carry the mismatch forward as a loud label rather than a silent no-op.
            extra_labels = (*extra_labels, reference_field_frequency_mismatch(upload_frequency_hz, body["frequency_hz"]))
        x0_m, y0_m, feed_label = _resolve_feed_xy_m(body, raw_reference_field, raw_nan_mask)
        feed_xy_m = (x0_m, y0_m)
        feeds_mm = [(x0_m * 1e3, y0_m * 1e3)]
        if feed_label is not None:
            extra_labels = (*extra_labels, feed_label)

    beams = [
        Beam(
            theta_deg=b["theta_deg"],
            phi_deg=b["phi_deg"],
            weight=complex(b["weight_real"], b["weight_imag"]),
            initial_phase_deg=b["initial_phase_deg"],
            polarization=b["polarization"],
            desired_polarization=b["desired_polarization"],
            oam_charge=b["oam_charge"],
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
            feed_xy_m=feed_xy_m,
            cp_phase_offset_deg=body["cp_phase_offset_deg"],
            multiplexing=body["multiplexing"],
        )
        layout_png_url = _layout_preview_url(
            job, design.x_m, design.y_m, design.control_map, "reactance (ohm)", "Impedance modulation"
        )
        aperture_field = design.aperture_field
        # Surface-wave designs are square patches of the dispersion table's geometry parameter:
        # carry the same per-cell L/W assignment shape Formulation B uses, so Export and the
        # patch-layout drawing work for this method too.
        assignments_a = None
        sizes_mm = [_geometry_param_to_mm(c.geometry_param, body["period_mm"]) for c in design.cells]
        if all(v is not None for v in sizes_mm):
            assignments_a = tuple(
                CellAssignment(
                    i=c.i, j=c.j, x_m=c.x_m, y_m=c.y_m,
                    required_phase_deg=c.control_value,
                    chosen_l_mm=v, chosen_w_mm=v,
                    achieved_phase_deg=c.achieved_value, phase_error_deg=c.phase_error_deg,
                    library_id=body["dispersion_id"], phase_coverage_ok=c.reachable,
                )
                for c, v in zip(design.cells, sizes_mm, strict=True)
            )
        patch_layout_png_url = (
            _patch_layout_url(job, assignments_a, design.x_m, body["period_mm"]) if assignments_a else None
        )
        return {
            "formulation": design.formulation,
            "x_m": design.x_m,
            "y_m": design.y_m,
            "control_map": design.control_map,
            "cells": design.cells,
            "assignments": assignments_a,
            "patch_layout_png_url": patch_layout_png_url,
            "labels": design.labels + extra_labels,
            "unreachable_fraction": design.unreachable_fraction,
            "phase_error_stats": design.phase_error_stats,
            "masked_cells": masked_cells,
            # PHASE14 Part C: Formulation A now carries a (documented-approximate) aperture field
            # too, so Pattern prediction works for both methods -- previously only B did.
            "aperture_field_real": aperture_field.real,
            "aperture_field_imag": aperture_field.imag,
            "period_mm": body["period_mm"],
            "layout_png_url": layout_png_url,
            "feeds_mm": feeds_mm,  # PHASE15 Part C: feed marker(s), for Export/Pattern
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
        design = synthesize_formulation_b(
            reference_field, beams, body["frequency_hz"], candidates, period_m,
            feed_xy_m=feed_xy_m, cp_phase_offset_deg=body["cp_phase_offset_deg"], multiplexing=body["multiplexing"],
        )

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
        patch_layout_png_url = _patch_layout_url(job, assignments, design.x_m, body["period_mm"])
        return {
            "patch_layout_png_url": patch_layout_png_url,
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
            "feeds_mm": feeds_mm,  # PHASE15 Part C: feed marker(s), for Export/Pattern
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
    accuracy: Literal["fast", "accurate"] = "fast"  # PHASE13D item 2 ("Words"): the Advanced
    # accuracy control the brief asked for.


def _auto_library_id(body_dict):
    mesh_cells = AUTO_LIBRARY_MESH_CELLS_BY_ACCURACY[body_dict["accuracy"]]
    payload = {
        "kind": "hmsa_auto_library",
        "l_fractions": _l_fractions_of_period(mesh_cells),
        "mesh_cells": mesh_cells,
        **body_dict,
    }
    return "auto_" + config_hash(payload)[:20]


def _run_ensure_library_job(job, body_dict, library_id):
    mesh_cells = AUTO_LIBRARY_MESH_CELLS_BY_ACCURACY[body_dict["accuracy"]]
    l_fractions = _l_fractions_of_period(mesh_cells)
    l_values_mm = tuple(round(body_dict["period_mm"] * f, 6) for f in l_fractions)
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
            "cells_x": mesh_cells,
            "cells_y": mesh_cells,
            "floquet_order": mesh_cells,
        },
        "reference_planes": {"top_air_thickness": 5.0e-3},
    }

    def _progress(event):
        job.progress = event["done"] / max(event["total"], 1)
        job.progress_message = f"{event['done']}/{event['total']} cell sizes solved"

    cells = build_square_library_cells(base_config, l_values_mm, body_dict["frequency_hz"], progress=_progress)
    out_dir = job_output_dir(job.id)
    path = write_square_library_json(cells, out_dir / "unit_cell_library.json")

    # PHASE13D item 1: report the achievable phase-lookup range as soon as the library exists, so
    # Step 2 can show it before the user commits to a synthesis run (and Step 3 can explain a poor
    # unreachable-fraction result in terms of it).
    table = library_table_at_frequency({"cells": cells}, body_dict["frequency_hz"] / 1e9)
    density = phase_density(table)
    coverage_meta = {"phase_span_deg": density.phase_span_deg, "resolution_deg": density.median_gap_deg}

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
            **coverage_meta,
        },
    )
    return {
        "library_id": library_id,
        "n_geometries": len(cells),
        "l_values_mm": l_values_mm,
        "already_built": False,
        **coverage_meta,
    }


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
        meta = library_registry.list_libraries()[library_id]["meta"]
        result = {
            "library_id": library_id,
            "already_built": True,
            "phase_span_deg": meta.get("phase_span_deg"),
            "resolution_deg": meta.get("resolution_deg"),
        }
        job = request.app.state.jobs.submit("hmsa_ensure_library", lambda job: result)
        return job.to_dict()
    except KeyError:
        pass
    key = config_hash({"kind": "hmsa_ensure_library", **body.model_dump()})
    job = request.app.state.jobs.submit(
        "hmsa_ensure_library", _run_ensure_library_job, body.model_dump(), library_id, cache_key=key
    )
    return job.to_dict()


class EnsureDispersionRequest(BaseModel):
    period_mm: float
    eps_r: float
    tan_delta: float
    thickness_mm: float
    frequency_hz: float
    accuracy: Literal["fast", "accurate"] = "fast"


def _auto_dispersion_id(body_dict):
    mesh_cells = AUTO_LIBRARY_MESH_CELLS_BY_ACCURACY[body_dict["accuracy"]]
    payload = {
        "kind": "hmsa_auto_dispersion",
        "l_fractions": _l_fractions_of_period(mesh_cells),
        "mesh_cells": mesh_cells,
        **body_dict,
    }
    return "autodisp_" + config_hash(payload)[:20]


def _write_dispersion_csv(table, path):
    import csv as _csv

    with open(path, "w", newline="") as fh:
        writer = _csv.writer(fh)
        writer.writerow(["geometry_param", "beta_rad_per_m", "x_s_ohm"])
        for s in table.samples:
            writer.writerow([s.geometry_param, s.beta_rad_per_m, s.x_s_ohm])


def _dispersion_x_s_range(table):
    values = [s.x_s_ohm for s in table.samples if s.x_s_ohm is not None]
    return (min(values), max(values)) if values else (None, None)


def _dispersion_summary(table):
    """PHASE14 Part B: X_avg/modulation depth set automatically from the table (centre of the
    usable range, largest symmetric depth around that centre) -- the fixed "100 ohm / 0.3"
    defaults this replaces were never validated against any real board."""
    x_min, x_max = _dispersion_x_s_range(table)
    if x_min is None or x_min <= 0 or x_max <= x_min:
        return {"x_s_min_ohm": x_min, "x_s_max_ohm": x_max, "x_avg_ohm": None, "modulation_index": None}
    x_avg = (x_min + x_max) / 2.0
    modulation_index = (x_max - x_min) / (x_max + x_min)
    return {
        "x_s_min_ohm": x_min,
        "x_s_max_ohm": x_max,
        "x_avg_ohm": x_avg,
        "modulation_index": min(modulation_index, 0.99),  # reactance_from_modulation requires < 1
    }


def _run_ensure_dispersion_job(job, body_dict, dispersion_id):
    mesh_cells = AUTO_LIBRARY_MESH_CELLS_BY_ACCURACY[body_dict["accuracy"]]
    l_fractions = _l_fractions_of_period(mesh_cells)
    period_m = body_dict["period_mm"] * 1e-3
    l_values_m = [period_m * f for f in l_fractions]
    substrate = Substrate("user", body_dict["eps_r"], body_dict["tan_delta"], body_dict["thickness_mm"] * 1e-3)
    metal = Metal()
    mesh = MeshSettings(cells_x=mesh_cells, cells_y=mesh_cells, floquet_order=mesh_cells)

    done = [0]

    def unit_cell_at(l_m):
        done[0] += 1
        job.progress = done[0] / max(len(l_values_m), 1)
        job.progress_message = f"{done[0]}/{len(l_values_m)} cell sizes solved"
        return square_patch_unit_cell(period_m, substrate, metal, l_m)

    table = build_dispersion_table_from_sheet_model(unit_cell_at, l_values_m, body_dict["frequency_hz"], mesh, period_m=period_m)

    csv_path = upload_dir(dispersion_id) / "dispersion.csv"
    _write_dispersion_csv(table, csv_path)
    summary = _dispersion_summary(table)
    return {"dispersion_id": dispersion_id, "n_geometries": len(table.samples), "already_built": False, **summary}


@router.post("/ensure_dispersion")
def ensure_dispersion(body: EnsureDispersionRequest, request: Request):
    """PHASE14 Part B: the "Surface-wave (eigenmode)" design method's own auto-built unit-cell
    data, the impedance-modulated-method counterpart to `ensure_library`. Builds (or reuses, via
    the job cache and a deterministic on-disk path) a `dispersion.sheet_model` table for the
    chosen board/frequency, written as a CSV in the same schema an uploaded HFSS table would use
    (`dispersion.hfss_table.read_dispersion_csv`) -- so Formulation A's synthesis code needs no
    changes at all to consume it: it's just another `dispersion_id`."""
    dispersion_id = _auto_dispersion_id(body.model_dump())
    csv_path = upload_dir(dispersion_id) / "dispersion.csv"
    if csv_path.is_file():
        table = read_dispersion_csv(csv_path, body.frequency_hz)
        summary = _dispersion_summary(table)
        result = {"dispersion_id": dispersion_id, "already_built": True, "n_geometries": len(table.samples), **summary}
        job = request.app.state.jobs.submit("hmsa_ensure_dispersion", lambda job: result)
        return job.to_dict()
    key = config_hash({"kind": "hmsa_ensure_dispersion", **body.model_dump()})
    job = request.app.state.jobs.submit(
        "hmsa_ensure_dispersion", _run_ensure_dispersion_job, body.model_dump(), dispersion_id, cache_key=key
    )
    return job.to_dict()


# ============================================================================================
# PHASE15 Part C (gap-fill) / Part D: multi-port sub-apertures
# ============================================================================================


class PortBeamSpec(BaseModel):
    """A beam within ONE port's own beam list -- no `feed_index` (unlike `ObjectBeamSpec`): a
    port's beams are always radial about that SAME port's own feed, by definition of "sub-
    aperture"."""

    theta_deg: float
    phi_deg: float
    weight_real: float = 1.0
    weight_imag: float = 0.0
    initial_phase_deg: float = 0.0
    desired_polarization: str = "natural"
    oam_charge: int = 0


class PortSpec(BaseModel):
    beams: "list[PortBeamSpec]"


class PortLayoutRequest(BaseModel):
    """`preset`: `"2x2_plus_1_large"` (Mishra's own layout, needs `small_side_mm`),
    `"1xn"` (needs `n`/`cell_side_mm`), or `"custom"` (needs `regions_mm`, each either
    `[x_min,x_max,y_min,y_max]` or `[x_min,x_max,y_min,y_max,feed_x,feed_y]`, millimetres)."""

    preset: str
    small_side_mm: "float | None" = None
    n: "int | None" = None
    cell_side_mm: "float | None" = None
    gap_mm: float = 0.0
    regions_mm: "list[list[float]] | None" = None


class HmsaSynthesizePortsRequest(BaseModel):
    """PHASE15 Part C gap-fill: N independent sub-aperture syntheses (`hologram.
    multi_port_synthesis.synthesize_ports`) sharing one formulation/library (`MULTI_PORT_SHARED_METHOD`
    -- a documented scope reduction, see the phase doc). Every port's feed is analytic (TM0 about
    that port's own centre, `port_eps_r`/`port_thickness_mm`) -- no per-port `.fld` upload in this
    first API version."""

    formulation: str  # "A" | "B"
    frequency_hz: float
    period_mm: float
    layout: PortLayoutRequest
    ports: "list[PortSpec]"
    library_id: "str | None" = None
    dispersion_id: "str | None" = None
    x_avg_ohm: "float | None" = None
    modulation_index: "float | None" = None
    rescale_to_beta_bar: bool = True
    cp_phase_offset_deg: float = 0.0
    port_eps_r: float
    port_thickness_mm: float
    min_lambda_sw: float = 1.0


def _build_port_regions(layout: dict):
    preset = layout["preset"]
    if preset == "2x2_plus_1_large":
        if layout["small_side_mm"] is None:
            raise config_invalid(message="layout preset '2x2_plus_1_large' needs small_side_mm")
        return layout_preset_2x2_plus_1_large(layout["small_side_mm"] * 1e-3, gap_m=layout["gap_mm"] * 1e-3)
    if preset == "1xn":
        if layout["n"] is None or layout["cell_side_mm"] is None:
            raise config_invalid(message="layout preset '1xn' needs n and cell_side_mm")
        return layout_preset_1xn(layout["n"], layout["cell_side_mm"] * 1e-3, gap_m=layout["gap_mm"] * 1e-3)
    if preset == "custom":
        if not layout["regions_mm"]:
            raise config_invalid(message="layout preset 'custom' needs regions_mm")
        rects_m = [tuple(v * 1e-3 for v in r) for r in layout["regions_mm"]]
        return custom_layout(rects_m)
    raise config_invalid(message=f"unknown layout preset {preset!r}; expected '2x2_plus_1_large', '1xn', or 'custom'")


def _port_cell_to_m(formulation, period_mm):
    """Returns a `GeometryCell -> (length_m, width_m) | None` function for
    `multi_port_synthesis.combined_board_lw_grid_m`, matching Formulation A/B's own geometry-
    parameter convention (`_geometry_param_to_mm` for A, plain L_mm/W_mm for B)."""

    def cell_m(cell):
        if formulation == "B":
            if cell.geometry_param is None:
                return None
            w_mm = cell.geometry_param_w if cell.geometry_param_w is not None else cell.geometry_param
            return cell.geometry_param * 1e-3, w_mm * 1e-3
        mm = _geometry_param_to_mm(cell.geometry_param, period_mm)
        return (mm * 1e-3, mm * 1e-3) if mm is not None else None

    return cell_m


def _run_synthesize_ports_job(job, body: dict):
    regions = _build_port_regions(body["layout"])
    if len(body["ports"]) != len(regions):
        raise config_invalid(
            message=f"layout '{body['layout']['preset']}' has {len(regions)} port(s) but "
            f"{len(body['ports'])} port beam list(s) were given"
        )
    port_beams = {}
    for region, port in zip(regions, body["ports"], strict=True):
        port_beams[region.port_id] = [
            Beam(
                theta_deg=b["theta_deg"], phi_deg=b["phi_deg"],
                weight=complex(b["weight_real"], b["weight_imag"]), initial_phase_deg=b["initial_phase_deg"],
                desired_polarization=b["desired_polarization"], oam_charge=b["oam_charge"],
            )
            for b in port["beams"]
        ]

    period_m = body["period_mm"] * 1e-3
    k_swg = tm0_beta(body["frequency_hz"], body["port_eps_r"], body["port_thickness_mm"] * 1e-3)

    if body["formulation"] == "A":
        if body["dispersion_id"] is None or body["x_avg_ohm"] is None or body["modulation_index"] is None:
            raise config_invalid(message="formulation A needs dispersion_id, x_avg_ohm and modulation_index")
        dispersion_table = read_dispersion_csv(upload_dir(body["dispersion_id"]) / "dispersion.csv", body["frequency_hz"])

        def per_port(ref, beams, feed_xy_m):
            return synthesize_formulation_a(
                ref, beams, body["frequency_hz"], body["x_avg_ohm"], body["modulation_index"],
                dispersion_table, period_m, rescale_to_beta_bar=body["rescale_to_beta_bar"],
                feed_xy_m=feed_xy_m, cp_phase_offset_deg=body["cp_phase_offset_deg"],
            )
    elif body["formulation"] == "B":
        if body["library_id"] is None:
            raise config_invalid(message="formulation B needs library_id")
        lib = UnitCellLibrary.from_json(library_registry.resolve_path(body["library_id"]))
        candidates = library_candidates_at_frequency(lib, body["library_id"], body["frequency_hz"])
        if not candidates:
            raise library_frequency_unavailable(
                f"library {body['library_id']!r} has no record at {body['frequency_hz'] / 1e9:g} GHz"
            )

        def per_port(ref, beams, feed_xy_m):
            return synthesize_formulation_b(
                ref, beams, body["frequency_hz"], candidates, period_m,
                feed_xy_m=feed_xy_m, cp_phase_offset_deg=body["cp_phase_offset_deg"],
            )
    else:
        raise config_invalid(message=f"formulation must be 'A' or 'B', got {body['formulation']!r}")

    multi = synthesize_ports(
        regions, port_beams, per_port, period_m, k_swg=k_swg, min_lambda_sw=body["min_lambda_sw"]
    )

    isolation = []
    for flag in multi.isolation_flags:
        entry = {
            "port_a": flag.port_a, "port_b": flag.port_b, "gap_mm": flag.gap_m * 1e3,
            "gap_lambda_sw": flag.gap_lambda_sw, "flagged": flag.flagged,
            "guard_row_suggestion_mm": flag.guard_row_suggestion_m * 1e3 if flag.guard_row_suggestion_m else None,
        }
        isolation.append(entry)
    isolation_labels = tuple(
        isolation_edge_distance(f.port_a, f.port_b, f.gap_lambda_sw) for f in multi.isolation_flags if f.flagged
    )

    board = combined_board_lw_grid_m(multi, _port_cell_to_m(body["formulation"], body["period_mm"]))
    patch_layout_png_url = None
    try:
        cells_mm = []
        for i, row in enumerate(board.lw_grid_m):
            for j, (l_m, w_m) in enumerate(row):
                if l_m <= 0.0 or w_m <= 0.0:
                    continue
                cx, cy = board.cell_center_m(i, j)
                cells_mm.append((cx * 1e3, cy * 1e3, l_m * 1e3, w_m * 1e3))
        out_dir = job_output_dir(job.id)
        n_y, n_x = len(board.lw_grid_m), len(board.lw_grid_m[0])
        path = patch_layout_png(
            cells_mm, body["period_mm"], n_x, n_y, out_dir / "patch_layout.png",
            title=f"Multi-port board — {len(multi.ports)} ports, {n_x} × {n_y} cells", feed_xy_mm=None,
        )
        patch_layout_png_url = f"/api/files/{job.id}/{path.name}"
    except Exception:
        patch_layout_png_url = None

    ports_result = []
    for p in multi.ports:
        ports_result.append(
            {
                "port_id": p.port_id,
                "feed_x_mm": p.feed_xy_m[0] * 1e3,
                "feed_y_mm": p.feed_xy_m[1] * 1e3,
                "x_min_mm": p.region.x_min_m * 1e3, "x_max_mm": p.region.x_max_m * 1e3,
                "y_min_mm": p.region.y_min_m * 1e3, "y_max_mm": p.region.y_max_m * 1e3,
                "unreachable_fraction": p.design.unreachable_fraction,
                "labels": p.design.labels,
                "n_cells": len(p.design.cells),
                "x_m": p.design.x_m,
                "y_m": p.design.y_m,
                "aperture_field_real": p.design.aperture_field.real if p.design.aperture_field is not None else None,
                "aperture_field_imag": p.design.aperture_field.imag if p.design.aperture_field is not None else None,
            }
        )

    return {
        "ports": ports_result,
        "isolation": isolation,
        "labels": (MULTI_PORT_SHARED_METHOD, PORT_ISOLATION_NOT_FULL_WAVE, *isolation_labels),
        "formulation": body["formulation"],
        "period_mm": body["period_mm"],
        "patch_layout_png_url": patch_layout_png_url,
        "board_width_mm": board.width_m * 1e3,
        "board_height_mm": board.height_m * 1e3,
        "feeds_mm": [(p.feed_xy_m[0] * 1e3, p.feed_xy_m[1] * 1e3) for p in multi.ports],
        "lw_grid_m": board.lw_grid_m,
    }


@router.post("/synthesize_ports")
def synthesize_ports_endpoint(body: HmsaSynthesizePortsRequest, request: Request):
    if body.formulation == "B" and body.library_id is not None:
        try:
            library_registry.resolve_path(body.library_id)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc
    job = request.app.state.jobs.submit("hmsa_synthesize_ports", _run_synthesize_ports_job, body.model_dump())
    return job.to_dict()


# ============================================================================================
# PHASE15 Part E / Part D: the Coded (1-bit / N-bit) design method
# ============================================================================================


class EnsureAnnularSlotLibraryRequest(BaseModel):
    """PHASE15 Part E: a small real R-sweep (`library_sweep.build_annular_slot_library_cells`),
    `a_mm`/`b_mm`/`r_mm` fixed, `R` swept over a handful of values between `r_mm` and the largest
    radius that still fits the patch -- so the Coded method's "best-separated pair" search has a
    real, solved family of annular-slot candidates to search, the same way the square-patch
    library already does for plain patches."""

    period_mm: float
    eps_r: float
    tan_delta: float
    thickness_mm: float
    frequency_hz: float
    a_mm: float
    b_mm: float
    r_mm: float
    accuracy: Literal["fast", "accurate"] = "fast"


_ANNULAR_SLOT_R_SWEEP_POINTS = 8


def _annular_slot_r_values_mm(a_mm, b_mm, r_mm):
    r_hard_max = min(a_mm, b_mm) / 2.0
    r_max = r_hard_max * 0.95  # a small safety margin inside AnnularSlotPatch's own R<=min(a,b)/2
    r_min = min(r_mm * 1.3, r_max - 1e-6)
    if r_min >= r_max:
        raise config_invalid(message=f"r_mm={r_mm:g} leaves no room for an outer radius R below {r_hard_max:g}")
    n = _ANNULAR_SLOT_R_SWEEP_POINTS
    return [round(r_min + i * (r_max - r_min) / (n - 1), 6) for i in range(n)]


def _auto_annular_slot_library_id(body_dict):
    payload = {"kind": "hmsa_auto_annular_slot_library", **body_dict}
    return "auto_annular_" + config_hash(payload)[:20]


def _run_ensure_annular_slot_library_job(job, body_dict, library_id):
    mesh_cells = AUTO_LIBRARY_MESH_CELLS_BY_ACCURACY[body_dict["accuracy"]]
    r_values_mm = _annular_slot_r_values_mm(body_dict["a_mm"], body_dict["b_mm"], body_dict["r_mm"])
    base_config = {
        "period_x": body_dict["period_mm"] * 1e-3,
        "period_y": body_dict["period_mm"] * 1e-3,
        "substrate": {
            "name": "user", "eps_r": body_dict["eps_r"], "tan_delta": body_dict["tan_delta"],
            "thickness": body_dict["thickness_mm"] * 1e-3, "backing": "pec",
        },
        "geometry": {
            "type": "patch_annular_slot", "a": body_dict["a_mm"] * 1e-3, "b": body_dict["b_mm"] * 1e-3,
            "R": r_values_mm[0] * 1e-3, "r": body_dict["r_mm"] * 1e-3,
        },  # placeholder, overwritten per R
        "excitation": {"frequency": body_dict["frequency_hz"]},  # placeholder
        "mesh": {"cells_x": mesh_cells, "cells_y": mesh_cells, "floquet_order": mesh_cells},
        "reference_planes": {"top_air_thickness": 5.0e-3},
    }

    def _progress(event):
        job.progress = event["done"] / max(event["total"], 1)
        job.progress_message = f"{event['done']}/{event['total']} slot radii solved"

    cells = build_annular_slot_library_cells(
        base_config, body_dict["a_mm"], body_dict["b_mm"], body_dict["r_mm"], r_values_mm,
        body_dict["frequency_hz"], progress=_progress,
    )
    out_dir = job_output_dir(job.id)
    path = write_annular_slot_library_json(cells, out_dir / "annular_slot_library.json")
    library_registry.register(
        library_id, path,
        meta={
            "kind": "annular_slot", "period_mm": body_dict["period_mm"], "backing": "pec",
            "a_mm": body_dict["a_mm"], "b_mm": body_dict["b_mm"], "r_mm": body_dict["r_mm"],
            "n_geometries": len(cells), "auto_built": True, "frequency_hz": body_dict["frequency_hz"],
        },
    )
    return {"annular_slot_library_id": library_id, "n_geometries": len(cells), "r_values_mm": r_values_mm, "already_built": False}


@router.post("/ensure_annular_slot_library")
def ensure_annular_slot_library(body: EnsureAnnularSlotLibraryRequest, request: Request):
    library_id = _auto_annular_slot_library_id(body.model_dump())
    try:
        library_registry.resolve_path(library_id)
        job = request.app.state.jobs.submit(
            "hmsa_ensure_annular_slot_library", lambda job: {"annular_slot_library_id": library_id, "already_built": True}
        )
        return job.to_dict()
    except KeyError:
        pass
    key = config_hash({"kind": "hmsa_ensure_annular_slot_library", **body.model_dump()})
    job = request.app.state.jobs.submit(
        "hmsa_ensure_annular_slot_library", _run_ensure_annular_slot_library_job, body.model_dump(),
        library_id, cache_key=key,
    )
    return job.to_dict()


class TwoStateSourceRequest(BaseModel):
    """PHASE15 Part E's three two-state-cell options: (a) `"square_only"` -- the best-separated
    pair within one square-patch library (`library_id`); (b) `"square_and_annular"` -- the
    best-separated pair across a square-patch library AND an annular-slot one
    (`annular_slot_library_id`, from `/ensure_annular_slot_library`) -- Mishra's own full-patch-vs-
    annular-slot pairing; (c) `"uploaded"` -- two HFSS-measured complex values given directly."""

    source: str = "square_and_annular"
    annular_slot_library_id: "str | None" = None
    uploaded_label_0: str = "state 0"
    uploaded_value_0_real: "float | None" = None
    uploaded_value_0_imag: "float | None" = None
    uploaded_label_1: str = "state 1"
    uploaded_value_1_real: "float | None" = None
    uploaded_value_1_imag: "float | None" = None


def _resolve_two_state_choice(two_state: dict, library_id, frequency_hz):
    source = two_state["source"]
    if source == "uploaded":
        for key in ("uploaded_value_0_real", "uploaded_value_0_imag", "uploaded_value_1_real", "uploaded_value_1_imag"):
            if two_state[key] is None:
                raise config_invalid(message=f"two_state source 'uploaded' needs {key}")
        return two_state_from_values(
            two_state["uploaded_label_0"], complex(two_state["uploaded_value_0_real"], two_state["uploaded_value_0_imag"]),
            two_state["uploaded_label_1"], complex(two_state["uploaded_value_1_real"], two_state["uploaded_value_1_imag"]),
        )
    if library_id is None:
        raise config_invalid(message=f"two_state source {source!r} needs library_id (a square-patch library)")
    lib = UnitCellLibrary.from_json(library_registry.resolve_path(library_id))
    candidates = candidates_from_square_library(library_candidates_at_frequency(lib, library_id, frequency_hz))
    if source == "square_only":
        return best_separated_pair(candidates)
    if source == "square_and_annular":
        if two_state["annular_slot_library_id"] is None:
            raise config_invalid(message="two_state source 'square_and_annular' needs annular_slot_library_id")
        annular_path = library_registry.resolve_path(two_state["annular_slot_library_id"])
        annular_cells = json.loads(annular_path.read_text())["cells"]
        annular_candidates = candidates_from_annular_slot_library(annular_cells, frequency_hz / 1e9)
        return best_separated_pair(candidates + annular_candidates)
    raise config_invalid(message=f"unknown two_state source {source!r}")


class HmsaSynthesizeCodedRequest(BaseModel):
    """PHASE15 Part E: the Coded design method. Reuses the single-aperture reference-field
    resolution (`upload_id` or `feeds`) and beam spec exactly as `/synthesize` does -- only the
    control-map/geometry-mapping step differs (`hologram.coded_synthesis`)."""

    upload_id: "str | None" = None
    feeds: "list[FeedSpecRequest] | None" = None
    multi_feed_strategy: str = "superposition"
    multi_feed_eps_r: "float | None" = None
    multi_feed_thickness_mm: "float | None" = None
    frequency_hz: float
    period_mm: float
    object_beams: "list[ObjectBeamSpec]"
    n_bits: int = 1
    library_id: "str | None" = None  # required for n_bits > 1; also used by two_state sources
    two_state: TwoStateSourceRequest = TwoStateSourceRequest()
    cp_phase_offset_deg: float = 0.0
    feed_x_mm: "float | None" = None
    feed_y_mm: "float | None" = None
    # PHASE15 Part E: substrate params for an immediate quantization-loss measurement (Part B's own
    # vector model) against the continuous (unquantized) design -- optional; omitted means "not
    # computed" (never a fabricated/assumed loss number).
    eps_r: "float | None" = None
    tan_delta: "float | None" = None
    thickness_mm: "float | None" = None


def _two_state_geometry_to_cell_m(geometry):
    """A `two_state_cells.TwoStateCell.geometry` dict (millimetres) -> the `(length_m, width_m)` /
    annular-slot-dict cell entry `export.array_export`'s mixed-cell writers expect (metres), or
    `None` if this state carries no exportable geometry at all (the "uploaded" two-state source,
    `two_state_cells.two_state_from_values` with no `geometry_0`/`geometry_1` given -- an
    HFSS-measured reflection coefficient with no physical cell shape recorded in this app; Export
    then correctly refuses with "no geometry" rather than this function fabricating one)."""
    if geometry.get("kind") == "annular_slot":
        return {
            "kind": "annular_slot", "a_m": geometry["a_mm"] * 1e-3, "b_m": geometry["b_mm"] * 1e-3,
            "R_m": geometry["R_mm"] * 1e-3, "r_m": geometry["r_mm"] * 1e-3,
        }
    if "L_mm" in geometry and "W_mm" in geometry:
        return (geometry["L_mm"] * 1e-3, geometry["W_mm"] * 1e-3)
    return None


def _run_synthesize_coded_job(job, body: dict):
    extra_labels = ("CODED_DESIGN_METHOD",)
    feeds_mm = None
    if body["feeds"] is not None:
        reference_field, feed_xy_m, feeds_mm, multi_feed_labels = _resolve_multi_feed(body, body["frequency_hz"])
        extra_labels = (*extra_labels, *multi_feed_labels)
    else:
        if body["upload_id"] is None:
            raise config_invalid(message="synthesize_coded needs either upload_id or feeds")
        raw_reference_field = _load_reference_field(body["upload_id"])
        raw_nan_mask = np.isnan(raw_reference_field.phase())
        reference_field, _masked = _mask_feed_region(raw_reference_field)
        x0_m, y0_m, feed_label = _resolve_feed_xy_m(body, raw_reference_field, raw_nan_mask)
        feed_xy_m = (x0_m, y0_m)
        feeds_mm = [(x0_m * 1e3, y0_m * 1e3)]
        if feed_label is not None:
            extra_labels = (*extra_labels, feed_label)

    beams = [
        Beam(
            theta_deg=b["theta_deg"], phi_deg=b["phi_deg"], weight=complex(b["weight_real"], b["weight_imag"]),
            initial_phase_deg=b["initial_phase_deg"], desired_polarization=b["desired_polarization"],
            oam_charge=b["oam_charge"],
        )
        for b in body["object_beams"]
    ]
    period_m = body["period_mm"] * 1e-3
    n_bits = body["n_bits"]

    if n_bits == 1:
        choice = _resolve_two_state_choice(body["two_state"], body["library_id"], body["frequency_hz"])
        design = synthesize_coded_one_bit(
            reference_field, beams, body["frequency_hz"], choice, period_m,
            feed_xy_m=feed_xy_m, cp_phase_offset_deg=body["cp_phase_offset_deg"],
        )
        warning = separation_warning(choice)
        state_labels = (two_state_phase_separation(choice.achieved_separation_deg, n_bits=1),)
        two_state_info = {
            "state_0": {"label": choice.state_0.label, "geometry": choice.state_0.geometry,
                        "phase_deg": float(np.degrees(np.angle(choice.state_0.value)))},
            "state_1": {"label": choice.state_1.label, "geometry": choice.state_1.geometry,
                        "phase_deg": float(np.degrees(np.angle(choice.state_1.value)))},
            "achieved_separation_deg": choice.achieved_separation_deg,
            "warning": warning,
        }
        ny, nx = design.x_m.shape
        lw_grid_m = [[None] * nx for _ in range(ny)]
        for c in design.cells:
            state = choice.state_0 if int(c.geometry_param) == 0 else choice.state_1
            lw_grid_m[c.i][c.j] = _two_state_geometry_to_cell_m(state.geometry)
    else:
        if body["library_id"] is None:
            raise config_invalid(message="n_bits > 1 needs library_id (a continuous phase library)")
        lib = UnitCellLibrary.from_json(library_registry.resolve_path(body["library_id"]))
        candidates = library_candidates_at_frequency(lib, body["library_id"], body["frequency_hz"])
        if not candidates:
            raise library_frequency_unavailable(
                f"library {body['library_id']!r} has no record at {body['frequency_hz'] / 1e9:g} GHz"
            )
        design = synthesize_coded_n_bit(
            reference_field, beams, body["frequency_hz"], n_bits, candidates, period_m,
            feed_xy_m=feed_xy_m, cp_phase_offset_deg=body["cp_phase_offset_deg"],
        )
        state_labels = ()
        two_state_info = None
        lw_grid_m = None  # Export for N>1 reuses the existing Formulation-B-shaped `assignments`

    assignments = None
    if n_bits > 1:
        assignments = tuple(
            CellAssignment(
                i=c.i, j=c.j, x_m=c.x_m, y_m=c.y_m, required_phase_deg=c.control_value,
                chosen_l_mm=c.geometry_param, chosen_w_mm=c.geometry_param_w,
                achieved_phase_deg=c.achieved_value, phase_error_deg=c.phase_error_deg,
                library_id=body["library_id"], phase_coverage_ok=c.reachable,
            )
            for c in design.cells
        )

    cell_field = resample_to_cell_grid(reference_field, period_m)
    continuous_field = continuous_reference_aperture_field(cell_field, design.continuous_phase_deg)

    quantization_loss_label = None
    if body["eps_r"] is not None and body["tan_delta"] is not None and body["thickness_mm"] is not None:
        # Part E: "report quantization loss against the continuous design" -- a REAL measurement
        # via the same vector pattern model Part B already validated, not a textbook bound.
        try:
            substrate = Substrate("user", body["eps_r"], body["tan_delta"], body["thickness_mm"] * 1e-3)
            green = LayeredGreenFunction(substrate, body["frequency_hz"], backing="pec")
            dx = float(design.x_m[0, 1] - design.x_m[0, 0])
            dy = float(design.y_m[1, 0] - design.y_m[0, 0])
            jx_q, jy_q = aperture_current_from_scalar(design.x_m, design.y_m, design.aperture_field, feed_xy_m=feed_xy_m)
            jx_c, jy_c = aperture_current_from_scalar(design.x_m, design.y_m, continuous_field, feed_xy_m=feed_xy_m)
            pat_q = compute_vector_pattern(design.x_m, design.y_m, jx_q, jy_q, body["frequency_hz"], green, dx, dy)
            pat_c = compute_vector_pattern(design.x_m, design.y_m, jx_c, jy_c, body["frequency_hz"], green, dx, dy)
            loss_db = pat_c.directivity_estimate_dbi - pat_q.directivity_estimate_dbi
            quantization_loss_label = quantization_loss_db(loss_db, n_bits)
        except Exception:
            quantization_loss_label = None

    patch_layout_png_url = None
    try:
        out_dir = job_output_dir(job.id)
        ny, nx = design.x_m.shape
        if n_bits == 1:
            cells_mm = []
            for c in design.cells:
                state = choice.state_0 if int(c.geometry_param) == 0 else choice.state_1
                geom = state.geometry
                if geom.get("kind") == "annular_slot":
                    cells_mm.append((c.x_m * 1e3, c.y_m * 1e3, geom["a_mm"], geom["b_mm"], geom["R_mm"], geom["r_mm"]))
                else:
                    cells_mm.append((c.x_m * 1e3, c.y_m * 1e3, geom.get("L_mm"), geom.get("W_mm")))
        else:
            cells_mm = [
                (c.x_m * 1e3, c.y_m * 1e3, c.geometry_param, c.geometry_param_w or c.geometry_param)
                for c in design.cells
            ]
        path = patch_layout_png(
            cells_mm, body["period_mm"], nx, ny, out_dir / "patch_layout.png",
            title=f"Coded design — {nx} × {ny} cells, {n_bits}-bit", feed_xy_mm=(feed_xy_m[0] * 1e3, feed_xy_m[1] * 1e3),
        )
        patch_layout_png_url = f"/api/files/{job.id}/{path.name}"
    except Exception:
        patch_layout_png_url = None

    return {
        "n_bits": n_bits,
        "x_m": design.x_m,
        "y_m": design.y_m,
        "period_mm": body["period_mm"],
        "continuous_phase_deg": design.continuous_phase_deg,
        "quantized_phase_deg": design.quantized_phase_deg,
        "level_map": design.level_map,
        "cells": design.cells,
        "assignments": assignments,
        "two_state": two_state_info,
        "labels": extra_labels + state_labels + ((quantization_loss_label,) if quantization_loss_label else ()),
        "aperture_field_real": design.aperture_field.real,
        "aperture_field_imag": design.aperture_field.imag,
        "continuous_aperture_field_real": continuous_field.real if continuous_field is not None else None,
        "continuous_aperture_field_imag": continuous_field.imag if continuous_field is not None else None,
        "patch_layout_png_url": patch_layout_png_url,
        "feeds_mm": feeds_mm,
        "lw_grid_m": lw_grid_m,  # n_bits==1 only -- the two-state export path reads this directly
        # (1-bit states may be a rectangle OR an annular-slot dict; `assignments` above (n_bits>1
        # only) is rectangle-only, the existing reflectarray-shaped export path)
    }


@router.post("/synthesize_coded")
def synthesize_coded_endpoint(body: HmsaSynthesizeCodedRequest, request: Request):
    if body.library_id is not None:
        try:
            library_registry.resolve_path(body.library_id)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc
    job = request.app.state.jobs.submit("hmsa_synthesize_coded", _run_synthesize_coded_job, body.model_dump())
    return job.to_dict()
