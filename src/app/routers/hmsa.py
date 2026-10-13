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
from ...hologram.auto_states import auto_states
from ...hologram.coded_synthesis import (
    apply_feed_keepout_one_bit,
    aperture_coupling_loss_db,
    continuous_reference_aperture_field,
    synthesize_coded_n_bit,
    synthesize_coded_one_bit,
    textbook_quantization_loss_db,
)
from ...hologram.hmsa_synthesis import synthesize_formulation_a, synthesize_formulation_b
from ...hologram.multi_feed import (
    FeedSpec,
    feed_xy_m_for_beams,
    segmented_reference_field_and_regions,
    superposed_reference_field,
)
from ...hologram.multi_port_synthesis import (
    combined_board_lw_grid_m,
    export_frame_feeds_and_rects_mm,
    export_frame_offset_m,
    mask_design_to_region,
    region_cell_mask,
    synthesize_ports,
)
from ...hologram.pointing_correction import measure_beam_aims, refine_design
from ...hologram.ports import (
    check_no_overlap,
    custom_layout,
    layout_preset_1xn,
    layout_preset_2x2_plus_1_large,
    layout_preset_grid,
    layout_preset_hex4,
    snap_side_to_cells,
)
from ...hologram.two_state_cells import (
    BETA_SPREAD_LIMIT_FRACTION,
    EDGE_PHASE_WARNING_DEG,
    LIBRARY_REFERENCE_OFFSET_M,
    MAGNITUDE_MISMATCH_WARNING_DB,
    TwoStateCell,
    best_valid_pair,
    candidates_from_annular_slot_library,
    candidates_from_square_library,
    lambda_sw_over_cell,
    mixed_surface_mode,
    patch_plane_gamma,
    resonance_proximity_patch_plane_deg,
    screen_candidates,
    separation_warning,
)
from ...labels import (
    EXPERIMENTAL,
    FEED_REGION_MASKED,
    MULTI_PORT_SHARED_METHOD,
    PORT_ISOLATION_NOT_FULL_WAVE,
    feed_position_auto_detected,
    isolation_edge_distance,
    quantization_loss_db,
    reference_field_frequency_mismatch,
    two_state_phase_separation,
    two_state_surface_wave_check,
)
from ...library_sweep import (
    build_annular_slot_library_cells,
    build_square_library_cells,
    write_annular_slot_library_json,
    write_square_library_json,
)
from ...materials import Metal, Substrate
from ...pattern.vector_aperture_pattern import (
    ORTHOGONAL_POLARIZATION,
    aperture_current_from_scalar,
    aperture_efficiency_estimate,
    beam_local_peak_fine,
    beam_polarization_metrics,
    compute_vector_pattern,
    copol_component,
    first_null_half_width_deg,
    mirror_beam_level_db,
    strongest_secondary_beam,
)
from ...phase import principal_phase_deg, wrapped_difference_deg
from ...phase8b.object_wave import Beam
from ...reflectarray.designer import CellAssignment
from ...unit_cell_library import UnitCellLibrary
from .. import library_registry
from ..jobs import config_hash
from ..jsonutil import to_jsonable
from ..library_query import available_frequencies_hz, library_candidates_at_frequency
from ..plots import layout_map_png, patch_layout_png, reference_phase_map_png
from ..storage import app_data_dir, job_output_dir, upload_dir

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
    # The substrate of the pattern check used by the optional stages below (needed for both).
    eps_r: "float | None" = None
    tan_delta: "float | None" = None
    thickness_mm: "float | None" = None
    # Optional stages (off / on by default as noted): pointing correction, and matching the requested
    # relative beam levels for 2+ beams (on by default).
    pointing_correction: bool = False
    pointing_iterations: int = 2
    balance_levels: bool = True
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
            FeedSpec(
                x_m=f["x_mm"] * 1e-3,
                y_m=f["y_mm"] * 1e-3,
                amplitude=complex(f["amplitude_real"], f["amplitude_imag"]),
                reference_field=rf,
            )
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
        raise config_invalid(
            message=f"multi_feed_strategy must be 'superposition' or 'segmentation', got {body['multi_feed_strategy']!r}"
        )

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


def _refined_design(body, beams, synthesize, feed_xy_m):
    """One single-aperture design through the optional stages (`hologram.pointing_correction`).
    Returns `(design, refinement)`; `refinement` is None when no stage applies."""
    needs_pattern = body["pointing_correction"] or (body["balance_levels"] and len(beams) > 1)
    if not needs_pattern:
        return synthesize(beams), None
    if body["eps_r"] is None or body["tan_delta"] is None or body["thickness_mm"] is None:
        if body["pointing_correction"]:
            raise config_invalid(
                message="pointing_correction needs eps_r, tan_delta and thickness_mm (the pattern check needs them)"
            )
        return synthesize(beams), {"balance": {"skipped": "needs eps_r, tan_delta and thickness_mm"}}
    substrate = Substrate("user", body["eps_r"], body["tan_delta"], body["thickness_mm"] * 1e-3)
    green = LayeredGreenFunction(substrate, body["frequency_hz"], backing="pec")
    period_m = body["period_mm"] * 1e-3

    def measure(design, beams_now):
        return measure_beam_aims(design, beams_now, green, period_m, body["frequency_hz"], feed_xy_m)

    design, _used, report = refine_design(
        beams,
        synthesize,
        measure,
        pointing=body["pointing_correction"],
        balance=body["balance_levels"],
        pointing_iterations=body["pointing_iterations"],
    )
    return design, report


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
            extra_labels = (
                *extra_labels,
                reference_field_frequency_mismatch(upload_frequency_hz, body["frequency_hz"]),
            )
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
        design, refinement = _refined_design(
            body,
            beams,
            lambda bs: synthesize_formulation_a(
                reference_field,
                bs,
                body["frequency_hz"],
                body["x_avg_ohm"],
                body["modulation_index"],
                dispersion_table,
                period_m,
                rescale_to_beta_bar=body["rescale_to_beta_bar"],
                feed_xy_m=feed_xy_m,
                cp_phase_offset_deg=body["cp_phase_offset_deg"],
                multiplexing=body["multiplexing"],
            ),
            feed_xy_m,
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
                    i=c.i,
                    j=c.j,
                    x_m=c.x_m,
                    y_m=c.y_m,
                    required_phase_deg=c.control_value,
                    chosen_l_mm=v,
                    chosen_w_mm=v,
                    achieved_phase_deg=c.achieved_value,
                    phase_error_deg=c.phase_error_deg,
                    library_id=body["dispersion_id"],
                    phase_coverage_ok=c.reachable,
                )
                for c, v in zip(design.cells, sizes_mm, strict=True)
            )
        patch_layout_png_url = (
            _patch_layout_url(job, assignments_a, design.x_m, body["period_mm"]) if assignments_a else None
        )
        return {
            "formulation": design.formulation,
            "refinement": refinement,
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
        design, refinement = _refined_design(
            body,
            beams,
            lambda bs: synthesize_formulation_b(
                reference_field,
                bs,
                body["frequency_hz"],
                candidates,
                period_m,
                feed_xy_m=feed_xy_m,
                cp_phase_offset_deg=body["cp_phase_offset_deg"],
                multiplexing=body["multiplexing"],
            ),
            feed_xy_m,
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
            "refinement": refinement,
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


def _square_library_base_config(body):
    mesh_cells = AUTO_LIBRARY_MESH_CELLS_BY_ACCURACY[body["accuracy"]]
    l_values_mm = tuple(round(body["period_mm"] * f, 6) for f in _l_fractions_of_period(mesh_cells))
    return {
        "period_x": body["period_mm"] * 1e-3,
        "period_y": body["period_mm"] * 1e-3,
        "substrate": {
            "name": "user",
            "eps_r": body["eps_r"],
            "tan_delta": body["tan_delta"],
            "thickness": body["thickness_mm"] * 1e-3,
            "backing": "pec",
        },
        "geometry": {"type": "rectangle", "length": l_values_mm[0], "width": l_values_mm[0]},  # placeholder
        "excitation": {"frequency": body["frequency_hz"]},  # placeholder
        "mesh": {
            "cells_x": mesh_cells,
            "cells_y": mesh_cells,
            "floquet_order": mesh_cells,
        },
        "reference_planes": {"top_air_thickness": 5.0e-3},
    }


def _annular_library_base_config(body, placeholder_R_mm):
    mesh_cells = AUTO_LIBRARY_MESH_CELLS_BY_ACCURACY[body["accuracy"]]
    return {
        "period_x": body["period_mm"] * 1e-3,
        "period_y": body["period_mm"] * 1e-3,
        "substrate": {
            "name": "user",
            "eps_r": body["eps_r"],
            "tan_delta": body["tan_delta"],
            "thickness": body["thickness_mm"] * 1e-3,
            "backing": "pec",
        },
        "geometry": {
            "type": "patch_annular_slot",
            "a": body["a_mm"] * 1e-3,
            "b": body["b_mm"] * 1e-3,
            "R": placeholder_R_mm * 1e-3,
            "r": body["r_mm"] * 1e-3,
        },
        "excitation": {"frequency": body["frequency_hz"]},  # placeholder
        "mesh": {"cells_x": mesh_cells, "cells_y": mesh_cells, "floquet_order": mesh_cells},
        "reference_planes": {"top_air_thickness": 5.0e-3},
    }


def _run_ensure_library_job(job, body_dict, library_id):
    mesh_cells = AUTO_LIBRARY_MESH_CELLS_BY_ACCURACY[body_dict["accuracy"]]
    l_fractions = _l_fractions_of_period(mesh_cells)
    l_values_mm = tuple(round(body_dict["period_mm"] * f, 6) for f in l_fractions)
    base_config = _square_library_base_config(body_dict)

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

    table = build_dispersion_table_from_sheet_model(
        unit_cell_at, l_values_m, body_dict["frequency_hz"], mesh, period_m=period_m
    )

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


class TwoStateSourceRequest(BaseModel):
    """PHASE15 Part E's three two-state-cell options: (a) `"square_only"` -- the best valid pair
    within one square-patch library (`library_id`); (b) `"square_and_annular"` -- the best valid
    pair across a square-patch library AND an annular-slot one (`annular_slot_library_id`, from
    `/ensure_annular_slot_library`) -- Mishra's own full-patch-vs-annular-slot pairing; (c)
    `"uploaded"` -- two HFSS-measured complex values given directly (taken at the patch plane)."""

    source: str = "square_and_annular"
    annular_slot_library_id: "str | None" = None
    uploaded_label_0: str = "state 0"
    uploaded_value_0_real: "float | None" = None
    uploaded_value_0_imag: "float | None" = None
    uploaded_label_1: str = "state 1"
    uploaded_value_1_real: "float | None" = None
    uploaded_value_1_imag: "float | None" = None
    uploaded_geometry_0: "dict | None" = None  # optional cell geometry (two_state_cells' own dict format, mm)
    uploaded_geometry_1: "dict | None" = None
    # Phase 22: source "auto" finds the cell states itself (empty + squares, MoM, mesh-corrected, snapped to
    # a 0.05 mm grid) from the substrate, cell size, frequency and number of bits. `auto_fixed_L_mm` (1-bit)
    # computes that square's Gamma instead of searching.
    auto_fixed_L_mm: "float | None" = None
    # smallest gap left between neighbouring patches (a fabrication limit): a square is at most P - gap wide
    auto_min_gap_mm: float = 0.2


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
    `"1xn"` (needs `n`/`cell_side_mm`), `"grid"` (needs `rows`, `cols`, `cell_side_mm`: an N x M grid of
    equal sub-apertures, serpentine port order from the top-right), or `"custom"` (needs `regions_mm`,
    each either `[x_min,x_max,y_min,y_max]` or `[x_min,x_max,y_min,y_max,feed_x,feed_y]`, millimetres).
    `snap_side_to_cells` rounds the preset's side to a whole number of cells (off by default here so
    existing requests are unchanged; the app sends it on)."""

    preset: str
    small_side_mm: "float | None" = None
    n: "int | None" = None
    rows: "int | None" = None
    cols: "int | None" = None
    cell_side_mm: "float | None" = None
    gap_mm: float = 0.0
    regions_mm: "list[list[float]] | None" = None
    snap_side_to_cells: bool = False
    # Phase 23, preset "hex4": four regular hexagons (two offset rows). `side_mm` = hexagon side =
    # circumradius, or `flat_to_flat_mm` instead; `gap_mm` is shared with the other presets;
    # `orientation` "pointy" (two rows) or "flat" (two columns); `trim_board` cuts the exported
    # substrate/ground to the hexagon outline instead of the cluster's bounding rectangle.
    side_mm: "float | None" = None
    flat_to_flat_mm: "float | None" = None
    orientation: str = "pointy"
    trim_board: bool = False


class HmsaSynthesizePortsRequest(BaseModel):
    """PHASE15 Part C gap-fill: N independent sub-aperture syntheses (`hologram.
    multi_port_synthesis.synthesize_ports`) sharing one formulation/library (`MULTI_PORT_SHARED_METHOD`
    -- a documented scope reduction, see the phase doc). Every port's feed is analytic (TM0 about
    that port's own centre, `port_eps_r`/`port_thickness_mm`) -- no per-port `.fld` upload in this
    first API version."""

    formulation: str  # "A" | "B" | "coded"
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
    # "coded" only: the substrate's loss tangent for the surface-wave screen, and the two-state source
    # (shared by every port -- each port is single-feed, so the one-feed rule of /synthesize_coded holds)
    port_tan_delta: "float | None" = None
    two_state: TwoStateSourceRequest = TwoStateSourceRequest()
    # "coded" only: the surface-wave beta the phase is built with. "state": each state's own solved
    # beta (pair must agree within BETA_SPREAD_LIMIT_FRACTION). "model": beta of the mixed 1-bit surface.
    # ("paper", Mishra eq. 1, was removed in Phase 22.)
    beta_mode: str = "state"
    min_lambda_sw: float = 1.0
    # Optional pointing correction (off by default): re-aim each beam by its measured miss, 1-2 passes.
    pointing_correction: bool = False
    pointing_iterations: int = 2
    # Match the requested relative beam levels (2+ beams per port), on by default.
    balance_levels: bool = True
    # "coded" only: 1 = two states (`two_state`), >1 = 2**n_bits phase levels matched against the
    # continuous library `library_id` (as /synthesize_coded does for a single aperture).
    n_bits: int = 1
    # "coded" 1-bit only: force the cell(s) under each feed to state 0 (empty) -- the coax passes
    # through the board there. The number of cells changed is reported per port.
    feed_keepout: bool = True


def _snap_layout_sides(layout: dict, period_mm):
    """Returns `(layout, notes)`: with `snap_side_to_cells`, the preset's side rounded to the nearest
    whole number of cells (`notes`: one `{requested_mm, snapped_mm, n_cells}` per snapped side)."""
    if not layout.get("snap_side_to_cells") or period_mm is None or layout["preset"] == "custom":
        return layout, []
    key = "small_side_mm" if layout["preset"] == "2x2_plus_1_large" else "cell_side_mm"
    if layout.get(key) is None:
        return layout, []
    snapped_m, n_cells = snap_side_to_cells(layout[key] * 1e-3, period_mm * 1e-3)
    note = {"requested_mm": layout[key], "snapped_mm": snapped_m * 1e3, "n_cells": n_cells}
    return {**layout, key: snapped_m * 1e3}, [note]


def _build_port_regions(layout: dict, period_mm=None):
    layout, _notes = _snap_layout_sides(layout, period_mm)
    preset = layout["preset"]
    if preset == "2x2_plus_1_large":
        if layout["small_side_mm"] is None:
            raise config_invalid(message="layout preset '2x2_plus_1_large' needs small_side_mm")
        return layout_preset_2x2_plus_1_large(layout["small_side_mm"] * 1e-3, gap_m=layout["gap_mm"] * 1e-3)
    if preset == "1xn":
        if layout["n"] is None or layout["cell_side_mm"] is None:
            raise config_invalid(message="layout preset '1xn' needs n and cell_side_mm")
        return layout_preset_1xn(layout["n"], layout["cell_side_mm"] * 1e-3, gap_m=layout["gap_mm"] * 1e-3)
    if preset == "grid":
        if layout["rows"] is None or layout["cols"] is None or layout["cell_side_mm"] is None:
            raise config_invalid(message="layout preset 'grid' needs rows, cols and cell_side_mm")
        try:
            return layout_preset_grid(
                layout["rows"], layout["cols"], layout["cell_side_mm"] * 1e-3, gap_m=layout["gap_mm"] * 1e-3
            )
        except ValueError as exc:
            raise config_invalid(message=str(exc)) from exc
    if preset == "hex4":
        if (layout["side_mm"] is None) == (layout["flat_to_flat_mm"] is None):
            raise config_invalid(message="layout preset 'hex4' needs exactly one of side_mm and flat_to_flat_mm")
        try:
            return layout_preset_hex4(
                side_m=None if layout["side_mm"] is None else layout["side_mm"] * 1e-3,
                flat_to_flat_m=None if layout["flat_to_flat_mm"] is None else layout["flat_to_flat_mm"] * 1e-3,
                gap_m=layout["gap_mm"] * 1e-3,
                orientation=layout["orientation"],
                period_m=None if period_mm is None else period_mm * 1e-3,
            )
        except ValueError as exc:
            raise config_invalid(message=str(exc)) from exc
    if preset == "custom":
        if not layout["regions_mm"]:
            raise config_invalid(message="layout preset 'custom' needs regions_mm")
        rects_m = [tuple(v * 1e-3 for v in r) for r in layout["regions_mm"]]
        try:
            regions = custom_layout(rects_m)
            check_no_overlap(regions)
        except ValueError as exc:
            raise config_invalid(message=f"custom layout: {exc}") from exc
        return regions
    raise config_invalid(
        message=f"unknown layout preset {preset!r}; expected '2x2_plus_1_large', '1xn', 'grid', 'hex4', or 'custom'"
    )


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


def _geometry_mm(geometry_m):
    """One `cell_geometry_m` entry (metres) -> mm for the result: `[L, W]`, an annular-slot dict, or
    `None` for "no patch"."""
    if geometry_m is None:
        return None
    if isinstance(geometry_m, dict):
        return {k: (v if k == "kind" else round(v * 1e3, 9)) for k, v in geometry_m.items()}
    length_m, width_m = geometry_m
    if length_m <= 0.0 or width_m <= 0.0:
        return None
    return [round(length_m * 1e3, 9), round(width_m * 1e3, 9)]  # m -> mm without float residue


def _port_patches_mm(port, cell_geometry_m, period_mm):
    """Phase 21 C3: this port's own cells at their OWN exact positions, for export -- no shared
    lattice. `cells` lists only cells that carry a patch, as `{i, j, geometry}` (geometry in mm);
    the cell centre is `x_mid + (j - (nx-1)/2)*period`, `y_mid + (i - (ny-1)/2)*period`, with
    `(x_mid, y_mid)` the sub-aperture's centre (the design grid is centred in its rectangle)."""
    ny, nx = port.design.x_m.shape
    cells = []
    for c in port.design.cells:
        geometry = _geometry_mm(cell_geometry_m(c))
        if geometry is not None:
            cells.append({"i": c.i, "j": c.j, "geometry": geometry})
    region = port.region
    layout = {
        "nx": nx,
        "ny": ny,
        "period_mm": period_mm,
        "x_mid_mm": (region.x_min_m + region.x_max_m) / 2.0 * 1e3,
        "y_mid_mm": (region.y_min_m + region.y_max_m) / 2.0 * 1e3,
        "cells": cells,
    }
    return layout, len(cells)


def _ports_reports(multi, port_beams, body, cell_geometry_m):
    """Phase 21 B3: per port, per requested beam: where the beam actually points, the pointing error,
    its directivity and its cross-pol, from the port's own vector pattern. `{}` when the loss tangent
    is unknown (the pattern needs the substrate)."""
    if body["port_tan_delta"] is None:
        return {}
    f = body["frequency_hz"]
    period_m = body["period_mm"] * 1e-3
    substrate = Substrate("user", body["port_eps_r"], body["port_tan_delta"], body["port_thickness_mm"] * 1e-3)
    green = LayeredGreenFunction(substrate, f, backing="pec")
    wavelength_m = 299792458.0 / f
    reports = {}
    for p in multi.ports:
        design = p.design
        x_m, y_m = design.x_m, design.y_m
        jx, jy = aperture_current_from_scalar(x_m, y_m, design.aperture_field, feed_xy_m=p.feed_xy_m)
        dx = float(x_m[0, 1] - x_m[0, 0])
        dy = float(y_m[1, 0] - y_m[0, 0])
        pat = compute_vector_pattern(x_m, y_m, jx, jy, f, green, dx, dy)
        total = np.abs(pat.e_theta) ** 2 + np.abs(pat.e_phi) ** 2
        peak_power = float(total.max())
        extent_m = float(x_m.max() - x_m.min()) + dx
        beams = port_beams[p.port_id]
        aims = [
            beam_local_peak_fine(x_m, y_m, jx, jy, f, green, dx, dy, b.theta_deg, b.phi_deg) for b in beams
        ]
        fine_peak = max((a[2] for a in aims), default=0.0)
        per_beam = []
        for idx, (b, (theta_a, phi_a, power)) in enumerate(zip(beams, aims, strict=True)):
            req = _unit_vector(b.theta_deg, b.phi_deg)
            got = _unit_vector(theta_a, phi_a)
            pointing_error = float(np.degrees(np.arccos(np.clip(float(req @ got), -1.0, 1.0))))
            pol = b.desired_polarization
            cross_level_db = xpd_db = None
            if pol in ORTHOGONAL_POLARIZATION:
                gi = int(np.argmin(np.abs(pat.theta_deg[:, 0] - theta_a)))
                gj = int(np.argmin(np.abs(((pat.phi_deg[0, :] - phi_a + 180.0) % 360.0) - 180.0)))
                cross = copol_component(pat.e_theta, pat.e_phi, ORTHOGONAL_POLARIZATION[pol])
                cross_level_db = float(20.0 * np.log10(max(abs(cross[gi, gj]), 1e-12) / np.sqrt(peak_power)))
                others = [(o.theta_deg, o.phi_deg) for k, o in enumerate(beams) if k != idx]
                xpd_db = beam_polarization_metrics(pat, pol, b.theta_deg, b.phi_deg, extent_m, wavelength_m, others).point_xpd_db
            level_db = float(10.0 * np.log10(max(power, 1e-30) / fine_peak)) if fine_peak > 0.0 else None
            per_beam.append(
                {
                    "requested_theta_deg": b.theta_deg,
                    "requested_phi_deg": b.phi_deg,
                    "theta_deg": float(theta_a),
                    "phi_deg": float(phi_a),
                    "pointing_error_deg": pointing_error,
                    "level_db": level_db,
                    "directivity_dbi": float(pat.directivity_estimate_dbi + 10.0 * np.log10(max(power, 1e-30) / peak_power)),
                    "polarization": pol,
                    "cross_pol_level_db": cross_level_db,
                    "xpd_db": xpd_db,
                }
            )
        eff, _dmax = aperture_efficiency_estimate(pat.directivity_estimate_dbi, p.region.area_m2, f)
        reports[p.port_id] = {
            "peak_directivity_dbi": float(pat.directivity_estimate_dbi),
            "aperture_efficiency_upper_bound": eff,
            "beams": per_beam,
        }
    return reports


def _unit_vector(theta_deg, phi_deg):
    t, p = np.radians(theta_deg), np.radians(phi_deg)
    return np.array([np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)])


def _state_table(multi, body):
    """The phase table of a coded multi-port design, reported ONCE (every port shares it): for 1 bit
    the two states' own phases (see `two_state`); for N bits one row per level -- its ideal phase,
    how many cells (all ports) use it, and the mean achieved phase / error of those cells."""
    n_bits = body["n_bits"]
    if n_bits == 1:
        return None  # `two_state.states` carries the pair
    step = 360.0 / 2**n_bits
    rows = []
    for level in range(2**n_bits):
        errors, achieved = [], []
        for p in multi.ports:
            lm = np.asarray(p.design.level_map)
            for c in p.design.cells:
                if int(lm[c.i, c.j]) == level and c.achieved_value is not None:
                    achieved.append(c.achieved_value)
                    errors.append(c.phase_error_deg if c.phase_error_deg is not None else 0.0)
        rows.append(
            {
                "level": level,
                "target_phase_deg": level * step,
                "n_cells": len(achieved),
                "mean_achieved_phase_deg": float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(achieved))))))
                if achieved
                else None,
                "mean_phase_error_deg": float(np.mean(errors)) if errors else None,
            }
        )
    return rows


def _run_synthesize_ports_job(job, body: dict):
    regions = _build_port_regions(body["layout"], body["period_mm"])
    if len(body["ports"]) != len(regions):
        raise config_invalid(
            message=f"layout '{body['layout']['preset']}' has {len(regions)} port(s) but "
            f"{len(body['ports'])} port beam list(s) were given"
        )
    port_beams = {}
    for region, port in zip(regions, body["ports"], strict=True):
        port_beams[region.port_id] = [
            Beam(
                theta_deg=b["theta_deg"],
                phi_deg=b["phi_deg"],
                weight=complex(b["weight_real"], b["weight_imag"]),
                initial_phase_deg=b["initial_phase_deg"],
                desired_polarization=b["desired_polarization"],
                oam_charge=b["oam_charge"],
            )
            for b in port["beams"]
        ]

    period_m = body["period_mm"] * 1e-3
    k_swg = tm0_beta(body["frequency_hz"], body["port_eps_r"], body["port_thickness_mm"] * 1e-3)

    if body["formulation"] == "A":
        if body["dispersion_id"] is None or body["x_avg_ohm"] is None or body["modulation_index"] is None:
            raise config_invalid(message="formulation A needs dispersion_id, x_avg_ohm and modulation_index")
        dispersion_table = read_dispersion_csv(
            upload_dir(body["dispersion_id"]) / "dispersion.csv", body["frequency_hz"]
        )

        def synthesize_one(ref, beams, feed_xy_m):
            return synthesize_formulation_a(
                ref,
                beams,
                body["frequency_hz"],
                body["x_avg_ohm"],
                body["modulation_index"],
                dispersion_table,
                period_m,
                rescale_to_beta_bar=body["rescale_to_beta_bar"],
                feed_xy_m=feed_xy_m,
                cp_phase_offset_deg=body["cp_phase_offset_deg"],
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

        def synthesize_one(ref, beams, feed_xy_m):
            return synthesize_formulation_b(
                ref,
                beams,
                body["frequency_hz"],
                candidates,
                period_m,
                feed_xy_m=feed_xy_m,
                cp_phase_offset_deg=body["cp_phase_offset_deg"],
            )
    elif body["formulation"] == "coded":
        n_bits = body["n_bits"]
        auto_n_bit_info = {}
        if body["beta_mode"] == "paper":
            raise config_invalid(message="beta_mode 'paper' was removed; use 'model'")
        if n_bits < 1:
            raise config_invalid(message=f"n_bits must be >= 1, got {n_bits!r}")
        if n_bits == 1:
            if body["port_tan_delta"] is None:
                raise config_invalid(
                    message="formulation 'coded' needs port_tan_delta (the surface-wave check needs the loss tangent)"
                )
            if body["beta_mode"] not in ("state", "model"):
                raise config_invalid(message=f"beta_mode must be 'state' or 'model', got {body['beta_mode']!r}")
            substrate = Substrate("user", body["port_eps_r"], body["port_tan_delta"], body["port_thickness_mm"] * 1e-3)
            # model: one beta for the whole mixed surface, so the states' own betas need not agree
            spread_limit = BETA_SPREAD_LIMIT_FRACTION if body["beta_mode"] == "state" else float("inf")
            auto_info = {}
            choice, refusals = _resolve_two_state_choice(
                body["two_state"],
                body["library_id"],
                body["frequency_hz"],
                substrate,
                period_m,
                beta_spread_limit=spread_limit,
                beta_mode=body["beta_mode"],
                job=job,
                auto_info_out=auto_info,
            )
            beta_model = mixed_surface_mode(choice, substrate, body["frequency_hz"], period_m).beta_rad_per_m
            if body["beta_mode"] == "state":
                betas = (choice.surface_wave_0.beta_rad_per_m, choice.surface_wave_1.beta_rad_per_m)
            else:
                betas = (beta_model, beta_model)

            def synthesize_one(ref, beams, feed_xy_m):
                design = synthesize_coded_one_bit(
                    ref,
                    beams,
                    body["frequency_hz"],
                    choice,
                    period_m,
                    feed_xy_m=feed_xy_m,
                    cp_phase_offset_deg=body["cp_phase_offset_deg"],
                    surface_wave_betas_rad_per_m=betas,
                )
                if body["feed_keepout"]:
                    design, _n_changed = apply_feed_keepout_one_bit(design, feed_xy_m, period_m)
                return design
        else:
            if body["two_state"]["source"] == "auto":
                # Phase 22: the 2**n_bits states are found automatically (no library needed)
                if body["port_tan_delta"] is None:
                    raise config_invalid(message="automatic cell states need port_tan_delta (the MoM solve needs it)")
                auto_substrate = Substrate(
                    "user", body["port_eps_r"], body["port_tan_delta"], body["port_thickness_mm"] * 1e-3
                )
                auto = _compute_auto_states(n_bits, body["two_state"], auto_substrate, body["frequency_hz"], period_m, job)
                auto_n_bit_info.update(auto_states_to_dict(auto))
                candidates = auto.library_candidates()
            else:
                # same validation messages as the single-aperture /synthesize_coded N-bit path
                if body["library_id"] is None:
                    raise config_invalid(message="n_bits > 1 needs library_id (a continuous phase library)")
                lib = UnitCellLibrary.from_json(library_registry.resolve_path(body["library_id"]))
                candidates = library_candidates_at_frequency(lib, body["library_id"], body["frequency_hz"])
                if not candidates:
                    raise library_frequency_unavailable(
                        f"library {body['library_id']!r} has no record at {body['frequency_hz'] / 1e9:g} GHz"
                    )

            def synthesize_one(ref, beams, feed_xy_m):
                return synthesize_coded_n_bit(
                    ref,
                    beams,
                    body["frequency_hz"],
                    n_bits,
                    candidates,
                    period_m,
                    feed_xy_m=feed_xy_m,
                    cp_phase_offset_deg=body["cp_phase_offset_deg"],
                )

    else:
        raise config_invalid(message=f"formulation must be 'A', 'B' or 'coded', got {body['formulation']!r}")

    # Phase 23: a polygon port keeps only the cells inside its outline (`synthesize_ports` passes
    # `region=` for polygons only, so a rectangle layout's calls are exactly as before).
    synthesize_one_unmasked = synthesize_one

    def synthesize_one(ref, beams, feed_xy_m, region=None):
        design = synthesize_one_unmasked(ref, beams, feed_xy_m)
        if region is None:
            return design
        earlier = [r for r in regions if r.port_id < region.port_id]
        return mask_design_to_region(design, region_cell_mask(region, design.x_m, design.y_m, earlier))

    refinement_reports = []
    per_port = synthesize_one
    region_ids = [r.port_id for r in regions]
    refine = body["pointing_correction"] or (
        body["balance_levels"] and any(len(pb["beams"]) > 1 for pb in body["ports"])
    )
    if body["pointing_correction"] and body["port_tan_delta"] is None:
        raise config_invalid(
            message="pointing_correction needs port_tan_delta (the pattern check needs the loss tangent)"
        )
    if refine and body["port_tan_delta"] is not None:
        port_substrate = Substrate("user", body["port_eps_r"], body["port_tan_delta"], body["port_thickness_mm"] * 1e-3)
        port_green = LayeredGreenFunction(port_substrate, body["frequency_hz"], backing="pec")

        def per_port(ref, beams, feed_xy_m, region=None):
            def measure(design, beams_now):
                return measure_beam_aims(design, beams_now, port_green, period_m, body["frequency_hz"], feed_xy_m)

            design, _used, report = refine_design(
                beams,
                lambda bs: synthesize_one(ref, bs, feed_xy_m, region),
                measure,
                pointing=body["pointing_correction"],
                balance=body["balance_levels"],
                pointing_iterations=body["pointing_iterations"],
            )
            refinement_reports.append({"port_id": region_ids[len(refinement_reports)], **report})
            return design

    multi = synthesize_ports(regions, port_beams, per_port, period_m, k_swg=k_swg, min_lambda_sw=body["min_lambda_sw"])

    isolation = []
    for flag in multi.isolation_flags:
        entry = {
            "port_a": flag.port_a,
            "port_b": flag.port_b,
            "gap_mm": flag.gap_m * 1e3,
            "gap_lambda_sw": flag.gap_lambda_sw,
            "flagged": flag.flagged,
            "guard_row_suggestion_mm": flag.guard_row_suggestion_m * 1e3 if flag.guard_row_suggestion_m else None,
        }
        isolation.append(entry)
    isolation_labels = tuple(
        isolation_edge_distance(f.port_a, f.port_b, f.gap_lambda_sw) for f in multi.isolation_flags if f.flagged
    )

    if body["formulation"] == "coded" and body["n_bits"] > 1:
        cell_geometry_m = _port_cell_to_m("B", body["period_mm"])
        board = combined_board_lw_grid_m(multi, cell_geometry_m)
        two_state_info = None
    elif body["formulation"] == "coded":

        def coded_cell_m(cell):
            state = choice.state_0 if int(cell.geometry_param) == 0 else choice.state_1
            return _two_state_geometry_to_cell_m(state.geometry)

        cell_geometry_m = coded_cell_m
        board = combined_board_lw_grid_m(multi, coded_cell_m)
        two_state_info = {
            "achieved_separation_deg": choice.achieved_separation_deg,
            "beta_bar_rad_per_m": choice.beta_bar_rad_per_m,
            "beta_spread_fraction": choice.beta_spread_fraction,
            "refused_states": refusals,
            "beta_mode": body["beta_mode"],
            "beta_used_rad_per_m": betas[0] if body["beta_mode"] != "state" else None,
            "beta_model_rad_per_m": beta_model,
            "state_betas_rad_per_m": [choice.surface_wave_0.beta_rad_per_m, choice.surface_wave_1.beta_rad_per_m],
            "auto": auto_info or None,
            # the exact phases the aperture field was built with (state.value, not the patch-plane shift);
            # Phase 23's sidelobe optimiser rebuilds each cell's two possible aperture values from them
            "achieved_phase_deg": [
                float(principal_phase_deg(choice.state_0.value)),
                float(principal_phase_deg(choice.state_1.value)),
            ],
            # the phase table, once for the whole board (every port shares the pair)
            "states": [
                {
                    "label": st.label,
                    "geometry": st.geometry,
                    "magnitude_db": float(20.0 * np.log10(abs(st.value))),
                    "phase_deg_patch_plane": float(np.degrees(np.angle(st.patch_plane_value(body["frequency_hz"])))),
                }
                for st in (choice.state_0, choice.state_1)
            ],
        }
    else:
        cell_geometry_m = _port_cell_to_m(body["formulation"], body["period_mm"])
        board = combined_board_lw_grid_m(multi, cell_geometry_m)
        two_state_info = None
    patch_layout_png_url = None
    try:
        cells_mm = []
        for i, row in enumerate(board.lw_grid_m):
            for j, cell_geometry in enumerate(row):
                cx, cy = board.cell_center_m(i, j)
                if isinstance(cell_geometry, dict):
                    cells_mm.append(
                        (
                            cx * 1e3,
                            cy * 1e3,
                            cell_geometry["a_m"] * 1e3,
                            cell_geometry["b_m"] * 1e3,
                            cell_geometry["R_m"] * 1e3,
                            cell_geometry["r_m"] * 1e3,
                        )
                    )
                    continue
                l_m, w_m = cell_geometry
                if l_m <= 0.0 or w_m <= 0.0:
                    continue
                cells_mm.append((cx * 1e3, cy * 1e3, l_m * 1e3, w_m * 1e3))
        out_dir = job_output_dir(job.id)
        n_y, n_x = len(board.lw_grid_m), len(board.lw_grid_m[0])
        path = patch_layout_png(
            cells_mm,
            body["period_mm"],
            n_x,
            n_y,
            out_dir / "patch_layout.png",
            title=f"Multi-port board — {len(multi.ports)} ports, {n_x} × {n_y} cells",
            feed_xy_mm=[(p.feed_xy_m[0] * 1e3, p.feed_xy_m[1] * 1e3) for p in multi.ports],
            outlines_mm=[[(x * 1e3, y * 1e3) for x, y in p.region.vertices_m] for p in multi.ports if p.region.is_polygon]
            or None,
        )
        patch_layout_png_url = f"/api/files/{job.id}/{path.name}"
    except Exception:
        patch_layout_png_url = None

    reports = _ports_reports(multi, port_beams, body, cell_geometry_m)
    ports_result = []
    for p in multi.ports:
        patches, n_patches = _port_patches_mm(p, cell_geometry_m, body["period_mm"])
        ports_result.append(
            {
                "report": reports.get(p.port_id),
                "patch_count": n_patches,
                "patches_mm": patches,
                "keepout_cells_changed": getattr(p.design, "keepout_cells_changed", None)
                if body["formulation"] == "coded" and body["n_bits"] == 1 and body["feed_keepout"]
                else None,
                "port_id": p.port_id,
                "feed_x_mm": p.feed_xy_m[0] * 1e3,
                "feed_y_mm": p.feed_xy_m[1] * 1e3,
                "x_min_mm": p.region.x_min_m * 1e3,
                "x_max_mm": p.region.x_max_m * 1e3,
                "y_min_mm": p.region.y_min_m * 1e3,
                "y_max_mm": p.region.y_max_m * 1e3,
                **(
                    {"outline_mm": [[x * 1e3, y * 1e3] for x, y in p.region.vertices_m]}
                    if p.region.is_polygon
                    else {}
                ),
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
        "n_bits": body["n_bits"] if body["formulation"] == "coded" else None,
        "two_state": two_state_info,
        "state_table": _state_table(multi, body) if body["formulation"] == "coded" else None,
        "auto_states": (auto_n_bit_info or None) if body["formulation"] == "coded" and body["n_bits"] > 1 else None,
        "refinement": refinement_reports,
        "period_mm": body["period_mm"],
        "patch_layout_png_url": patch_layout_png_url,
        "board_width_mm": board.width_m * 1e3,
        "board_height_mm": board.height_m * 1e3,
        "feeds_mm": [(p.feed_xy_m[0] * 1e3, p.feed_xy_m[1] * 1e3) for p in multi.ports],
        "export_feeds_mm": export_frame_feeds_and_rects_mm(multi, board)[0],
        "export_offset_mm": [v * 1e3 for v in export_frame_offset_m(multi, board)],
        "lw_grid_m": board.lw_grid_m,
        **({"layout_preset": "hex4", "trim_board": bool(body["layout"]["trim_board"])} if any(r.is_polygon for r in regions) else {}),
    }


@router.post("/synthesize_ports")
def synthesize_ports_endpoint(body: HmsaSynthesizePortsRequest, request: Request):
    if body.formulation in ("B", "coded") and body.library_id is not None:
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
    r_values_mm = _annular_slot_r_values_mm(body_dict["a_mm"], body_dict["b_mm"], body_dict["r_mm"])
    base_config = _annular_library_base_config(body_dict, placeholder_R_mm=r_values_mm[0])

    def _progress(event):
        job.progress = event["done"] / max(event["total"], 1)
        job.progress_message = f"{event['done']}/{event['total']} slot radii solved"

    cells = build_annular_slot_library_cells(
        base_config,
        body_dict["a_mm"],
        body_dict["b_mm"],
        body_dict["r_mm"],
        r_values_mm,
        body_dict["frequency_hz"],
        progress=_progress,
    )
    out_dir = job_output_dir(job.id)
    path = write_annular_slot_library_json(cells, out_dir / "annular_slot_library.json")
    library_registry.register(
        library_id,
        path,
        meta={
            "kind": "annular_slot",
            "period_mm": body_dict["period_mm"],
            "backing": "pec",
            "a_mm": body_dict["a_mm"],
            "b_mm": body_dict["b_mm"],
            "r_mm": body_dict["r_mm"],
            "n_geometries": len(cells),
            "auto_built": True,
            "frequency_hz": body_dict["frequency_hz"],
        },
    )
    return {
        "annular_slot_library_id": library_id,
        "n_geometries": len(cells),
        "r_values_mm": r_values_mm,
        "already_built": False,
    }


@router.post("/ensure_annular_slot_library")
def ensure_annular_slot_library(body: EnsureAnnularSlotLibraryRequest, request: Request):
    library_id = _auto_annular_slot_library_id(body.model_dump())
    try:
        library_registry.resolve_path(library_id)
        job = request.app.state.jobs.submit(
            "hmsa_ensure_annular_slot_library",
            lambda job: {"annular_slot_library_id": library_id, "already_built": True},
        )
        return job.to_dict()
    except KeyError:
        pass
    key = config_hash({"kind": "hmsa_ensure_annular_slot_library", **body.model_dump()})
    job = request.app.state.jobs.submit(
        "hmsa_ensure_annular_slot_library",
        _run_ensure_annular_slot_library_job,
        body.model_dump(),
        library_id,
        cache_key=key,
    )
    return job.to_dict()


def _compute_auto_states(n_bits, two_state, substrate, frequency_hz, period_m, job=None):
    """Phase 22: the automatic cell states, cached on disk under the app data dir (a re-run is instant)."""

    def progress(message):
        if job is not None:
            job.progress_message = message

    return auto_states(
        n_bits,
        substrate,
        frequency_hz,
        period_m,
        fixed_L_mm=two_state.get("auto_fixed_L_mm") if n_bits == 1 else None,
        min_gap_mm=two_state.get("auto_min_gap_mm", 0.2),
        cache_dir=app_data_dir() / "auto_states_cache",
        progress=progress,
    )


WEAK_PAIR_SEPARATION_DEG = 150.0  # a fixed-side 1-bit pair below this separation is flagged


def auto_states_to_dict(auto):
    """JSON form of an `AutoStates`: the table the UI shows (state, geometry, |Gamma|, patch-plane phase,
    separation, sensitivity, labels) plus timings."""
    states = auto.states
    rows = []
    for s in states:
        sep = None
        if len(states) == 2 and s.level == 1:
            sep = abs(wrapped_difference_deg(states[0].phase_deg, s.phase_deg))
        rows.append(
            {
                "level": s.level,
                "label": s.label,
                "geometry": s.geometry,
                "gamma_real": s.gamma.real,
                "gamma_imag": s.gamma.imag,
                "magnitude_db": s.magnitude_db,
                "phase_deg_patch_plane": s.phase_deg,
                "target_phase_deg": float(((s.target_phase_deg + 180.0) % 360.0) - 180.0),
                "phase_error_deg": s.phase_error_deg,
                "separation_deg": sep,
                "sensitivity_deg_per_0p01mm": s.sensitivity_deg_per_0p01mm,
                "mesh_uncertainty_deg": s.mesh_uncertainty_deg,
                "labels": list(s.labels),
                "message": s.message,
            }
        )
    warnings = []
    if auto.fixed_L_mm is not None and auto.n_bits == 1 and rows[1]["separation_deg"] < WEAK_PAIR_SEPARATION_DEG:
        warnings.append(
            f"Weak pair: {auto.fixed_L_mm:.2f} mm gives {rows[1]['separation_deg']:.0f}° separation "
            f"(below {WEAK_PAIR_SEPARATION_DEG:.0f}°); leave the fixed side blank to search."
        )
    return {
        "n_bits": auto.n_bits,
        "fixed_L_mm": auto.fixed_L_mm,
        "warnings": warnings,
        "states": rows,
        "solve_seconds": auto.solve_seconds,
        "n_solves": auto.n_solves,
        "n_cached": auto.n_cached,
        "min_gap_mm": auto.min_gap_mm,
        "reference_plane": "patch plane (z = 0)",
    }


def _resolve_two_state_choice(
    two_state: dict,
    library_id,
    frequency_hz,
    substrate,
    period_m,
    beta_spread_limit=BETA_SPREAD_LIMIT_FRACTION,
    beta_mode="state",
    job=None,
    auto_info_out=None,
):
    """Every candidate state is screened first (`screen_candidates`: TM0 root, clear of
    patch-plane resonance, inside the sheet model's range); the pair is then
    the best one whose BOTH states pass with agreeing betas (`best_valid_pair`), never the raw
    best-separated pair. Returns `(choice, refusals)`, `refusals` naming every refused state."""
    source = two_state["source"]
    if source == "auto":
        auto = _compute_auto_states(1, two_state, substrate, frequency_hz, period_m, job)
        if auto_info_out is not None:
            auto_info_out.update(auto_states_to_dict(auto))
        candidates = auto.cells()
    elif source == "uploaded":
        for key in ("uploaded_value_0_real", "uploaded_value_0_imag", "uploaded_value_1_real", "uploaded_value_1_imag"):
            if two_state[key] is None:
                raise config_invalid(message=f"two_state source 'uploaded' needs {key}")
        candidates = [
            TwoStateCell(
                two_state["uploaded_label_0"],
                "uploaded_state_0",
                complex(two_state["uploaded_value_0_real"], two_state["uploaded_value_0_imag"]),
                two_state["uploaded_geometry_0"] or {},
            ),
            TwoStateCell(
                two_state["uploaded_label_1"],
                "uploaded_state_1",
                complex(two_state["uploaded_value_1_real"], two_state["uploaded_value_1_imag"]),
                two_state["uploaded_geometry_1"] or {},
            ),
        ]
    else:
        if library_id is None:
            raise config_invalid(message=f"two_state source {source!r} needs library_id (a square-patch library)")
        lib = UnitCellLibrary.from_json(library_registry.resolve_path(library_id))
        candidates = candidates_from_square_library(library_candidates_at_frequency(lib, library_id, frequency_hz))
        if source == "square_only":
            pass
        elif source == "square_and_annular":
            if two_state["annular_slot_library_id"] is None:
                raise config_invalid(message="two_state source 'square_and_annular' needs annular_slot_library_id")
            annular_path = library_registry.resolve_path(two_state["annular_slot_library_id"])
            annular_cells = json.loads(annular_path.read_text())["cells"]
            candidates = candidates + candidates_from_annular_slot_library(annular_cells, frequency_hz / 1e9)
        else:
            raise config_invalid(message=f"unknown two_state source {source!r}")
    # an automatic pair with one beta for the whole mixed surface (model/paper) may sit on the patch resonance
    refuse_resonance = not (source == "auto" and beta_mode != "state")
    modes, reasons = screen_candidates(candidates, frequency_hz, substrate, period_m, refuse_resonance)
    choice = best_valid_pair(
        candidates,
        modes,
        frequency_hz,
        substrate,
        period_m=period_m,
        reasons=reasons,
        beta_spread_limit=beta_spread_limit,
    )
    refusals = [
        {"label": c.label, "record_id": c.record_id, "reason": r}
        for c, r in zip(candidates, reasons, strict=True)
        if r is not None
    ]
    return choice, refusals


def _state_s11_at_frequency(state, body, frequency_hz):
    """The library-referenced (5 mm above the patch) S11 of one two-state candidate RE-SOLVED at
    `frequency_hz` with the same solver settings its library used -- for the +-5% robustness check
    on the chosen pair only. `None` for a state with no re-solvable geometry (an uploaded value) or
    a solve that did not converge."""
    geom = state.geometry
    if geom.get("kind") == "annular_slot":
        base = _annular_library_base_config({**body, "frequency_hz": frequency_hz}, placeholder_R_mm=geom["R_mm"])
        docs = build_annular_slot_library_cells(
            base, geom["a_mm"], geom["b_mm"], geom["r_mm"], [geom["R_mm"]], frequency_hz
        )
    elif "L_mm" in geom:
        docs = build_square_library_cells(
            _square_library_base_config({**body, "frequency_hz": frequency_hz}), [geom["L_mm"]], frequency_hz
        )
    else:
        return None
    doc = docs[0]
    if doc["solver_status"] != "OK" or not doc["S11"]:
        return None
    entry = doc["S11"][0]
    return complex(entry["real"], entry["imag"])


def _frequency_checks(choice, body):
    """The state's |S11| and phase difference at f*0.95, f and f*1.05 (patch-plane phases). Returns
    `(rows, note)`: `rows` is None when a state cannot be re-solved, with the reason in `note`."""
    f = body["frequency_hz"]
    rows = []
    for scale in (0.95, 1.05):
        f_k = f * scale
        values = []
        for state in (choice.state_0, choice.state_1):
            gamma_top = _state_s11_at_frequency(state, body, f_k)
            if gamma_top is None:
                return None, f"{state.label} could not be re-solved at {f_k / 1e9:.3f} GHz (uploaded or no solution)."
            values.append(patch_plane_gamma(gamma_top, f_k, LIBRARY_REFERENCE_OFFSET_M))
        rows.append(_frequency_row(f_k, values[0], values[1]))
    rows.insert(1, _frequency_row(f, choice.state_0.patch_plane_value(f), choice.state_1.patch_plane_value(f)))
    return rows, None


def _frequency_row(f_hz, gamma_0, gamma_1):
    sep = abs(wrapped_difference_deg(principal_phase_deg(gamma_0), principal_phase_deg(gamma_1)))
    return {
        "frequency_GHz": f_hz / 1e9,
        "state_0_magnitude_db": float(20.0 * np.log10(abs(gamma_0))),
        "state_1_magnitude_db": float(20.0 * np.log10(abs(gamma_1))),
        "state_0_phase_patch_deg": float(np.degrees(np.angle(gamma_0))),
        "state_1_phase_patch_deg": float(np.degrees(np.angle(gamma_1))),
        "separation_deg": float(sep),
    }


class AutoStatesRequest(BaseModel):
    """Phase 22: find the coded-cell states for a substrate/cell/frequency/bits without designing anything
    (the Coded panel's "Compute cell states" button; results are cached, so the later design run reuses them)."""

    eps_r: float
    tan_delta: float
    thickness_mm: float
    period_mm: float
    frequency_hz: float
    n_bits: int = 1
    fixed_L_mm: "float | None" = None
    min_gap_mm: float = 0.2


def _run_auto_states_job(job, body: dict):
    substrate = Substrate("user", body["eps_r"], body["tan_delta"], body["thickness_mm"] * 1e-3)
    two_state = {"auto_fixed_L_mm": body["fixed_L_mm"], "auto_min_gap_mm": body["min_gap_mm"]}
    auto = _compute_auto_states(body["n_bits"], two_state, substrate, body["frequency_hz"], body["period_mm"] * 1e-3, job)
    return {**auto_states_to_dict(auto), "labels": (EXPERIMENTAL,)}


@router.post("/auto_states")
def auto_states_endpoint(body: AutoStatesRequest, request: Request):
    if body.n_bits < 1 or body.n_bits > 3:
        raise HTTPException(422, f"n_bits must be 1, 2 or 3, got {body.n_bits}")
    if body.fixed_L_mm is not None and body.n_bits != 1:
        raise HTTPException(422, "fixed_L_mm is for 1-bit designs")
    if not (0.0 <= body.min_gap_mm < body.period_mm / 2.0):
        raise HTTPException(422, f"min_gap_mm must be in [0, period/2), got {body.min_gap_mm}")
    job = request.app.state.jobs.submit("hmsa_auto_states", _run_auto_states_job, body.model_dump())
    return job.to_dict()


class HmsaSynthesizeCodedRequest(BaseModel):
    """PHASE15 Part E: the Coded design method. Reuses the single-aperture reference-field
    resolution (`upload_id` or `feeds`) and beam spec exactly as `/synthesize` does -- only the
    control-map/geometry-mapping step differs (`hologram.coded_synthesis`)."""

    pointing_correction: bool = False
    pointing_iterations: int = 2
    balance_levels: bool = True

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
    # Substrate of the library the two-state cells were solved on. Required for the 1-bit method
    # (the surface-wave check and the +-5% re-solves need it); optional for N-bit (then the
    # quantization-loss measurement is not computed).
    eps_r: "float | None" = None
    tan_delta: "float | None" = None
    thickness_mm: "float | None" = None
    accuracy: Literal["fast", "accurate"] = "fast"


def _two_state_geometry_to_cell_m(geometry):
    """A `two_state_cells.TwoStateCell.geometry` dict (millimetres) -> the `(length_m, width_m)` /
    annular-slot-dict cell entry `export.array_export`'s mixed-cell writers expect (metres), or
    `None` if this state carries no exportable geometry at all (the "uploaded" two-state source,
    `two_state_cells.two_state_from_values` with no `geometry_0`/`geometry_1` given -- an
    HFSS-measured reflection coefficient with no physical cell shape recorded in this app; Export
    then correctly refuses with "no geometry" rather than this function fabricating one)."""
    if geometry.get("kind") == "empty":
        return (0.0, 0.0)
    if geometry.get("kind") == "annular_slot":
        return {
            "kind": "annular_slot",
            "a_m": geometry["a_mm"] * 1e-3,
            "b_m": geometry["b_mm"] * 1e-3,
            "R_m": geometry["R_mm"] * 1e-3,
            "r_m": geometry["r_mm"] * 1e-3,
        }
    if "L_mm" in geometry and "W_mm" in geometry:
        return (geometry["L_mm"] * 1e-3, geometry["W_mm"] * 1e-3)
    return None


def _ref_complex(reference_field):
    return (
        reference_field.Ex if reference_field.amplitude_available else np.exp(1j * np.radians(reference_field.phase()))
    )


def _pattern_power(x_m, y_m, aperture_field, feed_xy_m, body, green, dx, dy):
    jx, jy = aperture_current_from_scalar(x_m, y_m, aperture_field, feed_xy_m=feed_xy_m)
    pat = compute_vector_pattern(x_m, y_m, jx, jy, body["frequency_hz"], green, dx, dy)
    power = np.abs(pat.e_theta) ** 2 + np.abs(pat.e_phi) ** 2
    theta_rad = np.radians(pat.theta_deg)
    dtheta = np.radians(pat.theta_deg[1, 0] - pat.theta_deg[0, 0])
    dphi = np.radians(360.0 / pat.phi_deg.shape[1])
    radiated = float(np.sum(power * np.sin(theta_rad) * dtheta * dphi))
    return pat, float(power.max()), radiated


def _quantization_report(design, cell_field, feed_xy_m, body, green, dx, dy, n_bits):
    """Quantization loss, defined explicitly (never a single number of unknown origin):
    the CONTINUOUS design `E_ref*exp(j*psi_cont)` is compared with the IDEAL quantized one
    `E_ref*exp(j*level)` (levels 0/180deg for 1 bit) -- same aperture amplitude, so equal input
    power. Reported three ways: (1) the peak far-field power ratio at each design's own main beam
    (the textbook quantity, the main-beam gain loss); (2) the aperture coupling
    `-10log10(|<E_q,E_c>|^2/(|E_q|^2|E_c|^2))` (the same quantity computed without a pattern grid);
    (3) the state-error loss, ideal levels vs the ACHIEVED measured states, computed as an aperture
    coupling so a constant phase offset common to both states (a reference-plane shift) does not
    count as error -- it moves the whole aperture's phase, not the pattern. The directivity
    difference is also returned, for reference only: directivity normalises by the power the
    design RADIATES, and a 1-bit design radiates less, so the directivity loss is much smaller than
    the peak-gain loss (see docs/history/PHASE15_POLARIZATION_VORTEX.md Part E)."""
    ref_complex = _ref_complex(cell_field)
    e_cont = ref_complex * np.exp(1j * np.radians(design.continuous_phase_deg))
    e_ideal = ref_complex * np.exp(1j * np.radians(design.quantized_phase_deg))
    e_ach = design.aperture_field
    pat_c, peak_c, rad_c = _pattern_power(design.x_m, design.y_m, e_cont, feed_xy_m, body, green, dx, dy)
    pat_i, peak_i, rad_i = _pattern_power(design.x_m, design.y_m, e_ideal, feed_xy_m, body, green, dx, dy)
    pat_a, peak_a, rad_a = _pattern_power(design.x_m, design.y_m, e_ach, feed_xy_m, body, green, dx, dy)
    return {
        "definition": "continuous design E_ref*exp(j*psi_cont) vs ideal quantized levels E_ref*exp(j*level); equal aperture power",
        "n_bits": n_bits,
        "textbook_db": textbook_quantization_loss_db(n_bits),
        "peak_loss_db": float(10.0 * np.log10(peak_c / peak_i)),
        "coupling_loss_db": aperture_coupling_loss_db(e_ideal, e_cont),
        "state_error_loss_db": aperture_coupling_loss_db(e_ideal, e_ach),
        "radiated_power_ratio_db": float(10.0 * np.log10(rad_i / rad_c)),
        "directivity_difference_db": float(pat_c.directivity_estimate_dbi - pat_i.directivity_estimate_dbi),
        "main_beam_continuous": (float(pat_c.main_beam_theta_deg), float(pat_c.main_beam_phi_deg)),
        "main_beam_ideal": (float(pat_i.main_beam_theta_deg), float(pat_i.main_beam_phi_deg)),
        "pattern_achieved": pat_a,
    }


def _pattern_checks(pat_a, design, body):
    wavelength_m = 299792458.0 / body["frequency_hz"]
    extent_m = float(design.x_m.max() - design.x_m.min()) + body["period_mm"] * 1e-3
    exclusion = first_null_half_width_deg(extent_m, wavelength_m, pat_a.main_beam_theta_deg)
    secondary = strongest_secondary_beam(pat_a, exclusion)
    mirror_db, mirror_theta, mirror_phi = mirror_beam_level_db(pat_a)
    return {
        "main_beam_theta_deg": float(pat_a.main_beam_theta_deg),
        "main_beam_phi_deg": float(pat_a.main_beam_phi_deg),
        "main_lobe_exclusion_deg": exclusion,
        "strongest_secondary": None
        if secondary is None
        else {
            "theta_deg": secondary[0],
            "phi_deg": secondary[1],
            "level_db": secondary[2],
        },
        "mirror_beam": {"theta_deg": mirror_theta, "phi_deg": mirror_phi, "level_db": mirror_db},
    }


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
            theta_deg=b["theta_deg"],
            phi_deg=b["phi_deg"],
            weight=complex(b["weight_real"], b["weight_imag"]),
            initial_phase_deg=b["initial_phase_deg"],
            desired_polarization=b["desired_polarization"],
            oam_charge=b["oam_charge"],
        )
        for b in body["object_beams"]
    ]
    period_m = body["period_mm"] * 1e-3
    n_bits = body["n_bits"]
    has_substrate = body["eps_r"] is not None and body["tan_delta"] is not None and body["thickness_mm"] is not None
    substrate = (
        Substrate("user", body["eps_r"], body["tan_delta"], body["thickness_mm"] * 1e-3) if has_substrate else None
    )
    warnings = []
    extra_checks = {}

    state_labels = ()
    auto_n_bit_info = {}
    if n_bits == 1:
        if substrate is None:
            raise config_invalid(
                message="the 1-bit coded method needs the substrate (eps_r, tan_delta, thickness_mm) for its surface-wave check"
            )
        if body["feeds"] is not None:
            raise config_invalid(
                message="the 1-bit coded method with multiple feeds is not supported yet (one feed defines the surface-wave reference)"
            )
        auto_info = {}
        choice, refusals = _resolve_two_state_choice(
            body["two_state"],
            body["library_id"],
            body["frequency_hz"],
            substrate,
            period_m,
            job=job,
            auto_info_out=auto_info,
        )
        beta_bar = choice.beta_bar_rad_per_m
        design, refinement = _refined_design(
            body,
            beams,
            lambda bs: synthesize_coded_one_bit(
                reference_field,
                bs,
                body["frequency_hz"],
                choice,
                period_m,
                feed_xy_m=feed_xy_m,
                cp_phase_offset_deg=body["cp_phase_offset_deg"],
                surface_wave_betas_rad_per_m=(
                    choice.surface_wave_0.beta_rad_per_m,
                    choice.surface_wave_1.beta_rad_per_m,
                ),
            ),
            feed_xy_m,
        )
        if not design.path_converged:
            warnings.append(
                f"The path-integrated beta did not settle within {design.path_iterations} iterations; "
                "the last state map is reported."
            )
        sep_warning = separation_warning(choice)
        if sep_warning:
            warnings.append(sep_warning)
        state_labels = (
            two_state_phase_separation(choice.achieved_separation_deg, n_bits=1),
            two_state_surface_wave_check(beta_bar, choice.beta_spread_fraction, BETA_SPREAD_LIMIT_FRACTION),
        )
        states_report = []
        for state, sw in ((choice.state_0, choice.surface_wave_0), (choice.state_1, choice.surface_wave_1)):
            gamma_patch = state.patch_plane_value(body["frequency_hz"])
            states_report.append(
                {
                    "label": state.label,
                    "record_id": state.record_id,
                    "geometry": state.geometry,
                    "magnitude_db": float(20.0 * np.log10(abs(state.value))),
                    "phase_deg_library_reference": float(np.degrees(np.angle(state.value))),
                    "phase_deg_patch_plane": float(np.degrees(np.angle(gamma_patch))),
                    "resonance_proximity_deg": resonance_proximity_patch_plane_deg(gamma_patch),
                    "beta_rad_per_m": sw.beta_rad_per_m,
                    "x_s_ohm": sw.x_s_ohm,
                    "stop_band_warning": bool(sw.stop_band_warning),
                }
            )
        mag_gap_db = abs(states_report[0]["magnitude_db"] - states_report[1]["magnitude_db"])
        if mag_gap_db > MAGNITUDE_MISMATCH_WARNING_DB:
            warnings.append(
                f"The two states' magnitudes differ by {mag_gap_db:.2f} dB (limit {MAGNITUDE_MISMATCH_WARNING_DB:g} dB)."
            )
        frequency_rows, frequency_note = _frequency_checks(choice, body)
        if frequency_rows is not None:
            for row in frequency_rows:
                gap = abs(row["state_0_magnitude_db"] - row["state_1_magnitude_db"])
                if gap > MAGNITUDE_MISMATCH_WARNING_DB:
                    warnings.append(f"At {row['frequency_GHz']:.3f} GHz the two magnitudes differ by {gap:.2f} dB.")
        delta_beta = abs(choice.surface_wave_0.beta_rad_per_m - choice.surface_wave_1.beta_rad_per_m)
        rho_max = float(np.max(np.hypot(design.x_m - feed_xy_m[0], design.y_m - feed_xy_m[1])) + period_m / 2)
        edge_error_deg = float(np.degrees(0.5 * delta_beta * rho_max))
        if edge_error_deg > EDGE_PHASE_WARNING_DEG:
            warnings.append(
                f"The beta mismatch between the states leaves up to {edge_error_deg:.0f}deg of phase error at the "
                f"aperture edge (limit {EDGE_PHASE_WARNING_DEG:g}deg)."
            )
        two_state_info = {
            "reference_plane": "patch plane (z=0); library S11 de-embedded from its 5 mm reference",
            "states": states_report,
            "achieved_separation_deg": choice.achieved_separation_deg,
            "beta_bar_rad_per_m": beta_bar,
            "beta_spread_fraction": choice.beta_spread_fraction,
            "beta_spread_limit_fraction": BETA_SPREAD_LIMIT_FRACTION,
            "beta_edge_phase_error_deg": edge_error_deg,
            "refused_states": refusals,
            "auto": auto_info or None,
            "path_iterations": design.path_iterations,
            "path_converged": design.path_converged,
            "cell_size_mm": body["period_mm"],
            "lambda_sw_mm": 2.0 * np.pi / beta_bar * 1e3,
            "lambda_sw_over_cell": lambda_sw_over_cell(beta_bar, period_m),
            "lambda0_mm": 299792458.0 / body["frequency_hz"] * 1e3,
            "frequency_checks": frequency_rows,
            "frequency_checks_note": frequency_note,
            "warnings": warnings,
        }
        ny, nx = design.x_m.shape
        lw_grid_m = [[None] * nx for _ in range(ny)]
        for c in design.cells:
            state = choice.state_0 if int(c.geometry_param) == 0 else choice.state_1
            lw_grid_m[c.i][c.j] = _two_state_geometry_to_cell_m(state.geometry)
    else:
        if body["two_state"]["source"] == "auto":
            # Phase 22: the 2**n_bits states are found automatically (no library needed)
            if substrate is None:
                raise config_invalid(message="automatic cell states need the substrate (eps_r, tan_delta, thickness_mm)")
            auto = _compute_auto_states(n_bits, body["two_state"], substrate, body["frequency_hz"], period_m, job)
            auto_n_bit_info.update(auto_states_to_dict(auto))
            candidates = auto.library_candidates()
        else:
            if body["library_id"] is None:
                raise config_invalid(message="n_bits > 1 needs library_id (a continuous phase library)")
            lib = UnitCellLibrary.from_json(library_registry.resolve_path(body["library_id"]))
            candidates = library_candidates_at_frequency(lib, body["library_id"], body["frequency_hz"])
            if not candidates:
                raise library_frequency_unavailable(
                    f"library {body['library_id']!r} has no record at {body['frequency_hz'] / 1e9:g} GHz"
                )
        design, refinement = _refined_design(
            body,
            beams,
            lambda bs: synthesize_coded_n_bit(
                reference_field,
                bs,
                body["frequency_hz"],
                n_bits,
                candidates,
                period_m,
                feed_xy_m=feed_xy_m,
                cp_phase_offset_deg=body["cp_phase_offset_deg"],
            ),
            feed_xy_m,
        )
        two_state_info = None
        lw_grid_m = None  # Export for N>1 reuses the existing Formulation-B-shaped `assignments`

    assignments = None
    if n_bits > 1:
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

    cell_field = resample_to_cell_grid(reference_field, period_m)
    continuous_field = continuous_reference_aperture_field(cell_field, design.continuous_phase_deg)

    quantization = None
    quantization_label = None
    if substrate is not None:
        green = LayeredGreenFunction(substrate, body["frequency_hz"], backing="pec")
        dx = float(design.x_m[0, 1] - design.x_m[0, 0])
        dy = float(design.y_m[1, 0] - design.y_m[0, 0])
        report = _quantization_report(design, cell_field, feed_xy_m, body, green, dx, dy, n_bits)
        pat_a = report.pop("pattern_achieved")
        quantization = report
        quantization_label = quantization_loss_db(report["peak_loss_db"], n_bits, report["textbook_db"])
        extra_checks = _pattern_checks(pat_a, design, body)

    layout_url = _layout_preview_url(
        job,
        design.x_m,
        design.y_m,
        design.quantized_phase_deg,
        "target phase (deg)",
        f"Coded design target phase, {n_bits}-bit",
    )
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
            cells_mm,
            body["period_mm"],
            nx,
            ny,
            out_dir / "patch_layout.png",
            title=f"Coded design — {nx} × {ny} cells, {n_bits}-bit",
            feed_xy_mm=(feed_xy_m[0] * 1e3, feed_xy_m[1] * 1e3),
        )
        patch_layout_png_url = f"/api/files/{job.id}/{path.name}"
    except Exception:
        patch_layout_png_url = None

    labels = extra_labels + state_labels + ((quantization_label,) if quantization_label else ())
    return {
        "n_bits": n_bits,
        "refinement": refinement,
        "x_m": design.x_m,
        "y_m": design.y_m,
        "period_mm": body["period_mm"],
        "continuous_phase_deg": design.continuous_phase_deg,
        "quantized_phase_deg": design.quantized_phase_deg,
        "level_map": design.level_map,
        "cells": design.cells,
        "assignments": assignments,
        "two_state": two_state_info,
        "auto_states": auto_n_bit_info or None,
        "quantization": quantization,
        "pattern_checks": extra_checks or None,
        "labels": labels,
        "aperture_field_real": design.aperture_field.real,
        "aperture_field_imag": design.aperture_field.imag,
        "continuous_aperture_field_real": continuous_field.real if continuous_field is not None else None,
        "continuous_aperture_field_imag": continuous_field.imag if continuous_field is not None else None,
        "layout_png_url": layout_url,
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
