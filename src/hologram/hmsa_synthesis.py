"""HMSA hologram synthesis: wires the imported reference field, an object beam (or beams), and
Phase 8B's two hologram formulations into a per-cell control-quantity map, then maps that onto
geometry via a dispersion library (`src/dispersion/`) or, for the phase-difference formulation, a
normal-incidence-phase library (the same machinery the reflectarray designer uses).

- **Formulation A (Fong et al. 2010, impedance-modulated)**: needs `X_avg`/`modulation_index`
  (fabrication calibration constants Phase 8B already found this project cannot derive from its
  own data -- `phase8b.control_variable`) supplied directly by the caller/app user. Labelled
  `EXPERIMENTAL` for that reason; maps X(x,y) to geometry via a `dispersion.hfss_table.
  DispersionTable`'s `geometry_for_x_s`.
- **Formulation B (phase-difference)**: needs only `arg(E_ref)`/`arg(E_obj)` (works with a
  phase-only imported field). Maps required_phase(x,y) to geometry via the SAME normal-incidence
  S11-phase library machinery the reflectarray designer uses
  (`phase8d_rev.response_matching.nearest_response_match`).

Both paths carry `LOCAL_PERIODICITY_APPROXIMATION` (a per-cell lookup-table control quantity is
being used to represent a coupled, continuously-varying surface) and a cells-per-fringe check
(`fields.reference_field_grid.cells_per_fringe`/`aliasing_warning`) against the imported field's
own measured local wavenumber.
"""

from dataclasses import dataclass, replace

import numpy as np

from ..fields.reference_field_grid import aliasing_warning, cells_per_fringe, resample_to_cell_grid
from ..labels import (
    EXPERIMENTAL,
    LOCAL_PERIODICITY_APPROXIMATION,
    REFERENCE_AMPLITUDE_ASSUMED_UNIFORM,
    object_wave_normalised,
    phase_coverage_shortfall,
    reference_wave_not_rescaled,
    reference_wave_rescaled_to_beta_bar,
)
from ..phase8b.holographic_modulation import (
    interference_modulation,
    phase_difference_hologram_from_phases_deg,
    reactance_from_modulation,
)
from ..phase8d_rev.response_matching import nearest_response_match
from ..reflectarray.designer import PHASE_ERROR_UNREACHABLE_THRESHOLD_DEG
from .polarization_compensation import combine_beams_with_polarization
from .sector_multiplexing import sector_index_map


@dataclass(frozen=True)
class GeometryCell:
    i: int
    j: int
    x_m: float
    y_m: float
    control_value: float  # X(x,y) [ohm] for Formulation A, or required_phase [deg] for B
    geometry_param: "float | None"  # None if unreachable; L_mm for Formulation B, the dispersion
    # table's geometry parameter for A
    achieved_value: "float | None"
    phase_error_deg: "float | None"  # Formulation B only (phase-quantization error); None for A
    reachable: bool
    geometry_param_w: "float | None" = None  # Formulation B only: the matched cell's W_mm (the
    # library search is over the full L,W grid, not just square cells -- PHASE13B B2 needs both
    # to build an exportable per-cell L/W assignment)


@dataclass(frozen=True)
class HmsaDesign:
    formulation: str  # "A_impedance_modulation" | "B_phase_difference"
    x_m: np.ndarray
    y_m: np.ndarray
    control_map: np.ndarray
    cells: "tuple[GeometryCell, ...]"
    labels: tuple
    unreachable_fraction: float
    phase_error_stats: "dict | None"  # Formulation B only: {mean,rms,max}_deg, n_cells, n_unreachable
    aperture_field: "np.ndarray | None"  # Formulation B only: the RADIATED field for Pattern
    # prediction -- E_ref(x,y) * exp(j*achieved_modulation_phase), NOT the modulation phasor
    # alone (see `synthesize_formulation_b`'s docstring: using the modulation phasor by itself
    # silently drops the reference wave's own phase, which must cancel it back out)


def _combined_object_field(
    x_m, y_m, object_beams, frequency_hz, feed_xy_m=(0.0, 0.0), cp_phase_offset_deg=0.0, multiplexing="full_aperture"
):
    """PHASE15 Part A: the object-wave combination both formulations feed into their interference
    kernel, now polarization-compensated (`hologram.polarization_compensation.
    combine_beams_with_polarization`, replacing the old bare `phase8b.object_wave.combine_beams_scalar`
    call) and optionally sector-multiplexed (paper eq. 17: `multiplexing="sectors"` restricts each
    cell to its own nearest-azimuth beam, `sector_multiplexing.sector_index_map`; default
    `"full_aperture"` is eq. 16, every beam superposed everywhere -- unchanged from pre-Phase-15).
    `multiplexing="full_aperture"` with every beam left at its `desired_polarization="natural"`,
    `oam_charge=0` defaults makes this call identical to the old `combine_beams_scalar` sum (see
    `combine_beams_with_polarization`'s own docstring for why) -- the bit-identical-reproduction
    invariant existing designs rely on."""
    if multiplexing not in ("full_aperture", "sectors"):
        raise ValueError(f"multiplexing must be 'full_aperture' or 'sectors', got {multiplexing!r}")
    sector_owner_index = sector_index_map(x_m, y_m, feed_xy_m, object_beams) if multiplexing == "sectors" else None
    return combine_beams_with_polarization(
        x_m, y_m, object_beams, frequency_hz, feed_xy_m, cp_phase_offset_deg, sector_owner_index
    )


def object_wave_scale(e_obj):
    """Divisor that stops a summed object wave from exceeding unit amplitude. N beams in phase reach
    |E_obj| = N, which pushes Psi past its [-1, 1] range; dividing by the peak magnitude keeps the
    relative beam weights. A sum that already peaks at or below 1 is returned unchanged (divisor 1.0),
    so single beams at their own weight keep their pre-fix modulation exactly."""
    return max(1.0, float(np.max(np.abs(e_obj))))


def formulation_a_control_map_and_scale(
    reference_field,
    object_beams,
    frequency_hz,
    x_avg_ohm,
    modulation_index,
    feed_xy_m=(0.0, 0.0),
    cp_phase_offset_deg=0.0,
    multiplexing="full_aperture",
):
    """`(control_map, object_wave_divisor)` for `formulation_a_control_map`; the divisor is 1.0 unless
    the multi-beam sum saturated (see `object_wave_scale`)."""
    if not reference_field.amplitude_available:
        raise ValueError("Formulation A needs a complex reference field (amplitude available)")
    e_obj = _combined_object_field(
        reference_field.x_m,
        reference_field.y_m,
        object_beams,
        frequency_hz,
        feed_xy_m,
        cp_phase_offset_deg,
        multiplexing,
    )
    scale = object_wave_scale(e_obj)
    psi = interference_modulation(reference_field.Ex, e_obj / scale)
    return reactance_from_modulation(psi, x_avg_ohm, modulation_index), scale


def formulation_a_control_map(
    reference_field,
    object_beams,
    frequency_hz,
    x_avg_ohm,
    modulation_index,
    feed_xy_m=(0.0, 0.0),
    cp_phase_offset_deg=0.0,
    multiplexing="full_aperture",
):
    """Fong et al. 2010: X(x,y) = X_avg*(1 + M*Re{E_ref*conj(E_obj)}). Requires a COMPLEX
    reference field (amplitude available) -- see `_ensure_complex_reference` (PHASE13B B4) for the
    phase-only fallback callers should apply before this. `feed_xy_m`/`cp_phase_offset_deg`/
    `multiplexing`: PHASE15 Part A, see `_combined_object_field`. The summed object wave is kept at
    unit peak amplitude (see `object_wave_scale`)."""
    control, _scale = formulation_a_control_map_and_scale(
        reference_field,
        object_beams,
        frequency_hz,
        x_avg_ohm,
        modulation_index,
        feed_xy_m,
        cp_phase_offset_deg,
        multiplexing,
    )
    return control


def _ensure_complex_reference(reference_field):
    """PHASE13B B4: the live review found Formulation A refused to run at all without a complex
    (magnitude+phase) reference field, even though the project's only real reference-field export
    is phase-only. The brief permits running phase-only here too, with UNIFORM amplitude assumed
    (never silently -- the caller must apply `REFERENCE_AMPLITUDE_ASSUMED_UNIFORM`). Returns
    `(field_with_Ex, amplitude_was_assumed)`."""
    if reference_field.amplitude_available:
        return reference_field, False
    ex = np.exp(1j * np.radians(reference_field.phase()))
    return replace(reference_field, Ex=ex, Ey=np.zeros_like(ex), amplitude_available=True), True


def _rescale_reference_to_beta_bar(cell_field, raw_reference_field, dispersion_table, x_avg_ohm):
    """PHASE13B B3: the imported field is the bare-slab TM0 wave; once patches are present the
    guided wavenumber becomes beta_bar, set by the hologram's own mean surface reactance X_avg
    (looked up from the SAME dispersion table Formulation A already needs). Using the bare-slab
    phase as-is accumulates a phase error ~ (beta_bar - beta_slab)*r across the aperture. Rather
    than the brief's idealized analytic ψ_ref ~ A(r)*exp(-j*beta_bar*r) (which would discard every
    measured deviation from a perfectly radial wave), this applies the additive correction
    `phase_new(r) = phase_measured(r) + (beta_slab - beta_bar)*r` -- replaces only the mean radial
    slope, keeps the field's own measured detail, and reduces to the brief's formula exactly when
    the measured field already IS purely radial.

    Returns `(rescaled_cell_field, beta_bar, beta_slab, worst_case_error_deg)`. Raises `ValueError`
    if the dispersion table has no `beta_rad_per_m` column (a `x_s_ohm`-only table, e.g. from a
    reactance-only HFSS export, cannot supply beta_bar at all -- the caller falls back to
    reporting the accumulated-error estimate without applying a correction)."""
    from ..fields.reference_field_grid import local_wavenumber_field

    beta_field = local_wavenumber_field(_as_fld_like(raw_reference_field))
    beta_slab = float(np.nanmedian(beta_field))
    geometry_param, _achieved_xs, _within = dispersion_table.geometry_for_x_s(x_avg_ohm)
    beta_bar = dispersion_table.beta_at(geometry_param)  # raises ValueError if no beta column

    r_m = np.hypot(cell_field.x_m, cell_field.y_m)
    correction_rad = (beta_slab - beta_bar) * r_m
    ex_new = cell_field.Ex * np.exp(1j * correction_rad)
    worst_case_error_deg = float(np.degrees(abs(beta_bar - beta_slab) * np.max(r_m)))
    rescaled = replace(cell_field, Ex=ex_new, Ey=np.zeros_like(ex_new), phase_deg=np.degrees(np.angle(ex_new)))
    return rescaled, beta_bar, beta_slab, worst_case_error_deg


def formulation_b_control_map(
    reference_field,
    object_beams,
    frequency_hz,
    feed_xy_m=(0.0, 0.0),
    cp_phase_offset_deg=0.0,
    multiplexing="full_aperture",
):
    """Phase-difference hologram: required_phase(x,y) = wrap(arg(E_obj) - arg(E_ref)). Works with
    a phase-only reference field (the brief's own "the one case where a phase-only reference is
    sufficient"). `feed_xy_m`/`cp_phase_offset_deg`/`multiplexing`: PHASE15 Part A, see
    `_combined_object_field`."""
    e_obj = _combined_object_field(
        reference_field.x_m,
        reference_field.y_m,
        object_beams,
        frequency_hz,
        feed_xy_m,
        cp_phase_offset_deg,
        multiplexing,
    )
    obj_phase_deg = np.degrees(np.angle(e_obj))
    return phase_difference_hologram_from_phases_deg(reference_field.phase(), obj_phase_deg)


def map_control_to_geometry_dispersion(x_m, y_m, control_map, dispersion_table):
    """Formulation A path: X(x,y) -> geometry via `DispersionTable.geometry_for_x_s`."""
    ny, nx = np.shape(x_m)
    cells = []
    for i in range(ny):
        for j in range(nx):
            target = float(control_map[i, j])
            gp, achieved, within = dispersion_table.geometry_for_x_s(target)
            cells.append(GeometryCell(i, j, float(x_m[i, j]), float(y_m[i, j]), target, gp, achieved, None, within))
    return tuple(cells)


def map_control_to_geometry_library(
    x_m, y_m, required_phase_deg_map, library_candidates, threshold_deg=PHASE_ERROR_UNREACHABLE_THRESHOLD_DEG
):
    """Formulation B path: required_phase(x,y) -> geometry via the normal-incidence-phase
    library, reusing `nearest_response_match` unmodified (phase-only match: unit-magnitude target,
    same convention `reflectarray.designer.assign_cells_from_library` uses). `reachable` requires
    BOTH a candidate to exist AND its phase error to be within `threshold_deg` (PHASE13B A4 --
    previously any candidate at all, however large its error, counted as reachable; phase
    libraries are periodic so there is no true out-of-range case the way a monotonic dispersion
    table has, but a large quantization error is still a real design problem)."""
    candidates = list(library_candidates)
    ny, nx = np.shape(x_m)
    cells = []
    for i in range(ny):
        for j in range(nx):
            target_deg = float(required_phase_deg_map[i, j])
            target_complex = np.exp(1j * np.radians(target_deg))
            matches = nearest_response_match(target_complex, candidates)
            best = matches[0] if matches else None
            achieved = float(np.degrees(np.angle(best.complex_value))) if best else None
            phase_error = best.phase_error_deg if best else None
            reachable = best is not None and phase_error <= threshold_deg
            cells.append(
                GeometryCell(
                    i,
                    j,
                    float(x_m[i, j]),
                    float(y_m[i, j]),
                    target_deg,
                    best.L_mm if best else None,
                    achieved,
                    phase_error,
                    reachable,
                    geometry_param_w=best.W_mm if best else None,
                )
            )
    return tuple(cells)


def synthesize_formulation_a(
    reference_field,
    object_beams,
    frequency_hz,
    x_avg_ohm,
    modulation_index,
    dispersion_table,
    period_m,
    rescale_to_beta_bar=True,
    feed_xy_m=(0.0, 0.0),
    cp_phase_offset_deg=0.0,
    multiplexing="full_aperture",
):
    # PHASE13B A2: synthesize on the unit-cell-centre grid of the CHOSEN period, not one "cell"
    # per raw field sample (a 151x151, 2mm-pitch import would otherwise give 22801 "cells" instead
    # of the real e.g. 37x37 the period implies).
    cell_field = resample_to_cell_grid(reference_field, period_m)
    cell_field, amplitude_assumed = _ensure_complex_reference(cell_field)  # PHASE13B B4
    extra_labels = [EXPERIMENTAL]
    if amplitude_assumed:
        extra_labels.append(REFERENCE_AMPLITUDE_ASSUMED_UNIFORM)

    try:  # PHASE13B B3: rescale to beta_bar (on by default), or at least report the risk of not doing so
        rescaled_field, beta_bar, _beta_slab, worst_case_error_deg = _rescale_reference_to_beta_bar(
            cell_field, reference_field, dispersion_table, x_avg_ohm
        )
        if rescale_to_beta_bar:
            cell_field = rescaled_field
            extra_labels.append(reference_wave_rescaled_to_beta_bar(beta_bar))
        else:
            extra_labels.append(reference_wave_not_rescaled(worst_case_error_deg))
    except ValueError:
        pass  # the dispersion table has no beta_rad_per_m column -- beta_bar is not computable at all

    control_map, object_scale = formulation_a_control_map_and_scale(
        cell_field,
        object_beams,
        frequency_hz,
        x_avg_ohm,
        modulation_index,
        feed_xy_m,
        cp_phase_offset_deg,
        multiplexing,
    )
    if object_scale > 1.0:
        extra_labels.append(object_wave_normalised(object_scale))
    cells = map_control_to_geometry_dispersion(cell_field.x_m, cell_field.y_m, control_map, dispersion_table)

    # PHASE14 Part C: an equivalent aperture field for pattern prediction -- NOT the same model as
    # Formulation B's (a reflectarray-style phase-front reconstruction), since this surface
    # radiates by a different mechanism (leaky-wave scattering off the reactance MODULATION, not a
    # re-reflected phase front). Fong et al.'s own derivation is a first-order (small modulation
    # depth) perturbation result: only the AC part of the surface admittance -- proportional to
    # `X_avg*M*Psi`, i.e. to `X(x,y)/X_avg - 1` -- couples the bound surface wave to radiation; the
    # DC part (X_avg itself) stays purely reactive and non-radiating by design. To that same first
    # order, the locally radiated field is proportional to the exciting surface wave times that
    # admittance perturbation: `E_radiated(x,y) ~ E_ref(x,y) * (X_achieved(x,y)/X_avg - 1)`, using
    # the ACHIEVED (dispersion-table-matched, quantization-error-included) reactance per cell --
    # the real design's expected radiated field, not the idealized continuous target (mirroring
    # exactly how Formulation B uses the achieved, not the ideal, phase). This is a documented
    # approximation (a rigorous leaky-wave/full-wave radiation calculation is out of scope here),
    # not a claim of the same fidelity as a full leaky-wave analysis.
    achieved_x_ohm = np.array([c.achieved_value if c.achieved_value is not None else x_avg_ohm for c in cells]).reshape(
        cell_field.x_m.shape
    )
    psi_achieved = achieved_x_ohm / x_avg_ohm - 1.0
    ref_complex = cell_field.Ex if cell_field.amplitude_available else np.exp(1j * np.radians(cell_field.phase()))
    aperture_field = ref_complex * psi_achieved

    return _finish_design(
        "A_impedance_modulation",
        reference_field,
        cell_field.x_m,
        cell_field.y_m,
        control_map,
        cells,
        period_m,
        extra_labels=tuple(extra_labels),
        aperture_field=aperture_field,
    )


def synthesize_formulation_b(
    reference_field,
    object_beams,
    frequency_hz,
    library_candidates,
    period_m,
    feed_xy_m=(0.0, 0.0),
    cp_phase_offset_deg=0.0,
    multiplexing="full_aperture",
):
    # PHASE13B B4: B uses a normal-incidence S11-phase library -- a reflectarray control variable,
    # not the surface-wave one (dispersion/impedance) a monopole-fed hologram actually needs.
    # Recommended path is A; B is offered as a usable-today approximation, labelled accordingly.
    cell_field = resample_to_cell_grid(reference_field, period_m)
    control_map = formulation_b_control_map(
        cell_field, object_beams, frequency_hz, feed_xy_m, cp_phase_offset_deg, multiplexing
    )
    cells = map_control_to_geometry_library(cell_field.x_m, cell_field.y_m, control_map, library_candidates)

    # The RADIATED field a modulated hologram actually produces is E_ref(x,y) * Modulation(x,y),
    # not the modulation phasor alone: if Modulation = exp(j*(obj_phase - ref_phase)) exactly
    # cancelled the reference wave's own phase, E_ref*Modulation = |E_ref|*exp(j*obj_phase) --
    # recovering the pure desired beam. Using the modulation phasor by itself (as an earlier
    # version of this code did) silently drops that cancellation and predicts a systematically
    # wrong beam direction (empirically: the example file, a 20deg/0deg object beam, predicted a
    # peak at 51deg/130deg instead of ~20deg/0deg -- fixed by this multiplication, verified against
    # an idealized zero-quantization-error case landing at 21deg/0deg). Uses the achieved
    # (library-matched, quantization-error-included) phase per cell -- the REAL design's expected
    # radiated field, not the idealized continuous target.
    achieved_deg = np.array([c.achieved_value if c.achieved_value is not None else 0.0 for c in cells]).reshape(
        cell_field.x_m.shape
    )
    ref_complex = cell_field.Ex if cell_field.amplitude_available else np.exp(1j * np.radians(cell_field.phase()))
    aperture_field = ref_complex * np.exp(1j * np.radians(achieved_deg))

    return _finish_design(
        "B_phase_difference",
        reference_field,
        cell_field.x_m,
        cell_field.y_m,
        control_map,
        cells,
        period_m,
        extra_labels=(EXPERIMENTAL,),
        aperture_field=aperture_field,
    )


def _finish_design(
    formulation, raw_reference_field, x_m, y_m, control_map, cells, period_m, extra_labels, aperture_field
):
    """`raw_reference_field` is the RAW (un-resampled) imported field, used only for the
    cells-per-fringe sanity check (a property of the imported wave itself, best measured at its
    native resolution); `x_m`/`y_m`/`control_map`/`cells` are already on the resampled cell grid."""
    labels = [LOCAL_PERIODICITY_APPROXIMATION, *extra_labels]
    if raw_reference_field.dx_m is not None:
        from ..fields.reference_field_grid import local_wavenumber_field

        try:
            beta_field = local_wavenumber_field(_as_fld_like(raw_reference_field))
            beta_typical = float(np.nanmedian(beta_field))
            cpf = cells_per_fringe(beta_typical, period_m)
            warning = aliasing_warning(cpf)
            if warning:
                labels.append(f"CELLS_PER_FRINGE({cpf:.1f})")
        except Exception:
            pass  # best-effort sanity check only; never blocks synthesis
    unreachable = sum(1 for c in cells if not c.reachable)
    unreachable_fraction = unreachable / len(cells) if cells else 0.0
    phase_errors = [c.phase_error_deg for c in cells if c.phase_error_deg is not None]
    stats = None
    if phase_errors:
        errors = np.array(phase_errors)
        stats = {
            "mean_deg": float(np.mean(errors)),
            "rms_deg": float(np.sqrt(np.mean(errors**2))),
            "max_deg": float(np.max(errors)),
            "n_cells": len(cells),
            "n_unreachable": unreachable,
        }
        if unreachable > 0:
            worst = max(c.phase_error_deg for c in cells if not c.reachable)
            labels.append(phase_coverage_shortfall(worst))
    return HmsaDesign(
        formulation=formulation,
        x_m=x_m,
        y_m=y_m,
        control_map=control_map,
        cells=cells,
        labels=tuple(labels),
        unreachable_fraction=unreachable_fraction,
        phase_error_stats=stats,
        aperture_field=aperture_field,
    )


def _as_fld_like(reference_field):
    """A tiny duck-typed shim so `local_wavenumber_field` (which reads `.value`/.x_m/.y_m) can be
    reused on a `ReferenceFieldGrid` (which stores phase in `.phase_deg`, degrees) without
    duplicating the unwrap/gradient logic."""
    from dataclasses import dataclass as _dc

    @_dc(frozen=True)
    class _Shim:
        x_m: object
        y_m: object
        value: object

    phase_rad = np.radians(reference_field.phase())
    return _Shim(reference_field.x_m, reference_field.y_m, phase_rad)
