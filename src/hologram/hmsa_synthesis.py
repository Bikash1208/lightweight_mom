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

from dataclasses import dataclass

import numpy as np

from ..errors import impedance_unreachable
from ..fields.reference_field_grid import aliasing_warning, cells_per_fringe
from ..labels import EXPERIMENTAL, LOCAL_PERIODICITY_APPROXIMATION
from ..phase8b.holographic_modulation import (
    interference_modulation,
    phase_difference_hologram_from_phases_deg,
    reactance_from_modulation,
)
from ..phase8b.object_wave import combine_beams_scalar
from ..phase8d_rev.response_matching import nearest_response_match


@dataclass(frozen=True)
class GeometryCell:
    i: int
    j: int
    x_m: float
    y_m: float
    control_value: float            # X(x,y) [ohm] for Formulation A, or required_phase [deg] for B
    geometry_param: "float | None"  # None if unreachable
    achieved_value: "float | None"
    reachable: bool


@dataclass(frozen=True)
class HmsaDesign:
    formulation: str   # "A_impedance_modulation" | "B_phase_difference"
    x_m: np.ndarray
    y_m: np.ndarray
    control_map: np.ndarray
    cells: "tuple[GeometryCell, ...]"
    labels: tuple
    unreachable_fraction: float


def formulation_a_control_map(reference_field, object_beams, frequency_hz, x_avg_ohm, modulation_index):
    """Fong et al. 2010: X(x,y) = X_avg*(1 + M*Re{E_ref*conj(E_obj)}). Requires a COMPLEX
    reference field (amplitude available), since Formulation A's kernel needs it."""
    if not reference_field.amplitude_available:
        raise ValueError("Formulation A needs a complex reference field (amplitude available)")
    e_obj = combine_beams_scalar(reference_field.x_m, reference_field.y_m, object_beams, frequency_hz)
    psi = interference_modulation(reference_field.Ex, e_obj)
    return reactance_from_modulation(psi, x_avg_ohm, modulation_index)


def formulation_b_control_map(reference_field, object_beams, frequency_hz):
    """Phase-difference hologram: required_phase(x,y) = wrap(arg(E_obj) - arg(E_ref)). Works with
    a phase-only reference field (the brief's own "the one case where a phase-only reference is
    sufficient")."""
    e_obj = combine_beams_scalar(reference_field.x_m, reference_field.y_m, object_beams, frequency_hz)
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
            cells.append(GeometryCell(i, j, float(x_m[i, j]), float(y_m[i, j]), target, gp, achieved, within))
    return tuple(cells)


def map_control_to_geometry_library(x_m, y_m, required_phase_deg_map, library_candidates):
    """Formulation B path: required_phase(x,y) -> geometry via the normal-incidence-phase
    library, reusing `nearest_response_match` unmodified (phase-only match: unit-magnitude target,
    same convention `reflectarray.designer.assign_cells_from_library` uses)."""
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
            # "reachable" here means the library has ANY candidate; phase libraries are periodic
            # (any phase is nominally reachable to within quantization) so there is no true
            # out-of-range case the way a monotonic dispersion table has -- quantization error is
            # reported via achieved vs target instead.
            cells.append(GeometryCell(i, j, float(x_m[i, j]), float(y_m[i, j]), target_deg,
                                       best.L_mm if best else None, achieved, best is not None))
    return tuple(cells)


def synthesize_formulation_a(reference_field, object_beams, frequency_hz, x_avg_ohm, modulation_index,
                              dispersion_table, period_m):
    control_map = formulation_a_control_map(reference_field, object_beams, frequency_hz, x_avg_ohm, modulation_index)
    cells = map_control_to_geometry_dispersion(reference_field.x_m, reference_field.y_m, control_map, dispersion_table)
    return _finish_design("A_impedance_modulation", reference_field, control_map, cells, period_m, extra_labels=(EXPERIMENTAL,))


def synthesize_formulation_b(reference_field, object_beams, frequency_hz, library_candidates, period_m):
    control_map = formulation_b_control_map(reference_field, object_beams, frequency_hz)
    cells = map_control_to_geometry_library(reference_field.x_m, reference_field.y_m, control_map, library_candidates)
    return _finish_design("B_phase_difference", reference_field, control_map, cells, period_m, extra_labels=())


def _finish_design(formulation, reference_field, control_map, cells, period_m, extra_labels):
    labels = [LOCAL_PERIODICITY_APPROXIMATION, *extra_labels]
    if reference_field.dx_m is not None:
        # Estimate local wavenumber from the ORIGINAL reference field's own phase, not from
        # control_map (which for Formulation A is a reactance, not a phase).
        from ..fields.reference_field_grid import local_wavenumber_field

        try:
            beta_field = local_wavenumber_field(_as_fld_like(reference_field))
            beta_typical = float(np.nanmedian(beta_field))
            cpf = cells_per_fringe(beta_typical, period_m)
            warning = aliasing_warning(cpf)
            if warning:
                labels.append(f"CELLS_PER_FRINGE({cpf:.1f})")
        except Exception:
            pass   # best-effort sanity check only; never blocks synthesis
    unreachable = sum(1 for c in cells if not c.reachable)
    unreachable_fraction = unreachable / len(cells) if cells else 0.0
    return HmsaDesign(
        formulation=formulation, x_m=reference_field.x_m, y_m=reference_field.y_m,
        control_map=control_map, cells=cells, labels=tuple(labels), unreachable_fraction=unreachable_fraction,
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
