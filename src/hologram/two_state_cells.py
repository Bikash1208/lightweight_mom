"""Phase 15 Part E: two-state (1-bit) unit-cell selection for the Coded design method (Mishra,
Barad, Sarkar, IEEE TVT 2026). A 1-bit coded hologram needs two reflective cell states with
reflection phases ~180deg apart and comparable magnitude -- Mishra's own pair is a full square
patch vs. a patch with an annular slot (`src.geometry.AnnularSlotPatch`), but the achieved
separation is a MEASUREMENT of whatever two real candidates are chosen, never assumed to be
exactly 180deg (`separation_warning` reports the shortfall honestly, per this project's own
scientific-honesty convention -- see CLAUDE.md).

Reference planes: the auto-built libraries store S11 referenced 5 mm ABOVE the patch
(`routers/hmsa.py`'s library base configs use `top_air_thickness=5e-3`), while the Phase 14
surface-wave model needs the patch-plane value. `TwoStateCell.reference_offset_m` records that
offset, and `patch_plane_value` de-embeds it with `sparameters.Sparameters.at_planes`'s own
phase rule (`S11_patch = S11_top * exp(+2j*k0*offset)`).

Surface-wave check (the pair-validity rule): a state is usable only if the Phase 14 sheet model
finds a TM0 root at the design frequency (`sheet_model.find_beta_from_gamma`). A pair is valid
only if BOTH states carry the wave and their beta values agree within `BETA_SPREAD_LIMIT_FRACTION`
of their mean. That limit is this project's own choice (not from the paper) and is reported with
every design rather than hidden.
"""

from dataclasses import dataclass

import numpy as np

from ..dispersion.sheet_model import find_beta_from_sheet_admittance, sheet_admittance_from_gamma
from ..errors import LwmomError
from ..frequency import free_space_wavenumber
from ..phase import principal_phase_deg, wrapped_difference_deg

IDEAL_SEPARATION_DEG = 180.0
SEPARATION_WARNING_THRESHOLD_DEG = 30.0  # how far from 180deg before this is flagged, not silently used
RESONANCE_PHASE_PATCH_PLANE_DEG = 0.0  # open-circuit patch resonance, at the patch plane (Huang & Encinar)
RESONANCE_REFUSAL_DEG = 30.0  # a state this close to the patch-plane resonance is refused, not used
REFUSED_NO_SURFACE_WAVE = "NO_SURFACE_WAVE"
REFUSED_NEAR_RESONANCE = "NEAR_PATCH_PLANE_RESONANCE"
REFUSED_OUTSIDE_MODEL_RANGE = "OUTSIDE_SHEET_MODEL_RANGE"
MAGNITUDE_MISMATCH_WARNING_DB = 1.0
BETA_SPREAD_LIMIT_FRACTION = 0.05  # |beta_a - beta_b| / mean(beta) above this -> the pair is invalid
EDGE_PHASE_WARNING_DEG = 30.0  # beta mismatch (dbeta/2)*rho_max at the aperture edge above this -> flagged
LIBRARY_REFERENCE_OFFSET_M = 5.0e-3  # the auto-built libraries' S11 reference plane, above the patch


@dataclass(frozen=True)
class TwoStateCell:
    """One candidate reflective state: a human label, its own record id (traceable back to a
    library entry or an uploaded row), the achieved complex reflection coefficient (referenced
    `reference_offset_m` above the patch plane, 0 meaning the patch plane itself), and whatever
    geometry parameters identify it for export (`{"kind": "rectangle", "L_mm", "W_mm"}` or
    `{"kind": "annular_slot", "a_mm", "b_mm", "R_mm", "r_mm"}`)."""

    label: str
    record_id: str
    value: complex
    geometry: dict
    reference_offset_m: float = 0.0

    def patch_plane_value(self, frequency_hz):
        return patch_plane_gamma(self.value, frequency_hz, self.reference_offset_m)


@dataclass(frozen=True)
class TwoStateChoice:
    state_0: TwoStateCell
    state_1: TwoStateCell
    achieved_separation_deg: float
    surface_wave_0: "object | None" = None  # sheet_model.SheetModelMode, when checked
    surface_wave_1: "object | None" = None
    beta_bar_rad_per_m: "float | None" = None  # mean of the two states' beta, when checked
    beta_spread_fraction: "float | None" = None


def patch_plane_gamma(value, frequency_hz, reference_offset_m):
    """A reflection coefficient referenced `reference_offset_m` above the patch, moved to z=0."""
    if reference_offset_m == 0.0:
        return complex(value)
    k0 = free_space_wavenumber(frequency_hz)
    return complex(value) * np.exp(2j * k0 * reference_offset_m)


def candidates_from_square_library(candidates):
    """`candidates`: `(record_id, L_mm, W_mm, library_id, complex_value)` tuples (the format
    `phase8d_rev.response_matching.nearest_response_match`/`hmsa_synthesis.map_control_to_geometry_library`
    already use) -- wraps each as a `TwoStateCell` labelled by its own square side length."""
    return [
        TwoStateCell(
            label=f"square patch L={l_mm:g}mm",
            record_id=record_id,
            value=complex(value),
            geometry={"kind": "rectangle", "L_mm": l_mm, "W_mm": w_mm},
            reference_offset_m=LIBRARY_REFERENCE_OFFSET_M,
        )
        for record_id, l_mm, w_mm, _library_id, value in candidates
    ]


def candidates_from_annular_slot_library(cell_docs, frequency_ghz):
    """`cell_docs`: `library_sweep.build_annular_slot_library_cells`'s own return list -- wraps
    each successfully-solved cell, at the frequency sample nearest `frequency_ghz`, as a
    `TwoStateCell`."""
    out = []
    for doc in cell_docs:
        if doc["solver_status"] != "OK" or not doc["frequency_GHz"]:
            continue
        idx = min(range(len(doc["frequency_GHz"])), key=lambda i: abs(doc["frequency_GHz"][i] - frequency_ghz))
        entry = doc["S11"][idx]
        value = complex(entry["real"], entry["imag"])
        out.append(
            TwoStateCell(
                label=f"annular slot R={doc['R_mm']:g}mm",
                record_id=doc["cell_id"],
                value=value,
                geometry={
                    "kind": "annular_slot",
                    "a_mm": doc["a_mm"],
                    "b_mm": doc["b_mm"],
                    "R_mm": doc["R_mm"],
                    "r_mm": doc["r_mm"],
                },
                reference_offset_m=LIBRARY_REFERENCE_OFFSET_M,
            )
        )
    return out


def resonance_proximity_patch_plane_deg(gamma):
    """Distance in degrees from the patch-plane resonance phase (0deg: the open-circuit patch
    resonance of a reflectarray element, Huang & Encinar). The library's own stored reference
    plane is 5 mm above the patch, where this same resonance appears at a different phase, so
    the check must run on `patch_plane_value`, never on the library-referenced value."""
    return abs(wrapped_difference_deg(principal_phase_deg(gamma), RESONANCE_PHASE_PATCH_PLANE_DEG))


def screen_gamma(gamma, substrate, frequency_hz, period_m=None, refuse_near_resonance=True):
    """The pair-validity screen for one patch-plane reflection coefficient. Returns
    `(sheet_model_mode, None)` for a usable state, or `(None, reason)` with `reason` one of the
    `REFUSED_*` names. Checks run in order: a TM0 root must exist (Phase 14 sheet model), the
    state must not sit within `RESONANCE_REFUSAL_DEG` of the patch-plane resonance, and the sheet
    model's own range must hold (below the stop band, and the period under half a free-space
    wavelength). No reactance-sign test: for a lossless cell Z = j*eta0*cot(theta/2), so two states
    180deg apart always have opposite-sign reactance and such a test would refuse every pair."""
    from ..dispersion.sheet_model import find_beta_from_gamma

    if not np.isfinite(gamma) or abs(1.0 + gamma) < 1e-9:
        return None, REFUSED_NO_SURFACE_WAVE
    try:
        mode = find_beta_from_gamma(gamma, substrate, frequency_hz, period_m=period_m)
    except LwmomError:
        return None, REFUSED_NO_SURFACE_WAVE
    if refuse_near_resonance and resonance_proximity_patch_plane_deg(gamma) <= RESONANCE_REFUSAL_DEG:
        return None, REFUSED_NEAR_RESONANCE
    if period_m is not None:
        wavelength_m = 299792458.0 / frequency_hz
        if mode.stop_band_warning or period_m >= wavelength_m / 2.0:
            return None, REFUSED_OUTSIDE_MODEL_RANGE
    return mode, None


def screen_state(cell, frequency_hz, substrate, period_m=None, refuse_near_resonance=True):
    """`screen_gamma` applied to one `TwoStateCell`'s patch-plane value."""
    return screen_gamma(
        cell.patch_plane_value(frequency_hz), substrate, frequency_hz, period_m, refuse_near_resonance
    )


def screen_candidates(candidates, frequency_hz, substrate, period_m=None, refuse_near_resonance=True):
    """`screen_state` over every candidate: `(modes, reasons)`, both aligned with `candidates`
    (`modes[i]` None and `reasons[i]` a `REFUSED_*` name for a refused state).

    `refuse_near_resonance=False` (Phase 22, the automatic pair with beta_mode "model"/"paper"): the
    near-resonance refusal exists because a state's own beta is only trustworthy away from the patch
    resonance; with one beta for the whole mixed surface (model/paper) that reason does not apply, and
    the 180 degree partner of the empty cell is often exactly a near-resonant square."""
    modes, reasons = [], []
    for cell in candidates:
        mode, reason = screen_state(cell, frequency_hz, substrate, period_m, refuse_near_resonance)
        modes.append(mode)
        reasons.append(reason)
    return modes, reasons


def surface_wave_of_state(cell, frequency_hz, substrate, period_m=None):
    """The screened TM0 solution for one state (`screen_state`), or `None` if the state is refused."""
    return screen_state(cell, frequency_hz, substrate, period_m)[0]


def paper_surface_wave_beta(frequency_hz, eps_r):
    """Mishra eq. (1): beta_sw = k0 * sqrt((eps_r + 1) / 2), one value for the whole mixed surface."""
    return free_space_wavenumber(frequency_hz) * np.sqrt((eps_r + 1.0) / 2.0)


def mixed_surface_mode(choice, substrate, frequency_hz, period_m):
    """Model mode: beta of the mixed 1-bit surface. The two states cover equal areas, so their
    sheet admittances combine as the mean of the two `Y_sheet` values (parallel loading), then
    `sheet_model.find_beta_from_sheet_admittance` solves the transverse resonance for that mean."""
    y_0 = sheet_admittance_from_gamma(choice.state_0.patch_plane_value(frequency_hz), substrate, frequency_hz)
    y_1 = sheet_admittance_from_gamma(choice.state_1.patch_plane_value(frequency_hz), substrate, frequency_hz)
    return find_beta_from_sheet_admittance(0.5 * (y_0 + y_1), substrate, frequency_hz, period_m=period_m)


def best_separated_pair(candidates):
    """Brute-force O(n^2) search (library sizes here are at most a few dozen to ~100 candidates --
    this project never guesses at a smarter search without measuring that it is needed) for the
    pair of `candidates` whose reflection PHASES are closest to 180deg apart -- real, measured
    values only. `candidates` may mix cell families (e.g. square-patch AND annular-slot
    `TwoStateCell`s together), which is exactly how Mishra's own two-state pair (a plain patch vs.
    one WITH an annular slot) is found in practice. No surface-wave check -- see `best_valid_pair`."""
    candidates = list(candidates)
    if len(candidates) < 2:
        raise ValueError(f"best_separated_pair needs at least 2 candidates, got {len(candidates)}")
    best = None
    for i in range(len(candidates)):
        phase_i = principal_phase_deg(candidates[i].value)
        for j in range(i + 1, len(candidates)):
            phase_j = principal_phase_deg(candidates[j].value)
            sep = abs(wrapped_difference_deg(phase_i, phase_j))
            if best is None or sep > best[0]:
                best = (sep, candidates[i], candidates[j])
    sep, a, b = best
    return TwoStateChoice(state_0=a, state_1=b, achieved_separation_deg=float(sep))


def pair_beta_spread_fraction(beta_a, beta_b):
    mean = 0.5 * (beta_a + beta_b)
    return abs(beta_a - beta_b) / mean


def best_valid_pair(
    candidates,
    surface_waves,
    frequency_hz,
    substrate,
    period_m=None,
    beta_spread_limit=BETA_SPREAD_LIMIT_FRACTION,
    reasons=None,
):
    """The Part E pair rule: the largest measured phase separation among pairs where BOTH states
    carry the TM0 surface wave (`surface_waves[i]` not None) and whose beta values agree within
    `beta_spread_limit` of their mean. `surface_waves` is aligned with `candidates` (from
    `surface_wave_of_state`). Returns a `TwoStateChoice` carrying both states' beta/X_s and the
    mean beta used in the coded phase; raises `two_state_no_valid_pair` if no pair qualifies."""
    from ..errors import two_state_no_valid_pair

    candidates = list(candidates)
    surface_waves = list(surface_waves)
    if len(candidates) != len(surface_waves):
        raise ValueError("best_valid_pair: candidates and surface_waves must be the same length")
    best = None
    for i in range(len(candidates)):
        sw_i = surface_waves[i]
        if sw_i is None:
            continue
        phase_i = principal_phase_deg(candidates[i].patch_plane_value(frequency_hz))
        for j in range(i + 1, len(candidates)):
            sw_j = surface_waves[j]
            if sw_j is None:
                continue
            if pair_beta_spread_fraction(sw_i.beta_rad_per_m, sw_j.beta_rad_per_m) > beta_spread_limit:
                continue
            phase_j = principal_phase_deg(candidates[j].patch_plane_value(frequency_hz))
            sep = abs(wrapped_difference_deg(phase_i, phase_j))
            if best is None or sep > best[0]:
                best = (sep, i, j)
    if best is None:
        counts = ""
        if reasons is not None:
            tally = {}
            for reason in reasons:
                if reason is not None:
                    tally[reason] = tally.get(reason, 0) + 1
            counts = " Refused: " + ", ".join(f"{k} x{v}" for k, v in sorted(tally.items())) + "." if tally else ""
        raise two_state_no_valid_pair(
            message=f"none of {len(candidates)} candidate states forms a pair that carries the TM0 surface wave "
            f"with beta values within {beta_spread_limit:.0%} of each other.{counts}"
        )
    sep, i, j = best
    beta_a = surface_waves[i].beta_rad_per_m
    beta_b = surface_waves[j].beta_rad_per_m
    return TwoStateChoice(
        state_0=candidates[i],
        state_1=candidates[j],
        achieved_separation_deg=float(sep),
        surface_wave_0=surface_waves[i],
        surface_wave_1=surface_waves[j],
        beta_bar_rad_per_m=0.5 * (beta_a + beta_b),
        beta_spread_fraction=pair_beta_spread_fraction(beta_a, beta_b),
    )


def two_state_from_values(label_0, value_0, label_1, value_1, geometry_0=None, geometry_1=None):
    """Part E option (c): an HFSS-measured two-state table supplied directly by the user, rather
    than searched from one of this project's own libraries. The values are taken to be referenced
    at the patch plane (the user's own HFSS reference plane, stated in the report)."""
    a = TwoStateCell(label_0, "uploaded_state_0", complex(value_0), geometry_0 or {})
    b = TwoStateCell(label_1, "uploaded_state_1", complex(value_1), geometry_1 or {})
    sep = abs(wrapped_difference_deg(principal_phase_deg(a.value), principal_phase_deg(b.value)))
    return TwoStateChoice(state_0=a, state_1=b, achieved_separation_deg=float(sep))


def separation_warning(choice: TwoStateChoice, threshold_deg=SEPARATION_WARNING_THRESHOLD_DEG):
    """`None` if the achieved separation is within `threshold_deg` of the ideal 180deg; otherwise a
    plain-language warning string naming the actual shortfall -- never silently accepted as "close
    enough" (this project's own scientific-honesty convention, CLAUDE.md)."""
    shortfall = abs(IDEAL_SEPARATION_DEG - choice.achieved_separation_deg)
    if shortfall <= threshold_deg:
        return None
    return (
        f"These two cell states achieve only {choice.achieved_separation_deg:.1f}deg of phase "
        f"separation (ideal: 180deg) -- {shortfall:.1f}deg short. 1-bit coding will lose more "
        f"beam efficiency than an ideal pair; consider a different library or geometry pair."
    )


def lambda_sw_over_cell(beta_rad_per_m, period_m):
    """`lambda_sw / cell` -- the surface-wave wavelength in units of the cell pitch (the same
    `lambda_sw` `hologram/ports.py`'s isolation rule uses)."""
    return float((2.0 * np.pi / beta_rad_per_m) / period_m)
