"""Sec. 11-13: the requirement object, the (honest) candidate-period heuristic, and the
CANDIDATE-vs-VALIDATED distinction (Sec. 12's central rule: "A period is VALIDATED only when a
corresponding unit-cell electromagnetic dataset exists and satisfies the required validity
conditions"). Nothing here invents a universal "p = 0.2*lambda"-style answer (Sec. 12, 30): the
one heuristic this module offers is the Nyquist-Shannon sampling theorem applied to the
holographic modulation's own transverse spatial frequency (already the quantity
`phase8b.holographic_modulation.interference_modulation` computes as
``exp(-j*(k_ref - k_obj).r)``) -- a textbook sampling argument, not a project-specific number, and
it is always labelled `CANDIDATE_HEURISTIC`, never `VALIDATED`.
"""

import math
from dataclasses import dataclass

from ..feasibility_grating import grating_check
from ..phase8b.object_wave import transverse_wavenumber_deg
from .query import UnitCellQuery, search_unit_cells


@dataclass(frozen=True)
class HMSARequirements:
    """Sec. 11: a future design's actual electromagnetic requirements. Every field is optional
    (``None`` = unconstrained/not yet specified); ``fixed_period_mm`` in particular MUST stay
    `None` unless the caller explicitly wants a constrained-period design (Sec. 11's own rule) --
    normally period is an *output* of `candidate_periods_from_sampling`/`validate_period`, not an
    input here."""

    frequency_hz: float | None = None
    frequency_band_hz: tuple | None = None            # (min, max)
    substrate_name: str | None = None
    substrate_thickness_mm: float | None = None
    desired_polarization: str | None = None
    incidence_condition: str | None = None
    reference_field: object | None = None             # e.g. a phase8b.reference_wave.ReferenceField
    object_wave: object | None = None                 # e.g. a phase8b.object_wave.Beam or list thereof
    beam_directions_deg: tuple = ()                    # ((theta_deg, phi_deg), ...)
    required_complex_response: complex | None = None
    spatial_sampling_requirement: str | None = None
    fabrication_constraints: dict | None = None
    fixed_period_mm: tuple | None = None               # (p_x_mm, p_y_mm); None unless explicitly requested


@dataclass(frozen=True)
class PeriodAssessment:
    p_x_mm: float
    p_y_mm: float
    status: str                       # "VALIDATED" | "CANDIDATE_ONLY"
    grating_safe: bool
    grating_onset_hz: float | None
    matching_library_ids: tuple
    note: str


def required_transverse_spatial_frequency_per_m(beam_directions_deg, frequency_hz):
    """Sec. 13: the maximum transverse spatial frequency the holographic modulation must
    represent, from the actual k-vector formula (`phase8b.object_wave.transverse_wavenumber_deg`,
    itself cross-checked against `excitation.PlaneWaveExcitation`) -- the modulation pattern
    between any two beams varies as ``exp(-j*(k_a - k_b).r)`` (the same interference term
    `phase8b.holographic_modulation.interference_modulation` already computes), so its spatial
    frequency is the transverse k-vector *difference*, not either beam's k-vector alone. Raises
    ValueError with fewer than two beam directions (no interference term, no spatial frequency to
    sample)."""
    if len(beam_directions_deg) < 2:
        raise ValueError(
            "required_transverse_spatial_frequency_per_m: need at least two beam directions "
            "(e.g. reference + one object beam) to define an interference spatial frequency"
        )
    kvecs = [transverse_wavenumber_deg(frequency_hz, theta, phi) for theta, phi in beam_directions_deg]
    max_k = 0.0
    for i in range(len(kvecs)):
        for j in range(i + 1, len(kvecs)):
            kx_a, ky_a = kvecs[i]
            kx_b, ky_b = kvecs[j]
            max_k = max(max_k, math.hypot(kx_a - kx_b, ky_a - ky_b))
    return max_k


def candidate_periods_from_sampling(k_required_max_per_m, n_candidates=5):
    """Sec. 12-13: candidate (not validated) periods from the Nyquist-Shannon sampling theorem
    applied to ``k_required_max_per_m`` -- at least two unit cells per fringe of the fastest
    spatial variation the hologram needs, i.e. `P <= pi / k_required_max_per_m`. This is a
    textbook sampling bound (cited here, not invented as a project-specific number); it produces
    *candidates* to be checked against `validate_period`, not an answer by itself (Sec. 12's
    explicit CANDIDATE-vs-VALIDATED distinction). Raises ValueError if ``k_required_max_per_m`` is
    not positive (a zero/negative requirement has no finite Nyquist bound)."""
    if k_required_max_per_m <= 0.0:
        raise ValueError(f"k_required_max_per_m must be > 0, got {k_required_max_per_m!r}")
    p_max_m = math.pi / k_required_max_per_m
    fractions = [1.0, 0.8, 0.6, 0.4, 0.2][:max(1, n_candidates)]
    return tuple(
        {
            "p_mm": p_max_m * frac * 1e3,
            "status": "CANDIDATE_HEURISTIC",
            "rationale": (
                f"Nyquist-Shannon sampling of the required transverse spatial frequency "
                f"({k_required_max_per_m:.6g} rad/m): P <= pi/k_max = {p_max_m*1e3:.6g} mm; "
                f"this candidate uses {frac:g} of that bound. Not a validated period -- see "
                "validate_period."
            ),
        }
        for frac in fractions
    )


def validate_period(universal_db, topology, p_x_mm, p_y_mm, frequency_hz=None, substrate_name=None):
    """Sec. 12: a period is `VALIDATED` only if `universal_db` actually has a registered
    `LibraryEntry` at exactly this (topology, period[, frequency, substrate]) -- otherwise
    `CANDIDATE_ONLY`, regardless of how well it satisfies the sampling heuristic above. Also
    reports the *unit cell's own* grating-lobe safety at ``frequency_hz`` (normal incidence,
    `feasibility_grating.grating_check`, unmodified) as separate, additional information -- a
    period can be `CANDIDATE_ONLY` for lack of data yet still be grating-safe, or vice versa; the
    two are never conflated into one verdict."""
    query = UnitCellQuery(topology=topology, p_x_mm=p_x_mm, p_y_mm=p_y_mm, frequency_hz=frequency_hz,
                           substrate_name=substrate_name)
    result = search_unit_cells(universal_db, query)
    status = "VALIDATED" if result.status == "MATCH" else "CANDIDATE_ONLY"

    grating_safe, onset_hz = None, None
    if frequency_hz is not None:
        gc = grating_check(p_x_mm * 1e-3, p_y_mm * 1e-3, frequency_hz)
        grating_safe, onset_hz = (not gc.is_grating_risk), gc.grating_onset_hz

    note = (
        f"period ({p_x_mm:g}, {p_y_mm:g}) mm is {status} for topology={topology!r}"
        + (f", frequency={frequency_hz/1e9:g} GHz" if frequency_hz is not None else "")
        + f" (query filter trail: {result.filter_trail})."
    )
    return PeriodAssessment(
        p_x_mm=p_x_mm, p_y_mm=p_y_mm, status=status, grating_safe=grating_safe,
        grating_onset_hz=onset_hz, matching_library_ids=tuple(e.library_id for e in result.matching_entries),
        note=note,
    )
