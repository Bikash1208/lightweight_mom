"""Sec. 13-14: complex-plane (never phase-coverage-alone) controllability of the actual unit-cell
library against a required response set R_required. Reuses `phase8c.scattering_analysis` (complex-
state/coverage statistics) and `phase8d_rev.response_matching` (per-point nearest-neighbour
complex-distance search) -- both unmodified.
"""

from dataclasses import dataclass

import numpy as np

from ..phase import wrapped_difference_deg
from ..phase8c.scattering_analysis import candidate_coverage, distinct_complex_and_phase_states
from ..phase8d_rev.response_matching import nearest_response_match


@dataclass(frozen=True)
class LibraryComplexPlaneSummary:
    component: str
    n_geometries: int
    n_distinct_complex_states: int
    n_distinct_phase_states: int
    amplitude_min: float
    amplitude_max: float
    phase_coverage_deg: float
    phase_largest_gap_deg: float
    phase_magnitude_correlation: float | None


def library_complex_plane_summary(rows, component_name):
    """Sec. 13: Re/Im/|C|/angle(C) distribution over the available geometry states, plus 2-D
    complex-state density -- computed from `phase8c.scattering_analysis`, never phase alone."""
    states = distinct_complex_and_phase_states(rows)
    cov = candidate_coverage(rows, component_name)
    return LibraryComplexPlaneSummary(
        component=component_name,
        n_geometries=states.n_geometries,
        n_distinct_complex_states=states.n_distinct_complex_states,
        n_distinct_phase_states=states.n_distinct_phase_states,
        amplitude_min=cov.magnitude_min,
        amplitude_max=cov.magnitude_max,
        phase_coverage_deg=cov.coverage_deg,
        phase_largest_gap_deg=cov.largest_gap_deg,
        phase_magnitude_correlation=cov.phase_magnitude_correlation,
    )


@dataclass(frozen=True)
class PhaseOnlyControllability:
    n_required_points: int
    mean_phase_error_deg: float
    max_phase_error_deg: float
    uncontrolled_amplitude_min: float
    uncontrolled_amplitude_max: float
    note: str


def phase_only_controllability(required_phase_deg, rows):
    """Sec. 13-14 for a phase-only requirement (Candidates B/C): the nearest achievable phase for
    every required point, and the amplitude range that comes along *uncontrolled* with whichever
    geometry realizes that phase (Sec. 14's explicit warning: 'if the formulation requires
    amplitude as well as phase, phase-only coverage is insufficient' -- reported here as a fact
    about the library, not hidden by only reporting phase error)."""
    required = np.atleast_1d(np.asarray(required_phase_deg, dtype=float)).ravel()
    lib_phases = np.array([r.phase_deg for r in rows])
    lib_mags = np.array([r.mag for r in rows])
    errors = np.array([np.min(np.abs(wrapped_difference_deg(lib_phases, target))) for target in required])
    return PhaseOnlyControllability(
        n_required_points=int(required.size),
        mean_phase_error_deg=float(np.mean(errors)),
        max_phase_error_deg=float(np.max(errors)),
        uncontrolled_amplitude_min=float(lib_mags.min()),
        uncontrolled_amplitude_max=float(lib_mags.max()),
        note="phase-only formulation: amplitude at each required point is whatever the "
        "phase-matched geometry happens to have, not separately controlled.",
    )


@dataclass(frozen=True)
class ComplexPlaneControllability:
    n_required_points: int
    mean_complex_error: float
    max_complex_error: float
    mean_amplitude_error: float
    max_amplitude_error: float
    mean_phase_error_deg: float
    max_phase_error_deg: float


def complex_plane_controllability(required_complex, rows):
    """Sec. 13-14 for a full complex requirement (Candidate D's per-point transfer-operator
    entries): the nearest-neighbour complex-plane match for every required point
    (`phase8d_rev.response_matching.nearest_response_match`, reused), reporting amplitude, phase
    and complex error separately (never collapsed, Sec. 14's own rule, `response_matching.py`'s
    module docstring)."""
    required = np.atleast_1d(np.asarray(required_complex, dtype=complex)).ravel()
    candidates = [(f"L{r.l_mm:g}_W{r.w_mm:g}", r.l_mm, r.w_mm, "lib", complex(r.real, r.imag)) for r in rows]
    complex_errors, amp_errors, phase_errors = [], [], []
    for target in required:
        best = nearest_response_match(target, candidates)[0]
        complex_errors.append(best.complex_error)
        amp_errors.append(best.amplitude_error)
        phase_errors.append(best.phase_error_deg)
    return ComplexPlaneControllability(
        n_required_points=int(required.size),
        mean_complex_error=float(np.mean(complex_errors)),
        max_complex_error=float(np.max(complex_errors)),
        mean_amplitude_error=float(np.mean(amp_errors)),
        max_amplitude_error=float(np.max(amp_errors)),
        mean_phase_error_deg=float(np.mean(phase_errors)),
        max_phase_error_deg=float(np.max(phase_errors)),
    )
