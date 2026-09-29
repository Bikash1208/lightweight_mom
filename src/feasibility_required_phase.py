"""Phase 8A Sec. 17-18: the required-holographic-phase computation step, and generic statistics
over whatever required-phase array is supplied.

Per the brief: "The gate must NOT implement a second holographic formulation. Use the
authoritative existing implementation." This project has none
(`feasibility_phase_source.py`'s module docstring lists the inspected evidence).
`compute_required_phase` therefore never derives a `required_phase(x, y)` array from a reference
field/beam set -- it only reports, honestly, that this step cannot be performed here, which is
itself the correct Sec. 38 failure ("missing holographic formulation"), not a workaround.

`required_phase_statistics`, in contrast, is a pure function of *any* already-available
required-phase array (real, once a generator exists, or a synthetic one supplied directly by a
test or a caller who computed it elsewhere) -- Sec. 18's statistics do not themselves require a
generator, only a phase array, so they are implemented in full and exercised by synthetic inputs
(brief Sec. 44 Tests 10-11).
"""

from dataclasses import dataclass

import numpy as np

from .feasibility_phase_source import HOLOGRAPHIC_GENERATOR_IMPLEMENTED
from .phase_extraction import phase_coverage


@dataclass(frozen=True)
class RequiredPhaseResult:
    available: bool
    status: str            # "HOLOGRAPHIC_FORMULATION_MISSING" (today) | "OK" (once implemented)
    note: str


def compute_required_phase(reference_field_validation, beam_validations, polarization):
    """Would call into this project's holographic generator, if one existed. It does not; this
    function is the explicit, documented stop (Sec. 17), never a substitute formulation."""
    if HOLOGRAPHIC_GENERATOR_IMPLEMENTED:
        raise NotImplementedError(
            "a holographic generator now exists in this project but "
            "feasibility_required_phase.compute_required_phase was not updated to call it"
        )
    return RequiredPhaseResult(
        available=False,
        status="HOLOGRAPHIC_FORMULATION_MISSING",
        note=(
            "This project has no implemented holographic-synthesis formulation (README.md, "
            "DATABASE.md, src/__init__.py). Phase 8A does not implement a second one to fill the "
            "gap (brief Sec. 17, Sec. 49). The required holographic phase distribution "
            "required_phase(x, y) therefore cannot be computed for this design, independent of "
            "whether the supplied reference field and beam definitions are themselves valid."
        ),
    )


@dataclass(frozen=True)
class RequiredPhaseStatistics:
    n_samples: int
    min_required_phase_deg: float
    max_required_phase_deg: float
    circular_coverage_deg: float
    largest_gap_deg: float
    mean_phase_deg: float          # circular mean (atan2 of mean sin/cos), not an arithmetic mean


def required_phase_statistics(required_phase_deg):
    """Sec. 18: statistics over an already-available required-phase array, phase treated as
    circular throughout (`phase_extraction.phase_coverage`, reused rather than re-derived)."""
    phases = np.asarray(required_phase_deg, dtype=float)
    if phases.size == 0:
        raise ValueError("required_phase_statistics: empty input")
    table = [{"phase_deg": float(p)} for p in phases]
    cov = phase_coverage(table)
    mean_deg = float(np.degrees(np.arctan2(
        np.mean(np.sin(np.radians(phases))), np.mean(np.cos(np.radians(phases))),
    )))
    return RequiredPhaseStatistics(
        n_samples=int(phases.size),
        min_required_phase_deg=cov.phase_min_deg,
        max_required_phase_deg=cov.phase_max_deg,
        circular_coverage_deg=cov.coverage_deg,
        largest_gap_deg=cov.largest_gap_deg,
        mean_phase_deg=mean_deg,
    )
