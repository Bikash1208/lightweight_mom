"""Phase 8C Sec. 21-22: phase-magnitude coupling report and the quantitative phase-only
justification test.

Sec. 22 is explicit: "This is not a subjective decision. Define quantitative criteria before
applying the classification. If no defensible threshold exists from the formulation: REPORT
UNRESOLVED rather than inventing a threshold." Neither Phase 8B candidate formulation
(`holographic_modulation.py`) specifies a numeric magnitude tolerance anywhere, so this module
reports `PHASE_ONLY_APPROXIMATION_UNRESOLVED` for both -- not because the magnitude data itself is
ambiguous (it is reported in full, Sec. 21), but because no formulation-derived threshold exists to
judge it against.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class PhaseOnlyTestResult:
    component: str
    magnitude_range: tuple
    magnitude_median: float
    magnitude_std: float
    threshold_available: bool
    status: str  # "PHASE_ONLY_APPROXIMATION_JUSTIFIED" | "..._NOT_JUSTIFIED" | "..._UNRESOLVED"
    note: str


def phase_only_justification_test(candidate_coverage):
    """``candidate_coverage``: a `scattering_analysis.CandidateCoverage`. Neither Formulation A nor
    B (Phase 8B) defines a magnitude tolerance for treating their control variable as phase-only,
    so this always returns UNRESOLVED for this project today -- an honest gap, not a computed
    verdict from an invented number."""
    note = (
        f"Magnitude range [{candidate_coverage.magnitude_min:.4f}, "
        f"{candidate_coverage.magnitude_max:.4f}] (median {candidate_coverage.magnitude_median:.4f}, "
        f"std {candidate_coverage.magnitude_std:.4f}) is reported in full, but neither Phase 8B "
        "candidate formulation (holographic_modulation.py) specifies a quantitative magnitude "
        "tolerance for a phase-only approximation to be considered acceptable -- Formulation A "
        "needs the full complex interference kernel (magnitude matters); Formulation B's "
        "phase-difference kernel uses only phase by construction, but that does not itself say how "
        "much the *library's* magnitude may vary and still be treated as a a good approximation of "
        "unit magnitude. Reporting UNRESOLVED rather than inventing a threshold (Sec. 22)."
    )
    return PhaseOnlyTestResult(
        component=candidate_coverage.component,
        magnitude_range=(candidate_coverage.magnitude_min, candidate_coverage.magnitude_max),
        magnitude_median=candidate_coverage.magnitude_median,
        magnitude_std=candidate_coverage.magnitude_std,
        threshold_available=False,
        status="PHASE_ONLY_APPROXIMATION_UNRESOLVED",
        note=note,
    )
