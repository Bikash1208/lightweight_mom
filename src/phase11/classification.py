"""Sec. 18 (Phase 11N): classify the metasurface's spatial-coupling behavior into one of the four
named cases, from the measured `convergence.ConvergenceResult`s (never assumed beforehand, Sec. 4).
"""

from dataclasses import dataclass

from .convergence import CONVERGED, NOT_CONVERGED, STILL_SUBSTANTIALLY_CHANGING

CASE_A_LOCAL = "CASE_A_LOCAL"
CASE_B_MODERATELY_NONLOCAL = "CASE_B_MODERATELY_NONLOCAL"
CASE_C_STRONGLY_NONLOCAL = "CASE_C_STRONGLY_NONLOCAL"
CASE_D_PATTERN_DEPENDENT = "CASE_D_PATTERN_DEPENDENT"


@dataclass(frozen=True)
class ClassificationResult:
    per_center_case: dict
    overall_case: str
    pattern_dependence_detected: bool
    rationale: str


def _case_for_one_center(convergence_result, converged_at_small_n_threshold=3):
    if convergence_result.status == CONVERGED:
        if (convergence_result.converged_at_n or 99) <= converged_at_small_n_threshold:
            return CASE_A_LOCAL
        return CASE_B_MODERATELY_NONLOCAL
    if convergence_result.status == STILL_SUBSTANTIALLY_CHANGING:
        return CASE_C_STRONGLY_NONLOCAL
    # NOT_CONVERGED (measurable but not yet "substantial") by the largest tested size: still
    # changing, but not dramatically -- the honest middle ground is "moderately nonlocal," not a
    # silent upgrade to "local."
    return CASE_B_MODERATELY_NONLOCAL


def classify_overall(convergence_results, pattern_dependence_detected, converged_at_small_n_threshold=3):
    """Sec. 18: one case per center, plus an overall project classification. Sec. 18 Case D
    (pattern-dependent) is reported additionally/independently of the per-center local/nonlocal
    finding -- a metasurface can be simultaneously "converges by 9x9 for a checkerboard pattern"
    AND "a different arrangement with the same local statistics gives a substantially different
    answer," which is exactly what Phase 11F/11H/11I are designed to detect."""
    per_center = {
        cr.center_name: _case_for_one_center(cr, converged_at_small_n_threshold) for cr in convergence_results
    }
    # Sec. 18: the overall case is the worst (least local) one measured across centers -- a
    # holographic surface must accommodate whichever center geometry behaves the least locally.
    order = (CASE_A_LOCAL, CASE_B_MODERATELY_NONLOCAL, CASE_C_STRONGLY_NONLOCAL)
    worst = max(per_center.values(), key=lambda v: order.index(v) if v in order else 0)

    if pattern_dependence_detected:
        overall = CASE_D_PATTERN_DEPENDENT
        rationale = (
            "different spatial arrangements with similar local statistics produced significantly "
            "different center responses (Sec. 18 Case D) -- local descriptors (including the "
            "coupling-radius series alone) are insufficient regardless of the per-center "
            f"convergence finding ({worst})."
        )
    else:
        overall = worst
        rationale = f"worst-case per-center finding across all tested centers: {worst}."

    return ClassificationResult(per_center, overall, pattern_dependence_detected, rationale)
