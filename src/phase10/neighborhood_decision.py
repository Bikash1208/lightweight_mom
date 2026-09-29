"""Sec. 16: 5_CELL_SUFFICIENT / 9_CELL_REQUIRED / LARGER_NEIGHBORHOOD_REQUIRED, decided from
measured predictive error -- never chosen because a smaller model is cheaper (Sec. 16's explicit
prohibition). The two criteria reused below (``NOISE_FLOOR_RELATIVE``, ``SUBSTANTIAL_DEVIATION_
RELATIVE``) are not invented for this module: they are ``phase9.local_periodicity``'s own
already-sourced constants (a measured solver self-consistency floor, and a self-referential "at
least half the signal's own magnitude" substantiality test), reused rather than duplicated so
Phase 10's decision is held to exactly the same non-arbitrary standard Phase 9's classification
was.
"""

from dataclasses import dataclass

from ..phase9.local_periodicity import NOISE_FLOOR_RELATIVE, SUBSTANTIAL_DEVIATION_RELATIVE

NEIGHBORHOOD_DECISIONS = ("5_CELL_SUFFICIENT", "9_CELL_REQUIRED", "LARGER_NEIGHBORHOOD_REQUIRED")


@dataclass(frozen=True)
class FiveCellPredictiveError:
    """What a cardinal-only (5-cell) model misses, at one contrast level: it would predict the
    center response using only the 4 cardinal values (``s_5cell_prediction``, i.e. Test 5D's
    "N=E=S=W != C" case, which by construction leaves the diagonals at the center geometry) while
    the true environment also perturbs the 4 diagonals (``s_9cell_true``, the genuine 9-cell
    uniform reference)."""

    contrast_label: str
    s0: complex
    s_5cell_prediction: complex
    s_9cell_true: complex
    prediction_error: float
    relative_to_s0: float


def five_cell_predictive_error(contrast_label, s0, s_5cell_prediction, s_9cell_true):
    err = abs(s_9cell_true - s_5cell_prediction)
    return FiveCellPredictiveError(
        contrast_label,
        s0,
        s_5cell_prediction,
        s_9cell_true,
        err,
        err / abs(s0) if abs(s0) > 0 else float("inf"),
    )


def classify_per_contrast(predictive_error: FiveCellPredictiveError):
    """Sec. 16, per contrast level. `"5_CELL_SUFFICIENT"` only when ignoring the diagonals is
    indistinguishable from the solver's own numerical floor; `"9_CELL_REQUIRED"` when it is at
    least as large as half the center cell's own baseline response (self-evidently not
    negligible); otherwise the ambiguous middle ground is resolved conservatively toward
    `"9_CELL_REQUIRED"` (never silently assumed sufficient, matching this project's established
    discipline of not claiming validation without clear supporting evidence)."""
    if predictive_error.relative_to_s0 <= NOISE_FLOOR_RELATIVE:
        return "5_CELL_SUFFICIENT"
    return "9_CELL_REQUIRED"


@dataclass(frozen=True)
class NeighborhoodDecision:
    per_contrast: dict  # contrast_label -> FiveCellPredictiveError
    per_contrast_classification: dict
    long_range_relative: dict  # direction -> relative-to-s0 deviation at the tested contrast
    long_range_significant: bool
    final_decision: str
    rationale: str


def decide_neighborhood_size(predictive_errors_by_contrast, long_range_cases, s0_5x5):
    """``predictive_errors_by_contrast``: {label: FiveCellPredictiveError}. ``long_range_cases``:
    {direction: long_range.LongRangeCase}. The final decision uses the worst (most demanding)
    tested contrast level, since a holographic surface can contain extreme environments
    (Sec. 26's own principle, applied here to the neighborhood-size decision as well) --
    consistent with, not weaker than, any individual contrast level's own classification."""
    per_classification = {label: classify_per_contrast(pe) for label, pe in predictive_errors_by_contrast.items()}
    worst_label = max(predictive_errors_by_contrast, key=lambda k: predictive_errors_by_contrast[k].relative_to_s0)
    worst = predictive_errors_by_contrast[worst_label]
    worst_classification = per_classification[worst_label]

    long_range_relative = {
        d: (c.delta_s_abs / abs(s0_5x5) if abs(s0_5x5) > 0 else float("inf")) for d, c in long_range_cases.items()
    }
    long_range_significant = any(v > NOISE_FLOOR_RELATIVE for v in long_range_relative.values())

    if worst_classification == "5_CELL_SUFFICIENT":
        final = "5_CELL_SUFFICIENT"
        rationale = (
            f"worst-case contrast ({worst_label}): ignoring the diagonals changes the predicted "
            f"center response by only {worst.relative_to_s0:.2e} of its own baseline magnitude, "
            f"at or below the solver's measured self-consistency floor ({NOISE_FLOOR_RELATIVE:.0e})."
        )
    elif long_range_significant and worst.relative_to_s0 >= SUBSTANTIAL_DEVIATION_RELATIVE:
        final = "LARGER_NEIGHBORHOOD_REQUIRED"
        rationale = (
            f"worst-case contrast ({worst_label}): the diagonal-ignoring error is substantial "
            f"({worst.relative_to_s0:.3f} of the baseline magnitude, >= the "
            f"{SUBSTANTIAL_DEVIATION_RELATIVE:.1f} substantiality criterion) AND distance-2 "
            f"(two-period) neighbors still produce a measurable effect above the solver's noise "
            f"floor -- a 3x3 (9-cell) neighborhood does not appear to capture the full coupling."
        )
    else:
        final = "9_CELL_REQUIRED"
        rationale = (
            f"worst-case contrast ({worst_label}): the diagonal-ignoring error "
            f"({worst.relative_to_s0:.3e} of the baseline magnitude) exceeds the solver's "
            f"self-consistency floor, so the cardinal-only (5-cell) model is not sufficient; "
            f"distance-2 neighbor effects "
            + (
                "were also measurable but the worst-case deviation did not reach the "
                f"{SUBSTANTIAL_DEVIATION_RELATIVE:.1f} substantiality criterion, so a 3x3 (9-cell) "
                "neighborhood is reported as the required (not merely convenient) size."
                if long_range_significant
                else "remained at or below the solver's noise floor, so a 3x3 (9-cell) neighborhood "
                "is reported as sufficient beyond the immediate ring."
            )
        )

    return NeighborhoodDecision(
        predictive_errors_by_contrast,
        per_classification,
        long_range_relative,
        long_range_significant,
        final,
        rationale,
    )
