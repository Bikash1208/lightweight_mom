"""Sec. 16, 47: the non-arbitrary (reused, sourced) 5-cell/9-cell/larger-neighborhood decision."""

from src.phase10.neighborhood_decision import (NOISE_FLOOR_RELATIVE, SUBSTANTIAL_DEVIATION_RELATIVE,
                                                classify_per_contrast, decide_neighborhood_size,
                                                five_cell_predictive_error)


def test_reuses_phase9s_own_constants_not_new_ones():
    from src.phase9.local_periodicity import NOISE_FLOOR_RELATIVE as p9_floor
    from src.phase9.local_periodicity import SUBSTANTIAL_DEVIATION_RELATIVE as p9_sub
    assert NOISE_FLOOR_RELATIVE == p9_floor
    assert SUBSTANTIAL_DEVIATION_RELATIVE == p9_sub


def test_five_cell_sufficient_when_prediction_matches_the_true_9cell_response():
    s0 = 1.0 + 0j
    pe = five_cell_predictive_error("weak", s0, 1.0000001 + 0j, 1.0000001 + 0j)
    assert classify_per_contrast(pe) == "5_CELL_SUFFICIENT"


def test_nine_cell_required_when_prediction_error_is_measurable_but_not_substantial():
    s0 = 1.0 + 0j
    pe = five_cell_predictive_error("moderate", s0, 1.0 + 0j, 1.1 + 0j)   # 10% relative error
    assert classify_per_contrast(pe) == "9_CELL_REQUIRED"


def test_decide_neighborhood_size_escalates_to_larger_when_long_range_also_significant():
    s0 = 1.0 + 0j
    predictive_errors = {
        "weak": five_cell_predictive_error("weak", s0, 1.0 + 0j, 1.0 + 0j),
        "strong_high": five_cell_predictive_error("strong_high", s0, 1.0 + 0j, 2.0 + 0j),  # 100% error
    }

    class _LR:
        def __init__(self, abs_):
            self.delta_s_abs = abs_

    long_range_cases = {"N2": _LR(0.1), "S2": _LR(0.1), "E2": _LR(0.1), "W2": _LR(0.1)}
    decision = decide_neighborhood_size(predictive_errors, long_range_cases, s0)
    assert decision.final_decision == "LARGER_NEIGHBORHOOD_REQUIRED"
    assert decision.long_range_significant is True


def test_decide_neighborhood_size_stays_9cell_when_long_range_is_negligible():
    s0 = 1.0 + 0j
    predictive_errors = {
        "weak": five_cell_predictive_error("weak", s0, 1.0 + 0j, 1.0 + 0j),
        "moderate": five_cell_predictive_error("moderate", s0, 1.0 + 0j, 1.05 + 0j),
    }

    class _LR:
        def __init__(self, abs_):
            self.delta_s_abs = abs_

    long_range_cases = {"N2": _LR(1e-9), "S2": _LR(1e-9), "E2": _LR(1e-9), "W2": _LR(1e-9)}
    decision = decide_neighborhood_size(predictive_errors, long_range_cases, s0)
    assert decision.final_decision == "9_CELL_REQUIRED"
    assert decision.long_range_significant is False
