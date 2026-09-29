"""Part 3: the Verify page's cost estimate/refusal gate (fast) and one real end-to-end check at
the cheapest measured size (N=6, ~7s)."""

import pytest

from src.verify import supercell_check as vc


def test_cost_estimate_at_measured_points_is_exact():
    for n, t in vc._MEASURED_RUNTIMES:
        cost = vc.estimate_cost(n)
        assert cost.is_measured
        assert cost.estimated_runtime_s == pytest.approx(t)
        assert not cost.refused


def test_cost_extrapolates_above_max_measured():
    cost = vc.estimate_cost(16)
    assert cost.is_extrapolated
    assert cost.estimated_runtime_s > vc._MEASURED_RUNTIMES[-1][1]
    assert not cost.refused


def test_cost_refuses_far_beyond_measured_range():
    cost = vc.estimate_cost(25)
    assert cost.refused
    assert "no measured data" in cost.refusal_reason


def test_run_verification_raises_on_refused_size_unless_forced():
    with pytest.raises(ValueError, match="no measured data"):
        vc.run_verification(25)


@pytest.mark.slow
def test_run_verification_end_to_end_smallest_size():
    result = vc.run_verification(6)
    assert result.n_cells == 6
    assert result.cost.is_measured
    assert result.predicted_vs_measured_intended_fraction_diff < 0.05
    # PHASE13B A10: the predicted intended-order fraction must be available directly (previously
    # the app's summary line showed a blank value here, only the |predicted - measured| diff).
    assert abs(result.predicted_intended_fraction - result.aperture_result.fraction_intended) == pytest.approx(
        result.predicted_vs_measured_intended_fraction_diff, abs=1e-9
    )
