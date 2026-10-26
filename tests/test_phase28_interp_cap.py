"""`_interp_cap` of the Phase 28 sweep script: the smallest cap (log-interpolated between grid points) whose loss is within the tolerance."""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "research"))

import phase28_feasibility_sweep as sw  # noqa: E402


def _curve(caps, losses):
    return [{"alpha_cap_np_per_m": c, "loss": v} for c, v in zip(caps, losses, strict=True)]


def test_returns_the_log_interpolated_crossing_on_a_monotone_curve():
    cv = _curve([1.0, 2.0, 4.0], [3.0, 1.5, 0.5])
    cap = sw._interp_cap(cv, "loss", 1.0)  # halfway in loss between 1.5 and 0.5 -> halfway in log(cap) between 2 and 4
    assert cap == math.exp(0.5 * (math.log(2.0) + math.log(4.0)))
    # the result lies between the neighbouring grid points and its linear-in-loss interpolant equals the tolerance
    assert 2.0 < cap < 4.0


def test_returns_the_first_grid_cap_when_even_it_is_within_tolerance():
    assert sw._interp_cap(_curve([1.0, 2.0], [0.4, 0.2]), "loss", 1.0) == 1.0


def test_returns_none_when_the_largest_cap_still_exceeds_the_tolerance():
    assert sw._interp_cap(_curve([1.0, 2.0, 4.0], [3.0, 2.0, 1.2]), "loss", 1.0) is None


def test_a_tolerance_equal_to_a_grid_loss_returns_that_grid_cap():
    cap = sw._interp_cap(_curve([1.0, 2.0, 4.0], [3.0, 1.5, 1.0]), "loss", 1.0)
    assert abs(cap - 4.0) < 1e-12


def test_non_monotone_curve_returns_the_first_crossing_documented_behaviour():
    cv = _curve([1.0, 2.0, 4.0, 8.0], [3.0, 0.8, 1.4, 0.9])
    cap = sw._interp_cap(cv, "loss", 1.0)
    assert 1.0 < cap < 2.0  # the first crossing, not the later one
