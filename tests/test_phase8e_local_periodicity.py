"""Sec. 15: local-periodic approximation is always reported unvalidated."""

import math

import pytest

from src.phase8e.local_periodicity import assess_local_periodicity


def test_always_reports_unvalidated():
    result = assess_local_periodicity(4.8e-3, 100.0)
    assert result.status == "LOCAL_PERIODICITY_UNVALIDATED"


def test_cells_per_fringe_is_computed_correctly():
    period_m = 4.8e-3
    k_t = 2.0 * math.pi / 0.024  # fringe period = 24 mm
    result = assess_local_periodicity(period_m, k_t)
    assert result.cells_per_fringe == pytest.approx(0.024 / period_m)


def test_none_k_t_gives_none_cells_per_fringe_but_still_unvalidated():
    result = assess_local_periodicity(4.8e-3, None)
    assert result.cells_per_fringe is None
    assert result.status == "LOCAL_PERIODICITY_UNVALIDATED"


def test_zero_k_t_does_not_raise_a_divide_by_zero():
    result = assess_local_periodicity(4.8e-3, 0.0)
    assert result.cells_per_fringe is None
