"""Sec. 16, 28: the convergence criterion, declared before evaluation."""

import pytest

from src.phase11.convergence import CONVERGED, NOT_CONVERGED, STILL_SUBSTANTIALLY_CHANGING, classify_convergence
from src.phase11.coupling_radius import ShellIncrement


def test_reuses_phase9s_own_sourced_constants():
    from src.phase9.local_periodicity import NOISE_FLOOR_RELATIVE as p9_floor
    from src.phase9.local_periodicity import SUBSTANTIAL_DEVIATION_RELATIVE as p9_sub
    from src.phase11.convergence import NOISE_FLOOR_RELATIVE, SUBSTANTIAL_DEVIATION_RELATIVE

    assert p9_floor == NOISE_FLOOR_RELATIVE
    assert p9_sub == SUBSTANTIAL_DEVIATION_RELATIVE


def test_converged_when_last_increment_at_noise_floor():
    incs = [ShellIncrement(1, 3, 1e-9, 0.1), ShellIncrement(3, 5, 1e-12, 1e-7)]
    result = classify_convergence(incs, "C_test")
    assert result.status == CONVERGED
    assert result.converged_at_n == 5


def test_still_substantially_changing_when_any_increment_is_large():
    incs = [ShellIncrement(1, 3, 1e-9, 0.6), ShellIncrement(3, 5, 1e-10, 0.05)]
    result = classify_convergence(incs, "C_test")
    assert result.status == STILL_SUBSTANTIALLY_CHANGING


def test_not_converged_middle_ground():
    incs = [ShellIncrement(1, 3, 1e-9, 0.1), ShellIncrement(3, 5, 1e-10, 0.05)]
    result = classify_convergence(incs, "C_test")
    assert result.status == NOT_CONVERGED


def test_raises_on_empty_increments():
    with pytest.raises(ValueError):
        classify_convergence([], "C_test")
