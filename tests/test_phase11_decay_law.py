"""Sec. 17, 28: no automatic decay-law fitting -- observe first, fit only when visibly monotonic,
and always report residuals/failures alongside any fitted law."""

import pytest

from src.phase11.decay_law import fit_exponential_if_monotonic, observe_decay


def test_observe_decay_sorts_and_detects_monotonic_trend():
    pairs = [(3, 0.5), (1, 2.0), (2, 1.0)]
    obs = observe_decay(pairs)
    assert obs.distances == (1, 2, 3)
    assert obs.magnitudes == (2.0, 1.0, 0.5)
    assert obs.monotonically_decreasing is True
    assert obs.candidate_trend == "APPROXIMATELY_MONOTONIC_DECAY"


def test_observe_decay_flags_non_monotonic_data_honestly():
    pairs = [(1, 1.0), (2, 3.0), (3, 0.5)]
    obs = observe_decay(pairs)
    assert obs.monotonically_decreasing is False
    assert obs.candidate_trend == "NON_MONOTONIC_NO_CLEAR_TREND"


def test_fit_refuses_when_trend_is_not_monotonic():
    pairs = [(1, 1.0), (2, 3.0), (3, 0.5)]
    obs = observe_decay(pairs)
    fit = fit_exponential_if_monotonic(obs)
    assert fit.law == "NONE_FITTED"
    assert "no decay law was fitted" in fit.note


def test_fit_reports_parameters_and_quality_for_a_clean_exponential():
    import numpy as np
    distances = [1, 2, 3, 4, 5]
    magnitudes = [10.0 * np.exp(-d / 2.0) for d in distances]
    obs = observe_decay(list(zip(distances, magnitudes)))
    fit = fit_exponential_if_monotonic(obs)
    assert fit.law == "EXPONENTIAL_CANDIDATE_FIT"
    assert fit.r_squared > 0.99
    assert fit.fitted_params["length_scale_cells"] == pytest.approx(2.0, rel=1e-2)
