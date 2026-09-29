"""Sec. 17 (Phase 11M): report the directly computed interaction magnitude vs. distance FIRST,
never automatically fitting a decay law. A candidate analytical trend is only investigated if the
plotted data visibly support one, and even then every fitted parameter, interval, residual and
failure case is reported alongside it -- never presented as an electromagnetic law without physical
justification.
"""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class DecayObservation:
    distances: tuple
    magnitudes: tuple
    monotonically_decreasing: bool
    candidate_trend: str   # "NONE_ASSESSED" | "APPROXIMATELY_MONOTONIC_DECAY" | "NON_MONOTONIC_NO_CLEAR_TREND"


def observe_decay(distance_magnitude_pairs):
    """Sec. 17: no fitting -- just the raw, sorted (distance, |interaction|) series and a plain
    monotonicity check (the only "trend" claim this function is willing to make without a human
    inspecting the plot first)."""
    pairs = sorted(distance_magnitude_pairs, key=lambda p: p[0])
    distances = tuple(p[0] for p in pairs)
    magnitudes = tuple(p[1] for p in pairs)
    is_monotonic = all(magnitudes[i] >= magnitudes[i + 1] for i in range(len(magnitudes) - 1))
    trend = "APPROXIMATELY_MONOTONIC_DECAY" if is_monotonic else "NON_MONOTONIC_NO_CLEAR_TREND"
    return DecayObservation(distances, magnitudes, is_monotonic, trend)


@dataclass(frozen=True)
class DecayFitReport:
    law: str
    fitted_params: dict
    fitting_interval: tuple
    residual_rms: float
    r_squared: float
    failure_cases: tuple
    note: str


def fit_exponential_if_monotonic(observation: DecayObservation):
    """Sec. 17: ONLY invoked when ``observation.candidate_trend ==
    'APPROXIMATELY_MONOTONIC_DECAY'`` (the caller's responsibility to check -- this function
    still re-checks and refuses otherwise). Fits |R(d)| = A * exp(-d / L) by ordinary least squares
    on log-magnitude, reports fit quality and residuals explicitly, and is never described as an
    established electromagnetic law -- only as a compact numerical summary of what was measured."""
    if observation.candidate_trend != "APPROXIMATELY_MONOTONIC_DECAY":
        return DecayFitReport(
            "NONE_FITTED", {}, (min(observation.distances), max(observation.distances)),
            float("nan"), float("nan"), observation.distances,
            "the data do not show a visibly monotonic trend; no decay law was fitted (Sec. 17).",
        )
    d = np.asarray(observation.distances, dtype=float)
    m = np.asarray(observation.magnitudes, dtype=float)
    positive = m > 0
    if positive.sum() < 2:
        return DecayFitReport("NONE_FITTED", {}, (float(d.min()), float(d.max())), float("nan"),
                               float("nan"), observation.distances,
                               "fewer than 2 positive magnitude points; no fit attempted.")
    log_m = np.log(m[positive])
    d_fit = d[positive]
    slope, intercept = np.polyfit(d_fit, log_m, 1)
    predicted = np.exp(intercept + slope * d_fit)
    residual_rms = float(np.sqrt(np.mean((m[positive] - predicted) ** 2)))
    ss_res = float(np.sum((m[positive] - predicted) ** 2))
    ss_tot = float(np.sum((m[positive] - np.mean(m[positive])) ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    length_scale = -1.0 / slope if slope != 0 else float("inf")
    failures = tuple(observation.distances[i] for i in range(len(m)) if not positive[i])
    return DecayFitReport(
        "EXPONENTIAL_CANDIDATE_FIT", {"amplitude": float(np.exp(intercept)), "length_scale_cells": length_scale},
        (float(d_fit.min()), float(d_fit.max())), residual_rms, r_squared, failures,
        "a compact numerical summary of the measured trend, NOT an asserted electromagnetic decay law.",
    )
