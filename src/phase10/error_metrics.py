"""Sec. 25-26: complex/relative/magnitude/phase error metrics plus the full percentile report
(mean, median, RMS, 95th percentile, max) Sec. 26 requires for every blind-test error
distribution -- "do not report only the mean." Built directly on ``phase9.error_metrics``'s
already-verified primitives (complex/magnitude/wrapped-phase error), not a re-implementation.
"""

from dataclasses import dataclass

import numpy as np

from ..phase9.error_metrics import complex_error, magnitude_error, phase_error_deg

RELATIVE_ERROR_EPSILON = 1e-15  # Sec. 25: explicitly documented denominator floor, used only to
# avoid a literal division by zero for a reference of exactly 0;
# far below any physically meaningful dipole-moment magnitude in
# this project (~1e-9 to 1e-8 A.m), so it never masks a real error.


def relative_complex_error(pred, ref):
    ref = np.asarray(ref, dtype=complex)
    return complex_error(pred, ref) / np.maximum(np.abs(ref), RELATIVE_ERROR_EPSILON)


@dataclass(frozen=True)
class ErrorDistributionReport:
    n: int
    mean_complex_error: float
    median_complex_error: float
    rms_complex_error: float
    p95_complex_error: float
    max_complex_error: float
    mean_relative_error: float
    p95_relative_error: float
    max_relative_error: float
    mean_magnitude_error: float
    p95_magnitude_error: float
    max_magnitude_error: float
    mean_phase_error_deg: float
    p95_phase_error_deg: float
    max_phase_error_deg: float


def _stats(values):
    v = np.asarray(values, dtype=float)
    return float(np.mean(v)), float(np.median(v)), float(np.percentile(v, 95)), float(np.max(v))


def compute_error_distribution(pred, ref):
    """Sec. 26's full report, computed once over an array of predicted/reference complex values
    (e.g. one per blind-test sample)."""
    pred = np.atleast_1d(np.asarray(pred, dtype=complex))
    ref = np.atleast_1d(np.asarray(ref, dtype=complex))
    if pred.shape != ref.shape:
        raise ValueError(f"compute_error_distribution: shape mismatch {pred.shape} vs {ref.shape}")

    ce = complex_error(pred, ref)
    re = relative_complex_error(pred, ref)
    me = magnitude_error(pred, ref)
    pe = phase_error_deg(pred, ref)

    ce_mean, ce_med, ce_p95, ce_max = _stats(ce)
    re_mean, _, re_p95, re_max = _stats(re)
    me_mean, _, me_p95, me_max = _stats(me)
    pe_mean, _, pe_p95, pe_max = _stats(pe)
    rms = float(np.sqrt(np.mean(ce**2)))

    return ErrorDistributionReport(
        n=int(pred.size),
        mean_complex_error=ce_mean,
        median_complex_error=ce_med,
        rms_complex_error=rms,
        p95_complex_error=ce_p95,
        max_complex_error=ce_max,
        mean_relative_error=re_mean,
        p95_relative_error=re_p95,
        max_relative_error=re_max,
        mean_magnitude_error=me_mean,
        p95_magnitude_error=me_p95,
        max_magnitude_error=me_max,
        mean_phase_error_deg=pe_mean,
        p95_phase_error_deg=pe_p95,
        max_phase_error_deg=pe_max,
    )
