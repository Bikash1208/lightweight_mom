"""Sec. 3, 9-12: complex/magnitude/phase error metrics shared by 9A (solver vs. HFSS) and usable
anywhere a complex-quantity comparison is needed. Phase is always compared wrapped
(`phase.wrapped_difference_deg`, Phase 1, unmodified) -- never a naive subtraction across the
+-180 deg boundary (Sec. 12's explicit rule).
"""

from dataclasses import dataclass

import numpy as np

from ..phase import principal_phase_deg, wrapped_difference_deg


@dataclass(frozen=True)
class ComplexErrorStats:
    n_points: int
    mean_complex_error: float
    rms_complex_error: float
    max_complex_error: float
    mean_magnitude_error: float
    rms_magnitude_error: float
    max_magnitude_error: float
    mean_abs_phase_error_deg: float
    max_abs_phase_error_deg: float


def complex_error(a, b):
    """|a - b|, element-wise. Rule 3/4: this is a physical amplitude-and-phase difference,
    never normalized to hide a magnitude mismatch."""
    return np.abs(np.asarray(a, dtype=complex) - np.asarray(b, dtype=complex))


def magnitude_error(a, b):
    return np.abs(np.abs(np.asarray(a, dtype=complex)) - np.abs(np.asarray(b, dtype=complex)))


def phase_error_deg(a, b):
    """Wrapped |arg(a) - arg(b)| [deg] (Sec. 12)."""
    phase_a = principal_phase_deg(np.asarray(a, dtype=complex))
    phase_b = principal_phase_deg(np.asarray(b, dtype=complex))
    return np.abs(wrapped_difference_deg(phase_a, phase_b))


def compute_error_stats(a, b):
    a, b = np.atleast_1d(np.asarray(a, dtype=complex)), np.atleast_1d(np.asarray(b, dtype=complex))
    if a.shape != b.shape:
        raise ValueError(f"compute_error_stats: shape mismatch {a.shape} vs {b.shape}")
    ce, me, pe = complex_error(a, b), magnitude_error(a, b), phase_error_deg(a, b)
    return ComplexErrorStats(
        n_points=int(a.size),
        mean_complex_error=float(np.mean(ce)),
        rms_complex_error=float(np.sqrt(np.mean(ce**2))),
        max_complex_error=float(np.max(ce)),
        mean_magnitude_error=float(np.mean(me)),
        rms_magnitude_error=float(np.sqrt(np.mean(me**2))),
        max_magnitude_error=float(np.max(me)),
        mean_abs_phase_error_deg=float(np.mean(pe)),
        max_abs_phase_error_deg=float(np.max(pe)),
    )
