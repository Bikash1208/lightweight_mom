"""Shared argument validation."""

import math
import numbers


def finite_real(name, value):
    """Return ``value`` as float; reject bool, complex, NaN, inf, non-numbers."""
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise ValueError(f"{name} must be a finite real number, got {value!r}")
    v = float(value)
    if not math.isfinite(v):
        raise ValueError(f"{name} must be a finite real number, got {value!r}")
    return v
