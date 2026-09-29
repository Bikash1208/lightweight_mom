"""Frequency validation and derived free-space quantities (scalar, SI)."""

import math

from ._validate import finite_real
from .constants import C0


def validate_frequency(f_hz):
    """Return ``f_hz`` as a float, requiring a finite, strictly positive real."""
    f = finite_real("frequency", f_hz)
    if f <= 0.0:
        raise ValueError(f"frequency must be > 0 Hz, got {f!r}")
    return f


def angular_frequency(f_hz):
    """omega = 2*pi*f  [rad/s]."""
    return 2.0 * math.pi * validate_frequency(f_hz)


def free_space_wavenumber(f_hz):
    """k0 = omega / c0 = 2*pi*f / c0  [rad/m]."""
    return angular_frequency(f_hz) / C0


def free_space_wavelength(f_hz):
    """lambda0 = c0 / f  [m]."""
    return C0 / validate_frequency(f_hz)
