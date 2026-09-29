"""Phase wrapping and unwrapping (degrees; this is an I/O-boundary module, see units.py).

Conventions
-----------
* Principal value: every stored phase lies in the half-open interval (-180, 180] degrees.
  ``+180`` is used for the negative real axis, never ``-180``: the two are the same
  physical phase and a database must not contain both spellings of it.
* Phase is ``arg(S)``, the angle of the complex phasor under exp(+j w t); a delay is a
  negative phase (CONVENTIONS.md).
* Unwrapping: given samples p_0, p_1, ... taken along an ordered axis (frequency, or a
  swept geometry parameter), the unwrapped phase is

      u_0 = p_0,       u_i = u_{i-1} + wrap(p_i - p_{i-1}).

  It is the continuous phase *only if* the true phase changes by less than 180 deg between
  neighbouring samples. Nothing in the data can reveal a violation: an undersampled sweep
  unwraps cleanly to a wrong answer. ``max_wrapped_step_deg`` reports the largest wrapped
  step so a caller can refuse suspicious data (a step near 180 deg means the sampling is
  at the limit; steps above ~90 deg deserve a finer grid).
* The unwrapped phase is anchored to the first sample; it is only defined up to a global
  multiple of 360 deg.
* A phase step of exactly +-180 deg is genuinely ambiguous; ``wrap`` resolves it to +180.
"""

import numpy as np

# Default warning threshold for the largest wrapped step, half of the 180 deg ambiguity limit.
STEP_WARN_DEG = 90.0


def _finite_array(name, values):
    a = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(a)):
        raise ValueError(f"{name} contains non-finite values")
    return a


def wrap_deg(phase_deg):
    """Map angles [deg] to the principal interval (-180, 180]. Scalars in, scalars out."""
    p = _finite_array("phase", phase_deg)
    w = p - 360.0 * np.ceil(p / 360.0 - 0.5)
    return float(w) if w.ndim == 0 else w


def principal_phase_deg(z):
    """arg(z) in degrees on (-180, 180]; ``z`` complex scalar or array (non-zero expected)."""
    z = np.asarray(z, dtype=complex)
    if not np.all(np.isfinite(z)):
        raise ValueError("phasor contains non-finite values")
    return wrap_deg(np.degrees(np.angle(z)))


def wrapped_difference_deg(a_deg, b_deg):
    """Shortest signed angular distance a - b on (-180, 180]: 179 - (-179) = -2, not 358."""
    return wrap_deg(_finite_array("a", a_deg) - _finite_array("b", b_deg))


def wrapped_steps_deg(phase_deg, axis=-1):
    """Wrapped differences between neighbouring samples along ``axis``."""
    p = _finite_array("phase", phase_deg)
    return wrap_deg(np.diff(p, axis=axis))


def max_wrapped_step_deg(phase_deg, axis=-1):
    """Largest |wrapped step| between neighbouring samples along ``axis`` (0 for < 2 samples)."""
    steps = wrapped_steps_deg(phase_deg, axis)
    return float(np.max(np.abs(steps))) if steps.size else 0.0


def unwrap_deg(phase_deg, axis=-1):
    """Continuous phase along ``axis`` (see the module docstring for when it is valid)."""
    p = _finite_array("phase", phase_deg)
    if p.ndim == 0:
        raise ValueError("unwrap needs an array with an axis to unwrap along")
    if p.shape[axis] == 0:
        return p.copy()
    first = np.take(p, [0], axis=axis)
    steps = wrapped_steps_deg(p, axis)
    return np.concatenate([first, first + np.cumsum(steps, axis=axis)], axis=axis)
