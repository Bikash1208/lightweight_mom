"""The TM0 surface-wave dispersion of a bare grounded dielectric slab (air | dielectric | PEC
ground, no patches) -- closed-form, standard result (e.g. Pozar, *Microwave Engineering*, grounded
dielectric slab waveguide). Used as the field-import page's sanity-check baseline
(PHASE13_WEBAPP_PROMPT.md: "compare with the analytic TM0 of the bare grounded slab... as a sanity
check, not a requirement").

**Dispersion equation** (TM0, the fundamental mode -- no cutoff frequency): for a slab of
thickness `h` and relative permittivity `eps_r`, propagation constant `beta` along the surface,

    kz0 = sqrt(eps_r*k0^2 - beta^2)        (transverse wavenumber inside the dielectric)
    alpha0 = sqrt(beta^2 - k0^2)            (transverse decay constant in air, beta > k0)
    kz0 * tan(kz0*h) = eps_r * alpha0

solved here for `beta` in `(k0, sqrt(eps_r)*k0)` by bisection (`scipy.optimize.brentq`) -- the
unique TM0 root in that interval for the electrically-thin slabs this project uses (verified
against PHASE13_WEBAPP_PROMPT.md's own cited numbers for RT/Duroid 5880 at 1.575 mm: this
function reproduces beta = 213/269/326 rad/m at 10/12.5/15 GHz to within 0.4 rad/m,
`tests/test_dispersion_analytic_tm0.py`).
"""

import numpy as np
from scipy.optimize import brentq

from ..frequency import free_space_wavenumber


def tm0_beta(frequency_hz, eps_r, thickness_m, bracket_shrink=1e-7):
    """The TM0 surface-wave propagation constant [rad/m] of a bare grounded dielectric slab.
    Raises `ValueError` if no root is found in `(k0, sqrt(eps_r)*k0)` (e.g. `eps_r <= 1`)."""
    if eps_r <= 1.0:
        raise ValueError(f"eps_r must be > 1 for a grounded-slab surface wave to exist, got {eps_r!r}")
    k0 = free_space_wavenumber(frequency_hz)

    def dispersion_equation(beta):
        kz0 = np.sqrt(eps_r * k0**2 - beta**2)
        alpha0 = np.sqrt(beta**2 - k0**2)
        return kz0 * np.tan(kz0 * thickness_m) - eps_r * alpha0

    lo = k0 * (1.0 + bracket_shrink)
    hi = k0 * np.sqrt(eps_r) * (1.0 - bracket_shrink)
    if lo >= hi:
        raise ValueError(f"no valid (k0, sqrt(eps_r)*k0) bracket for eps_r={eps_r!r}, thickness={thickness_m!r}")
    return float(brentq(dispersion_equation, lo, hi))


def tm0_n_eff(frequency_hz, eps_r, thickness_m):
    """beta / k0 -- the surface wave's effective index, > 1 by construction (bound mode)."""
    k0 = free_space_wavenumber(frequency_hz)
    return tm0_beta(frequency_hz, eps_r, thickness_m) / k0
