"""Independent analytic reference: uniform multilayer slab at oblique incidence (exp(+jwt)).

Characteristic-matrix (Macleod) method with N = n - j k.  Used only as a test oracle for
``rcwa.py``; it shares no code with it.
"""

from __future__ import annotations

import cmath
import math


def stack_transmission(layers, freq_hz: float, theta_deg: float, pol: str):
    """layers: iterable of (thickness_m, eps_r, tan_delta) from the incident side, in air | ... | air.

    Returns (r, t): E-field coefficients referenced to the first and last interface, with the
    tangential-E convention (r_p equals r_s at normal incidence).
    """
    c0 = 299_792_458.0
    lam = c0 / freq_hz
    s0 = math.sin(math.radians(theta_deg))
    cos0 = math.cos(math.radians(theta_deg))

    def adm(n_c, cos_t):
        return n_c * cos_t if pol == "TE" else n_c / cos_t

    y0 = adm(1.0, cos0)
    bb, cc = 1.0 + 0j, y0 + 0j  # placeholder, overwritten below
    m = [[1.0 + 0j, 0j], [0j, 1.0 + 0j]]
    for d, eps, tand in layers:
        n_c = cmath.sqrt(eps * (1 - 1j * tand))  # exp(+jwt): Im(eps) < 0
        cos_t = cmath.sqrt(1 - (s0 / n_c) ** 2)
        if cos_t.real < 0:
            cos_t = -cos_t
        delta = 2 * math.pi * n_c * d * cos_t / lam
        y = adm(n_c, cos_t)
        mi = [
            [cmath.cos(delta), 1j * cmath.sin(delta) / y],
            [1j * y * cmath.sin(delta), cmath.cos(delta)],
        ]
        m = [
            [m[0][0] * mi[0][0] + m[0][1] * mi[1][0], m[0][0] * mi[0][1] + m[0][1] * mi[1][1]],
            [m[1][0] * mi[0][0] + m[1][1] * mi[1][0], m[1][0] * mi[0][1] + m[1][1] * mi[1][1]],
        ]
    ys = y0  # air substrate
    bb = m[0][0] + m[0][1] * ys
    cc = m[1][0] + m[1][1] * ys
    r = (y0 * bb - cc) / (y0 * bb + cc)
    t = 2 * y0 / (y0 * bb + cc)
    return r, t
