"""Circular-polarisation view of the TE/TM reflection matrix (Phase 20 R3). Convention: docs/theory/CONVENTIONS.md, "Circular basis".

`ObliqueSParameters.s11` is the 2x2 power-wave reflection in the (TE, TM) basis, TM = tangential E along p, TE along q = z x p, for the
incident AND the reflected wave. With that convention the (TM, TE, k) triad is left-handed for the incident wave (k = -z at normal incidence)
and right-handed for the reflected wave, which is why a plain metal plate flips handedness. Handedness is defined relative to each wave's OWN
propagation direction (IEEE, exp(+j w t)): at normal incidence

    incident  (k = -z):  RHCP = (x + j y)/sqrt2 = (TM + j TE)/sqrt2,   LHCP = (TM - j TE)/sqrt2
    reflected (k = +z):  RHCP = (x - j y)/sqrt2 = (TM - j TE)/sqrt2,   LHCP = (TM + j TE)/sqrt2

and the same relations hold in power-wave amplitudes at oblique incidence (the TE/TM amplitudes carry the same sqrt(cos theta) factor).
`circular_reflection(s11)[a, b]` is the reflected `a`-handed amplitude per unit incident `b`-handed amplitude, a, b in (R, L).
"""

import numpy as np

R, L = 0, 1  # index of RHCP, LHCP in the circular matrix
_INC = {R: np.array([1.0, 1.0j]) / np.sqrt(2.0), L: np.array([1.0, -1.0j]) / np.sqrt(2.0)}  # (TM, TE) components
_REF = {R: np.array([1.0, -1.0j]) / np.sqrt(2.0), L: np.array([1.0, 1.0j]) / np.sqrt(2.0)}


def circular_reflection(s11):
    """2x2 complex matrix `[[r_RR, r_RL], [r_LR, r_LL]]` (reflected row, incident column) from `s11` in the (TE, TM) basis."""
    s = np.asarray(s11, dtype=complex)
    if s.shape != (2, 2):
        raise ValueError("s11 must be a 2x2 reflection matrix in the (TE, TM) basis")
    s_tm_te = s[::-1, ::-1]  # reorder to (TM, TE)
    out = np.zeros((2, 2), dtype=complex)
    for a in (R, L):
        for b in (R, L):
            out[a, b] = np.vdot(_REF[a], s_tm_te @ _INC[b])
    return out


def cartesian_reflection(s11, phi=0.0):
    """Reflection in the fixed Cartesian (x, y) frame at NORMAL incidence from `s11` (TE, TM) for the plane of incidence at azimuth `phi`:
    TM = p = (cos phi, sin phi), TE = q = (-sin phi, cos phi). Rows: reflected (x, y); columns: incident (x, y)."""
    s = np.asarray(s11, dtype=complex)[::-1, ::-1]  # (TM, TE)
    b = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])  # columns p, q
    return b @ s @ b.T


def jones_to_te_tm(r_cart, phi=0.0):
    """Inverse of `cartesian_reflection`: an (x, y) reflection matrix as the (TE, TM) s11 for azimuth `phi` (normal incidence)."""
    b = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])
    s_tm_te = b.T @ np.asarray(r_cart, dtype=complex) @ b
    return s_tm_te[::-1, ::-1]
