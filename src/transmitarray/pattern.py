"""Far field of the exit aperture of the lens from per-cell complex E (Ex, Ey) and cell geometry.

Model (documented approximation, not a full-wave result): each cell radiates a square patch of tangential
field (Ex, Ey) = A_c exp(-j k (r_c + u_inc . (r - r_c))) J_c [1, 0]^T -- the incident wave keeps its own tilt
inside the cell, only the cell's transmission J_c is piecewise constant; Huygens obliquity (1 + cos t)/2; Ludwig-3 co/cross for an
x-polarised aperture, co = K Fx, cross = K Fy with  F(u) = sum_c E_c p^2 sinc(k(ux-uix)p/2) sinc(k(uy-uiy)p/2) exp(+j k u . r_c)
(exp(+jwt), u = (sin t cos p, sin t sin p) toward the beam, t from the boresight / exit axis).
Normalised so the horn's total radiated power is 1: U = G / (4 pi) = realised gain / 4 pi.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .feed import field_amplitude


def exit_field(layout, q, x, y, theta_f_deg, r_m, k0, jones_list):
    """(Ex, Ey, tau) per cell: E0 normalised so the horn radiates unit power; tau = cell power transmission."""
    from .cells import basis_matrix

    amp = field_amplitude(theta_f_deg, r_m, q) * math.sqrt((2 * q + 1) / (2 * math.pi))
    ex = np.empty(len(x), complex)
    ey = np.empty(len(x), complex)
    tau = np.empty(len(x))
    for i, j in enumerate(jones_list):
        inc = amp[i] * np.exp(-1j * k0 * r_m[i])
        out = j @ np.array([1.0, 0.0])
        ex[i], ey[i] = inc * out
        m_inv = np.linalg.inv(basis_matrix(theta_f_deg[i], math.degrees(math.atan2(y[i], x[i]))))
        a = m_inv @ np.array([1.0, 0.0])
        b = m_inv @ out
        tau[i] = float(np.sum(np.abs(b) ** 2) / np.sum(np.abs(a) ** 2))
    return ex, ey, tau


def incident_power_per_cell(theta_f_deg, r_m, q, pitch_m):
    """Power through each cell footprint for unit horn power: (2q+1)/(2 pi) * A^2 cos(t_f) p^2."""
    amp2 = field_amplitude(theta_f_deg, r_m, q) ** 2
    return (2 * q + 1) / (2 * math.pi) * amp2 * np.cos(np.radians(theta_f_deg)) * pitch_m**2


@dataclass
class FarField:
    theta_deg: np.ndarray
    phi_deg: np.ndarray
    u_co: np.ndarray  # U (= gain / 4 pi), shape (n_theta, n_phi)
    u_cross: np.ndarray

    def gain_dbi(self, which="co"):
        u = self.u_co if which == "co" else self.u_cross
        return 10 * np.log10(np.maximum(4 * math.pi * u, 1e-30))

    def radiated_power(self):
        th = np.radians(self.theta_deg)
        dth, dph = th[1] - th[0], math.radians(self.phi_deg[1] - self.phi_deg[0])
        w = np.sin(th)[:, None]
        return float(np.sum((self.u_co + self.u_cross) * w) * dth * dph)


def far_field(x, y, ex, ey, k0, pitch_m, theta_deg, phi_deg, chunk=2000, u_inc=None) -> FarField:
    th = np.radians(theta_deg)[:, None]
    ph = np.radians(phi_deg)[None, :]
    ux = (np.sin(th) * np.cos(ph)).ravel()
    uy = (np.sin(th) * np.sin(ph)).ravel()
    fx = np.zeros(ux.size, complex)
    fy = np.zeros(ux.size, complex)
    uix, uiy = (np.zeros(len(x)), np.zeros(len(x))) if u_inc is None else u_inc
    for s in range(0, ux.size, chunk):
        sl = slice(s, s + chunk)
        ph_mat = np.exp(1j * k0 * (np.outer(ux[sl], x) + np.outer(uy[sl], y)))
        # square-patch integral of exp(+j k (u - u_inc) . (r - r_c)); u_inc is the cell's own incident direction
        ph_mat *= (
            np.sinc(k0 * (ux[sl, None] - uix[None, :]) * pitch_m / 2 / math.pi)
            * np.sinc(k0 * (uy[sl, None] - uiy[None, :]) * pitch_m / 2 / math.pi)
            * pitch_m**2
        )
        fx[sl] = ph_mat @ ex
        fy[sl] = ph_mat @ ey
    kfac = (k0 / (2 * math.pi)) ** 2 * (((1 + np.cos(th)) / 2) ** 2).repeat(len(phi_deg), axis=1).ravel()
    # eta cancels against E0^2 = 2 eta (2q+1)/(2 pi), already folded into (ex, ey) by exit_field
    shape = (len(theta_deg), len(phi_deg))
    return FarField(
        np.asarray(theta_deg),
        np.asarray(phi_deg),
        (kfac * np.abs(fx) ** 2).reshape(shape),
        (kfac * np.abs(fy) ** 2).reshape(shape),
    )


def beam_metrics(ff: FarField, beam_deg=30.0, search_deg=6.0, cone_deg=6.0):
    """Peaks near (+-beam, phi = 0/180), pointing, level balance, SLL and cross-pol."""
    th, ph = ff.theta_deg, ff.phi_deg
    g_co = ff.gain_dbi("co")
    out = {}
    for name, phi0 in (("beam_phi0", 0.0), ("beam_phi180", 180.0)):
        jphi = np.argmin(np.abs(((ph - phi0 + 180) % 360) - 180))
        sel = np.abs(th - beam_deg) <= search_deg
        # allow azimuth wander: search a small phi window
        jmask = np.abs(((ph - phi0 + 180) % 360) - 180) <= 6.0
        sub = g_co[np.ix_(sel, jmask)]
        i, j = np.unravel_index(np.argmax(sub), sub.shape)
        ti = np.nonzero(sel)[0][i]
        pj = np.nonzero(jmask)[0][j]
        out[name] = dict(
            peak_gain_dbi=float(sub[i, j]),
            theta_deg=float(th[ti]),
            phi_deg=float(ph[pj]),
            cross_at_peak_dbi=float(ff.gain_dbi("cross")[ti, pj]),
        )
        del jphi
    out["level_difference_db"] = out["beam_phi0"]["peak_gain_dbi"] - out["beam_phi180"]["peak_gain_dbi"]
    # SLL: max co-pol outside cones around both beams
    thg, phg = np.meshgrid(np.radians(th), np.radians(ph), indexing="ij")
    ux, uy = np.sin(thg) * np.cos(phg), np.sin(thg) * np.sin(phg)
    mask = np.ones(g_co.shape, bool)
    for s in (+1, -1):
        b = math.sin(math.radians(beam_deg))
        ang = np.degrees(
            np.arccos(np.clip((ux * s * b + uy * 0.0 + np.cos(thg) * math.cos(math.radians(beam_deg))), -1, 1))
        )
        mask &= ang > cone_deg
    lower_peak = min(out["beam_phi0"]["peak_gain_dbi"], out["beam_phi180"]["peak_gain_dbi"])
    out["sll_db"] = float(np.max(g_co[mask]) - lower_peak)
    out["max_cross_dbi"] = float(np.max(ff.gain_dbi("cross")))
    return out
