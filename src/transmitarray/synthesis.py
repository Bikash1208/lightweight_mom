"""Dual-beam aperture-phase synthesis (weighted Gerchberg-Saxton) for a fixed aperture amplitude.

Far field of the exit aperture field a_c = A_c exp(j psi_c):   F(u) = sum_c a_c exp(+j k u.r_c)
(exp(+jwt); beam toward unit-vector projection u = (sin t cos p, sin t sin p)).
"""

from __future__ import annotations

import math

import numpy as np


class ApertureFFT:
    """Zero-padded FFT far field on the lens grid:  u = m * lambda / (n_fft * pitch)."""

    def __init__(self, ix, iy, n_grid, wavelength_m, pitch_m, n_fft=256):
        self.ix, self.iy, self.n_fft = np.asarray(ix), np.asarray(iy), n_fft
        self.n_grid, self.lam, self.p = n_grid, wavelength_m, pitch_m
        self.du = wavelength_m / (n_fft * pitch_m)

    def forward(self, a_cells):
        g = np.zeros((self.n_fft, self.n_fft), complex)
        g[self.ix, self.iy] = a_cells
        return np.fft.ifft2(g) * self.n_fft**2  # sum a exp(+j 2 pi (m i + n j)/N)

    def backward(self, far):
        g = np.fft.fft2(far) / self.n_fft**2
        return g[self.ix, self.iy]

    def u_axes(self):
        m = np.fft.fftfreq(self.n_fft, d=1.0 / self.n_fft)  # signed bins in ifft2 ordering
        return m * self.du

    def window_mask(self, u0, v0, radius_u):
        ux = self.u_axes()
        UX, UY = np.meshgrid(ux, ux, indexing="ij")
        return (UX - u0) ** 2 + (UY - v0) ** 2 <= radius_u**2


def beam_windows(fft: ApertureFFT, theta_deg: float, radius_u: float):
    s = math.sin(math.radians(theta_deg))
    return [fft.window_mask(s, 0.0, radius_u), fft.window_mask(-s, 0.0, radius_u)]


def window_efficiency(far, masks):
    """Fraction of the far-field energy inside each window, and in total."""
    tot = float(np.sum(np.abs(far) ** 2))
    per = [float(np.sum(np.abs(far[m]) ** 2)) / tot for m in masks]
    return per, sum(per)


def binary_baseline_phase(x_m, k0, theta_deg):
    """arg(exp(+j k sin t x) + exp(-j k sin t x)) = 0 / pi."""
    return np.where(np.cos(k0 * math.sin(math.radians(theta_deg)) * x_m) >= 0, 0.0, math.pi)


def _group_mean(values, groups):
    """Complex mean of ``values`` over each group of equal integer ``groups`` (mirror-image cells), broadcast back."""
    n = groups.max() + 1
    cnt = np.maximum(np.bincount(groups, minlength=n), 1)
    re = np.bincount(groups, weights=values.real, minlength=n) / cnt
    im = np.bincount(groups, weights=values.imag, minlength=n) / cnt
    return (re + 1j * im)[groups]


def gerchberg_saxton(
    amp,
    x_m,
    fft: ApertureFFT,
    theta_deg,
    radius_u,
    n_iter=400,
    k0=None,
    seed_phase=None,
    out_floor=0.0,
    shaped=False,
    u_center=None,
    mirror_groups=None,
):
    """Maximise equal power into two windows around +-theta_deg (x-z plane).  Returns phase [rad] per cell.

    Projection: far-field magnitude is zeroed (times ``out_floor``) outside the windows, and the two
    windows are rescaled to carry equal energy; the aperture amplitude is restored to ``amp`` each step.
    ``mirror_groups`` (integer id per cell) forces the aperture field to be identical over each group of
    mirror-image cells every iteration, so the synthesised phase keeps the lens's D2 symmetry.
    ``shaped=True`` additionally pins each window's magnitude to the *profile of the ideal single beam* aimed at
    its centre (aperture phase -k u_b x), so the beam cannot wander inside the window (pointing is fixed).
    """
    if u_center is None:
        masks = beam_windows(fft, theta_deg, radius_u)
    else:
        masks = [fft.window_mask(u_center, 0.0, radius_u), fft.window_mask(-u_center, 0.0, radius_u)]
    union = masks[0] | masks[1]
    psi = binary_baseline_phase(x_m, k0, theta_deg) if seed_phase is None else np.asarray(seed_phase, float)
    profiles = None
    if shaped:
        sgn = (+1.0, -1.0)
        profiles = [
            np.abs(fft.forward(amp * np.exp(-1j * k0 * sg * math.sin(math.radians(theta_deg)) * x_m))) * m
            for sg, m in zip(sgn, masks, strict=True)
        ]
    best = (-1.0, psi.copy())
    for _ in range(n_iter):
        far = fft.forward(amp * np.exp(1j * psi))
        _per, tot = window_efficiency(far, masks)
        if tot > best[0]:
            best = (tot, psi.copy())
        e = [float(np.sum(np.abs(far[m]) ** 2)) for m in masks]
        target = 0.5 * (e[0] + e[1])
        new = far * out_floor
        if profiles is None:
            for m, em in zip(masks, e, strict=True):
                new[m] = far[m] * math.sqrt(target / em)
        else:
            for m, prof in zip(masks, profiles, strict=True):
                ph = far[m] / np.maximum(np.abs(far[m]), 1e-300)
                new[m] = prof[m] * ph * math.sqrt(target / float(np.sum(prof[m] ** 2)))
        new[~union] = far[~union] * out_floor
        back = fft.backward(new)
        if mirror_groups is not None:
            back = _group_mean(back, mirror_groups)
        psi = np.angle(back)
    far = fft.forward(amp * np.exp(1j * psi))
    _per, tot = window_efficiency(far, masks)
    if tot > best[0]:
        best = (tot, psi.copy())
    return best[1], best[0]


def peak_u_along_x(amp, psi, x_m, k0, u_target, half_span=0.05, n=201):
    """Fine-grid location of the |F| maximum along u_x (u_y = 0) near +u_target, and its |F|^2."""
    ug = u_target + np.linspace(-half_span, half_span, n)
    f = np.abs(np.exp(1j * k0 * np.outer(ug, x_m)) @ (amp * np.exp(1j * psi))) ** 2
    i = int(np.argmax(f))
    return float(ug[i]), float(f[i])


def gs_pinned(amp, x_m, fft, theta_deg, radius_u, n_iter, k0, n_outer=8, tol_u=0.0006):
    """GS with the window centre shifted until the beam peak sits on sin(theta_deg) (within ``tol_u`` in u)."""
    target = math.sin(math.radians(theta_deg))
    centre = target
    psi, eff = gerchberg_saxton(amp, x_m, fft, theta_deg, radius_u, n_iter, k0, u_center=centre)
    for _ in range(n_outer):
        up, _f = peak_u_along_x(amp, psi, x_m, k0, target)
        if abs(up - target) < tol_u:
            break
        centre += target - up
        psi, eff = gerchberg_saxton(amp, x_m, fft, theta_deg, radius_u, n_iter, k0, u_center=centre)
    return psi, eff, centre
