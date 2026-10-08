"""Fourier coefficients of a z-invariant (extruded) dielectric pattern.

Pattern: one convex polygon (the octagonal pillar) of permittivity ``eps_in``
in a background ``eps_out``, on a rectangular lattice (px, py), origin at the
cell centre.  Expansion convention (matches ``rcwa.py``)::

    f(r) = sum_G f_G exp(+j G.r),   f_G = (1/A) * integral_cell f exp(-j G.r) dA

The convention only fixes the sign of G; for the centro-symmetric patterns used
here every coefficient is real and the sign is immaterial.
"""

from __future__ import annotations

import math

import numpy as np


def harmonic_indices(n_harm: int) -> tuple[np.ndarray, np.ndarray]:
    """Square truncation m, n in [-n_harm, n_harm]; returns flattened (m, n), m outer."""
    if int(n_harm) != n_harm or n_harm < 0:
        raise ValueError(f"n_harm must be a non-negative integer, got {n_harm!r}")
    r = np.arange(-n_harm, n_harm + 1)
    mm, nn = np.meshgrid(r, r, indexing="ij")
    return mm.ravel(), nn.ravel()


def octagon_vertices(circumradius: float, rotation_deg: float) -> np.ndarray:
    """Regular octagon, CCW, vertex k at angle rotation_deg + 45*k degrees.

    ``rotation_deg = 22.5`` puts the flats on the x/y axes (extent R*cos22.5 in x, y);
    ``rotation_deg = 0`` puts vertices on the axes (extent R).
    """
    ang = np.deg2rad(rotation_deg + 45.0 * np.arange(8))
    return circumradius * np.stack([np.cos(ang), np.sin(ang)], axis=1)


def polygon_ft(vertices: np.ndarray, gx: np.ndarray, gy: np.ndarray) -> np.ndarray:
    """Unnormalised FT  integral_polygon exp(-j G.r) dA  for a CCW polygon (any G, incl. 0).

    Divergence theorem with a = G/|G|^2:  integral = j * oint (G.n/|G|^2) exp(-j G.r) ds, and
    the line integral of exp(-j G.r) along an edge is L * exp(-j G.r_mid) * sinc(G.d/2).
    """
    gx = np.asarray(gx, float)
    gy = np.asarray(gy, float)
    g2 = gx**2 + gy**2
    out = np.zeros(gx.shape, complex)
    nz = g2 > 0
    area = 0.5 * float(
        np.sum(vertices[:, 0] * np.roll(vertices[:, 1], -1) - np.roll(vertices[:, 0], -1) * vertices[:, 1])
    )
    out[~nz] = area
    gxn, gyn, g2n = gx[nz], gy[nz], g2[nz]
    acc = np.zeros(gxn.shape, complex)
    for i in range(len(vertices)):
        r1, r2 = vertices[i], vertices[(i + 1) % len(vertices)]
        d = r2 - r1
        length = math.hypot(d[0], d[1])
        nx, ny = d[1] / length, -d[0] / length  # outward normal of a CCW polygon
        rm = 0.5 * (r1 + r2)
        gd = gxn * d[0] + gyn * d[1]
        acc += (gxn * nx + gyn * ny) * length * np.sinc(gd / (2 * np.pi)) * np.exp(-1j * (gxn * rm[0] + gyn * rm[1]))
    out[nz] = 1j * acc / g2n
    return out


def toeplitz_from_coeffs(coeff: np.ndarray, n_harm: int, m: np.ndarray, n: np.ndarray) -> np.ndarray:
    """N x N convolution matrix T[p, q] = coeff[m_p - m_q, n_p - n_q].

    ``coeff`` has shape (4 n_harm + 1, 4 n_harm + 1), index (dm + 2 n_harm, dn + 2 n_harm).
    """
    o = 2 * n_harm
    return coeff[(m[:, None] - m[None, :]) + o, (n[:, None] - n[None, :]) + o]


def difference_grid(n_harm: int) -> tuple[np.ndarray, np.ndarray]:
    """(dm, dn) arrays of shape (4M+1, 4M+1) covering all harmonic differences."""
    d = np.arange(-2 * n_harm, 2 * n_harm + 1)
    return np.meshgrid(d, d, indexing="ij")


def polygon_pattern_coeffs(
    vertices: np.ndarray, px: float, py: float, n_harm: int, val_in: complex, val_out: complex
) -> np.ndarray:
    """Coefficients of  f = val_out + (val_in - val_out) * chi_polygon  on the difference grid."""
    dm, dn = difference_grid(n_harm)
    gx = 2 * np.pi * dm / px
    gy = 2 * np.pi * dn / py
    chi = polygon_ft(vertices, gx, gy) / (px * py)
    c = (val_in - val_out) * chi
    c[2 * n_harm, 2 * n_harm] += val_out
    return c


def _segment_nearest(px_: np.ndarray, py_: np.ndarray, a: np.ndarray, b: np.ndarray):
    ab = b - a
    t = ((px_ - a[0]) * ab[0] + (py_ - a[1]) * ab[1]) / float(ab @ ab)
    t = np.clip(t, 0.0, 1.0)
    cx, cy = a[0] + t * ab[0], a[1] + t * ab[1]
    return cx, cy, (px_ - cx) ** 2 + (py_ - cy) ** 2


def normal_field_coeffs(vertices: np.ndarray, px: float, py: float, n_harm: int, n_grid: int):
    """Fourier coefficients of n_a n_b (xx, xy, yy) for the smooth-ish normal field n = grad(distance).

    n = (r - nearest boundary point)/|.| over periodic images.  The sign of n is irrelevant (P = n n^T).
    The coefficients are obtained by an exact matrix DFT on an ``n_grid``^2 pixel-centre grid.
    """
    u = (np.arange(n_grid) + 0.5) / n_grid - 0.5
    X, Y = np.meshgrid(u * px, u * py, indexing="ij")
    best = np.full(X.shape, np.inf)
    nxv = np.zeros(X.shape)
    nyv = np.zeros(X.shape)
    for ix in (-1, 0, 1):
        for iy in (-1, 0, 1):
            Xs, Ys = X - ix * px, Y - iy * py
            for i in range(len(vertices)):
                a, b = vertices[i], vertices[(i + 1) % len(vertices)]
                cx, cy, d2 = _segment_nearest(Xs, Ys, a, b)
                upd = d2 < best
                best = np.where(upd, d2, best)
                nxv = np.where(upd, Xs - cx, nxv)
                nyv = np.where(upd, Ys - cy, nyv)
    norm = np.sqrt(nxv**2 + nyv**2)
    # points that sit exactly on the boundary: fall back to the outward normal of the nearest edge
    tiny = norm < 1e-12 * max(px, py)
    if np.any(tiny):
        nxv = np.where(tiny, 1.0, nxv)
        nyv = np.where(tiny, 0.0, nyv)
        norm = np.where(tiny, 1.0, norm)
    nxv, nyv = nxv / norm, nyv / norm
    d = np.arange(-2 * n_harm, 2 * n_harm + 1)
    ex = np.exp(-2j * np.pi * np.outer(d, u))  # (nd, n_grid), phase of x_j/px = u_j
    fxx, fxy, fyy = nxv * nxv, nxv * nyv, nyv * nyv
    if abs(px - py) < 1e-12 * px:
        # The pillar is symmetric under x<->y; exact ties between periodic images on the x=y line are
        # otherwise broken arbitrarily.  Averaging with the transposed grid restores the symmetry.
        fxx, fyy = 0.5 * (fxx + fyy.T), 0.5 * (fyy + fxx.T)
        fxy = 0.5 * (fxy + fxy.T)
    out = []
    for f in (fxx, fxy, fyy):
        out.append((ex @ f @ ex.T) / n_grid**2)
    return out
