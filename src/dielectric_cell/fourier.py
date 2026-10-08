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


def mxy(n_harm) -> tuple[int, int]:
    """Truncation orders (Mx, My): an int means the same order on both axes."""
    mx, my = (n_harm, n_harm) if np.isscalar(n_harm) else tuple(n_harm)
    if int(mx) != mx or int(my) != my or mx < 0 or my < 0:
        raise ValueError(f"n_harm must be non-negative integer(s), got {n_harm!r}")
    return int(mx), int(my)


def harmonic_indices(n_harm) -> tuple[np.ndarray, np.ndarray]:
    """Rectangular truncation m in [-Mx, Mx], n in [-My, My]; returns flattened (m, n), m outer."""
    mx, my = mxy(n_harm)
    mm, nn = np.meshgrid(np.arange(-mx, mx + 1), np.arange(-my, my + 1), indexing="ij")
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


def toeplitz_from_coeffs(coeff: np.ndarray, n_harm, m: np.ndarray, n: np.ndarray) -> np.ndarray:
    """N x N convolution matrix T[p, q] = coeff[m_p - m_q, n_p - n_q].

    ``coeff`` has shape (4 Mx + 1, 4 My + 1), index (dm + 2 Mx, dn + 2 My).
    """
    mx, my = mxy(n_harm)
    return coeff[(m[:, None] - m[None, :]) + 2 * mx, (n[:, None] - n[None, :]) + 2 * my]


def difference_grid(n_harm) -> tuple[np.ndarray, np.ndarray]:
    """(dm, dn) arrays of shape (4 Mx + 1, 4 My + 1) covering all harmonic differences."""
    mx, my = mxy(n_harm)
    return np.meshgrid(np.arange(-2 * mx, 2 * mx + 1), np.arange(-2 * my, 2 * my + 1), indexing="ij")


def polygon_pattern_coeffs(
    vertices: np.ndarray, px: float, py: float, n_harm, val_in: complex, val_out: complex
) -> np.ndarray:
    """Coefficients of  f = val_out + (val_in - val_out) * chi_polygon  on the difference grid."""
    return region_pattern_coeffs([(vertices, val_in - val_out)], px, py, n_harm, val_out)


def region_pattern_coeffs(regions, px: float, py: float, n_harm, background: complex) -> np.ndarray:
    """Coefficients of  f = background + sum_i delta_i * chi_{polygon_i}  (regions: [(vertices, delta)])."""
    mx, my = mxy(n_harm)
    dm, dn = difference_grid(n_harm)
    gx = 2 * np.pi * dm / px
    gy = 2 * np.pi * dn / py
    c = np.zeros(dm.shape, complex)
    for vertices, delta in regions:
        c += delta * polygon_ft(vertices, gx, gy) / (px * py)
    c[2 * mx, 2 * my] += background
    return c


def _segment_nearest(px_: np.ndarray, py_: np.ndarray, a: np.ndarray, b: np.ndarray):
    ab = b - a
    t = ((px_ - a[0]) * ab[0] + (py_ - a[1]) * ab[1]) / float(ab @ ab)
    t = np.clip(t, 0.0, 1.0)
    cx, cy = a[0] + t * ab[0], a[1] + t * ab[1]
    return cx, cy, (px_ - cx) ** 2 + (py_ - cy) ** 2


def polygon_edges(vertices: np.ndarray) -> np.ndarray:
    """(K, 2, 2) array of the closed polygon's edges."""
    v = np.asarray(vertices, float)
    return np.stack([v, np.roll(v, -1, axis=0)], axis=1)


def normal_field_coeffs(edges: np.ndarray, px: float, py: float, n_harm, n_grid, xy_symmetric: bool = False):
    """Fourier coefficients of n_a n_b (xx, xy, yy) for the normal field n = grad(distance to ``edges``).

    ``edges`` is a (K, 2, 2) array of boundary segments (material interfaces only).  n = (r - nearest
    boundary point)/|.| over periodic images; the sign of n is irrelevant (P = n n^T).  ``n_grid`` is an int
    or an (nx, ny) pair of pixel-centre sample counts; coefficients come from an exact matrix DFT.
    ``xy_symmetric`` (single C4 cell only) symmetrises the field under x<->y.
    """
    mx, my = mxy(n_harm)
    if len(edges) == 0:  # no material interface: the normal field is irrelevant
        z = np.zeros((4 * mx + 1, 4 * my + 1), complex)
        return [z, z, z]
    ngx, ngy = (n_grid, n_grid) if np.isscalar(n_grid) else tuple(n_grid)
    ux = (np.arange(ngx) + 0.5) / ngx - 0.5
    uy = (np.arange(ngy) + 0.5) / ngy - 0.5
    X, Y = np.meshgrid(ux * px, uy * py, indexing="ij")
    best = np.full(X.shape, np.inf)
    nxv = np.zeros(X.shape)
    nyv = np.zeros(X.shape)
    for ix in (-1, 0, 1):
        for iy in (-1, 0, 1):
            Xs, Ys = X - ix * px, Y - iy * py
            for a, b in edges:
                cx, cy, d2 = _segment_nearest(Xs, Ys, a, b)
                upd = d2 < best
                best = np.where(upd, d2, best)
                nxv = np.where(upd, Xs - cx, nxv)
                nyv = np.where(upd, Ys - cy, nyv)
    norm = np.sqrt(nxv**2 + nyv**2)
    # points that sit exactly on the boundary: fall back to a fixed direction
    tiny = norm < 1e-12 * max(px, py)
    if np.any(tiny):
        nxv = np.where(tiny, 1.0, nxv)
        nyv = np.where(tiny, 0.0, nyv)
        norm = np.where(tiny, 1.0, norm)
    nxv, nyv = nxv / norm, nyv / norm
    fxx, fxy, fyy = nxv * nxv, nxv * nyv, nyv * nyv
    if xy_symmetric and ngx == ngy and abs(px - py) < 1e-12 * px:
        # The pillar is symmetric under x<->y; exact ties between periodic images on the x=y line are
        # otherwise broken arbitrarily.  Averaging with the transposed grid restores the symmetry.
        fxx, fyy = 0.5 * (fxx + fyy.T), 0.5 * (fyy + fxx.T)
        fxy = 0.5 * (fxy + fxy.T)
    dx = np.arange(-2 * mx, 2 * mx + 1)
    dy = np.arange(-2 * my, 2 * my + 1)
    ex = np.exp(-2j * np.pi * np.outer(dx, ux))
    ey = np.exp(-2j * np.pi * np.outer(dy, uy))
    return [(ex @ f @ ey.T) / (ngx * ngy) for f in (fxx, fxy, fyy)]
