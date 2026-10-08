"""Minimal RCWA (Fourier modal method) for  air | z-invariant layers | air  on a rectangular lattice.

Internals use the physics-text convention exp(-j w t) with fields ~ exp(+j(kx x + ky y + kz z));
:func:`solve_stack` conjugates the result so everything it returns follows the project convention
exp(+j w t) (``docs/theory/CONVENTIONS.md``; lossy dielectric ``eps_r (1 - j tan_delta)``).

Formulation (derived in ``docs/theory/DIELECTRIC_CELL_RCWA.md``):
    d/dz' [Ex;Ey] = A [Hx;Hy],   d/dz' [Hx;Hy] = B [Ex;Ey],   z' = k0 z,  H~ = -j eta0 H
with the eigenproblem  A B W = -gamma^2 W  per layer, Redheffer star-product assembly
(unconditionally stable, Rumpf 2011 layer S-matrices embedded in the air reference medium).

Field factorisation for a patterned layer:
    "laurent"  plain Laurent rule (Fourier convolution matrix of eps) -- slow (1/N) for TM-like fields;
    "normal"   normal-vector method (Schuster et al. 2007):  E_eff = [eps] P_t + [1/eps]^-1 P_n,
               with P_n = n n^T from the signed-distance normal field of the pillar boundary.
eps_zz always uses the Laurent rule (E_z is tangential to every pillar wall, hence continuous).

Reference planes (documented, see ``cell.py``): S11 at the first interface (top of the first layer), S21
at the last interface (bottom of the last layer); transmitted phase is *not* de-embedded to any air plane.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from . import fourier as fr


@dataclass(frozen=True)
class UniformLayer:
    thickness_m: float
    eps_r: float
    tan_delta: float


@dataclass(frozen=True)
class PolygonLayer:
    """Convex polygon of ``eps_r_in`` (loss ``tan_delta_in``) in a background (default air), extruded.

    A hole is ``eps_r_in = 1`` in a nylon background (``eps_r_out``, ``tan_delta_out``).
    """

    thickness_m: float
    vertices_m: tuple  # ((x, y), ...) CCW, metres
    eps_r_in: float
    tan_delta_in: float
    eps_r_out: float = 1.0
    tan_delta_out: float = 0.0


@dataclass(frozen=True)
class Region:
    """One painted polygon: value ``eps_r_in`` replaces the local value ``eps_r_prev`` (for 1/eps and eps).

    ``edges_m`` are the material-interface segments used for the normal field (default: every polygon
    edge); pass only the true interfaces when the polygon edge lies on a periodic cell border.
    """

    vertices_m: tuple
    eps_r_in: float
    tan_delta_in: float
    eps_r_prev: float
    tan_delta_prev: float
    edges_m: tuple | None = None


@dataclass(frozen=True)
class RegionLayer:
    """Layer = background + painted regions (nested/disjoint polygons); supercells and plate/hole layouts."""

    thickness_m: float
    eps_r_bg: float
    tan_delta_bg: float
    regions: tuple


def _branch(gamma2_neg: np.ndarray) -> np.ndarray:
    """gamma = sqrt(-lambda) on the decaying/forward branch for exp(-jwt) with exp(+j gamma z)."""
    g = np.sqrt(gamma2_neg.astype(complex))
    flip = g.imag < -1e-12 * np.maximum(np.abs(g), 1e-300)
    return np.where(flip, -g, g)


def _layer_matrices(layer, ctx):
    """Returns (E_eff 2Nx2N, Te_inv NxN) for the layer (exp(-jwt): eps = eps_r (1 + j tan_delta))."""
    n = ctx["N"]
    if isinstance(layer, UniformLayer):
        eps = layer.eps_r * (1 + 1j * layer.tan_delta)
        return eps * np.eye(2 * n), np.eye(n) / eps
    single = isinstance(layer, PolygonLayer)
    if single:
        verts = np.asarray(layer.vertices_m, float)
        reg = [Region(layer.vertices_m, layer.eps_r_in, layer.tan_delta_in, layer.eps_r_out, layer.tan_delta_out)]
        bg = layer.eps_r_out * (1 + 1j * layer.tan_delta_out)
        edges = fr.polygon_edges(verts)
    else:
        reg = list(layer.regions)
        bg = layer.eps_r_bg * (1 + 1j * layer.tan_delta_bg)
        edges = np.concatenate(
            [
                np.asarray(r.edges_m, float).reshape(-1, 2, 2)
                if r.edges_m is not None
                else fr.polygon_edges(np.asarray(r.vertices_m, float))
                for r in reg
            ]
        )
    M = ctx["n_harm"]
    cplx = lambda e, t: e * (1 + 1j * t)  # noqa: E731
    eps_regions, inv_regions = [], []
    for r in reg:
        v = np.asarray(r.vertices_m, float)
        ein, eprev = cplx(r.eps_r_in, r.tan_delta_in), cplx(r.eps_r_prev, r.tan_delta_prev)
        eps_regions.append((v, ein - eprev))
        inv_regions.append((v, 1.0 / ein - 1.0 / eprev))
    c_eps = fr.region_pattern_coeffs(eps_regions, ctx["px"], ctx["py"], M, bg)
    te = fr.toeplitz_from_coeffs(c_eps, M, ctx["m"], ctx["n"])
    te_inv = np.linalg.inv(te)
    if ctx["factorization"] == "laurent":
        z = np.zeros_like(te)
        return np.block([[te, z], [z, te]]), te_inv
    c_inv = fr.region_pattern_coeffs(inv_regions, ctx["px"], ctx["py"], M, 1.0 / bg)
    ti = fr.toeplitz_from_coeffs(c_inv, M, ctx["m"], ctx["n"])
    ti_inv = np.linalg.inv(ti)
    cxx, cxy, cyy = fr.normal_field_coeffs(edges, ctx["px"], ctx["py"], M, ctx["n_grid"], xy_symmetric=single)
    pxx = fr.toeplitz_from_coeffs(cxx, M, ctx["m"], ctx["n"])
    pxy = fr.toeplitz_from_coeffs(cxy, M, ctx["m"], ctx["n"])
    pyy = fr.toeplitz_from_coeffs(cyy, M, ctx["m"], ctx["n"])
    pn = np.block([[pxx, pxy], [pxy, pyy]])
    pt = np.eye(2 * n) - pn
    z = np.zeros_like(te)
    e_lau = np.block([[te, z], [z, te]])
    e_inv = np.block([[ti_inv, z], [z, ti_inv]])
    return e_lau @ pt + e_inv @ pn, te_inv


def _AB(e_eff, te_inv, kx, ky):
    n = kx.size
    KX, KY = np.diag(kx), np.diag(ky)
    exx, exy = e_eff[:n, :n], e_eff[:n, n:]
    eyx, eyy = e_eff[n:, :n], e_eff[n:, n:]
    eye = np.eye(n)
    a = np.block(
        [
            [-KX @ te_inv @ KY, KX @ te_inv @ KX - eye],
            [eye - KY @ te_inv @ KY, KY @ te_inv @ KX],
        ]
    )
    b = np.block(
        [
            [-KX @ KY - eyx, KX @ KX - eyy],
            [exx - KY @ KY, exy + KY @ KX],
        ]
    )
    return a, b


def _star(a, b):
    """Redheffer product a (upper) * b (lower); blocks keyed 11, 12, 21, 22."""
    eye = np.eye(a["11"].shape[0])
    t1 = np.linalg.solve(eye - b["11"] @ a["22"], np.hstack([b["11"] @ a["21"], b["12"]]))
    n = a["11"].shape[0]
    u, v = t1[:, :n], t1[:, n:]
    s11 = a["11"] + a["12"] @ u
    s12 = a["12"] @ v
    t2 = np.linalg.solve(eye - a["22"] @ b["11"], np.hstack([a["21"], a["22"] @ b["12"]]))
    s21 = b["21"] @ t2[:, :n]
    s22 = b["22"] + b["21"] @ t2[:, n:]
    return {"11": s11, "12": s12, "21": s21, "22": s22}


def solve_stack(
    layers,
    period_x_m: float,
    period_y_m: float,
    freq_hz: float,
    theta_deg: float,
    phi_deg: float,
    n_harm: int,
    factorization: str = "normal",
    n_grid: int | tuple | None = None,
    orders: tuple = (),
):
    """Total 4-block S-matrix of the stack in the (0,0)-order TE/TM basis, project convention exp(+jwt).

    Returns ``dict`` with ``r`` and ``t`` 2x2 complex matrices (rows = outgoing [s, p], columns =
    incoming [s, p]), plus ``power_other_orders_ok`` bookkeeping fields.  s = E perpendicular to the
    plane of incidence; p = in the plane.  Transmission uses the same unit vectors in and out;
    reflection uses the tangential-E convention (r_pp = r_ss at normal incidence).
    """
    if factorization not in ("laurent", "normal"):
        raise ValueError(f"factorization must be 'laurent' or 'normal', got {factorization!r}")
    c0 = 299_792_458.0
    lam0 = c0 / freq_hz
    th, ph = math.radians(theta_deg), math.radians(phi_deg)
    m, n = fr.harmonic_indices(n_harm)
    N = m.size
    kx = math.sin(th) * math.cos(ph) + m * lam0 / period_x_m
    ky = math.sin(th) * math.sin(ph) + n * lam0 / period_y_m
    mx_, my_ = fr.mxy(n_harm)
    if n_grid is None:
        n_grid = (max(128, 8 * (4 * mx_ + 1)), max(128, 8 * (4 * my_ + 1)))
        if mx_ == my_:
            n_grid = n_grid[0]
    ctx = dict(N=N, n_harm=n_harm, px=period_x_m, py=period_y_m, m=m, n=n, factorization=factorization, n_grid=n_grid)

    # air reference medium: W0 = I, gamma0 = sqrt(1 - kx^2 - ky^2)
    g0 = _branch(1 - kx**2 - ky**2)
    if np.any(np.abs(g0) < 1e-9):
        raise ValueError("a diffraction order is at grazing (Rayleigh anomaly); shift frequency/angle")
    KX, KY = np.diag(kx), np.diag(ky)
    eye = np.eye(N)
    b0 = np.block([[-KX @ KY, KX @ KX - eye], [eye - KY @ KY, KY @ KX]])
    g0_2 = np.concatenate([g0, g0])
    w0 = np.eye(2 * N)
    v0 = 1j * (b0 / g0_2[None, :])  # V_fwd = j B W Gamma^-1 with W = I

    total = None
    for layer in layers:
        e_eff, te_inv = _layer_matrices(layer, ctx)
        a, b = _AB(e_eff, te_inv, kx, ky)
        lam, w = np.linalg.eig(a @ b)
        gam = _branch(-lam)
        gam = np.where(np.abs(gam) < 1e-9, 1e-9, gam)
        v = 1j * ((b @ w) / gam[None, :])
        x = np.exp(1j * gam * (2 * math.pi / lam0) * layer.thickness_m)
        wi_w0 = np.linalg.solve(w, w0)
        vi_v0 = np.linalg.solve(v, v0)
        am, bm = wi_w0 + vi_v0, wi_w0 - vi_v0
        xb = x[:, None] * bm
        den = am - xb @ np.linalg.solve(am, xb)
        # Interface matching (derivation in docs/theory/DIELECTRIC_CELL_RCWA.md):
        #   S11 = den^-1 (X B A^-1 X A - B),   S21 = A^-1 X (A + B S11)
        s11 = np.linalg.solve(den, xb @ np.linalg.solve(am, x[:, None] * am) - bm)
        s12 = np.linalg.solve(am, x[:, None] * (am + bm @ s11))
        sl = {"11": s11, "12": s12, "21": s12, "22": s11}
        total = sl if total is None else _star(total, sl)

    # incident (0,0) order: transverse E = ex*xhat + ey*yhat
    i0 = int(np.flatnonzero((m == 0) & (n == 0))[0])
    ct, st = math.cos(th), math.sin(th)
    cphi, sphi = math.cos(ph), math.sin(ph)
    inc = {"s": (-sphi, cphi), "p": (ct * cphi, ct * sphi)}
    r_mat = np.zeros((2, 2), complex)
    t_mat = np.zeros((2, 2), complex)
    for j, pol in enumerate(("s", "p")):
        c_in = np.zeros(2 * N, complex)
        c_in[i0], c_in[N + i0] = inc[pol]
        for blk, mat in (("11", r_mat), ("21", t_mat)):
            out = total[blk] @ c_in
            ex, ey = out[i0], out[N + i0]
            mat[0, j] = -sphi * ex + cphi * ey
            mat[1, j] = (cphi * ex + sphi * ey) / ct
    # transmitted transverse E (ex, ey) of requested diffraction orders, per incident pol [s, p]
    orders_out = {}
    for mo, no in orders:
        idx = np.flatnonzero((m == mo) & (n == no))
        if idx.size == 0:
            raise ValueError(f"order {(mo, no)} outside the harmonic truncation")
        io = int(idx[0])
        res = np.zeros((2, 2), complex)  # [x/y, s/p]
        for j, pol in enumerate(("s", "p")):
            c_in = np.zeros(2 * N, complex)
            c_in[i0], c_in[N + i0] = inc[pol]
            out = total["21"] @ c_in
            res[:, j] = out[io], out[N + io]
        orders_out[(mo, no)] = np.conj(res)
    # exp(-jwt) -> exp(+jwt)
    return {
        "r": np.conj(r_mat),
        "t": np.conj(t_mat),
        "n_modes": N,
        "kx0": kx[i0],
        "ky0": ky[i0],
        "sin_theta": st,
        "orders": orders_out,
    }
