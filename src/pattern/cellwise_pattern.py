"""Phase 23 B: the radiated far field of a multi-port board as a LINEAR function of every cell's
aperture contribution, so that flipping one cell's state updates the pattern in O(directions).

`vector_aperture_pattern` evaluates `E = G(k) . FT[J] . cos(theta)` for the radial current
`J = rho_hat(x, y) * s(x, y)` about a feed, with `s` the scalar aperture field (`aperture_field`, the
reference wave times the cell's achieved reflection phasor). Everything after `s` is linear, so a cell
`n` at `(x_n, y_n)` contributes, per unit `s_n`,

    c[n, d] = dx * dy * exp(j (kx_d x_n + ky_d y_n)) * (a_x(d) cos(phi'_n) + a_y(d) sin(phi'_n))

to a scalar output (the co-polar component, or the cross-polar one) in direction `d`, where
`(cos phi'_n, sin phi'_n)` is the unit vector from the feed to the cell, and `a_x(d)`, `a_y(d)` are that
output's response to a unit x- / y-directed current at `d` (the spectral dyad, the `cos(theta)`
obliquity and the Ludwig-3 projection, all taken from `vector_aperture_pattern` unmodified). The
field of port `p` is `F_p(d) = sum_n s_n c[n, d]`; with several ports excited together,
`F(d) = sum_p w_p F_p(d)` with the complex excitations `w_p = a_p exp(j alpha_p)`.

`verify_against_vector_pattern` is the check the tests run: the cell-wise field equals
`_vector_far_field` (the function `compute_vector_pattern` is built on) to ~1e-10 relative.

Directions are an explicit list of `(theta, phi)` points, normally the visible disk of the u-v plane
(`uv_direction_set`), denser around the requested beams.
"""

from dataclasses import dataclass

import numpy as np

from ..frequency import free_space_wavenumber
from .vector_aperture_pattern import (
    _vector_far_field,
    aperture_current_from_scalar,
    copol_component,
)

DEFAULT_CHUNK_BYTES = 256 * 1024 * 1024


@dataclass(frozen=True)
class DirectionSet:
    """`theta_deg`, `phi_deg`, `u`, `v`: 1-D arrays of the same length (`u = sin(theta) cos(phi)`,
    `v = sin(theta) sin(phi)`). `tag`: int array, -1 for a background point, `k >= 0` for a point of
    the local patch around requested beam `k`."""

    theta_deg: np.ndarray
    phi_deg: np.ndarray
    u: np.ndarray
    v: np.ndarray
    tag: np.ndarray

    @property
    def size(self):
        return int(self.theta_deg.size)


def _from_uv(u, v, tag):
    u, v = np.asarray(u, dtype=float), np.asarray(v, dtype=float)
    s = np.clip(np.hypot(u, v), 0.0, 1.0)
    theta = np.degrees(np.arcsin(s))
    phi = np.degrees(np.arctan2(v, u)) % 360.0
    return DirectionSet(theta, phi, u, v, np.asarray(tag, dtype=int))


def _uv_lattice(step):
    n = int(np.ceil(1.0 / step))
    g = np.arange(-n, n + 1) * step
    uu, vv = np.meshgrid(g, g, indexing="xy")
    u, v = uu.ravel(), vv.ravel()
    ok = np.hypot(u, v) <= 1.0 + 1e-12
    return u[ok], v[ok]


def _nearest_beam(u, v, beams):
    """`(angle_deg, index)`: great-circle angle to the nearest beam and which beam it is."""
    ts = np.arcsin(np.clip(np.hypot(u, v), 0.0, 1.0))
    ps = np.arctan2(v, u)
    best = np.full(u.shape, np.inf)
    idx = np.zeros(u.shape, dtype=int)
    for k, (th, ph) in enumerate(beams):
        t0, p0 = np.radians(th), np.radians(ph)
        ang = np.degrees(np.arccos(np.clip(np.cos(ts) * np.cos(t0) + np.sin(ts) * np.sin(t0) * np.cos(ps - p0), -1, 1)))
        closer = ang < best
        best = np.where(closer, ang, best)
        idx = np.where(closer, k, idx)
    return best, idx


def uv_direction_set(
    du_far=0.025, beams_deg=(), du_near=0.012, near_radius_deg=20.0, patch_radius_deg=4.0, patch_step=0.004
):
    """The optimisation directions on the visible disk of the u-v plane, three disjoint lattices by
    angular distance `d` to the nearest requested beam (`beams_deg`: `[(theta, phi), ...]`):
    `d <= patch_radius_deg`: a fine `patch_step` lattice (the main-lobe peak and its pointing are
    measured there; these points carry `tag = beam index`); `d <= near_radius_deg`: a `du_near`
    lattice; elsewhere a `du_far` lattice. Without beams, only the `du_far` lattice."""
    beams = [(float(t), float(p)) for t, p in beams_deg]
    u_f, v_f = _uv_lattice(du_far)
    if not beams:
        return _from_uv(u_f, v_f, np.full(u_f.size, -1))
    d_f, _ = _nearest_beam(u_f, v_f, beams)
    u_n, v_n = _uv_lattice(du_near)
    d_n, _ = _nearest_beam(u_n, v_n, beams)
    u_p, v_p = _uv_lattice(patch_step)
    d_p, k_p = _nearest_beam(u_p, v_p, beams)
    far = d_f > near_radius_deg
    near = (d_n > patch_radius_deg) & (d_n <= near_radius_deg)
    patch = d_p <= patch_radius_deg
    u = np.concatenate([u_f[far], u_n[near], u_p[patch]])
    v = np.concatenate([v_f[far], v_n[near], v_p[patch]])
    tag = np.concatenate([np.full(int(far.sum()), -1), np.full(int(near.sum()), -1), k_p[patch]])
    return _from_uv(u, v, tag)


def theta_phi_direction_set(theta_step_deg=0.5, phi_step_deg=1.0, theta_max_deg=90.0):
    """A regular `(theta, phi)` grid, theta 0..theta_max inclusive, phi 0..360 exclusive -- row-major
    `(n_theta, n_phi)` flattened. The held-out grid (0.5 x 1 deg) is this, evaluated only, never optimised."""
    theta = np.arange(0.0, theta_max_deg + 1e-9, theta_step_deg)
    phi = np.arange(0.0, 360.0, phi_step_deg)
    tg, pg = np.meshgrid(theta, phi, indexing="ij")
    t, p = tg.ravel(), pg.ravel()
    u = np.sin(np.radians(t)) * np.cos(np.radians(p))
    v = np.sin(np.radians(t)) * np.sin(np.radians(p))
    return DirectionSet(t, p, u, v, np.full(t.size, -1)), (theta.size, phi.size)


def _unit_current_responses(directions: DirectionSet, frequency_hz, green_function, polarization):
    """`(a_x, a_y)`: the scalar output `polarization` ('rhcp', 'lhcp', 'hp', ... or 'theta'/'phi' for a
    raw component) for a unit x- and y-directed spectral current, at every direction -- computed with
    the SAME dyad / obliquity / projection `vector_aperture_pattern._vector_far_field` uses (a single
    1x1 cell at the origin with dx = dy = 1 gives `exp(0) * 1` as the transform)."""
    zero = np.zeros((1, 1))
    one = np.ones((1, 1))
    out = []
    for jx, jy in ((one, zero), (zero, one)):
        e_theta, e_phi = _vector_far_field(
            directions.theta_deg, directions.phi_deg, zero, zero, jx, jy, frequency_hz, green_function, 1.0, 1.0
        )
        if polarization == "theta":
            out.append(e_theta)
        elif polarization == "phi":
            out.append(e_phi)
        else:
            out.append(copol_component(e_theta, e_phi, polarization))
    return out[0], out[1]


def cell_contributions(
    x_cells_m,
    y_cells_m,
    feed_xy_m,
    directions: DirectionSet,
    frequency_hz,
    green_function,
    dx_m,
    dy_m,
    polarization,
    dtype=np.complex128,
    chunk_bytes=DEFAULT_CHUNK_BYTES,
):
    """`c[n, d]` (complex, shape `(N, D)`): cell `n`'s contribution to the `polarization` component in
    direction `d`, per unit scalar aperture field `s_n` (module docstring). Built chunked over directions
    so the transient `exp` matrix stays within `chunk_bytes`."""
    x = np.asarray(x_cells_m, dtype=float).ravel()
    y = np.asarray(y_cells_m, dtype=float).ravel()
    n = x.size
    jx_unit, jy_unit = aperture_current_from_scalar(x, y, np.ones(n), feed_xy_m=feed_xy_m)
    a_x, a_y = _unit_current_responses(directions, frequency_hz, green_function, polarization)
    k0 = free_space_wavenumber(frequency_hz)
    th, ph = np.radians(directions.theta_deg), np.radians(directions.phi_deg)
    kx, ky = k0 * np.sin(th) * np.cos(ph), k0 * np.sin(th) * np.sin(ph)
    out = np.empty((n, directions.size), dtype=dtype)
    per_col = n * 16 * 3  # exp matrix (complex128) + temporaries
    step = max(1, int(chunk_bytes // per_col))
    for lo in range(0, directions.size, step):
        hi = min(lo + step, directions.size)
        phase = np.outer(x, kx[lo:hi]) + np.outer(y, ky[lo:hi])
        block = np.exp(1j * phase)
        block *= (jx_unit[:, None] * a_x[None, lo:hi] + jy_unit[:, None] * a_y[None, lo:hi]) * (dx_m * dy_m)
        out[:, lo:hi] = block
    return out


def verify_against_vector_pattern(
    x_m, y_m, scalar_field, feed_xy_m, directions: DirectionSet, frequency_hz, green_function, dx_m, dy_m, polarization
):
    """`(cellwise, direct)` complex arrays of the `polarization` component: `sum_n s_n c[n, d]` against
    `_vector_far_field` on the full grid -- they agree to rounding (module docstring)."""
    s = np.asarray(scalar_field, dtype=complex)
    c = cell_contributions(
        x_m.ravel(), y_m.ravel(), feed_xy_m, directions, frequency_hz, green_function, dx_m, dy_m, polarization
    )
    cellwise = s.ravel() @ c
    jx, jy = aperture_current_from_scalar(x_m, y_m, s, feed_xy_m=feed_xy_m)
    e_theta, e_phi = _vector_far_field(
        directions.theta_deg, directions.phi_deg, x_m, y_m, jx, jy, frequency_hz, green_function, dx_m, dy_m
    )
    if polarization == "theta":
        direct = e_theta
    elif polarization == "phi":
        direct = e_phi
    else:
        direct = copol_component(e_theta, e_phi, polarization)
    return cellwise, direct


def combined_field(port_fields, excitations):
    """`sum_p excitations[p] * port_fields[p]` -- every port excited together. `excitations`: complex
    `a_p exp(j alpha_p)` per port."""
    total = None
    for f, w in zip(port_fields, excitations, strict=True):
        term = complex(w) * np.asarray(f)
        total = term if total is None else total + term
    return total
