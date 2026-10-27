"""Phase 29 design point B: ONE shared hexagonal aperture, 7 monopole ports, 7 holograms superimposed on the same cells.

Additive; reuses `tensor_cp_design` / `tensor_leaky_cp` unmodified. Labels of every result: MODEL_ONLY, NOT_HFSS_VALIDATED,
SOLVER_SUPERCELL, INTER_PORT_COUPLING_NOT_MODELLED (the port-alone fields superpose exactly under complex weights only if the
ports do not load each other: this is the model's assumption, HFSS is the check).

Model (the Phase 26 rotating-slot tensor cell, made multi-port). A cell has ONE slot (L, W, axis `a`). The surface wave of port j at
the cell is `s_j = sqrt(w_j / rho_j) exp(-j phi_j) (cos psi_j, sin psi_j)` (power per radian `w_j`, phase `phi_j`, ray azimuth
`psi_j` from port j), and the cell turns it into the current `J_j = kappa_j R_a s_j`, `R_a = [[cos 2a, sin 2a], [sin 2a, -cos 2a]]`,
`kappa_j = sqrt(2 alpha_j)` with `alpha_j` the table's leakage of that slot at port j's modulation wavenumber (this is exactly
`tensor_leaky_cp`'s `J = A exp(-j phi_ref)(cos chi, sin chi)`, `chi = 2a - psi`, `A = sqrt(2 alpha w / rho)`). Design: each port's own
ideal currents `t_j` (RHCP pencil at its (theta, phi), Gaussian conversion window `exp(-(r / w0)^2)` centred on the aperture so the
beam is broad) are fitted by ONE tensor `T = kappa R_a` per cell, least squares over the ports (closed form, real 2-parameter
`(p, q)`: p = b1 / S, q = b2 / S, S = sum |s_j|^2), then the slot size nearest in kappa; the surface-wave march is repeated on the
resulting cell map (analysis march) and the targets re-phased, a few times. Port-alone fields come from the final analysis march.
"""

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree

from ..pattern.vector_aperture_pattern import copol_component
from . import tensor_cp_design as t
from .leaky_forward import polygon_ray_distance
from .tensor_leaky_cp import modulation_k

BETA0 = 1.065
KEEPOUT_MM = 6.5
CELL_REACH_MM = 4.5


def hexagon_m(r_mm):
    return np.array(
        [
            [r_mm * 1e-3 * np.cos(np.radians(90 + 60 * i)), r_mm * 1e-3 * np.sin(np.radians(90 + 60 * i))]
            for i in range(6)
        ]
    )


def port_layout(ring_mm, tilt_deg=30.0, ring_side="opposite"):
    """Port 1 at the centre (broadside); ports 2..7: beam azimuth 0, 60, .. 300 deg, tilt `tilt_deg`; the monopole sits on a ring of
    radius `ring_mm` at the opposite azimuth (`ring_side='opposite'`: the surface wave sweeps the aperture towards the beam, a
    forward leaky wave) or at the beam azimuth ('same'). Returns `[(x_mm, y_mm, theta_deg, phi_deg)]`."""
    out = [(0.0, 0.0, 0.0, 0.0)]
    for k in range(6):
        ph = 60.0 * k
        pos = ph + (180.0 if ring_side == "opposite" else 0.0)
        out.append((ring_mm * np.cos(np.radians(pos)), ring_mm * np.sin(np.radians(pos)), tilt_deg, ph))
    return out


def aperture_cells(r_mm, ports):
    """Shared 6.4 mm lattice (through the origin) inside the hexagon of circumradius `r_mm`, minus a keep-out around every monopole."""
    n = int(np.ceil(r_mm / t.P_MM)) + 2
    k = np.arange(-n, n + 1)
    xg, yg = np.meshgrid(t.P_MM * k, t.P_MM * k)
    x, y = xg.ravel(), yg.ravel()
    ok = t.inside(hexagon_m(r_mm) * 1e3, x, y)
    for px, py, _, _ in ports:
        ok &= np.hypot(x - px, y - py) >= KEEPOUT_MM
    return x[ok], y[ok]


def _grid(poly_m, feed_m, n_psi, drho):
    psi = np.linspace(-np.pi, np.pi, n_psi, endpoint=False)
    r_edge = polygon_ray_distance(poly_m, feed_m, psi)
    n_rho = int(np.ceil(r_edge.max() / drho)) + 1
    return psi, r_edge, np.arange(n_rho) * drho


def _empty_march(psi, rho, r_edge, drho):
    shape = (psi.size, rho.size)
    return {k: np.zeros(shape) for k in ("alpha", "beta", "w", "phi_ref", "K")} | {
        "slot": np.full(shape, -1, dtype=int),
        "psi": psi,
        "rho": rho,
        "r_edge": r_edge,
        "drho": drho,
    }


def _finish(m, k0):
    w_edge = np.array(
        [m["w"][i, min(int(np.floor(m["r_edge"][i] / m["drho"])), m["rho"].size - 1)] for i in range(m["psi"].size)]
    )
    m["w_edge"], m["eta_conv"], m["eta_conv_by_ray"] = w_edge, float(np.mean(1.0 - w_edge)), 1.0 - w_edge
    return m


def march_design(table, poly_m, feed_m, k0, th0, ph0, q0, window, n_psi=360, drho=0.5e-3, rho_min=0.007):
    """Design march with a GLOBAL conversion density `q = q0 window(x, y)^2` (not normalised per ray, so rays that miss the window
    convert almost nothing). Same output keys as `tensor_leaky_cp.march_rays` (so `cell_design` consumes it)."""
    psi, r_edge, rho = _grid(poly_m, feed_m, n_psi, drho)
    m = _empty_march(psi, rho, r_edge, drho)
    inside = rho[None, :] <= r_edge[:, None]
    amax = np.nanmax(table.alpha, axis=0)
    cur_w, cur_phi, cur_beta = np.ones(psi.size), np.zeros(psi.size), np.full(psi.size, BETA0)
    for k, r in enumerate(rho):
        kvec = modulation_k(cur_beta, psi, th0, ph0)
        col = table.column(kvec)
        g = window(feed_m[0] + r * np.cos(psi), feed_m[1] + r * np.sin(psi))
        a_t = np.where(cur_w > 1e-12, q0 * g**2 * max(r, 1e-9) / (2.0 * np.maximum(cur_w, 1e-12)), amax[col])
        a_col = table.alpha[:, col]
        pick = np.nanargmin(np.abs(np.nan_to_num(a_col, nan=1e9) - a_t[None, :]), axis=0)
        a_used = np.where(inside[:, k] & (r >= rho_min), a_col[pick, np.arange(psi.size)], 0.0)
        b_used = table.beta[pick, col]
        m["alpha"][:, k], m["beta"][:, k], m["w"][:, k], m["phi_ref"][:, k], m["slot"][:, k], m["K"][:, k] = (
            a_used,
            b_used,
            cur_w,
            cur_phi,
            pick,
            kvec,
        )
        cur_w = cur_w * np.exp(-2.0 * a_used * drho)
        cur_phi = cur_phi + b_used * k0 * drho
        cur_beta = b_used
    return _finish(m, k0)


def march_analysis(table, poly_m, feed_m, k0, th0, ph0, cell_xy_m, cell_slot, n_psi=360, drho=0.5e-3, rho_min=0.007):
    """The same march over a GIVEN cell map: at every step the slot of the nearest cell (within `CELL_REACH_MM`; none -> no leakage)
    gives `alpha`, `beta` at the local modulation wavenumber of this port's beam."""
    psi, r_edge, rho = _grid(poly_m, feed_m, n_psi, drho)
    m = _empty_march(psi, rho, r_edge, drho)
    px = feed_m[0] + rho[None, :] * np.cos(psi)[:, None]
    py = feed_m[1] + rho[None, :] * np.sin(psi)[:, None]
    dist, idx = cKDTree(cell_xy_m).query(np.c_[px.ravel(), py.ravel()])
    node_slot = np.where(
        dist.reshape(px.shape) <= CELL_REACH_MM * 1e-3, np.asarray(cell_slot)[idx].reshape(px.shape), -1
    )
    inside = rho[None, :] <= r_edge[:, None]
    cur_w, cur_phi, cur_beta = np.ones(psi.size), np.zeros(psi.size), np.full(psi.size, BETA0)
    for k, r in enumerate(rho):
        kvec = modulation_k(cur_beta, psi, th0, ph0)
        col = table.column(kvec)
        s = node_slot[:, k]
        has = (s >= 0) & inside[:, k] & (r >= rho_min)
        a_used = np.where(has, np.nan_to_num(table.alpha[np.maximum(s, 0), col]), 0.0)
        b_used = np.where(has, table.beta[np.maximum(s, 0), col], BETA0)
        m["alpha"][:, k], m["beta"][:, k], m["w"][:, k], m["phi_ref"][:, k], m["slot"][:, k], m["K"][:, k] = (
            a_used,
            b_used,
            cur_w,
            cur_phi,
            s,
            kvec,
        )
        cur_w = cur_w * np.exp(-2.0 * a_used * drho)
        cur_phi = cur_phi + b_used * k0 * drho
        cur_beta = b_used
    return _finish(m, k0)


def incident(march, feed_m, x_m, y_m):
    """Per cell: surface-wave vector of one port `s = sqrt(w/rho) exp(-j phi) (cos psi, sin psi)` (n, 2) complex, the cell's `alpha`, `rho`, `psi`."""
    dx, dy = x_m - feed_m[0], y_m - feed_m[1]
    rho, psi_c = np.hypot(dx, dy), np.arctan2(dy, dx)
    psi, rg = march["psi"], march["rho"]
    pp = np.concatenate([[psi[-1] - 2 * np.pi], psi, [psi[0] + 2 * np.pi]])

    def interp(f):
        f = np.vstack([f[-1:], f, f[:1]])
        return RegularGridInterpolator((pp, rg), f, bounds_error=False, fill_value=None)(
            (psi_c, np.minimum(rho, rg[-1]))
        )

    w, ph = np.maximum(interp(march["w"]), 0.0), interp(march["phi_ref"])
    ip = np.argmin(np.abs(((psi_c[:, None] - psi[None, :] + np.pi) % (2 * np.pi)) - np.pi), axis=1)
    ir = np.minimum(np.round(rho / march["drho"]).astype(int), len(rg) - 1)
    amp = np.sqrt(w / np.maximum(rho, 1e-9)) * np.exp(-1j * ph)
    return np.stack([amp * np.cos(psi_c), amp * np.sin(psi_c)], axis=-1), march["alpha"][ip, ir], rho, psi_c, ph


def fit_tensor(s_list, t_list, omega=None):
    """Closed-form (weighted) least squares of `T = [[p, q], [q, -p]]` (real) over the ports: min sum_j omega_j |T s_j - t_j|^2.
    `s_list`, `t_list`: `(n_ports, n_cells, 2)` complex; `omega` `(n_ports, n_cells)` (default 1). Returns `p, q, S`."""
    sx, sy, tx, ty = s_list[..., 0], s_list[..., 1], t_list[..., 0], t_list[..., 1]
    om = 1.0 if omega is None else omega
    s2 = np.sum(om * (np.abs(sx) ** 2 + np.abs(sy) ** 2), axis=0)
    b1 = np.sum(om * np.real(np.conj(sx) * tx - np.conj(sy) * ty), axis=0)
    b2 = np.sum(om * np.real(np.conj(sy) * tx + np.conj(sx) * ty), axis=0)
    s2 = np.maximum(s2, 1e-30)
    return b1 / s2, b2 / s2, s2


def tensor_current(s, kappa, axis):
    """`J = kappa R_a s` for one port: `s` (n, 2) complex, `kappa` (n,), `axis` (n,) rad."""
    c2, s2 = np.cos(2 * axis), np.sin(2 * axis)
    return np.stack([kappa * (c2 * s[:, 0] + s2 * s[:, 1]), kappa * (s2 * s[:, 0] - c2 * s[:, 1])], axis=-1)


def design_shared(
    table,
    r_mm,
    ring_mm,
    w0_mm,
    q0,
    tilt_deg=30.0,
    ring_side="opposite",
    n_iter=4,
    n_psi=360,
    gamma=0.0,
    scheme="ls",
    hpbw_target_deg=32.0,
    ff_iters=(300, 150),
    log=None,
):
    """The shared-aperture design. `gamma`: port weight exponent of the fit, `omega_j = (|s_j|^2 / max_k |s_k|^2)^gamma` per cell
    (0: plain least squares over the ports; large: each cell follows the port whose wave is strongest there = spatial multiplexing).
    `scheme='ff'`: after the least-squares start, the whole cell map is synthesised in the PATTERN domain (`synthesise_farfield`: the far
    field is linear in the per-cell `(p, q)`, so the seven RHCP target beams of HPBW `hpbw_target_deg` are fitted directly, cross-talk
    and sidelobes included), then snapped to the discrete slot levels. Returns a dict: ports, cells (x/y mm, slot index/L/W, axis deg), per-port `J` (n, 2) complex,
    per-port `eta_conv`, the per-iteration fit residuals."""
    ports = port_layout(ring_mm, tilt_deg, ring_side)
    xm, ym = (a * 1e-3 for a in aperture_cells(r_mm, ports))
    poly = hexagon_m(r_mm)
    k0 = t.K0

    def window(x, y):
        return np.exp(-((x**2 + y**2) / (w0_mm * 1e-3) ** 2))

    n_p = len(ports)
    feeds = [(p[0] * 1e-3, p[1] * 1e-3) for p in ports]
    beams = [(np.radians(p[2]), np.radians(p[3])) for p in ports]
    des = []
    for j in range(n_p):
        m = march_design(table, poly, feeds[j], k0, beams[j][0], beams[j][1], q0, window, n_psi)
        d = t.cell_design(m, table, xm, ym, feeds[j], k0, beams[j][0], beams[j][1], hand=+1)
        des.append((m, d))
    kappa_des = [np.sqrt(2.0 * d["alpha"]) for _, d in des]
    s_list = np.zeros((n_p, xm.size, 2), dtype=complex)
    t_list = np.zeros_like(s_list)
    for j, (m, d) in enumerate(des):
        s_list[j] = incident(m, feeds[j], xm, ym)[0]
        t_list[j] = np.stack([d["jx"], d["jy"]], axis=-1)
    hist = []
    psi_pj = np.stack([np.arctan2(ym - f[1], xm - f[0]) for f in feeds])
    col = np.array([table.column(modulation_k(BETA0, psi_pj[j], *beams[j])) for j in range(n_p)])  # (n_ports, n_cells)
    for it in range(n_iter + 1):
        w2 = np.sum(np.abs(s_list) ** 2, axis=-1)
        omega = (w2 / np.maximum(w2.max(axis=0, keepdims=True), 1e-30)) ** gamma
        p, q, s2 = fit_tensor(s_list, t_list, omega)
        kappa_t, axis = np.hypot(p, q), np.mod(0.5 * np.arctan2(q, p), np.pi)
        # slot nearest in kappa, with alpha averaged over the ports weighted by |s_j|^2
        wgt = np.sum(np.abs(s_list) ** 2, axis=-1)  # (n_ports, n_cells)
        a_slot = np.array(
            [
                np.sum(wgt * np.nan_to_num(table.alpha[i][col]), axis=0) / np.maximum(np.sum(wgt, axis=0), 1e-30)
                for i in range(len(table.slots))
            ]
        )  # (n_slots, n_cells)
        slot = np.argmin(np.abs(np.sqrt(2.0 * a_slot) - kappa_t[None, :]), axis=0)
        ana = [
            march_analysis(table, poly, feeds[j], k0, beams[j][0], beams[j][1], np.c_[xm, ym], slot, n_psi)
            for j in range(n_p)
        ]
        out_j = []
        for j in range(n_p):
            s, alpha, rho, psi_c, ph = incident(ana[j], feeds[j], xm, ym)
            out_j.append((s, alpha, ph))
        err = 0.0
        for j in range(n_p):
            s, alpha, ph = out_j[j]
            kap = np.sqrt(2.0 * alpha)
            jj = tensor_current(s, kap, axis)
            tgt = _target(kappa_des[j], s, ph, xm - feeds[j][0], ym - feeds[j][1], k0, beams[j])
            err += float(np.sum(np.abs(jj - tgt) ** 2) / max(np.sum(np.abs(tgt) ** 2), 1e-30))
            s_list[j], t_list[j] = s, tgt
        hist.append(err / n_p)
        if log:
            log(
                f"iter {it}: mean relative fit error {err / n_p:.3f}, eta_conv "
                + " ".join(f"{m_['eta_conv']:.2f}" for m_ in ana)
            )
    if scheme == "ff":
        from ..green_function import LayeredGreenFunction

        green = LayeredGreenFunction(t.SUB, t.F, backing="pec")
        setup = _ff_setup(green, xm, ym)
        lev_k = np.sqrt(2.0 * np.sort(np.mean(a_slot, axis=1)))
        sel = np.argsort(np.mean(a_slot, axis=1))
        x0 = np.concatenate(
            [
                np.sqrt(2.0 * np.array([np.mean(a_slot, axis=1)[i] for i in slot])) * np.cos(2 * axis),
                np.sqrt(2.0 * np.array([np.mean(a_slot, axis=1)[i] for i in slot])) * np.sin(2 * axis),
            ]
        )
        _ = sel
        for it in range(n_iter):
            kap, axis, ffh = synthesise_farfield(setup, s_list, beams, lev_k, x0, hpbw_target_deg, ff_iters)
            a_slot = np.array(
                [
                    np.sum(wgt * np.nan_to_num(table.alpha[i][col]), axis=0) / np.maximum(np.sum(wgt, axis=0), 1e-30)
                    for i in range(len(table.slots))
                ]
            )
            slot = np.argmin(np.abs(np.sqrt(2.0 * a_slot) - kap[None, :]), axis=0)
            ana = [
                march_analysis(table, poly, feeds[j], k0, beams[j][0], beams[j][1], np.c_[xm, ym], slot, n_psi)
                for j in range(n_p)
            ]
            out_j = []
            for j in range(n_p):
                s_, alpha, rho, psi_c, ph = incident(ana[j], feeds[j], xm, ym)
                out_j.append((s_, alpha, ph))
                s_list[j] = s_
            wgt = np.sum(np.abs(s_list) ** 2, axis=-1)
            kk = np.sqrt(2.0 * np.array([np.mean(a_slot, axis=1)[i] for i in slot]))
            x0 = np.concatenate([kk * np.cos(2 * axis), kk * np.sin(2 * axis)])
            hist.append(ffh[-1])
            if log:
                log(
                    f"ff outer {it}: objective {ffh[0]:.3f} -> {ffh[-1]:.3f}, eta_conv "
                    + " ".join(f"{m_['eta_conv']:.2f}" for m_ in ana)
                )
    cur = []
    for j in range(n_p):
        s, alpha, ph = out_j[j]
        cur.append(tensor_current(s, np.sqrt(2.0 * alpha), axis))
    L = np.array([table.slots[s][0] for s in slot])
    W = np.array([table.slots[s][1] for s in slot])
    return {
        "r_mm": r_mm,
        "ring_mm": ring_mm,
        "w0_mm": w0_mm,
        "q0": q0,
        "ports": ports,
        "x_m": xm,
        "y_m": ym,
        "slot": slot,
        "slot_L_mm": L,
        "slot_W_mm": W,
        "axis_rad": axis,
        "J": cur,
        "eta_conv": [a["eta_conv"] for a in ana],
        "march": ana,
        "fit_error_history": hist,
        "feeds_m": feeds,
        "beams": beams,
    }


def _ff_setup(green, xm, ym, theta_step=4.0, phi_step=8.0):
    """Pattern-domain operator pieces on a hemisphere grid: per-direction phase matrix and the RHCP / LHCP field of a unit x / y current at the origin."""
    from ..pattern.vector_aperture_pattern import _vector_far_field

    th, ph = np.meshgrid(np.arange(0.0, 90.0 + 1e-9, theta_step), np.arange(0.0, 360.0, phi_step), indexing="ij")
    th, ph = th.ravel(), ph.ravel()
    kx, ky = (
        t.K0 * np.sin(np.radians(th)) * np.cos(np.radians(ph)),
        t.K0 * np.sin(np.radians(th)) * np.sin(np.radians(ph)),
    )
    bR, bL = np.zeros((th.size, 2), dtype=complex), np.zeros((th.size, 2), dtype=complex)
    for m, (jx, jy) in enumerate(((1.0, 0.0), (0.0, 1.0))):
        et, ep = _vector_far_field(th, ph, [0.0], [0.0], [jx], [jy], t.F, green, t.PITCH, t.PITCH)
        bR[:, m], bL[:, m] = rhcp_field(et, ep), copol_component(et, ep, "lhcp")
    return {"th": th, "ph": ph, "phase": np.exp(1j * (np.outer(kx, xm) + np.outer(ky, ym))), "bR": bR, "bL": bL}


def _ff_operator(setup, s):
    """`M` (2 n_dir, 2 n_cells) complex: rows RHCP then LHCP field, columns `p` then `q` of `T = [[p, q], [q, -p]]` acting on `s` (n, 2)."""
    ph, bR, bL = setup["phase"], setup["bR"], setup["bL"]
    sx, sy = s[:, 0][None, :], s[:, 1][None, :]

    def rows(b):
        return np.hstack([ph * (b[:, 0:1] * sx - b[:, 1:2] * sy), ph * (b[:, 0:1] * sy + b[:, 1:2] * sx)])

    return np.vstack([rows(bR), rows(bL)])


def _beam_target(setup, beam, hpbw_deg):
    th, ph = setup["th"], setup["ph"]
    c = np.cos(np.radians(th)) * np.cos(beam[0]) + np.sin(np.radians(th)) * np.sin(beam[0]) * np.cos(
        np.radians(ph) - beam[1]
    )
    ang = np.degrees(np.arccos(np.clip(c, -1, 1)))
    return np.exp(-2.0 * np.log(2.0) * (ang / hpbw_deg) ** 2)


def synthesise_farfield(setup, s_list, beams, kappa_levels, x0, hpbw_deg, iters=(300, 150)):
    """Projected-gradient synthesis of the per-cell `(p, q)` against the seven RHCP target beams (shape only: a free complex scale per
    port is eliminated in closed form, so the amplitude is set by the physical `kappa` range). `kappa_levels`: ascending allowed
    `kappa` values; first `iters[0]` iterations project onto the annulus [min, max], the last `iters[1]` onto the discrete levels.
    Returns `(kappa, axis_rad, objective_history)`."""
    n_p, n = s_list.shape[0], s_list.shape[1]
    ops = [_ff_operator(setup, s_list[j]) for j in range(n_p)]
    nd = setup["th"].size
    tgt = [np.concatenate([_beam_target(setup, b, hpbw_deg), np.zeros(nd)]).astype(complex) for b in beams]
    wdir = np.maximum(np.sin(np.radians(setup["th"])), 0.05)
    w = np.concatenate([wdir, wdir])
    tw = [float(np.real(np.vdot(tj, w * tj))) for tj in tgt]
    lo, hi = float(kappa_levels[0]), float(kappa_levels[-1])
    lev = np.asarray(kappa_levels, dtype=float)

    def project(x, hard):
        z = x[:n] + 1j * x[n:]
        r = np.abs(z)
        ang = np.where(r > 1e-12, z / np.maximum(r, 1e-30), 1.0)
        r2 = lev[np.argmin(np.abs(r[:, None] - lev[None, :]), axis=1)] if hard else np.clip(r, lo, hi)
        z = r2 * ang
        return np.concatenate([z.real, z.imag])

    def objective(x):
        f, g = 0.0, np.zeros_like(x)
        for j in range(n_p):
            m = ops[j] @ x
            lam = np.vdot(tgt[j], w * m) / tw[j]
            r = m - lam * tgt[j]
            f += float(np.real(np.vdot(r, w * r))) / tw[j]
            g += 2.0 * np.real(ops[j].conj().T @ (w * r)) / tw[j]
        return f, g

    v = np.random.default_rng(0).normal(size=2 * n)
    for _ in range(15):  # power iteration on the (shape-normalised) Hessian bound
        hv = sum(2.0 * np.real(o.conj().T @ (w * (o @ v))) / tw[j] for j, o in enumerate(ops))
        lip = float(np.linalg.norm(hv) / max(np.linalg.norm(v), 1e-30))
        v = hv / max(np.linalg.norm(hv), 1e-30)
    step = 1.0 / max(lip, 1e-12)
    x, hist = project(x0, False), []
    for k in range(sum(iters)):
        f, g = objective(x)
        hist.append(f / n_p)
        x = project(x - step * g, k >= iters[0])
    z = x[:n] + 1j * x[n:]
    return np.abs(z), np.mod(0.5 * np.angle(z), np.pi), hist


def _target(kappa_des, s, ph, dx, dy, k0, beam):
    """Port-j ideal current on the analysis march: `kappa_des |s| exp(-j phi)(cos chi, sin chi)`, `chi = phi - k_t . r` (RHCP)."""
    kt = k0 * np.sin(beam[0])
    chi = ph - kt * (dx * np.cos(beam[1]) + dy * np.sin(beam[1]))
    amp = kappa_des * np.sqrt(np.sum(np.abs(s) ** 2, axis=-1)) * np.exp(-1j * ph)
    return np.stack([amp * np.cos(chi), amp * np.sin(chi)], axis=-1)


def port_fields(design, green, theta, phi):
    """Port-alone fields `[(E_theta, E_phi)]` at the real cell positions (origin = aperture centre), each scaled so that its radiated
    power is `LAUNCH * eta_conv_j` of 1 W accepted at that port (no neighbour pickup, no port coupling)."""
    out = []
    for j, J in enumerate(design["J"]):
        d = {"jx": J[:, 0], "jy": J[:, 1]}
        em = t.pattern(d, design["x_m"], design["y_m"], green, theta, phi)
        s = np.sqrt(t.LAUNCH * design["eta_conv"][j] / t.prad(*em, theta, phi))
        out.append((em[0] * s, em[1] * s))
    return out


def rhcp_field(et, ep):
    return copol_component(et, ep, "rhcp")


def gain_db(et, ep, p_acc=1.0):
    g = np.abs(rhcp_field(et, ep)) ** 2 * 4 * np.pi * 10 ** (t.OFFSET_DB / 10) / p_acc
    return 10 * np.log10(np.maximum(g, 1e-30))


def _hpbw_equiv(g_db, theta, phi):
    tg, _ = np.meshgrid(theta, phi, indexing="ij")
    w = np.sin(np.radians(np.maximum(tg, 0.25))) * np.radians(theta[1] - theta[0]) * np.radians(phi[1] - phi[0])
    return float(np.degrees(2 * np.sqrt(np.sum(w[g_db >= g_db.max() - 3.0]) / np.pi)))


def scan_directions(theta_max=45.0, theta_step=3.0, phi_step=15.0):
    """Direction grid of the weight table: broadside + theta = step .. theta_max (all phi)."""
    out = [(0.0, 0.0)]
    for th in np.arange(theta_step, theta_max + 1e-9, theta_step):
        out += [(float(th), float(ph)) for ph in np.arange(0.0, 360.0, phi_step)]
    return out


def mrc_weights(fields, theta, phi, direction, n_active=3):
    """Weights `(7,)` complex, `sum |w|^2 = 1`, non-zero on the `n_active` ports with the largest RHCP field at `direction`
    (maximum-ratio combining: the best peak gain at that direction for those ports, ports taken as uncoupled)."""
    it, ip = int(round(direction[0] / (theta[1] - theta[0]))), int(round(direction[1] / (phi[1] - phi[0]))) % phi.size
    r = np.array([rhcp_field(*f)[it, ip] for f in fields])
    keep = np.argsort(-np.abs(r))[:n_active]
    w = np.zeros(len(fields), dtype=complex)
    w[keep] = np.conj(r[keep])
    return w / np.linalg.norm(w)


def combine(fields, w):
    et = sum(wj * f[0] for wj, f in zip(w, fields, strict=True) if wj != 0)
    ep = sum(wj * f[1] for wj, f in zip(w, fields, strict=True) if wj != 0)
    return et, ep


def scan_table(fields, theta, phi, directions=None, n_active=3):
    """For every direction: MRC weights, the combined RHCP pattern, its peak (gain, direction), pointing error, equivalent HPBW, SLL
    (RHCP maximum outside a cone of 1.5 x HPBW about the peak, relative to the peak). Gain: 1 W accepted in total (sum |w|^2 = 1),
    anchored by `OFFSET_DB`. Returns a list of dicts (weights as `[re, im]` pairs)."""
    from .. import hologram  # noqa: F401  (package context only)
    from .seven_hex_cluster import _unit, angle_between, lobes

    directions = scan_directions() if directions is None else directions
    tg, pg = np.meshgrid(theta, phi, indexing="ij")
    unit = _unit(tg, pg)
    rows = []
    for d in directions:
        w = mrc_weights(fields, theta, phi, d, n_active)
        g = gain_db(*combine(fields, w))
        lb = lobes(g, theta, phi, within_db=3.0)
        pk = (lb[0][1], lb[0][2])
        hp = _hpbw_equiv(g, theta, phi)
        ang = np.degrees(np.arccos(np.clip(np.sum(unit * _unit(*pk), axis=-1), -1, 1)))
        out = ang > 1.5 * hp
        rows.append(
            {
                "theta_deg": d[0],
                "phi_deg": d[1],
                "weights": [[float(x.real), float(x.imag)] for x in w],
                "peak_gain_dbi": float(lb[0][0]),
                "peak_theta_phi_deg": [float(pk[0]), float(pk[1])],
                "pointing_error_deg": angle_between(pk, d),
                "hpbw_deg": hp,
                "sll_db": float(g[out].max() - g.max()) if out.any() else float("nan"),
                "n_lobes_3db": len(lb),
            }
        )
    return rows


def scan_summary(rows, theta_split=30.0):
    """Pass test of the brief: pointing within +-3 deg and gain dip <= 3 dB (dip = best peak gain of the table minus this beam's) for theta <= 30;
    the 30-45 deg band reported as is."""
    gmax = max(r["peak_gain_dbi"] for r in rows)
    out = {"gain_ref_dbi": gmax}
    for name, sel in (
        ("inner_le_30", lambda r: r["theta_deg"] <= theta_split),
        ("band_30_45", lambda r: r["theta_deg"] > theta_split),
    ):
        rs = [r for r in rows if sel(r)]
        out[name] = {
            "n": len(rs),
            "max_pointing_error_deg": max(r["pointing_error_deg"] for r in rs),
            "max_gain_dip_db": max(gmax - r["peak_gain_dbi"] for r in rs),
            "min_gain_dbi": min(r["peak_gain_dbi"] for r in rs),
            "worst_sll_db": max(r["sll_db"] for r in rs),
            "mean_sll_db": float(np.mean([r["sll_db"] for r in rs])),
        }
    i = out["inner_le_30"]
    out["pass_inner"] = bool(i["max_pointing_error_deg"] <= 3.0 and i["max_gain_dip_db"] <= 3.0)
    return out


def optimise_weights(
    fields, theta, phi, direction, n_active=3, inner_deg=12.0, point_deg=3.0, margin=0.63, mu=(10.0, 100.0, 1000.0)
):
    """Weights `(7,)` for `direction` on the `n_active` ports with the strongest field there, by penalised maximisation of the power at the
    direction with: every point > `point_deg` away not above it (so the peak is at the direction within `point_deg`), and every point
    > `inner_deg` away below `margin` x it (a 2 dB lobe margin). Start: maximum-ratio weights. `sum |w|^2 = 1`."""
    from scipy.optimize import minimize

    w0 = mrc_weights(fields, theta, phi, direction, n_active)
    keep = np.nonzero(w0)[0]
    tg, pg = np.meshgrid(theta, phi, indexing="ij")
    from .seven_hex_cluster import _unit

    ang = np.degrees(np.arccos(np.clip(np.sum(_unit(tg, pg) * _unit(*direction), axis=-1), -1, 1))).ravel()
    cap = np.where(ang > inner_deg, margin, np.where(ang > point_deg, 1.0, np.inf))
    act = np.isfinite(cap)
    rk = np.stack([rhcp_field(*fields[j]).ravel() for j in keep])  # (k, npts)
    it, ip = int(round(direction[0] / (theta[1] - theta[0]))), int(round(direction[1] / (phi[1] - phi[0]))) % phi.size
    rd = np.array([rhcp_field(*fields[j])[it, ip] for j in keep])

    def unpack(x):
        z = x[: keep.size] + 1j * x[keep.size :]
        return z / max(np.linalg.norm(z), 1e-12)

    def f(x, m):
        z = unpack(x)
        ud = np.abs(z @ rd) ** 2
        u = np.abs(z @ rk[:, act]) ** 2
        pen = np.maximum(u - cap[act] * ud, 0.0)
        return -ud + m * float(np.sum(pen**2)) / max(ud, 1e-12) / act.sum() * 100.0

    x = np.concatenate([w0[keep].real, w0[keep].imag])
    for m in mu:
        x = minimize(f, x, args=(m,), method="L-BFGS-B", options={"maxiter": 60}).x
    w = np.zeros(len(fields), dtype=complex)
    w[keep] = unpack(x)
    return w
