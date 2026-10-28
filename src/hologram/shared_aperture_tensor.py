"""Phase 29: angle-dependent, full-tensor forward model of the B3 cell family and the Frey-style synthesis on it (MODEL_ONLY, additive; frozen modules untouched).

Replaces the radial reflection law of `shared_aperture_multiport` (`J = kappa R_a s`: unit-strength reflection of the incident vector about the slot axis, independent of the incidence angle) by the
Bloch sheet tensor of the cell from periodic MoM (`results/phase29/angle_lut_v1.npz`, `scripts/research/phase29_angle_lut.py`): `Y_s(slot, a, psi)`, `a` the slot orientation, `psi` the direction of the
incident surface wave AT THAT CELL (each feed's own local direction, so any feed position works, off-centre ones included).

Cell response. With `k = (cos psi, sin psi)`, `s_hat = z x k`, the part of the sheet that varies with the orientation (the part that radiates; its orientation average sets the phase constant) is
`dY(a, psi) = Y_s(a, psi) - <Y_s>_a(psi)`; the current driven by the TM0 surface wave `amp * k` (amp = sqrt(w / rho) exp(-j phi)) is
    `J = amp c_slot g_K (dY_kk k + dY_sk s_hat)`  (complex: an elliptical current, so cross-polar content is MODELLED, not assumed away),
`dY_kk = k^T dY k`, `dY_sk = s^T dY k`. `c_slot` makes the orientation-averaged strength of each slot equal to the validated collinear supercell-table strength `sqrt(2 alpha)` (so the angle dependence is
the LUT's, the absolute level the table's), `g_K = sqrt(alpha[slot, K] / mean_K alpha[slot])` keeps the table's modulation-wavenumber dependence. Energy: `alpha = |c g dY|^2 / 2`, radiated density
`|J|^2 = 2 alpha w / rho` (same bookkeeping as before). Phase constant: `beta = beta_table + beta_uniform(slot, a, psi) - <beta_uniform>_a`.
Synthesis (Frey et al., IEEE OJAP 2023, Section II-C/E generalised): per cell the slot (4 table slots) and orientation (2 deg grid) minimise the summed tensor error
`sum_n |s_n|^2 || f_n(slot, a) - kappa_des,n (cos(chi_n - psi_n), sin(chi_n - psi_n)) ||^2` over the feeds, each feed with ITS OWN local incidence angle `psi_n` and the full complex tensor (the
eigenvector-type assignment on all tensor components, by direct minimisation instead of an analytic average); `chi_n` the ideal objective direction (hologram offset `psi0_n` included);
the surface-wave march is then repeated on the resulting cell map (phases `phi_ref` update) `n_iter` times."""

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree

from . import shared_aperture_multiport as sa
from . import tensor_cp_design as t
from .tensor_leaky_cp import march_rays, modulation_k

N_DENSE = 90  # 2 deg grid in a and psi (period pi)


class AngleLUT:
    """Dense tables derived from `angle_lut_v1.npz` for a `SupercellTable`."""

    def __init__(self, npz_path, table):
        d = np.load(npz_path)
        self.slots = [tuple(x) for x in d["slots_mm"]]
        assert self.slots == [tuple(x) for x in table.slots], "LUT slots differ from the supercell table"
        y12, b12 = d["Y"], d["beta_over_k0"]
        self.step = np.pi / N_DENSE
        grid = np.arange(12) * (np.pi / 12)
        pad = np.concatenate([grid[-3:] - np.pi, grid, grid[:3] + np.pi])
        dense = np.arange(N_DENSE) * self.step
        aa, pp = np.meshgrid(dense, dense, indexing="ij")
        pts = np.stack([aa.ravel(), pp.ravel()], axis=-1)

        def up(arr):  # (12, 12) real -> (90, 90), periodic cubic
            p = np.concatenate([arr[-3:], arr, arr[:3]], axis=0)
            p = np.concatenate([p[:, -3:], p, p[:, :3]], axis=1)
            return RegularGridInterpolator((pad, pad), p, method="cubic")(pts).reshape(N_DENSE, N_DENSE)

        ns = len(self.slots)
        self.Y = np.zeros((ns, N_DENSE, N_DENSE, 2, 2), dtype=complex)
        self.beta = np.zeros((ns, N_DENSE, N_DENSE))
        for s in range(ns):
            for i in range(2):
                for j in range(2):
                    self.Y[s, :, :, i, j] = up(y12[s, :, :, i, j].real) + 1j * up(y12[s, :, :, i, j].imag)
            b = b12[s]
            self.beta[s] = up(np.where(np.isnan(b), np.nanmean(b), b))
        psi = dense
        k = np.stack([np.cos(psi), np.sin(psi)], axis=-1)  # (psi, 2)
        sh = np.stack([-np.sin(psi), np.cos(psi)], axis=-1)
        ybar = self.Y.mean(axis=1, keepdims=True)  # orientation average at fixed psi: (slot, 1, psi, 2, 2)
        dY = self.Y - ybar
        dkk = np.einsum("sapij,pi,pj->sap", dY, k, k)
        dsk = np.einsum("sapij,pi,pj->sap", dY, sh, k)
        g0 = -0.5 * np.angle(
            np.sum(dkk**2)
        )  # one global phase for all slots / angles / feeds (irrelevant for the far field)
        self.dkk, self.dsk = dkk * np.exp(1j * g0), dsk * np.exp(1j * g0)
        kappa_old = np.nanmean(np.sqrt(2.0 * np.nan_to_num(table.alpha, nan=np.nan)), axis=1)
        mag = np.sqrt(np.abs(self.dkk) ** 2 + np.abs(self.dsk) ** 2).mean(axis=(1, 2))
        self.c = kappa_old / mag  # per slot
        mean_alpha = np.nanmean(table.alpha, axis=1)
        self.g_col = np.sqrt(np.nan_to_num(table.alpha, nan=0.0) / mean_alpha[:, None])  # (slot, K column)
        self.dbeta = self.beta - self.beta.mean(axis=1, keepdims=True)
        self.table = table

    def _ip(self, psi):
        x = (np.asarray(psi) % np.pi) / self.step
        i0 = np.floor(x).astype(int) % N_DENSE
        return i0, (i0 + 1) % N_DENSE, x - np.floor(x)

    def tables_at(self, psi):
        """`dkk, dsk` (slot, axis, n) and `dbeta` at the directions `psi` (n,) (linear interpolation in psi)."""
        i0, i1, w = self._ip(psi)
        mix = lambda a: a[:, :, i0] * (1 - w) + a[:, :, i1] * w  # noqa: E731
        return mix(self.dkk), mix(self.dsk), mix(self.dbeta)

    def response(self, slot, ia, psi, col):
        """Complex `(f_kk, f_sk)` and `dbeta` for cells / nodes with slot index, axis index `ia` (0..89, 2 deg), incidence direction `psi`, K column."""
        i0, i1, w = self._ip(psi)
        pick = lambda a: a[slot, ia, i0] * (1 - w) + a[slot, ia, i1] * w  # noqa: E731
        sc = self.c[slot] * self.g_col[slot, col]
        return sc * pick(self.dkk), sc * pick(self.dsk), pick(self.dbeta)


def aperture_cells(r_mm, ports, centre_mm=(0.0, 0.0)):
    """Shared 6.4 mm lattice (through the ORIGIN) inside the hexagon of circumradius `r_mm` centred at `centre_mm`, minus a keep-out around every feed."""
    n = int(np.ceil((r_mm + np.hypot(*centre_mm)) / t.P_MM)) + 2
    k = np.arange(-n, n + 1)
    xg, yg = np.meshgrid(t.P_MM * k, t.P_MM * k)
    x, y = xg.ravel(), yg.ravel()
    ok = t.inside(sa.hexagon_m(r_mm) * 1e3 + np.asarray(centre_mm), x, y)
    for px, py, _, _ in ports:
        ok &= np.hypot(x - px, y - py) >= sa.KEEPOUT_MM
    return x[ok], y[ok]


def march_analysis_tensor(
    table, lut, poly_m, feed_m, k0, th0, ph0, cell_xy_m, cell_slot, cell_ia, n_psi=180, drho=0.5e-3, rho_min=0.007
):
    """`shared_aperture_multiport.march_analysis` with the angle-dependent cell: at each step `alpha = |f|^2 / 2` and `beta = beta_table + dbeta` for the local slot, orientation and ray direction."""
    from .tensor_leaky_cp import modulation_k as mk

    psi, r_edge, rho = sa._grid(poly_m, feed_m, n_psi, drho)
    m = sa._empty_march(psi, rho, r_edge, drho)
    px = feed_m[0] + rho[None, :] * np.cos(psi)[:, None]
    py = feed_m[1] + rho[None, :] * np.sin(psi)[:, None]
    dist, idx = cKDTree(cell_xy_m).query(np.c_[px.ravel(), py.ravel()])
    ok = dist.reshape(px.shape) <= sa.CELL_REACH_MM * 1e-3
    node_slot = np.where(ok, np.asarray(cell_slot)[idx].reshape(px.shape), -1)
    node_ia = np.asarray(cell_ia)[idx].reshape(px.shape)
    inside = rho[None, :] <= r_edge[:, None]
    cur_w, cur_phi, cur_beta = np.ones(psi.size), np.zeros(psi.size), np.full(psi.size, sa.BETA0)
    for k, r in enumerate(rho):
        col = table.column(mk(cur_beta, psi, th0, ph0))
        s = node_slot[:, k]
        has = (s >= 0) & inside[:, k] & (r >= rho_min)
        sn = np.maximum(s, 0)
        fkk, fsk, db = lut.response(sn, node_ia[:, k], psi, col)
        a_used = np.where(has, 0.5 * (np.abs(fkk) ** 2 + np.abs(fsk) ** 2), 0.0)
        b_used = np.where(has, np.nan_to_num(table.beta[sn, col]) + db, sa.BETA0)
        m["alpha"][:, k], m["beta"][:, k], m["w"][:, k], m["phi_ref"][:, k], m["slot"][:, k], m["K"][:, k] = (
            a_used,
            b_used,
            cur_w,
            cur_phi,
            s,
            mk(cur_beta, psi, th0, ph0),
        )
        cur_w = cur_w * np.exp(-2.0 * a_used * drho)
        cur_phi = cur_phi + b_used * k0 * drho
        cur_beta = b_used
    return sa._finish(m, k0)


def cell_current(march, lut, table, feed_m, x_m, y_m, slot, ia, th0, ph0):
    """Per-cell current of one feed (n, 2) complex (global x / y) from the analysis march and the angle-dependent response."""
    s, _, _, psi_c, _ = sa.incident(march, feed_m, x_m, y_m)
    kh = np.stack([np.cos(psi_c), np.sin(psi_c)], axis=-1)
    sh = np.stack([-np.sin(psi_c), np.cos(psi_c)], axis=-1)
    amp = np.sum(s * kh, axis=-1)
    col = table.column(modulation_k(sa.BETA0, psi_c, th0, ph0))
    fkk, fsk, _ = lut.response(slot, ia, psi_c, col)
    return amp[:, None] * (fkk[:, None] * kh + fsk[:, None] * sh)


def design_tensor(
    table,
    lut,
    r_mm,
    ports,
    offsets_deg,
    centre_mm=(0.0, 0.0),
    n_iter=2,
    n_psi=240,
    defocus_rad=0.0,
    eta_design=0.9,
    taper_c=0.5,
):
    """Frey-style shared-aperture design on the angle-dependent model (module docstring). `ports`: `[(x_mm, y_mm, theta_deg, phi_deg)]`, `centre_mm`: hexagon centre. Returns the dict
    `shared_aperture_multiport.port_fields` / `scan_table` consume (`x_m, y_m, J, eta_conv, ports, beams, feeds_m, march`, plus `slot`, `axis_rad`, `slot_L_mm`, `slot_W_mm`, `kappa_des`)."""
    k0 = t.K0
    xm, ym = (a * 1e-3 for a in aperture_cells(r_mm, ports, centre_mm))
    poly = sa.hexagon_m(r_mm) + np.asarray(centre_mm) * 1e-3
    n_p = len(ports)
    feeds = [(p[0] * 1e-3, p[1] * 1e-3) for p in ports]
    beams = [(np.radians(p[2]), np.radians(p[3])) for p in ports]
    off = np.radians(np.asarray(offsets_deg, dtype=float))
    defoc = np.broadcast_to(np.asarray(defocus_rad, dtype=float), (n_p,))
    cxy = np.c_[xm, ym]
    ideal = [
        march_rays(table, poly, feeds[j], k0, beams[j][0], beams[j][1], eta_design, taper_c=taper_c, n_psi=n_psi)
        for j in range(n_p)
    ]
    inc = [sa.incident(ideal[j], feeds[j], xm, ym) for j in range(n_p)]
    kappa_des = np.array([np.sqrt(2.0 * i[1]) for i in inc])
    psi_c = np.array([i[3] for i in inc])
    s_amp = np.array([np.sqrt(np.sum(np.abs(i[0]) ** 2, axis=-1)) for i in inc])
    ph_cur = np.array([i[4] for i in inc])
    dxy = [(xm - f[0], ym - f[1]) for f in feeds]
    dk, ds = [], []
    for j in range(n_p):  # tables at each feed's local psi
        a, b, _ = lut.tables_at(psi_c[j])
        col = table.column(modulation_k(sa.BETA0, psi_c[j], beams[j][0], beams[j][1]))
        sc = lut.c[:, None] * lut.g_col[:, col]  # (slot, n)
        dk.append(a * sc[:, None, :])
        ds.append(b * sc[:, None, :])
    slot = ia = None
    ana = None
    for _ in range(n_iter + 1):
        chi = np.array(
            [
                ph_cur[j]
                - k0 * np.sin(beams[j][0]) * (dxy[j][0] * np.cos(beams[j][1]) + dxy[j][1] * np.sin(beams[j][1]))
                + off[j]
                + defoc[j] * (xm**2 + ym**2) / (r_mm * 1e-3) ** 2
                for j in range(n_p)
            ]
        )
        tk, ts_ = (
            kappa_des * np.cos(chi - psi_c),
            kappa_des * np.sin(chi - psi_c),
        )  # target components along k_hat, s_hat, (n_p, n)
        err = np.zeros((len(table.slots), N_DENSE, xm.size))
        for j in range(n_p):
            wgt = s_amp[j] ** 2
            err += wgt * (np.abs(dk[j] - tk[j]) ** 2 + np.abs(ds[j] - ts_[j]) ** 2)
        flat = err.reshape(-1, xm.size)
        best = np.argmin(flat, axis=0)
        slot, ia = best // N_DENSE, best % N_DENSE
        ana = [
            march_analysis_tensor(table, lut, poly, feeds[j], k0, beams[j][0], beams[j][1], cxy, slot, ia, n_psi)
            for j in range(n_p)
        ]
        inc2 = [sa.incident(ana[j], feeds[j], xm, ym) for j in range(n_p)]
        s_amp = np.array([np.sqrt(np.sum(np.abs(i[0]) ** 2, axis=-1)) for i in inc2])
        ph_cur = np.array([i[4] for i in inc2])
    cur = [cell_current(ana[j], lut, table, feeds[j], xm, ym, slot, ia, beams[j][0], beams[j][1]) for j in range(n_p)]
    return {
        "r_mm": r_mm,
        "ports": ports,
        "x_m": xm,
        "y_m": ym,
        "slot": slot,
        "slot_L_mm": np.array([table.slots[s][0] for s in slot]),
        "slot_W_mm": np.array([table.slots[s][1] for s in slot]),
        "axis_rad": ia * lut.step,
        "J": cur,
        "eta_conv": [a["eta_conv"] for a in ana],
        "march": ana,
        "feeds_m": feeds,
        "beams": beams,
        "kappa_des": kappa_des,
        "centre_mm": tuple(centre_mm),
        "mode": "tensor",
    }
