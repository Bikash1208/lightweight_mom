"""Efficiency improvements I1-I7 of the K = 2 coded holographic board at EQUAL OR LOWER diode / line count (MODEL_ONLY).
Model (same as sidequest_coded_pol.py): hexagon (circumradius 150 mm), centre probe, crossed two-diode rim cell, state A RHCP (25, 0), state B LHCP (25, 180). Each cell has two elements (axes x, y);
element value m (relative to the baseline +-1 swing of the diode cell) in each state; cell current  J = sum_e m_e (pi/2) amp (u_e . rho_hat) u_e e^{-j beta rho}; amp from the radial leakage profile
alpha_i = A0 * mean_k sum_e (m_e,k u_e.rho_hat)^2 (power in = power out), ohmic loss AA0 * (ON-diode weight), edge power from the same profile.
Design = per-cell least squares to the complex target tau_k = g amp e^{-j eps_k} (co-pol projection incl. aperture phase), over option sets: static values (no diode), switched pairs (one diode), +lambda per diode.
Option sets:  bit: static +-1, switched (+1,-1)/(-1,+1).  I1: static values on a grid in [-Ms, Ms].  I2: switched pairs (m_OFF, m_ON) from the MoM sub-library (x ls/pos variants).
Library values are mapped m = (Y_xx - Y_ref)/H_base with the baseline cell (ls 2.8, pos 1.6) defining Y_ref (mean of OFF, ON) and H_base (half swing)."""
import json
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np

import sidequest_arch_study as st

THETA = np.radians(25.0)
PHIS = np.radians([0.0, 180.0])
POLS = ("RHCP", "LHCP")
BETA_K0 = 1.07
A0 = 0.70  # Np/m, baseline cell swing +-1 (converged estimate 0.65-0.75)
AA0 = 0.095  # Np/m, ohmic absorption for the baseline cell with half its diodes ON
LAUNCH0 = 0.627
AA_BG = 0.0  # Np/m background (dielectric + conductor) loss of the stack; 0 in the original RT5880 study
BETA = BETA_K0 * st.K0
AREA = 58.457e-3  # m^2, hexagon of side 150 mm
LAM0 = 299792458.0 / 10.54e9
DMAX_DBI = 10 * np.log10(4 * np.pi * AREA / LAM0**2)
U = {"x": np.array([1.0, 0.0]), "y": np.array([0.0, 1.0])}


class Lib:
    """Sub-library from coded_lib.json (mesh 64)."""

    def __init__(self, path="results/sidequest_hybrid/coded_lib.json"):
        rows = json.loads(Path(path).read_text(encoding="utf-8"))
        self.rows = rows
        sw = {}
        for r in rows:
            if r["kind"] == "sw":
                sw.setdefault((r["ls"], r["pos"]), {})[r["mode"]] = complex(*r["Y_mS"][0])
        self.base = sw[(2.8, 1.6)]
        self.yref = 0.5 * (self.base["OFF"].imag + self.base["ON"].imag)
        self.h = 0.5 * (self.base["OFF"].imag - self.base["ON"].imag)
        self.gon_base = self.base["ON"].real
        self.sw = {k: {"m_off": (v["OFF"].imag - self.yref) / self.h, "m_on": (v["ON"].imag - self.yref) / self.h, "g_on": v["ON"].real / self.gon_base} for k, v in sw.items()}
        self.static = {(r["ls"], r["mode"]): (complex(*r["Y_mS"][0]).imag - self.yref) / self.h for r in rows if r["kind"] == "st"}


def hexgeom():
    poly, x, y, _ = st.aperture("hex")
    rho, psi = np.hypot(x, y), np.arctan2(y, x)
    redge = st.redge(poly, None, psi)
    return x, y, rho, psi, redge


class Board:
    def __init__(self):
        self.x, self.y, self.rho, self.psi, self.redge = hexgeom()
        self.n = self.x.size
        self.rhat = np.array([np.cos(self.psi), np.sin(self.psi)])
        self.az = ((self.psi + np.pi) / (2 * np.pi) * 12).astype(int) % 12  # 12 azimuth sectors
        self.arr = [np.exp(1j * st.K0 * np.sin(THETA) * (self.x * np.cos(p) + self.y * np.sin(p))) for p in PHIS]

    def profile(self, alpha_i, aa_i):
        """Cumulative loss integral I_i = int (alpha + alpha_abs) drho from RHO_KEEP to rho_i along the cell's azimuth sector (radial bins of 6.4 mm, sector means)."""
        edges = np.arange(st.RHO_KEEP, 0.16 + 1e-9, 6.4e-3)
        out = np.zeros(self.n)
        e_edge = np.zeros(self.n)
        for s in range(12):
            sel = np.flatnonzero(self.az == s)
            b = np.clip(np.digitize(self.rho[sel], edges) - 1, 0, len(edges) - 2)
            tot = np.zeros(len(edges) - 1)
            cnt = np.zeros(len(edges) - 1)
            np.add.at(tot, b, alpha_i[sel] + aa_i[sel])
            np.add.at(cnt, b, 1.0)
            mean = np.where(cnt > 0, tot / np.maximum(cnt, 1), np.nan)
            # fill empty bins from the nearest filled one
            idx = np.arange(len(mean))
            ok = ~np.isnan(mean)
            mean = np.interp(idx, idx[ok], mean[ok])
            cum = np.concatenate([[0.0], np.cumsum(mean * np.diff(edges))])
            out[sel] = np.interp(self.rho[sel], edges, cum)
            e_edge[sel] = np.interp(self.redge[sel], edges, cum)
        return out, e_edge


def unit_currents(bd, amp):
    """J_unit for element e in {x, y, xy}: (pi/2) amp v_e e^{-j beta rho}, v_x = (rho_x, 0), v_y = (0, rho_y), v_xy = (rho_y, rho_x) (= Delta Y rho_hat for Delta Y = x x, y y, x y + y x);
    returns {e: (jx, jy)} and co-pol contributions c[e][k]."""
    j, c = {}, {}
    ph = np.exp(-1j * BETA * bd.rho)
    rx, ry = bd.rhat[0], bd.rhat[1]
    vec = {"x": (rx, 0 * rx), "y": (0 * ry, ry), "xy": (ry, rx)}
    for e, (v0, v1) in vec.items():
        f = (np.pi / 2) * amp * ph
        j[e] = (f * v0, f * v1)
        c[e] = []
        for k in range(2):
            t = st.JONES[POLS[k]]
            et = np.cos(THETA) * (np.cos(PHIS[k]) * j[e][0] + np.sin(PHIS[k]) * j[e][1])
            ep = -np.sin(PHIS[k]) * j[e][0] + np.cos(PHIS[k]) * j[e][1]
            c[e].append((np.conj(t[0]) * et + np.conj(t[1]) * ep) * bd.arr[k])
    return j, c


def combos(vs, vw):
    """Cell option list: each option = (mx(2,), my(2,), n_diodes, w_x ,w_y flags) with vs = static values (array), vw = switched pairs (array of (mA, mB), with the ON state index)."""
    opts = []
    for a in vs:
        for b in vs:
            opts.append(((a, a), (b, b), 0, 0, 0))
    for p in vw:
        for b in vs:
            opts.append(((p[0], p[1]), (b, b), 1, p[2], 0))
            opts.append(((b, b), (p[0], p[1]), 1, 0, p[2]))
    for p in vw:
        for q in vw:
            opts.append(((p[0], p[1]), (q[0], q[1]), 2, p[2], q[2]))
    return opts


def option_arrays(vs, pairs, tens=None):
    """vs: static values; pairs: list of (mA, mB, onA, onB, gon) switched pairs; tens: static tensor options (mx, my, mxy) (E1-type cells, both elements static, no diode).
    Returns mx (n,2), my (n,2), mxy (n,2), nd (n,), loss (n,2,2), flag_tensor (n,)."""
    mx, my, mxy, nd, ld, ft = [], [], [], [], [], []

    def add(a, b, c3, n, l, f=False):
        mx.append(a); my.append(b); mxy.append(c3); nd.append(n); ld.append(l); ft.append(f)

    z2 = (0.0, 0.0)
    for a in vs:
        for b in vs:
            add((a, a), (b, b), z2, 0, (z2, z2))
    for p in pairs:
        for b in vs:
            add((p[0], p[1]), (b, b), z2, 1, ((p[2] * p[4], p[3] * p[4]), z2))
            add((b, b), (p[0], p[1]), z2, 1, (z2, (p[2] * p[4], p[3] * p[4])))
    for p in pairs:
        for q in pairs:
            add((p[0], p[1]), (q[0], q[1]), z2, 2, ((p[2] * p[4], p[3] * p[4]), (q[2] * q[4], q[3] * q[4])))
    for t in (tens or []):
        add((t[0], t[0]), (t[1], t[1]), (t[2], t[2]), 0, (z2, z2), True)
    return np.array(mx, float), np.array(my, float), np.array(mxy, float), np.array(nd, float), np.array(ld, float), np.array(ft, bool)


def e1_options(m_max, levels=(0.35, 0.7, 1.0), n_angles=16):
    """Static tensor options of the E1-type slotted cell: Delta Y = delta [cos2a (xx - yy) + sin2a (xy + yx)], delta in levels * m_max, 2a in n_angles directions, plus the empty option."""
    out = [(0.0, 0.0, 0.0)]
    for lv in levels:
        for a in np.arange(n_angles) * 2 * np.pi / n_angles:
            out.append((lv * m_max * np.cos(a), -lv * m_max * np.cos(a), lv * m_max * np.sin(a)))
    return out


def bit_sets():
    return np.array([-1.0, 1.0]), [(1.0, -1.0, 0, 1, 1.0), (-1.0, 1.0, 1, 0, 1.0)]  # pair (mA, mB, A uses ON, B uses ON, g_on): m=+1 is OFF, m=-1 is ON


def design(bd, amp, vs, pairs, eps, lam, mu=0.0, g=1.0, w=(1.0, 1.0), scale=None, ring=None, g_ring=1.0, tens=None, e1_scale=None):
    """Per-cell least squares over the option product. Returns m (2 states, 3 elements x/y/xy, n), nd (n,), loss (2 states, n).
    scale (n,) multiplies all values (taper); e1_scale (n,) multiplies only the E1 tensor options."""
    _, c = unit_currents(bd, amp)
    mx, my, mxy, nd, ld, ft = option_arrays(vs, pairs, tens)
    sc = np.ones(bd.n) if scale is None else scale
    sce = np.ones(bd.n) if e1_scale is None else e1_scale
    fac = np.where(ft[None, :], (sc * sce)[:, None], sc[:, None])  # (n, nopt)
    cost = np.zeros((bd.n, mx.shape[0]))
    for k in range(2):
        gg = g if ring is None else np.where(ring, g_ring, g)
        tau = gg * amp * np.exp(-1j * eps[k])
        pred = (c["x"][k][:, None] * mx[None, :, k] + c["y"][k][:, None] * my[None, :, k] + c["xy"][k][:, None] * mxy[None, :, k]) * fac
        cost += w[k] * np.abs(tau[:, None] - pred) ** 2
    scale_c = np.mean(np.abs(g * amp) ** 2)
    cost = cost / scale_c + lam * nd[None, :]
    if ring is not None:
        cost = cost + 1e6 * (nd[None, :] > 0) * ring[:, None]  # edge-ring cells are static (no diode)
    wx = np.array([bd.rhat[0] ** 2, bd.rhat[1] ** 2])  # (u_e . rhat)^2
    if mu:
        lossc = np.zeros((bd.n, mx.shape[0]))
        for k in range(2):
            lossc += (ld[None, :, 0, k] * wx[0][:, None] + ld[None, :, 1, k] * wx[1][:, None])
        cost = cost + mu * lossc
    best = np.argmin(cost, axis=1)
    f_best = fac[np.arange(bd.n), best]
    m = np.zeros((2, 3, bd.n))
    for k in range(2):
        m[k, 0] = mx[best, k] * f_best
        m[k, 1] = my[best, k] * f_best
        m[k, 2] = mxy[best, k] * f_best
    loss = np.zeros((2, bd.n))
    for k in range(2):
        loss[k] = ld[best, 0, k] * wx[0] + ld[best, 1, k] * wx[1]
    return m, nd[best], loss


def realise(bd, amp, m):
    j, _ = unit_currents(bd, amp)
    cur = []
    for k in range(2):
        cur.append((sum(m[k, i] * j[e][0] for i, e in enumerate(("x", "y", "xy"))), sum(m[k, i] * j[e][1] for i, e in enumerate(("x", "y", "xy")))))
    return cur


def _vec2(bd, m):
    """|Delta Y rho_hat|^2 per state: (m_x rho_x + m_xy rho_y)^2 + (m_xy rho_x + m_y rho_y)^2."""
    rx, ry = bd.rhat
    return [(m[k, 0] * rx + m[k, 2] * ry) ** 2 + (m[k, 2] * rx + m[k, 1] * ry) ** 2 for k in range(2)]


def amplitude_profile(bd, m, loss):
    """alpha_i from the actual values, loss integral -> amp_i (radiated amplitude per unit m) and the converted fraction / edge power."""
    alpha = A0 * np.mean(_vec2(bd, m), axis=0)
    aa = AA_BG + 2 * AA0 * np.mean(loss, axis=0)
    integ, ie = bd.profile(alpha, aa)
    amp = np.sqrt(2 * A0 / bd.rho) * np.exp(-integ)
    w = 2 * alpha * np.exp(-2 * integ) / bd.rho * st.PITCH**2 / (2 * np.pi)
    return amp, float(w.sum()), float(np.mean(np.exp(-2 * ie))), alpha, aa


def evaluate(bd, amp, m, launch, conv):
    cur = realise(bd, amp, m)
    out = []
    a_id = np.sqrt(np.mean(_vec2(bd, m), axis=0)) * amp
    for k in range(2):
        et, ep = st.far_field(bd.x, bd.y, *cur[k])
        t = st.JONES[POLS[k]]
        j = st.jones_to_current(t, THETA, PHIS[k])
        idl = a_id * np.exp(-1j * st.K0 * np.sin(THETA) * (bd.x * np.cos(PHIS[k]) + bd.y * np.sin(PHIS[k])))
        ei, epi = st.far_field(bd.x, bd.y, j[0] * idl, j[1] * idl)
        mi = st.metrics(ei, epi, THETA, PHIS[k], POLS[k])
        mb = st.metrics(et, ep, THETA, PHIS[k], POLS[k])
        i0 = mi["_i"]
        ratio = (np.abs(et) ** 2 + np.abs(ep) ** 2)[i0] / (np.abs(ei) ** 2 + np.abs(epi) ** 2)[i0]
        gain = mi["D_dbi"] + 10 * np.log10(launch * conv) + 10 * np.log10(ratio)
        out.append({"gain_dbi": float(gain), "SLL_db": mb["SLL_db"], "AR_db": mb["AR_db"], "XPD_db": mb["XPD_db"], "peak": (mb["peak_theta_deg"], mb["peak_phi_deg"])})
    return out


def phase_resolution(bd, amp, m, eps, g=1.0):
    """rms phase error [deg] of the realised co-pol contribution against the target phase, amplitude weighted; effective uniform step sqrt(12) * rms."""
    _, c = unit_currents(bd, amp)
    e2, wsum = 0.0, 0.0
    for k in range(2):
        pred = m[k, 0] * c["x"][k] + m[k, 1] * c["y"][k] + m[k, 2] * c["xy"][k]
        tau = g * amp * np.exp(-1j * eps[k])
        d = np.angle(pred * np.conj(tau))
        wgt = np.abs(tau) ** 2
        e2 += np.sum(wgt * d**2)
        wsum += np.sum(wgt)
    rms = np.degrees(np.sqrt(e2 / wsum))
    return float(rms), float(np.sqrt(12) * rms)


def tune_lambda(bd, amp, vs, pairs, eps, target=766, **kw):  # kw: mu, g, w, scale, ring, g_ring, tens, e1_scale
    """Smallest lambda with diodes <= target (design only, amplitude profile frozen)."""
    lo, hi = 0.0, 3.0
    for _ in range(9):
        mid = 0.5 * (lo + hi)
        _, nd, _ = design(bd, amp, vs, pairs, eps, mid, **kw)
        if nd.sum() > target:
            lo = mid
        else:
            hi = mid
    return hi


def full_run(vs, pairs, eps, mu=0.0, g=1.0, launch=LAUNCH0, w=(1.0, 1.0), scale=None, target=766, bd=None, ring=None, g_ring=1.0, tens=None, e1_scale=None):
    """lambda tuned so that diodes <= target, 3 amplitude iterations."""
    bd = bd or Board()
    amp = np.sqrt(2 * A0 / bd.rho) * np.exp(-0.5 * (A0 + AA0) * (bd.rho - st.RHO_KEEP))
    kw = dict(mu=mu, g=g, w=w, scale=scale, ring=ring, g_ring=g_ring, tens=tens, e1_scale=e1_scale)
    for _ in range(3):
        lam = tune_lambda(bd, amp, vs, pairs, eps, target, **kw)
        m, nd, loss = design(bd, amp, vs, pairs, eps, lam, **kw)
        amp, conv, edge, alpha, aa = amplitude_profile(bd, m, loss)
    lam = tune_lambda(bd, amp, vs, pairs, eps, target, **kw)
    m, nd, loss = design(bd, amp, vs, pairs, eps, lam, **kw)
    amp2, conv, edge, alpha, aa = amplitude_profile(bd, m, loss)
    res = evaluate(bd, amp2, m, launch, conv)
    rms, step = phase_resolution(bd, amp2, m, eps, g)
    nsw = int((nd > 0).sum())
    words = {(1, 0)} if nsw else set()  # every switched diode is driven (ON,OFF) or (OFF,ON): one complement class per axis
    dmax = DMAX_DBI
    out = {"states": res, "lambda": float(lam), "static_frac": float(1 - nd.sum() / (2 * bd.n)), "diodes": int(nd.sum()), "lines_with_inverters": len(words), "SW_to_beam": conv, "edge": edge,
           "ap_eff": [10 ** ((r["gain_dbi"] - dmax) / 10) for r in res], "rms_phase_deg": rms, "eff_step_deg": step, "mean_abs_m": float(np.mean(np.abs(m))), "ON_loss": float(np.mean(loss) * 2),
           "alpha_mean": float(alpha.mean()), "aa_mean": float(aa.mean()), "e1_cell_frac": float(np.mean(np.abs(m[0, 2]) > 1e-9))}
    return out, (bd, amp2, m, nd, loss, alpha, aa)


def taper_scale(bd, alpha, aa, q_edge=0.15):
    """I3: per azimuth sector, the depth scale s(rho) <= 1 for (near) uniform radiated power density with edge power fraction q_edge, if the cells allow it (alpha_req <= alpha_cap = current mean alpha).
    alpha_req = C rho / (2 Q),  dQ/drho = -C rho - 2 aa Q,  Q(rho_keep) = 1. C = min(C for Q(R) = q_edge, largest C with alpha_req <= alpha_cap everywhere)."""
    scale = np.ones(bd.n)
    info = []
    for s in range(12):
        sel = np.flatnonzero(bd.az == s)
        rr = np.linspace(st.RHO_KEEP, float(bd.redge[sel].mean()), 120)
        edges = np.arange(st.RHO_KEEP, 0.16 + 1e-9, 6.4e-3)
        b = np.clip(np.digitize(bd.rho[sel], edges) - 1, 0, len(edges) - 2)
        tot, cnt, tota = np.zeros(len(edges) - 1), np.zeros(len(edges) - 1), np.zeros(len(edges) - 1)
        np.add.at(tot, b, alpha[sel]); np.add.at(cnt, b, 1.0); np.add.at(tota, b, aa[sel])
        ok = cnt > 0
        mid = 0.5 * (edges[1:] + edges[:-1])
        acap = np.interp(rr, mid[ok], (tot / np.maximum(cnt, 1))[ok])
        aav = np.interp(rr, mid[ok], (tota / np.maximum(cnt, 1))[ok])

        def solve(cc):
            q = np.ones(rr.size)
            for i in range(1, rr.size):
                dr = rr[i] - rr[i - 1]
                q[i] = max(q[i - 1] + dr * (-cc * rr[i - 1] - 2 * aav[i - 1] * q[i - 1]), 1e-6)
            return q, cc * rr / (2 * q)

        lo, hi = 0.0, 1e4
        for _ in range(40):  # largest C with alpha_req <= alpha_cap
            mid_c = 0.5 * (lo + hi)
            q, areq = solve(mid_c)
            if np.all(areq <= acap) and q[-1] > 1e-5:
                lo = mid_c
            else:
                hi = mid_c
        c_feas = lo
        lo, hi = 0.0, 1e4
        for _ in range(40):  # C giving Q(R) = q_edge
            mid_c = 0.5 * (lo + hi)
            q, _ = solve(mid_c)
            if q[-1] > q_edge:
                lo = mid_c
            else:
                hi = mid_c
        c_tgt = hi
        cc = min(c_feas, c_tgt)
        q, areq = solve(cc)
        sv = np.sqrt(np.clip(areq / acap, 0.0, 1.0))
        scale[sel] = np.interp(bd.rho[sel], rr, sv)
        info.append({"sector": s, "C_feasible": c_feas, "C_target": c_tgt, "edge_pred": float(q[-1]), "saturated": bool(c_feas < c_tgt)})
    return scale, info
