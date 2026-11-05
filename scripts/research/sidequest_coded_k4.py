"""K-state coded holographic board (K = 2 and 4), RT5880 chosen cell (rim disc, diode at 0.8 mm), MODEL_ONLY. Same physics as sidequest_coded_improve (per-cell least squares over option sets, power bookkeeping,
scalar element factor) generalised to K states: per element a switched option is a K-bit ON-pattern ("word", both complements included) of a diode cell with values (m_OFF, m_ON); static options are levels.
Unit: m = +1 / -1 = OFF / ON of the chosen cell (alpha_beam A0 for the 50/50 pattern). Library pairs/levels are the round-11 MoM rows renormalised to the (ls 2.8, pos 0.8) pair (normal-incidence-extracted shapes: flagged).
Lines per axis: distinct ON-patterns of the switched elements of that axis (without inverters) / distinct complement classes (with inverters)."""
import itertools
import json
import sys
import time
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np

import sidequest_arch_study as st
import sidequest_coded_improve as I

THETA = np.radians(25.0)
A0, AA0, AA_BG, BETA_K0, LAUNCH = 1.84, 0.094, 0.15, 1.079, 0.627
BETA = BETA_K0 * st.K0
LAM0 = 299792458.0 / 10.54e9
DMAX = 10 * np.log10(4 * np.pi * I.AREA / LAM0**2)
SETS = {
    "K2": ([0, 180], ["RHCP", "LHCP"]),
    "K4a": ([0, 90, 180, 270], ["RHCP", "LHCP", "LP-H", "LP-V"]),
    "K4b": ([0, 90, 180, 270], ["RHCP", "LHCP", "RHCP", "LHCP"]),
}


class KBoard:
    def __init__(self, phis_deg, pols):
        self.phis = np.radians(phis_deg)
        self.pols = pols
        self.K = len(pols)
        self.x, self.y, self.rho, self.psi, self.redge = I.hexgeom()
        self.n = self.x.size
        self.rhat = np.array([np.cos(self.psi), np.sin(self.psi)])
        self.az = ((self.psi + np.pi) / (2 * np.pi) * 12).astype(int) % 12
        self.arr = [np.exp(1j * st.K0 * np.sin(THETA) * (self.x * np.cos(p) + self.y * np.sin(p))) for p in self.phis]
        self.profile = I.Board.profile.__get__(self)

    def unit(self, amp):
        ph = np.exp(-1j * BETA * self.rho)
        rx, ry = self.rhat
        vec = {"x": (rx, 0 * rx), "y": (0 * ry, ry)}
        j, c = {}, {}
        for e, (v0, v1) in vec.items():
            f = (np.pi / 2) * amp * ph
            j[e] = (f * v0, f * v1)
            c[e] = []
            for k in range(self.K):
                t = st.JONES[self.pols[k]]
                et = np.cos(THETA) * (np.cos(self.phis[k]) * j[e][0] + np.sin(self.phis[k]) * j[e][1])
                ep = -np.sin(self.phis[k]) * j[e][0] + np.cos(self.phis[k]) * j[e][1]
                c[e].append((np.conj(t[0]) * et + np.conj(t[1]) * ep) * self.arr[k])
        return j, c


def library():
    """(statics, pairs): round-11 MoM library renormalised so that the (2.8, 0.8) pair is (+1, -1)."""
    lib = I.Lib()
    v0 = lib.sw[(2.8, 0.8)]
    c0, h0 = 0.5 * (v0["m_off"] + v0["m_on"]), 0.5 * (v0["m_off"] - v0["m_on"])
    norm = lambda m: (m - c0) / h0  # noqa: E731
    stat = np.array(sorted(norm(v) for v in lib.static.values()))
    pairs = {k: (norm(v["m_off"]), norm(v["m_on"]), v["g_on"] / v0["g_on"]) for k, v in lib.sw.items()}
    return stat, pairs


def words(K, allw=False):
    return [w for w in itertools.product((0, 1), repeat=K) if allw or 0 < sum(w) < K]


def canon(w):
    return min(tuple(w), tuple(1 - np.array(w)))


def options(K, vs, geoms, force_switched=False, allw=False, cls_x=None, cls_y=None):
    """Option arrays: mx, my (nopt, K); nd; loss (nopt, 2, K) ON-weight (g_on * word); wx, wy: word ids (or -1) for line counting."""
    W = words(K, allw)
    sw = []  # (values (K,), on (K,), id)
    for gi, (moff, mon, gon) in enumerate(geoms):
        for wi, w in enumerate(W):
            sw.append((np.array([mon if b else moff for b in w]), np.array(w) * gon, wi))
    swx = sw if cls_x is None else [t for t in sw if canon(W[t[2]]) in cls_x]
    swy = sw if cls_y is None else [t for t in sw if canon(W[t[2]]) in cls_y]
    mx, my, nd, ld, wx, wy = [], [], [], [], [], []
    z = np.zeros(K)

    def add(a, b, n, la, lb, ia, ib):
        mx.append(a); my.append(b); nd.append(n); ld.append((la, lb)); wx.append(ia); wy.append(ib)

    if not force_switched:
        for a in vs:
            for b in vs:
                add(np.full(K, a), np.full(K, b), 0, z, z, -1, -1)
        for s_ in swx:
            for b in vs:
                add(s_[0], np.full(K, b), 1, s_[1], z, s_[2], -1)
        for s_ in swy:
            for b in vs:
                add(np.full(K, b), s_[0], 1, z, s_[1], -1, s_[2])
    for s_ in swx:
        for q in swy:
            add(s_[0], q[0], 2, s_[1], q[1], s_[2], q[2])
    return (np.array(mx), np.array(my), np.array(nd, float), np.array(ld), np.array(wx), np.array(wy), W)


def design(bd, amp, opt, eps, lam, mu=0.0, g=1.0, w=None):
    mx, my, nd, ld, wxi, wyi, _W = opt
    K = bd.K
    w = np.ones(K) if w is None else w
    _, c = bd.unit(amp)
    cost = np.zeros((bd.n, mx.shape[0]))
    for k in range(K):
        tau = g * amp * np.exp(-1j * eps[k])
        pred = c["x"][k][:, None] * mx[None, :, k] + c["y"][k][:, None] * my[None, :, k]
        cost += w[k] * np.abs(tau[:, None] - pred) ** 2
    cost = cost / np.mean(np.abs(g * amp) ** 2) + lam * nd[None, :]
    wx = np.array([bd.rhat[0] ** 2, bd.rhat[1] ** 2])
    if mu:
        lossc = np.zeros((bd.n, mx.shape[0]))
        for k in range(K):
            lossc += ld[None, :, 0, k] * wx[0][:, None] + ld[None, :, 1, k] * wx[1][:, None]
        cost = cost + mu * lossc
    best = np.argmin(cost, axis=1)
    m = np.stack([np.stack([mx[best, k], my[best, k]]) for k in range(K)])  # (K, 2, n)
    loss = np.stack([ld[best, 0, k] * wx[0] + ld[best, 1, k] * wx[1] for k in range(K)])
    return m, nd[best], loss, (wxi[best], wyi[best])


def vec2(bd, m):
    rx, ry = bd.rhat
    return [(m[k, 0] * rx) ** 2 + (m[k, 1] * ry) ** 2 for k in range(bd.K)]


def amp_profile(bd, m, loss):
    alpha = A0 * np.mean(vec2(bd, m), axis=0)
    aa = AA_BG + 2 * AA0 * np.mean(loss, axis=0)
    integ, ie = bd.profile(alpha, aa)
    amp = np.sqrt(2 * A0 / bd.rho) * np.exp(-integ)
    w = 2 * alpha * np.exp(-2 * integ) / bd.rho * st.PITCH**2 / (2 * np.pi)
    return amp, float(w.sum()), float(np.mean(np.exp(-2 * ie)))


def evaluate(bd, amp, m, conv, launch=LAUNCH):
    j, _ = bd.unit(amp)
    a_id = np.sqrt(np.mean(vec2(bd, m), axis=0)) * amp
    out = []
    for k in range(bd.K):
        cur = (m[k, 0] * j["x"][0] + m[k, 1] * j["y"][0], m[k, 0] * j["x"][1] + m[k, 1] * j["y"][1])
        et, ep = st.far_field(bd.x, bd.y, *cur)
        jj = st.jones_to_current(st.JONES[bd.pols[k]], THETA, bd.phis[k])
        idl = a_id * np.exp(-1j * st.K0 * np.sin(THETA) * (bd.x * np.cos(bd.phis[k]) + bd.y * np.sin(bd.phis[k])))
        ei, epi = st.far_field(bd.x, bd.y, jj[0] * idl, jj[1] * idl)
        mi = st.metrics(ei, epi, THETA, bd.phis[k], bd.pols[k])
        mb = st.metrics(et, ep, THETA, bd.phis[k], bd.pols[k])
        i0 = mi["_i"]
        ratio = (np.abs(et) ** 2 + np.abs(ep) ** 2)[i0] / (np.abs(ei) ** 2 + np.abs(epi) ** 2)[i0]
        gain = mi["D_dbi"] + 10 * np.log10(launch * conv) + 10 * np.log10(ratio)
        out.append({"gain": float(gain), "SLL": mb["SLL_db"], "AR": mb["AR_db"], "XPD": mb["XPD_db"], "pk": (mb["peak_theta_deg"], mb["peak_phi_deg"]), "ap": 10 ** ((gain - DMAX) / 10)})
    return out


def lines(bd, ids, nd, W):
    res = {}
    for ax, wid in zip("xy", ids, strict=True):
        used = sorted({int(i) for i in wid if i >= 0})
        classes = {min(tuple(W[i]), tuple(1 - np.array(W[i]))) for i in used}
        res[ax] = {"words": len(used), "classes": len(classes)}
    return res


def tune(bd, amp, opt, eps, target, **kw):
    lo, hi = 0.0, 4.0
    for _ in range(9):
        mid = 0.5 * (lo + hi)
        if design(bd, amp, opt, eps, mid, **kw)[1].sum() > target:
            lo = mid
        else:
            hi = mid
    return hi


def run(bd, opt, eps, lam=None, target=None, mu=0.0, g=1.0, iters=3):
    amp = np.sqrt(2 * A0 / bd.rho) * np.exp(-0.5 * (A0 + AA0) * (bd.rho - st.RHO_KEEP))
    for _ in range(iters):
        l = lam if lam is not None else tune(bd, amp, opt, eps, target, mu=mu, g=g)
        m, nd, loss, ids = design(bd, amp, opt, eps, l, mu, g)
        amp, conv, edge = amp_profile(bd, m, loss)
    res = evaluate(bd, amp, m, conv)
    nl = lines(bd, ids, nd, opt[6])
    n_el = 2 * bd.n
    return {"states": res, "lam": float(l), "static_pct": 100 * (1 - nd.sum() / n_el), "diodes": int(nd.sum()), "lines": nl, "SW_to_beam": conv, "edge": edge,
            "ON_loss": float(np.mean(loss) * 2), "g": g}, (m, nd)


def summary(r):
    s = r["states"]
    return {"gain": [round(x["gain"], 2) for x in s], "SLL": [round(x["SLL"], 1) for x in s], "AR": [round(x["AR"], 2) for x in s], "XPD": [round(x["XPD"], 1) for x in s],
            "ap_eff_pct": [round(100 * x["ap"], 1) for x in s], "edge_pct": round(100 * r["edge"], 1), "SW_to_beam_pct": round(100 * r["SW_to_beam"], 1), "static_pct": round(r["static_pct"], 1),
            "diodes": r["diodes"], "diodes_per_state": round(r["diodes"] / len(s), 1), "lines_x_y": [r["lines"]["x"]["classes"], r["lines"]["y"]["classes"]],
            "lines_noinv_x_y": [r["lines"]["x"]["words"], r["lines"]["y"]["words"]], "g": r["g"], "lam": round(r["lam"], 3)}
