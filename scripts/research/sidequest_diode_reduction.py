"""Fewer diodes: ideas D1-D6 on the coded board (K = 2 and K = 4 CP set B), RT5880 chosen cell, MODEL_ONLY. Engine = sidequest_coded_k4 (per-cell least squares over option sets) extended with
  D1  power-ranked sparsification (switched elements ranked by amp^2 x code disagreement, the lowest forced static)
  D2  switched superpixels (1x2 / 2x2 groups of lattice cells share ONE option: one diode pair), group chosen only where cheaper than independent cells
  D3  one-axis switching (no cell with both elements switched; 'strict' = only the x element is ever switched)
  D5  state overlap: per-state pointing offsets (+-2 deg in theta, phi) and reference phases free (coordinate descent on the design objective)
  D6  reconfigurable core radius + static rim
  D4  one-diode rotated cell (separate module section: rotated disc, 8 axis angles).
A spec is a dict; execute(spec) returns the summary row used by the Pareto study."""
import itertools
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_arch_study as st
import sidequest_coded_k4 as K

PITCH = st.PITCH
THETA = K.THETA
STAT, PAIRS = K.library()
VS_AN = np.linspace(STAT.min(), STAT.max(), 9)
G2 = [PAIRS[(2.8, 0.8)], PAIRS[(2.6, 0.8)]]
CLS3 = {(0, 0, 1, 1), (0, 1, 0, 1), (0, 1, 1, 0)}
EPS_B = {"K2": [0.0, 0.0], "K4b": [0.0, 0.0, 0.0, np.pi]}
LAM_BASE = {"K2": 0.266, "K4b": 0.555}
_BD = {}


class GBoard(K.KBoard):
    """KBoard with per-state theta offsets and lattice indices."""

    def __init__(self, name, dphi=None, dtheta=None):
        ph, pol = K.SETS[name]
        dphi = np.zeros(len(ph)) if dphi is None else np.asarray(dphi, float)
        dth = np.zeros(len(ph)) if dtheta is None else np.asarray(dtheta, float)
        super().__init__(list(np.asarray(ph, float) + dphi), pol)
        self.thetas = THETA + np.radians(dth)
        self.arr = [np.exp(1j * st.K0 * np.sin(self.thetas[k]) * (self.x * np.cos(self.phis[k]) + self.y * np.sin(self.phis[k]))) for k in range(self.K)]
        self.ix = np.round(self.x / PITCH).astype(int)
        self.iy = np.round(self.y / PITCH).astype(int)

    def unit(self, amp):
        ph = np.exp(-1j * K.BETA * self.rho)
        rx, ry = self.rhat
        vec = {"x": (rx, 0 * rx), "y": (0 * ry, ry)}
        j, c = {}, {}
        for e, (v0, v1) in vec.items():
            f = (np.pi / 2) * amp * ph
            j[e] = (f * v0, f * v1)
            c[e] = []
            for k in range(self.K):
                t = st.JONES[self.pols[k]]
                et = np.cos(self.thetas[k]) * (np.cos(self.phis[k]) * j[e][0] + np.sin(self.phis[k]) * j[e][1])
                ep = -np.sin(self.phis[k]) * j[e][0] + np.cos(self.phis[k]) * j[e][1]
                c[e].append((np.conj(t[0]) * et + np.conj(t[1]) * ep) * self.arr[k])
        return j, c


def board(name, dphi=None, dtheta=None):
    key = (name, None if dphi is None else tuple(np.round(dphi, 3)), None if dtheta is None else tuple(np.round(dtheta, 3)))
    if key not in _BD:
        _BD[key] = GBoard(name, dphi, dtheta)
    return _BD[key]


def group_ids(bd, mode):
    if mode is None:
        return None
    if mode == "1x2x":
        key = (bd.ix // 2, bd.iy)
    elif mode == "1x2y":
        key = (bd.ix, bd.iy // 2)
    elif mode == "2x2":
        key = (bd.ix // 2, bd.iy // 2)
    else:
        raise ValueError(mode)
    _, gid = np.unique(np.stack(key, axis=1), axis=0, return_inverse=True)
    return gid.ravel()


def design_g(bd, amp, opt, eps, lam, mu=0.0, g=1.0, gid=None, forbid_sw=None, forbid_x=None, forbid_y=None, no_both=False, only_x=False):
    mx, my, nd, ld, wxi, wyi, _W = opt
    Kn = bd.K
    _, c = bd.unit(amp)
    err = np.zeros((bd.n, mx.shape[0]))
    for k in range(Kn):
        tau = g * amp * np.exp(-1j * eps[k])
        pred = c["x"][k][:, None] * mx[None, :, k] + c["y"][k][:, None] * my[None, :, k]
        err += np.abs(tau[:, None] - pred) ** 2
    err /= np.mean(np.abs(g * amp) ** 2)
    wx = np.array([bd.rhat[0] ** 2, bd.rhat[1] ** 2])
    if mu:
        for k in range(Kn):
            err += mu * (ld[None, :, 0, k] * wx[0][:, None] + ld[None, :, 1, k] * wx[1][:, None])
    bad = np.zeros_like(err, bool)
    swx, swy = (wxi >= 0), (wyi >= 0)
    if forbid_sw is not None:
        bad |= forbid_sw[:, None] & (nd > 0)[None, :]
    if forbid_x is not None:
        bad |= forbid_x[:, None] & swx[None, :]
    if forbid_y is not None:
        bad |= forbid_y[:, None] & swy[None, :]
    if no_both:
        bad |= (nd == 2)[None, :]
    if only_x:
        bad |= swy[None, :]
    cost = np.where(bad, np.inf, err + lam * nd[None, :])
    best = np.argmin(cost, axis=1)
    ndcell = nd[best].copy()
    jtot = float(cost[np.arange(bd.n), best].sum())
    if gid is not None:
        ng = gid.max() + 1
        ind = np.zeros(ng)
        np.add.at(ind, gid, cost[np.arange(bd.n), best])
        s_err = np.zeros((ng, mx.shape[0]))
        np.add.at(s_err, gid, np.where(bad, np.inf, err))
        sh = s_err + lam * nd[None, :]
        sbest = np.argmin(sh, axis=1)
        shc = sh[np.arange(ng), sbest]
        use = shc < ind - 1e-12
        jtot = float(np.where(use, shc, ind).sum())
        cellshared = use[gid]
        best = np.where(cellshared, sbest[gid], best)
        first = np.zeros(bd.n, bool)
        _, idx0 = np.unique(gid, return_index=True)
        first[idx0] = True
        ndcell = np.where(cellshared, np.where(first, nd[best], 0.0), nd[best])
    m = np.stack([np.stack([mx[best, k], my[best, k]]) for k in range(Kn)])
    loss = np.stack([ld[best, 0, k] * wx[0] + ld[best, 1, k] * wx[1] for k in range(Kn)])
    return m, ndcell, loss, (wxi[best], wyi[best]), jtot, best


def evaluate_g(bd, amp, m, conv, launch=K.LAUNCH):
    j, _ = bd.unit(amp)
    a_id = np.sqrt(np.mean(K.vec2(bd, m), axis=0)) * amp
    out = []
    for k in range(bd.K):
        cur = (m[k, 0] * j["x"][0] + m[k, 1] * j["y"][0], m[k, 0] * j["x"][1] + m[k, 1] * j["y"][1])
        et, ep = st.far_field(bd.x, bd.y, *cur)
        jj = st.jones_to_current(st.JONES[bd.pols[k]], bd.thetas[k], bd.phis[k])
        idl = a_id * np.exp(-1j * st.K0 * np.sin(bd.thetas[k]) * (bd.x * np.cos(bd.phis[k]) + bd.y * np.sin(bd.phis[k])))
        ei, epi = st.far_field(bd.x, bd.y, jj[0] * idl, jj[1] * idl)
        mi = st.metrics(ei, epi, bd.thetas[k], bd.phis[k], bd.pols[k])
        mb = st.metrics(et, ep, bd.thetas[k], bd.phis[k], bd.pols[k])
        i0 = mi["_i"]
        ratio = (np.abs(et) ** 2 + np.abs(ep) ** 2)[i0] / (np.abs(ei) ** 2 + np.abs(epi) ** 2)[i0]
        gain = mi["D_dbi"] + 10 * np.log10(launch * conv) + 10 * np.log10(ratio)
        out.append({"gain": float(gain), "SLL": mb["SLL_db"], "XPD": mb["XPD_db"], "ap": 10 ** ((gain - K.DMAX) / 10), "pk": (mb["peak_theta_deg"], mb["peak_phi_deg"])})
    return out


def run_g(name, opt, eps, lam, mu=0.05, g=2.0, iters=3, bd=None, **kw):
    bd = bd or board(name)
    amp = np.sqrt(2 * K.A0 / bd.rho) * np.exp(-0.5 * (K.A0 + K.AA0) * (bd.rho - st.RHO_KEEP))
    for _ in range(iters):
        m, ndc, loss, ids, J, best = design_g(bd, amp, opt, eps, lam, mu, g, **kw)
        amp, conv, edge = K.amp_profile(bd, m, loss)
    res = evaluate_g(bd, amp, m, conv)
    nl = K.lines(bd, ids, ndc, opt[6])
    return {"states": res, "diodes": int(ndc.sum()), "static_pct": 100 * (1 - ndc.sum() / (2 * bd.n)), "lines": nl, "SW_to_beam": conv, "edge": edge, "J": J}, (m, ndc, ids, amp, best)


def row(r, spec):
    s = r["states"]
    return {"spec": spec, "K": len(s), "gain": [round(x["gain"], 2) for x in s], "gain_mean": round(float(np.mean([x["gain"] for x in s])), 3), "SLL_worst": round(max(x["SLL"] for x in s), 1),
            "XPD_min": round(min(x["XPD"] for x in s), 1), "ap_eff_pct": round(100 * float(np.mean([x["ap"] for x in s])), 2), "diodes": r["diodes"], "static_pct": round(r["static_pct"], 1),
            "lines_x": r["lines"]["x"]["classes"], "lines_y": r["lines"]["y"]["classes"], "lines_noinv_x": r["lines"]["x"]["words"], "lines_noinv_y": r["lines"]["y"]["words"],
            "edge_pct": round(100 * r["edge"], 1), "SW_to_beam_pct": round(100 * r["SW_to_beam"], 1)}


def d5_search(name, lam, base_opt, eps0, passes=2):
    """Coordinate descent on the design objective J (LS error + lambda * diodes) over per-state (dphi, dtheta) in {-2, 0, 2} deg and eps in 60-degree steps (state 0: eps fixed)."""
    Kn = len(K.SETS[name][0])
    dphi, dth, eps = np.zeros(Kn), np.zeros(Kn), np.array(eps0, float)
    bd0 = board(name)
    amp = np.sqrt(2 * K.A0 / bd0.rho) * np.exp(-0.5 * (K.A0 + K.AA0) * (bd0.rho - st.RHO_KEEP))

    def cost(dp, dt, ep):
        bd = board(name, dp, dt)
        return design_g(bd, amp, base_opt, ep, lam, 0.05, 2.0)[4]

    cur = cost(dphi, dth, eps)
    for _ in range(passes):
        for k in range(Kn):
            best = (cur, dphi[k], dth[k], eps[k])
            for a, b, e in itertools.product((-2.0, 0.0, 2.0), (-2.0, 0.0, 2.0), np.radians(np.arange(0, 360, 60)) if k else [eps[k]]):
                dp, dt, ep = dphi.copy(), dth.copy(), eps.copy()
                dp[k], dt[k], ep[k] = a, b, e
                c = cost(dp, dt, ep)
                if c < best[0] - 1e-9:
                    best = (c, a, b, e)
            cur, dphi[k], dth[k], eps[k] = best
    return dphi, dth, eps, cur


def execute(spec):
    name = spec["K"]
    lam = spec.get("lam", LAM_BASE[name])
    cls = CLS3 if (name == "K4b" and spec.get("cls", True)) else None
    opt = K.options(len(K.SETS[name][0]), VS_AN, G2, cls_x=cls, cls_y=cls)
    eps = np.array(spec.get("eps", EPS_B[name]), float)
    bd = board(name, spec.get("dphi"), spec.get("dtheta"))
    kw = {}
    if spec.get("D2"):
        kw["gid"] = group_ids(bd, spec["D2"])
    if spec.get("D3") == "one":
        kw["no_both"] = True
    if spec.get("D3") == "strict_x":
        kw["only_x"] = True
    if spec.get("D6"):
        kw["forbid_sw"] = bd.rho > spec["D6"] * 1e-3
    if spec.get("D1"):
        # pass 1: loose design, rank switched elements by amp^2 x disagreement, keep the top fraction
        r1, (m1, nd1, ids1, amp1, best1) = run_g(name, opt, eps, spec.get("lam1", lam * 0.4), bd=bd, iters=2, **kw)
        rank, who = [], []
        for ax in (0, 1):
            sw = ids1[ax] >= 0
            dis = np.mean([np.abs(m1[k, ax] - m1[l, ax]) for k in range(bd.K) for l in range(k + 1, bd.K)], axis=0)
            rk = amp1**2 * dis
            for i in np.flatnonzero(sw):
                rank.append(rk[i]); who.append((ax, i))
        order = np.argsort(rank)[::-1]
        keep = order[: int(round(spec["D1"] * len(order)))]
        allow = np.zeros((2, bd.n), bool)
        for t in keep:
            allow[who[t][0], who[t][1]] = True
        kw["forbid_x"], kw["forbid_y"] = ~allow[0], ~allow[1]
        lam = spec.get("lam2", 0.0)
    r, (m_, ndc_, ids_, amp_, best_) = run_g(name, opt, eps, lam, bd=bd, **kw)
    out = row(r, {k: (v if not isinstance(v, np.ndarray) else v.tolist()) for k, v in spec.items()})
    out["x_static_cell_pct"] = round(100 * float(np.mean(ids_[0] < 0)), 1)
    out["y_static_cell_pct"] = round(100 * float(np.mean(ids_[1] < 0)), 1)
    out["both_static_cell_pct"] = round(100 * float(np.mean((ids_[0] < 0) & (ids_[1] < 0))), 1)
    out["pointing_err_deg"] = [round(float(np.hypot(r["states"][k]["pk"][0] - 25.0, ((r["states"][k]["pk"][1] - K.SETS[name][0][k] + 180) % 360 - 180) * np.sin(np.radians(25.0)))), 2) for k in range(bd.K)]
    if spec.get("D2"):
        out["nyquist_violation_pct"] = nyquist(bd, spec["D2"], amp_)
    return out


def nyquist(bd, mode, amp):
    """Per state: amplitude-weighted fraction of cells where the hologram phase advance across one group pitch exceeds pi (code under-sampled)."""
    dirs = {"1x2x": [(1, 0)], "1x2y": [(0, 1)], "2x2": [(1, 0), (0, 1)]}[mode]
    res = []
    for k in range(bd.K):
        gx = K.BETA * bd.rhat[0] - st.K0 * np.sin(bd.thetas[k]) * np.cos(bd.phis[k])
        gy = K.BETA * bd.rhat[1] - st.K0 * np.sin(bd.thetas[k]) * np.sin(bd.phis[k])
        bad = np.zeros(bd.n, bool)
        for d in dirs:
            bad |= np.abs(gx * d[0] + gy * d[1]) * 2 * PITCH > np.pi
        res.append(round(100 * float(np.sum(amp**2 * bad) / np.sum(amp**2)), 1))
    return res
