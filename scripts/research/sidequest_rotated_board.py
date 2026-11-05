"""D4: one-diode rotated disc cell (RT5880, diode axis chi_i per cell from 8 angles, rotated-load fix fcc6691), K = 2 and K = 4 CP set B. Cell numbers from rotated_cell_bloch.json (Bloch tensor at kt = beta);
one element per cell: J = m (u . rho_hat) u e^{-j beta rho}, m = +1/-1 = OFF/ON swing of the disc cell (alpha A0_rot for the 50/50 pattern). Options per cell: 8 axes x (static levels | switched ON-patterns).
Static levels: the normalised round-11 library range (assumed to transfer to the disc). Output results/sidequest_hybrid/rotated_board.json."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_arch_study as st
import sidequest_coded_k4 as K
import sidequest_diode_reduction as D

CHIS = np.radians(np.arange(8) * 22.5)
J = json.loads(Path("results/sidequest_hybrid/rotated_cell_bloch.json").read_text(encoding="utf-8"))
import os

CH = os.environ.get("ROT_CHI", "chi0")  # chi0: axis-aligned (exact mesh), chi45: staircased 45 deg (alpha x2.1): sensitivity
A0R = float(J[CH]["alpha_beam"])
AA0R = float(J[CH]["alpha_abs"])
VS = np.concatenate([[0.0], D.VS_AN])


def axis_contrib(bd, amp):
    ph = np.exp(-1j * K.BETA * bd.rho)
    C = np.zeros((len(CHIS), bd.K, bd.n), complex)
    cur = []
    for a, ch in enumerate(CHIS):
        u = np.array([np.cos(ch), np.sin(ch)])
        d = u[0] * bd.rhat[0] + u[1] * bd.rhat[1]
        f = (np.pi / 2) * amp * d * ph
        jx, jy = f * u[0], f * u[1]
        cur.append((jx, jy))
        for k in range(bd.K):
            t = st.JONES[bd.pols[k]]
            et = np.cos(bd.thetas[k]) * (np.cos(bd.phis[k]) * jx + np.sin(bd.phis[k]) * jy)
            ep = -np.sin(bd.phis[k]) * jx + np.cos(bd.phis[k]) * jy
            C[a, k] = (np.conj(t[0]) * et + np.conj(t[1]) * ep) * bd.arr[k]
    return C, cur


def options(Kn, cls=None):
    W = K.words(Kn)
    m, ax, nd, ld, wid = [], [], [], [], []
    for a in range(len(CHIS)):
        for v in VS:
            m.append(np.full(Kn, v)); ax.append(a); nd.append(0); ld.append(np.zeros(Kn)); wid.append(-1)
        for wi, w in enumerate(W):
            if cls is not None and K.canon(w) not in cls:
                continue
            m.append(np.array([-1.0 if b else 1.0 for b in w])); ax.append(a); nd.append(1); ld.append(np.array(w, float)); wid.append(wi)
    return np.array(m), np.array(ax), np.array(nd, float), np.array(ld), np.array(wid), W


def design(bd, amp, opt, eps, lam, mu, g, forbid_sw=None):
    m, ax, nd, ld, wid, _ = opt
    C, _ = axis_contrib(bd, amp)
    err = np.zeros((bd.n, m.shape[0]))
    for k in range(bd.K):
        tau = g * amp * np.exp(-1j * eps[k])
        err += np.abs(tau[:, None] - m[None, :, k] * C[ax, k, :].T) ** 2
    err /= np.mean(np.abs(g * amp) ** 2)
    ud2 = np.stack([(np.cos(c) * bd.rhat[0] + np.sin(c) * bd.rhat[1]) ** 2 for c in CHIS])  # (8, n)
    if mu:
        err += mu * np.sum(ld, axis=1)[None, :] * ud2[ax].T
    cost = err + lam * nd[None, :]
    if forbid_sw is not None:
        cost = np.where(forbid_sw[:, None] & (nd > 0)[None, :], np.inf, cost)
    best = np.argmin(cost, axis=1)
    return m[best], ax[best], nd[best], ld[best], wid[best], float(cost[np.arange(bd.n), best].sum())


def run(name, lam, mu=0.05, g=2.0, cls=None, forbid_R=None, iters=3):
    bd = D.board(name)
    Kn = bd.K
    opt = options(Kn, cls)
    eps = np.array(D.EPS_B[name], float)
    amp = np.sqrt(2 * A0R / bd.rho) * np.exp(-0.5 * (A0R + AA0R) * (bd.rho - st.RHO_KEEP))
    forbid = None if forbid_R is None else bd.rho > forbid_R * 1e-3
    for _ in range(iters):
        mv, axv, ndv, ldv, widv, Jc = design(bd, amp, opt, eps, lam, mu, g, forbid)
        ud2 = np.stack([(np.cos(c) * bd.rhat[0] + np.sin(c) * bd.rhat[1]) ** 2 for c in CHIS])[axv, np.arange(bd.n)]
        alpha = A0R * np.mean(mv**2, axis=1) * ud2
        aa = K.AA_BG + 2 * AA0R * np.mean(ldv, axis=1) * ud2
        integ, ie = bd.profile(alpha, aa)
        amp = np.sqrt(2 * A0R / bd.rho) * np.exp(-integ)
        conv = float((2 * alpha * np.exp(-2 * integ) / bd.rho * st.PITCH**2 / (2 * np.pi)).sum())
        edge = float(np.mean(np.exp(-2 * ie)))
    C, cur = axis_contrib(bd, amp)
    a_id = np.sqrt(np.mean(mv**2, axis=1) * ud2) * amp
    res = []
    for k in range(Kn):
        jx = mv[:, k] * np.array([cur[a][0][i] for i, a in enumerate(axv)])
        jy = mv[:, k] * np.array([cur[a][1][i] for i, a in enumerate(axv)])
        et, ep = st.far_field(bd.x, bd.y, jx, jy)
        jj = st.jones_to_current(st.JONES[bd.pols[k]], bd.thetas[k], bd.phis[k])
        idl = a_id * np.exp(-1j * st.K0 * np.sin(bd.thetas[k]) * (bd.x * np.cos(bd.phis[k]) + bd.y * np.sin(bd.phis[k])))
        ei, epi = st.far_field(bd.x, bd.y, jj[0] * idl, jj[1] * idl)
        mi = st.metrics(ei, epi, bd.thetas[k], bd.phis[k], bd.pols[k])
        mb = st.metrics(et, ep, bd.thetas[k], bd.phis[k], bd.pols[k])
        i0 = mi["_i"]
        ratio = (np.abs(et) ** 2 + np.abs(ep) ** 2)[i0] / (np.abs(ei) ** 2 + np.abs(epi) ** 2)[i0]
        gain = mi["D_dbi"] + 10 * np.log10(K.LAUNCH * conv) + 10 * np.log10(ratio)
        res.append({"gain": float(gain), "SLL": mb["SLL_db"], "XPD": mb["XPD_db"], "ap": 10 ** ((gain - K.DMAX) / 10), "pk": (mb["peak_theta_deg"], mb["peak_phi_deg"])})
    words = {K.canon(opt[5][int(w)]) for w in widv if w >= 0}
    nd_total = int(ndv.sum())
    spec = {"K": name, "id": "D4", "lam": lam, "D6": forbid_R, "cls": None if cls is None else "L3"}
    return {"spec": spec, "K": Kn, "gain": [round(x["gain"], 2) for x in res], "gain_mean": round(float(np.mean([x["gain"] for x in res])), 3), "SLL_worst": round(max(x["SLL"] for x in res), 1),
            "XPD_min": round(min(x["XPD"] for x in res), 1), "ap_eff_pct": round(100 * float(np.mean([x["ap"] for x in res])), 2), "diodes": nd_total, "static_pct": round(100 * (1 - nd_total / bd.n), 1),
            "lines_x": len(words), "lines_y": 0, "lines_noinv_x": len({int(w) for w in widv if w >= 0}), "lines_noinv_y": 0, "edge_pct": round(100 * edge, 1), "SW_to_beam_pct": round(100 * conv, 1),
            "axis_hist": np.bincount(axv, minlength=8).tolist(), "id": "D4"}


if __name__ == "__main__":
    rows = []
    for name in ("K2", "K4b"):
        lams = [0.1, 0.2, 0.4, 0.7, 1.0, 1.5] if name == "K2" else [0.2, 0.35, 0.55, 0.8, 1.2, 1.8]
        for l in lams:
            for cls in (None, D.CLS3) if name == "K4b" else (None,):
                rows.append(run(name, l, cls=cls))
                print({k: rows[-1][k] for k in ("spec", "gain", "SLL_worst", "diodes", "lines_x")}, flush=True)
        for rc in (65.0, 95.0):
            for l in (lams[1], lams[3]):
                rows.append(run(name, l, forbid_R=rc, cls=D.CLS3 if name == "K4b" else None))
                print({k: rows[-1][k] for k in ("spec", "gain", "SLL_worst", "diodes")}, flush=True)
    Path(f"results/sidequest_hybrid/rotated_board_{CH}.json").write_text(json.dumps({"cell": {"A0": A0R, "AA0": AA0R, "bloch": J}, "rows": rows}, indent=1, default=float), encoding="utf-8", newline="\n")
