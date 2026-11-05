"""RO3006 r 2.30 board v2: STATE-DEPENDENT beta (each state's hologram designed with its own real beta/k0) and the 'add-a-beam' variant.  MODEL_ONLY.
Generic joint design: per cell option o (statics, switched pairs), cost = sum over (state s, beam b wanted in s) |tau_sb - C_sb m_s(o)|^2 / norm + lam nd(o) + kappa x disagreement x m^2 (statics, swap only) + mu |W|^2 (wrong beam in state 1 / swap),
coordinate descent with the exact change of the wrong-beam sums; C_sb is computed with the REAL beta of state s.  Calibration offsets as in `sidequest_ro3006_design` (state 1 -1.78 dB, state 2 -2.7 dB).
    python sidequest_ro3006_v2.py swap | add"""
import json
import math
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import sidequest_ro3006_design as RD
import sidequest_sll_def as SLD

F = RD.F
OFFSET = RD.OFFSET
MAXD = 514
BETAS = (1.33, 1.26)  # beta/k0 of state 1 (OFF) / state 2 (ON): mid-range of the r 2.20 HFSS-offset estimate
RES = Path("results/sidequest_hybrid")


def lib_for(decouple, betas):
    lib = RD.load_lib(96, decouple)
    import sidequest_a0_direct as A

    dy = lib["dY_cal_mS"] if decouple else lib["dY_cal_mS"]  # dY_cal is the full MoM dY when decoupled, 0.75 x when coupled (load_lib)
    a0s = []
    for bk in betas:
        yb = A.ybar_for_beta(bk, RD.ER, RD.H * 1e-3, F)
        a0s.append(A.a0_direct(yb, dy * 1e-3, RD.ER, RD.H * 1e-3, F, tand=RD.TD)[0])
    lib["A0_states"] = [float(x) for x in a0s]
    lib["A0"] = float(sum(a0s) / 2)
    lib["beta_states"] = list(betas)
    return lib


def couplings(np, K, bd, amp, betas_rad):
    C = {}
    base = K.BETA
    for s in range(2):
        K.BETA = betas_rad[s]
        _, c = bd.unit(amp)
        for b in range(2):
            C[(s, b)] = c["x"][b] + c["y"][b]
    K.BETA = base
    return C


def design_gen(np, K, bd, amp, opt, C, wanted, eps, lam, mu, g, kappa, sweeps, maxd, wrong, seed=0, wts=None):
    """wanted: {s: [beams]}; wrong: list of (s, b) beams to suppress; eps[b]: beam phase."""
    mx, _my, nd, ld, _a, _b, _W = opt
    n, nopt = bd.n, mx.shape[0]
    s_norm = np.mean(np.abs(g * amp) ** 2)
    own = np.zeros((n, nopt))
    tau = {}
    for s, bs in wanted.items():
        for b in bs:
            amp_sb = g * amp * np.exp(-1j * eps[b]) * ((wts or {}).get((s, b), 1.0) if wts is not None else 1.0 / math.sqrt(len(bs)))
            tau[(s, b)] = amp_sb
            own += np.abs(amp_sb[:, None] - C[(s, b)][:, None] * mx[None, :, s]) ** 2
    own = own / s_norm + lam * nd[None, :]
    # ideal (unconstrained real) m per state and the disagreement
    ideal = {}
    for s, bs in wanted.items():
        num = sum(np.real(np.conj(C[(s, b)]) * tau[(s, b)]) for b in bs)
        den = sum(np.abs(C[(s, b)]) ** 2 for b in bs)
        ideal[s] = np.where(den > 1e-9 * den.max(), num / np.maximum(den, 1e-300), 0.0)  # cells on the axes do not couple to the LP-V beam (den = 0)
    dis = np.abs(ideal[0] - ideal[1])
    if kappa:  # neutral-static weight as in the strong routine: |tau_1 - tau_2|^2 of the two states' targets (swap), normalised to its mean
        if (0, 0) in tau and (1, 1) in tau:
            d2 = np.abs(tau[(0, 0)] - tau[(1, 1)]) ** 2 / s_norm
            d2 = d2 / d2.mean()
        else:
            d2 = dis ** 2 / np.mean(dis ** 2)
        own += kappa * d2[:, None] * (mx[None, :, 0] ** 2 + mx[None, :, 1] ** 2) / 2 * (nd[None, :] == 0)
    stat = nd == 0
    if maxd is None:  # priced allocation (lam per diode), no hard cap: the caller bisects lam
        cur = np.argmin(own, axis=1)
        maxd = 10 ** 9
    else:  # start: best static everywhere, the maxd/2 cells with the largest benefit from switching get their best switched option
        best_stat = np.argmin(np.where(stat[None, :], own, np.inf), axis=1)
        best_sw = np.argmin(np.where(~stat[None, :], own, np.inf), axis=1)
        ar = np.arange(n)
        benefit = own[ar, best_stat] - own[ar, best_sw]
        order = np.argsort(-benefit)[: int(maxd // 2)]
        order = order[benefit[order] > 0]
        cur = best_stat.copy()
        cur[order] = best_sw[order]
    pen = mu / (n * s_norm)
    W = {sb: np.sum(C[sb] * mx[cur, sb[0]]) for sb in wrong}
    rng = np.random.default_rng(seed)
    total_nd = float(nd[cur].sum())
    for _ in range(sweeps):
        changed = 0
        for i in rng.permutation(n):
            cost = own[i].copy()
            Wi = {}
            for sb in wrong:
                Wi[sb] = W[sb] - C[sb][i] * mx[cur[i], sb[0]]
                cost += pen * np.abs(Wi[sb] + C[sb][i] * mx[:, sb[0]]) ** 2
            cost = np.where((total_nd - nd[cur[i]] + nd) <= maxd, cost, np.inf)
            o = int(np.argmin(cost))
            if o != cur[i]:
                changed += 1
                for sb in wrong:
                    W[sb] = Wi[sb] + C[sb][i] * mx[o, sb[0]]
                total_nd += nd[o] - nd[cur[i]]
                cur[i] = o
        if changed == 0:
            break
    m = np.stack([np.stack([mx[cur, k], mx[cur, k]]) for k in range(2)])
    wx = np.array([bd.rhat[0] ** 2, bd.rhat[1] ** 2])
    loss = np.stack([ld[cur, 0, k] * wx[0] + ld[cur, 1, k] * wx[1] for k in range(2)])
    return m, nd[cur], loss, cur, ideal, dis


def evaluate_gen(np, K, D, bd, m, nd, launch, betas_rad, wanted, ideal=None, dis=None, scan=None):
    """Far field per state with its own beta; levels at both beam directions (4 deg max), SLL outside 21 deg of the wanted beams, isolation, balance."""
    st = D.st
    sw = nd > 0
    loss = np.zeros((2, bd.n))
    loss[1, sw] = 1.0
    amp = np.sqrt(2 * K.A0 / bd.rho) * np.exp(-0.5 * (K.A0 + K.AA0) * (bd.rho - st.RHO_KEEP))
    for _ in range(4):
        amp, conv, edge = K.amp_profile(bd, m, loss)
    out = {"conv": float(conv), "edge": float(edge), "diodes": int(nd.sum()), "states": []}
    base = K.BETA
    nb = (1 - launch) + launch * edge
    a_id = np.sqrt(np.mean(K.vec2(bd, m), axis=0)) * amp
    for k in range(2):
        beams = wanted[k]
        K.BETA = betas_rad[k]
        j, _ = bd.unit(amp)
        cur = (m[k, 0] * j["x"][0] + m[k, 1] * j["y"][0], m[k, 0] * j["x"][1] + m[k, 1] * j["y"][1])
        et, ep = st.far_field(bd.x, bd.y, *cur)
        inten = np.abs(et) ** 2 + np.abs(ep) ** 2
        # ideal single-beam reference at the first wanted beam (same convention as the earlier boards)
        b0 = beams[0]
        jj = st.jones_to_current(st.JONES[bd.pols[b0]], K.THETA, bd.phis[b0])
        idl = a_id * np.exp(-1j * st.K0 * np.sin(K.THETA) * (bd.x * np.cos(bd.phis[b0]) + bd.y * np.sin(bd.phis[b0])))
        ei, epi = st.far_field(bd.x, bd.y, jj[0] * idl, jj[1] * idl)
        mi = st.metrics(ei, epi, K.THETA, bd.phis[b0], bd.pols[b0])
        inten_i = np.abs(ei) ** 2 + np.abs(epi) ** 2
        i0 = mi["_i"]
        g0 = mi["D_dbi"] + 10 * np.log10(launch * conv)
        lvl = {}
        for b in range(2):
            cosd = np.cos(st.TG) * np.cos(K.THETA) + np.sin(st.TG) * np.sin(K.THETA) * np.cos(st.PG - bd.phis[b])
            near = np.degrees(np.arccos(np.clip(cosd, -1, 1))) <= 4.0
            lvl[b] = float(g0 + 10 * np.log10(inten[near].max() / inten_i[i0] + 1e-30)) + OFFSET[k]
        # sidelobes: outside 21 deg of every wanted beam
        mask = np.ones(inten.shape, bool)
        for b in beams:
            cosd = np.cos(st.TG) * np.cos(K.THETA) + np.sin(st.TG) * np.sin(K.THETA) * np.cos(st.PG - bd.phis[b])
            mask &= np.degrees(np.arccos(np.clip(cosd, -1, 1))) > 21.0
        ref = min(lvl[b] for b in beams)
        sl = float(g0 + 10 * np.log10(inten[mask].max() / inten_i[i0] + 1e-30)) + OFFSET[k]
        ipk = np.unravel_index(np.argmax(inten), inten.shape)
        ring = {t: float(10 * np.log10(np.median(inten[int(np.argmin(np.abs(np.degrees(st.TG[:, 0]) - t)))]) / inten.max() + 1e-30)) for t in (12, 24)}
        th_g, ph_g = np.degrees(st.TG[:, 0]), np.degrees(st.PG[0])
        inside = np.zeros(inten.shape, bool)
        u_peaks = []
        for b in beams:
            ins, u0b = SLD.analyze(th_g, ph_g, inten, peak=(25.0, float(np.degrees(bd.phis[b]))), return_mask=True)
            inside |= ins
            u_peaks.append(u0b)
        hpbw = {}
        for b in beams:
            a_ = SLD.analyze(th_g, ph_g, inten, peak=(25.0, float(np.degrees(bd.phis[b]))))
            hpbw[b] = {"elevation": float(a_["HPBW_elevation_deg"]), "azimuth_angular": float(a_["HPBW_azimuth_angular_deg"])}
        outside = np.where(inside | (st.TG > np.radians(90.0)), -np.inf, inten)
        sl_fn = float(10 * np.log10(outside.max() / max(u_peaks)))
        sl_fn_weak = float(10 * np.log10(outside.max() / min(u_peaks)))
        jfn = np.unravel_index(np.argmax(outside), outside.shape)
        out["states"].append({"SLL_firstnull_db": sl_fn, "SLL_firstnull_vs_weaker_db": sl_fn_weak, "SLL_firstnull_dir": [float(np.degrees(st.TG[jfn])), float(np.degrees(st.PG[jfn]))], "hpbw": hpbw, "beam_gain_raw_model_dbi": {"b1_(25,0)": lvl[0] - OFFSET[k], "b2_(25,180)": lvl[1] - OFFSET[k]}, "beam_gain_dbi": {"b1_(25,0)": lvl[0], "b2_(25,180)": lvl[1]}, "SLL_db": sl - max(lvl[b] for b in beams), "SLL_vs_weaker_beam_db": sl - ref,
                              "peak_dir": [float(np.degrees(st.TG[ipk])), float(np.degrees(st.PG[ipk]))], "wanted": beams, "balance_db": (lvl[0] - lvl[1]) if len(beams) == 2 else None,
                              "isolation_db": (lvl[beams[0]] - lvl[1 - beams[0]]) if len(beams) == 1 else None, "ring_median_rel_peak_db": ring,
                              "floor_db": 10 * math.log10(nb * RD.DNB) - max(lvl[b] for b in beams)})
    K.BETA = base
    out["NB"] = nb
    return out


def pattern_state(np, K, D, bd, m, nd, k, beta_rad, amp):
    st = D.st
    base = K.BETA
    K.BETA = beta_rad
    j, _ = bd.unit(amp)
    cur = (m[k, 0] * j["x"][0] + m[k, 1] * j["y"][0], m[k, 0] * j["x"][1] + m[k, 1] * j["y"][1])
    et, ep = st.far_field(bd.x, bd.y, *cur)
    K.BETA = base
    return np.abs(et) ** 2 + np.abs(ep) ** 2


def run(args):
    name, cfg = args
    warnings.simplefilter("ignore")
    mode = cfg["mode"]
    betas = cfg.get("betas", BETAS)
    lib = lib_for(cfg.get("decouple", True), betas)
    np, K, D, bd = RD.engine(lib, 0.71)
    betas_rad = [b * D.st.K0 for b in betas]
    lv = lib["levels"]
    levels = sorted({round(float(v), 4) for v in lv.values()} | {1.0})
    pairs = [lib["pairs"][f] for f in sorted(lib["pairs"])]
    opt = D.options(levels, pairs, [(0, 1)])
    g0 = cfg.get("g", 2.0)
    x = np.clip(bd.rho / bd.redge, 0, 1)
    w = (1.0 + cfg.get("inv_taper", 0.5) * x ** 2) * np.where(bd.rho > bd.redge - 12.8e-3, cfg.get("ring", 1.5), 1.0)
    g = g0 * w
    eps = [0.0, math.radians(cfg.get("eps_B", 180.0))]
    wanted = {0: [0], 1: [1]} if mode == "swap" else {0: [0], 1: [0, 1]}
    wrong = [(0, 1), (1, 0)] if mode == "swap" else [(0, 1)]
    kappa = cfg.get("kappa", 2.0) if mode == "swap" else 0.0
    mu = cfg.get("mu", 30.0)
    amp = np.sqrt(2 * K.A0 / bd.rho) * np.exp(-0.5 * (K.A0 + K.AA0) * (bd.rho - D.st.RHO_KEEP))
    wts = {(0, 0): 1.0, (1, 0): 1.0, (1, 1): cfg.get("w2", 1.0)} if mode == "add" else None

    def design_at(lam_, amp_, sweeps_):
        C_ = couplings(np, K, bd, amp_, betas_rad)
        return design_gen(np, K, bd, amp_, opt, C_, wanted, eps, lam_, mu, g, kappa, sweeps_, None, wrong, wts=wts)

    lo, hi, best = 0.0, 2.0, None
    for _ in range(7):  # priced allocation: smallest price per diode with <= MAXD diodes
        lam = 0.5 * (lo + hi)
        out = design_at(lam, amp, cfg.get('sweeps_bis', 4))
        amp, _c, _e = K.amp_profile(bd, out[0], out[2])
        if int(out[1].sum()) > MAXD:
            lo = lam
        else:
            hi = lam
            best = (lam,) + out
    if best is None:
        lam = hi
        best = (lam,) + design_at(lam, amp, cfg.get('sweeps_bis', 4))
    lam, m, nd, loss, cur, ideal, dis = best
    for _ in range(cfg.get('outer', 2)):  # re-design with the converged amplitude profile; raise the price until <= MAXD
        amp, _c, _e = K.amp_profile(bd, m, loss)
        for _try in range(5):
            out = design_at(lam, amp, cfg.get('sweeps_fin', 5))
            if int(out[1].sum()) <= MAXD:
                break
            lam *= 1.08
        m, nd, loss, cur, ideal, dis = out
    ev = evaluate_gen(np, K, D, bd, m, nd, 0.71, betas_rad, wanted)
    amp, _c, _e = K.amp_profile(bd, m, loss)
    ev["beam2_on_off_contrast_db"] = ev["states"][1]["beam_gain_dbi"]["b2_(25,180)"] - ev["states"][0]["beam_gain_dbi"]["b2_(25,180)"]
    # compromise cells and the switching requirement
    sw = nd > 0
    must = float(np.mean(dis > 0.5))
    comp_static = int(np.sum((~sw) & (dis > 0.5)))
    dev = np.maximum(np.abs(m[0, 0] - ideal[0]), np.abs(m[1, 0] - ideal[1]))
    comp_any = int(np.sum(dev > 0.5))
    ev.update({"name": name, "cfg": cfg, "lam": float(lam), "A0": lib["A0"], "A0_states": lib["A0_states"], "switched_cells": int(sw.sum()), "n_cells": int(bd.n), "switched_fraction": float(sw.mean()),
               "cells_that_must_switch_fraction": must, "compromise_static_cells": comp_static, "cells_deviating_>0.5_from_ideal": comp_any})
    # margin to break-up: each state's real beta varied alone
    if cfg.get("margin"):
        marg = []
        k0 = D.st.K0
        for k in range(2):
            for d in (-0.20, -0.15, -0.10, -0.05, 0.0, 0.05, 0.10, 0.15):
                u = pattern_state(np, K, D, bd, m, nd, k, betas_rad[k] + d * k0, amp)
                # gain relative to the d = 0 peak: loss of the intensity at the wanted beam directions and peak position
                st = D.st
                ipk = np.unravel_index(np.argmax(u), u.shape)
                lvls = []
                for b in wanted[k]:
                    cosd = np.cos(st.TG) * np.cos(K.THETA) + np.sin(st.TG) * np.sin(K.THETA) * np.cos(st.PG - bd.phis[b])
                    near = np.degrees(np.arccos(np.clip(cosd, -1, 1))) <= 4.0
                    lvls.append(float(u[near].max()))
                marg.append({"state": k + 1, "dbeta": d, "beam_intensity": lvls, "peak_dir": [float(np.degrees(st.TG[ipk])), float(np.degrees(st.PG[ipk]))], "peak_intensity": float(u.max())})
        ev["margin_raw"] = marg
    ev["cells"] = None
    if cfg.get("keep"):
        vals = {round(float(v), 4): k for k, v in lv.items()}
        pv = {round(float(v[0]), 4): f for f, v in lib["pairs"].items()}
        cells = []
        for i in range(bd.n):
            ix, iy = int(np.rint(bd.x[i] / 6.4e-3)), int(np.rint(bd.y[i] / 6.4e-3))
            v = round(float(m[0, 0, i]), 4)
            if nd[i] > 0:
                cells.append([ix, iy, "SW", "open", pv[v]])
            elif abs(v - 1.0) < 1e-6:
                cells.append([ix, iy, "ST", "open", 1.0])
            else:
                nm = vals[v]
                cells.append([ix, iy, "ST"] + ({"solid": ["none", 0.0], "strap_0.79": ["strap", 0.786], "strap_1.0": ["strap", 1.0]}.get(nm) or ["open", float(nm.split("_")[1])]))
        ev["cells"] = cells
    return ev


def summarize(e):
    s = e["states"]
    b = lambda st_, k: st_["beam_gain_dbi"][k]  # noqa: E731
    return "%-30s diodes %3d | S1 b1 %.1f b2 %.1f SLLfn %.1f | S2 b1 %.1f b2 %.1f SLLfn %.1f (vs weaker %.1f) bal %s | contrast(b2) %.1f dB | edge %.0f%% A0 %.1f | [21deg SLL %.1f / %.1f]" % (
        e["name"], e["diodes"], b(s[0], "b1_(25,0)"), b(s[0], "b2_(25,180)"), s[0]["SLL_firstnull_db"], b(s[1], "b1_(25,0)"), b(s[1], "b2_(25,180)"), s[1]["SLL_firstnull_db"], s[1]["SLL_firstnull_vs_weaker_db"], None if s[1]["balance_db"] is None else round(s[1]["balance_db"], 1),
        e["beam2_on_off_contrast_db"], 100 * e["edge"], e["A0"], s[0]["SLL_db"], s[1]["SLL_db"])


if __name__ == "__main__":
    what = sys.argv[1]
    if what == "swap":
        cfgs = [(f"v2 swap dec g{g} mu{mu}", {"mode": "swap", "decouple": True, "g": g, "mu": mu}) for g in (1.5, 2.0) for mu in (30.0, 70.0)]
    elif what == "add":
        cfgs = [(f"add w2={w2} eps{e}", {"mode": "add", "decouple": True, "g": 2.0, "mu": 30.0, "eps_B": float(e), "w2": w2}) for w2 in (0.3, 0.7, 1.0, 1.4, 2.0, 3.0) for e in (0, 90, 180, 270)]
    else:
        cfgs = [("custom", json.loads(what))]
    with ProcessPoolExecutor(min(len(cfgs), 32)) as ex:
        res = list(ex.map(run, cfgs))
    for e in res:
        print(summarize(e), flush=True)
    Path(RES / f"ro3006_v2_{what if what in ('swap', 'add') else 'custom'}.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline=chr(10))
