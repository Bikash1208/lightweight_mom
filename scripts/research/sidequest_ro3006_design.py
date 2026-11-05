"""RO3006 1.27 mm board (r 2.22 cell): design (i) = one diode group, one control line, default-OFF two-beam (25, 0) / (25, 180), LP-V, 150 mm hexagon, <= 514 diodes, with the non-diode sidelobe fixes.  MODEL_ONLY.
Cell data: MoM mesh 96 (final cell only), HFSS calibration: ON fixed (HFSS ON matched), OFF contrast x CAL (0.75 from the r 2.20 cell check at 90 deg: HFSS shift 0.342 GHz against 0.457), applied to the whole library as an affine map of the susceptance
(normalised levels unchanged; Ybar and dY change) -> direct A0 (first-order leakage, `sidequest_a0_direct`) and beta/k0.
Fixes: I1 analog static levels, I2 per-cell switched pairs (slit-length family, OFF / ON levels per geometry), I3 radial amplitude taper on the target, I6 edge-ring boost (target amplitude raised in the outer rows to deplete the surface wave),
I7 joint state optimisation with a sidelobe-power term (coordinate descent over the cell options, exact change of the objective per option, far-field matrix on a coarse grid).  I4 (launcher): bare launch 0.71 for the numbers; the pad (0.96) is UNVALIDATED and reported separately.
Calibration to HFSS (board_K2_default, RT5880, lam 0.28): gain offsets state 1 -1.78 dB (9.1 dBi), state 2 -2.7 dB (8.3 dBi).  Non-beam floor: NB = (1 - launch) + launch x edge, absolute level 10 log(NB x 7.45) [dBi] (HFSS strongest sidelobe 8.3 / 8.2 dBi in both states).
    python sidequest_ro3006_design.py ablation | final ...      (see main)"""
import json
import math
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
F = 10.5e9
ER, TD, H = 6.15, 0.002, 1.27
R_CELL = float(os.environ.get("RO_R", "2.22"))
CAL = 0.75
OFFSET = {0: -1.78, 1: -2.7}
DNB = 2.0   # non-beam floor peaking: HFSS ffd of board_K2_default: ring 60-90 + back = 71-76 % of the power, azimuthally uniform, D about 1.6 dBi (the 8.3 dBi sidelobes are aperture lobes at theta 12-26, not the floor)
MARGIN = (R_CELL + 0.2) * 1e-3  # disc + 0.2 mm to the outline
MAXD = 514
RES = Path("results/sidequest_hybrid")


def load_lib(mesh=96, decouple=False):
    d = json.loads((RES / f"ro3006_r{R_CELL}_mesh{mesh}_lib.json").read_text())
    rows = d["rows"]

    def g(k, f):
        return complex(*next(x for x in rows if x["kind"] == k and abs(x["fr"] - f) < 1e-9)["Y_mS"])

    yo, yn = g("open", 1.0), g("ON", 1.0)
    mid, half = (yo.imag + yn.imag) / 2, (yo.imag - yn.imag) / 2

    def m(y):
        return (y.imag - mid) / half

    lev = {"solid": m(g("none", 0.0)), "strap_0.79": m(g("strap", 0.786)), "strap_1.0": m(g("strap", 1.0))}
    for f in (0.357, 0.643, 0.786, 0.929):
        lev[f"open_{f}"] = m(g("open", f))
    pairs = {1.0: (1.0, -1.0, 1.0)}
    for f in (0.643, 0.786, 0.929):
        yon = g("ON", f)
        pairs[f] = (m(g("open", f)), m(yon), yon.real / yn.real)
    if decouple:
        # statics at full MoM depth; switched pairs limited to the calibrated contrast (ON fixed): m_off_cal = m_on + CAL (m_off - m_on); A0 and beta from the full MoM pair / mean
        pairs = {f: (v[1] + CAL * (v[0] - v[1]), v[1], v[2]) for f, v in pairs.items()}
        import sidequest_a0_direct as A0m

        a0f, bkf = A0m.a0_direct(mid * 1e-3, abs(yo - yn) * 1e-3, ER, H * 1e-3, F, tand=TD)
        return {"levels": lev, "pairs": pairs, "A0": float(a0f), "beta_k0": float(bkf), "Ybar_cal_mS": mid, "dY_cal_mS": abs(yo - yn), "Ybar_mom_mS": mid, "dY_mom_mS": abs(yo - yn), "Y_off": [yo.real, yo.imag], "Y_on": [yn.real, yn.imag], "mesh": mesh, "decoupled": True}
    # calibration: Y_cal = Y_ON + CAL (Y - Y_ON) (affine, ON fixed)
    mid_c = mid - half * (1 - CAL)
    half_c = CAL * half
    dY = 2 * half_c  # imaginary part of the calibrated pair contrast [mS]
    dYc = abs(yo - yn) * CAL
    import sidequest_a0_direct as A

    a0, bk = A.a0_direct(mid_c * 1e-3, dYc * 1e-3, ER, H * 1e-3, F, tand=TD)
    return {"levels": lev, "pairs": pairs, "A0": float(a0), "beta_k0": float(bk), "Ybar_cal_mS": mid_c, "dY_cal_mS": dYc, "Ybar_mom_mS": mid, "dY_mom_mS": abs(yo - yn), "Y_off": [yo.real, yo.imag], "Y_on": [yn.real, yn.imag], "mesh": mesh}


def engine(lib, launch=0.71):
    os.environ["K4_A0"] = "%.4f" % lib["A0"]
    os.environ["K4_BETA_K0"] = "%.4f" % lib["beta_k0"]
    import numpy as np
    import sidequest_coded_k4 as K
    import sidequest_k2_default_design as D

    K.A0 = lib["A0"]
    K.BETA_K0 = lib["beta_k0"]
    K.LAUNCH = launch
    K.AA_BG = 0.15 * (TD / 0.0009) ** 0.5
    D.setup(F, 150e-3, MARGIN)
    bd = K.KBoard([0, 180], ["LP-V", "LP-V"])
    return np, K, D, bd


def element_matrix(np, K, D, bd, amp, th_deg, ph_deg):
    """E[i, comp, d]: far field (theta, phi components) of a unit scalar-cell m at cell i (current f_i * rho_hat_i) in the directions d."""
    st = D.st
    j, _ = bd.unit(amp)
    fx, fy = j["x"][0], j["y"][1]  # f*rx, f*ry
    th, ph = np.radians(th_deg), np.radians(ph_deg)
    arr = np.exp(1j * st.K0 * np.sin(th)[None, :] * (np.cos(ph)[None, :] * bd.x[:, None] + np.sin(ph)[None, :] * bd.y[:, None]))
    et = (np.cos(th) * (np.cos(ph)[None, :] * fx[:, None] + np.sin(ph)[None, :] * fy[:, None])) * arr
    ep = (-np.sin(ph)[None, :] * fx[:, None] + np.cos(ph)[None, :] * fy[:, None]) * arr
    return np.stack([et, ep], axis=1)


def cd_design(np, K, D, bd, amp, opt, eps, lam, mu, g, nu, sweeps, m0, seed=0, mu_sep=0.0, sl_mode="far", kappa=0.0):
    mx, _my, nd, ld, _w1, _w2, _W = opt
    n, nopt = bd.n, mx.shape[0]
    _, c = bd.unit(amp)
    C = [c["x"][q] + c["y"][q] for q in range(2)]
    s = np.mean(np.abs(g * amp) ** 2)
    own = np.zeros((n, nopt))
    for k in range(2):
        tau = g * amp * np.exp(-1j * eps[k])
        pred = c["x"][k][:, None] * mx[None, :, k] + c["y"][k][:, None] * mx[None, :, k]
        own += np.abs(tau[:, None] - pred) ** 2
    own = own / s + lam * nd[None, :]
    if kappa:  # neutral static: where the two states' targets disagree, penalise radiating static levels
        t1, t2 = g * amp * np.exp(-1j * eps[0]), g * amp * np.exp(-1j * eps[1])
        dis = np.abs(t1 - t2) ** 2 / s
        dis = dis / dis.mean()
        own += kappa * dis[:, None] * (mx[None, :, 0] ** 2 + mx[None, :, 1] ** 2) / 2 * (nd[None, :] == 0)
    # directions: coarse grid + the two design directions
    th = np.arange(0.0, 89.0, 4.0)
    ph = np.arange(0.0, 359.0, 6.0)
    TG, PG = np.meshgrid(th, ph, indexing="ij")
    thd = np.concatenate([TG.ravel(), [25.0, 25.0]])
    phd = np.concatenate([PG.ravel(), [0.0, 180.0]])
    ib = [len(thd) - 2, len(thd) - 1]  # beam direction indices
    E = element_matrix(np, K, D, bd, amp, thd, phd)
    # exclusion of the main beam region: 35 deg around the beam of the state
    sl = []
    for k in range(2):
        cosd = np.cos(np.radians(thd)) * math.cos(math.radians(25)) + np.sin(np.radians(thd)) * math.sin(math.radians(25)) * np.cos(np.radians(phd) - math.radians(0 if k == 0 else 180))
        dist = np.degrees(np.arccos(np.clip(cosd, -1, 1)))
        if sl_mode == "ring":  # theta 8-32 deg cone, outside 12 deg of the own beam (includes the other beam)
            sl.append((thd >= 8.0) & (thd <= 32.0) & (dist > 12.0))
        else:
            sl.append(dist > 35.0)
    cur = m0.copy()
    rng = np.random.default_rng(seed)
    pen = mu / (n * s)
    sw_i = nd > 0
    for _ in range(sweeps):
        Fm = [np.einsum("i,icd->cd", mx[cur, k], E) for k in range(2)]
        Pk = [float(np.sum(np.abs(Fm[k][:, ib[k]]) ** 2)) for k in range(2)]
        Wt = [np.sum(C[1 - k] * mx[cur, k]) for k in range(2)]
        Ws = [np.sum((C[1 - k] * mx[cur, k])[nd[cur] == 0]) for k in range(2)]
        Ww = [np.sum((C[1 - k] * mx[cur, k])[nd[cur] > 0]) for k in range(2)]
        total_nd = float(nd[cur].sum())
        changed = 0
        for i in rng.permutation(n):
            mo = mx[cur[i]]
            Wi = [Wt[k] - C[1 - k][i] * mo[k] for k in range(2)]
            cost = own[i].copy()
            for k in range(2):
                cost += pen * np.abs(Wi[k] + C[1 - k][i] * mx[:, k]) ** 2
            if mu_sep:
                wsi = [Ws[k] - (C[1 - k][i] * mo[k] if nd[cur[i]] == 0 else 0.0) for k in range(2)]
                wwi = [Ww[k] - (C[1 - k][i] * mo[k] if nd[cur[i]] > 0 else 0.0) for k in range(2)]
                for k in range(2):
                    cost += (mu_sep / (n * s)) * np.where(nd == 0, np.abs(wsi[k] + C[1 - k][i] * mx[:, k]) ** 2 + np.abs(wwi[k]) ** 2, np.abs(wwi[k] + C[1 - k][i] * mx[:, k]) ** 2 + np.abs(wsi[k]) ** 2)
            if nu:
                for k in range(2):
                    Fi = Fm[k] - E[i] * mo[k]  # (2, ND)
                    Fn = Fi[None, :, :] + E[i][None, :, :] * mx[:, k][:, None, None]  # (nopt, 2, ND)
                    inten = np.abs(Fn[:, 0]) ** 2 + np.abs(Fn[:, 1]) ** 2
                    cost += nu * np.sum((inten[:, sl[k]] / Pk[k]) ** 2, axis=1) / sl[k].sum()
            allowed = (total_nd - nd[cur[i]] + nd) <= MAXD
            cost = np.where(allowed, cost, np.inf)
            o = int(np.argmin(cost))
            if o != cur[i]:
                changed += 1
                for k in range(2):
                    Fm[k] = Fm[k] - E[i] * mx[cur[i], k] + E[i] * mx[o, k]
                    Wt[k] = Wi[k] + C[1 - k][i] * mx[o, k]
                    Ws[k] += (C[1 - k][i] * (mx[o, k] if nd[o] == 0 else 0) - C[1 - k][i] * (mo[k] if nd[cur[i]] == 0 else 0))
                    Ww[k] += (C[1 - k][i] * (mx[o, k] if nd[o] > 0 else 0) - C[1 - k][i] * (mo[k] if nd[cur[i]] > 0 else 0))
                total_nd += nd[o] - nd[cur[i]]
                cur[i] = o
        if changed == 0:
            break
    wx = np.array([bd.rhat[0] ** 2, bd.rhat[1] ** 2])
    m = np.stack([np.stack([mx[cur, k], mx[cur, k]]) for k in range(2)])
    loss = np.stack([ld[cur, 0, k] * wx[0] + ld[cur, 1, k] * wx[1] for k in range(2)])
    return m, nd[cur], loss, cur


def evaluate(np, K, D, bd, m, nd, launch):
    """Full-grid far field per state with the amplitude profile of this m; returns calibrated numbers."""
    st = D.st
    sw = nd > 0
    loss = np.zeros((2, bd.n))
    loss[1, sw] = 1.0
    amp = np.sqrt(2 * K.A0 / bd.rho) * np.exp(-0.5 * (K.A0 + K.AA0) * (bd.rho - st.RHO_KEEP))
    for _ in range(4):
        amp, conv, edge = K.amp_profile(bd, m, loss)
    j, _ = bd.unit(amp)
    a_id = np.sqrt(np.mean(K.vec2(bd, m), axis=0)) * amp
    out = {"conv": float(conv), "edge": float(edge), "diodes": int(nd.sum()), "states": []}
    nb = (1 - launch) + launch * edge
    for k in range(2):
        cur = (m[k, 0] * j["x"][0] + m[k, 1] * j["y"][0], m[k, 0] * j["x"][1] + m[k, 1] * j["y"][1])
        et, ep = st.far_field(bd.x, bd.y, *cur)
        jj = st.jones_to_current(st.JONES[bd.pols[k]], K.THETA, bd.phis[k])
        idl = a_id * np.exp(-1j * st.K0 * np.sin(K.THETA) * (bd.x * np.cos(bd.phis[k]) + bd.y * np.sin(bd.phis[k])))
        ei, epi = st.far_field(bd.x, bd.y, jj[0] * idl, jj[1] * idl)
        mi = st.metrics(ei, epi, K.THETA, bd.phis[k], bd.pols[k])
        mb = st.metrics(et, ep, K.THETA, bd.phis[k], bd.pols[k])
        i0 = mi["_i"]
        inten = np.abs(et) ** 2 + np.abs(ep) ** 2
        inten_i = np.abs(ei) ** 2 + np.abs(epi) ** 2
        g_des = mi["D_dbi"] + 10 * np.log10(launch * conv) + 10 * np.log10(inten[i0] / inten_i[i0])
        scale = 10 ** (g_des / 10) / inten[i0]
        ipk = np.unravel_index(np.argmax(inten), inten.shape)
        lv = {}
        for q in range(2):
            cosd = np.cos(st.TG) * np.cos(K.THETA) + np.sin(st.TG) * np.sin(K.THETA) * np.cos(st.PG - bd.phis[q])
            near = np.degrees(np.arccos(np.clip(cosd, -1, 1))) <= 4.0
            lv[q] = float(10 * np.log10(scale * inten[near].max() + 1e-30))
        oth = lv[1 - k]
        g_cal = float(g_des) + OFFSET[k]
        gpk = float(10 * np.log10(scale * inten[ipk])) + OFFSET[k]
        floor_abs = 10 * math.log10(nb * DNB)
        sll_ap = float(mb["SLL_db"])
        ring = {t: float(10 * np.log10(np.median(inten[int(np.argmin(np.abs(np.degrees(st.TG[:, 0]) - t)))]) / inten.max() + 1e-30)) for t in (12, 24)}
        out["states"].append({"ring_median_rel_peak_db": ring, "gain_design_dir": g_cal, "gain_peak": gpk, "peak_dir": [float(np.degrees(st.TG[ipk])), float(np.degrees(st.PG[ipk]))], "SLL_aperture_db": sll_ap, "SLL_floor_db": floor_abs - gpk, "SLL_total_db": max(sll_ap, floor_abs - gpk),
                              "other_beam_dbi": oth + OFFSET[k], "isolation_db": g_cal - (oth + OFFSET[k]), "XPD_db": float(mb["XPD_db"])})
    out["NB"] = nb
    return out


def run_config(args):
    name, cfg = args
    warnings.simplefilter("ignore")
    global MAXD
    MAXD = cfg.get("maxd", 514)
    lib = load_lib(cfg.get("mesh", 96), cfg.get("decouple", False))
    np, K, D, bd = engine(lib, cfg.get("launch", 0.71))
    if cfg.get("rim_db"):  # absorbing rim (extra part): amplitude loss a over the outer 12.8 mm, 2 a w = ln(10^(db/10))
        a_rim = math.log(10 ** (cfg["rim_db"] / 10)) / (2 * 12.8e-3)
        rim = np.where(bd.rho > bd.redge - 12.8e-3, a_rim, 0.0)

        def ap(bd_, m_, loss_):
            alpha = K.A0 * np.mean(K.vec2(bd_, m_), axis=0)
            aa = K.AA_BG + 2 * K.AA0 * np.mean(loss_, axis=0) + rim
            integ, ie = bd_.profile(alpha, aa)
            amp_ = np.sqrt(2 * K.A0 / bd_.rho) * np.exp(-integ)
            w_ = 2 * alpha * np.exp(-2 * integ) / bd_.rho * D.st.PITCH ** 2 / (2 * np.pi)
            return amp_, float(w_.sum()), float(np.mean(np.exp(-2 * ie)))

        K.amp_profile = ap
    st = D.st
    lv = lib["levels"]
    if cfg.get("binary_statics"):
        levels = sorted({round(float(lv["open_0.929"]) if False else 1.0, 4), round(float(lv["strap_1.0"]), 4)})
    else:
        levels = sorted({round(float(v), 4) for v in lv.values()} | {1.0})
    pairs = [lib["pairs"][1.0]] if not cfg.get("multi_pairs") else [lib["pairs"][f] for f in sorted(lib["pairs"])]
    opt = D.options(levels, pairs, [(0, 1)])
    eps = [0.0, np.radians(180)]
    # target amplitude shape g(rho)
    g0 = cfg.get("g", 2.0)
    w = np.ones(bd.n)
    if cfg.get("taper") is not None:
        we = cfg["taper"]
        w = we + (1 - we) * (1 - np.clip(bd.rho / bd.redge, 0, 1) ** 2)
    if cfg.get("inv_taper"):  # inverse taper: target amplitude rising towards the rim (stronger modulation at large radius)
        w = w * (1.0 + cfg["inv_taper"] * np.clip(bd.rho / bd.redge, 0, 1) ** 2)
    if cfg.get("ring"):
        w = np.where(bd.rho > bd.redge - 12.8e-3, w * cfg["ring"], w)
    g = g0 * w
    mu = cfg.get("mu", 30.0)
    # stage 1: design_null with lam bisected to <= MAXD diodes
    lo, hi, r = 0.01, 1.5, None
    for _ in range(7):
        mid = 0.5 * (lo + hi)
        rr = D.run_design(bd, opt, eps, mid, g, mu=mu)
        if int(rr["nd"].sum()) > MAXD:
            lo = mid
        else:
            hi = mid
            r = rr
    if r is None:
        r = D.run_design(bd, opt, eps, hi, g, mu=mu)
    lam = float(r["lam"])
    m, nd, amp = r["m"], r["nd"], r["amp"]
    if cfg.get("nu") or cfg.get("kappa"):
        # stage 2: sidelobe-aware coordinate descent (start = design_null assignment)
        mxo = opt[0]
        cur = np.array([int(np.argmin(np.abs(mxo[:, 0] - m[0, 0, i]) + np.abs(mxo[:, 1] - m[1, 0, i]) + 1.0 * (opt[2] != nd[i]))) for i in range(bd.n)])
        for _it in range(cfg.get("outer", 3)):
            m, nd, loss, cur = cd_design(np, K, D, bd, amp, opt, eps, lam, mu, g, cfg.get("nu", 0.0), cfg.get("sweeps", 4), cur, mu_sep=cfg.get("mu_sep", 0.0), sl_mode=cfg.get("sl_mode", "far"), kappa=cfg.get("kappa", 0.0))
            amp, _c, _e = K.amp_profile(bd, m, loss)
    ev = evaluate(np, K, D, bd, m, nd, cfg.get("launch", 0.71))
    ev.update({"name": name, "cfg": {k: v for k, v in cfg.items()}, "lam": lam, "A0": lib["A0"], "beta_k0": lib["beta_k0"]})
    if cfg.get("keep"):
        import numpy as np2
        vals = {round(float(v), 4): k for k, v in lv.items()}
        pv = {round(float(v[0]), 4): f for f, v in lib["pairs"].items()}
        cells = []
        for i in range(bd.n):
            ix, iy = int(np2.rint(bd.x[i] / 6.4e-3)), int(np2.rint(bd.y[i] / 6.4e-3))
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


def summarize(ev):
    s = ev["states"]
    return "%-26s diodes %3d | gain(design dir) %.1f / %.1f | peak dir (%.0f,%.0f)/(%.0f,%.0f) | SLL ap %.1f / %.1f floor %.1f / %.1f total %.1f / %.1f | iso %.1f / %.1f | edge %.0f%% conv %.0f%%" % (
        ev["name"], ev["diodes"], s[0]["gain_design_dir"], s[1]["gain_design_dir"], *s[0]["peak_dir"], *s[1]["peak_dir"], s[0]["SLL_aperture_db"], s[1]["SLL_aperture_db"], s[0]["SLL_floor_db"], s[1]["SLL_floor_db"],
        s[0]["SLL_total_db"], s[1]["SLL_total_db"], s[0]["isolation_db"], s[1]["isolation_db"], 100 * ev["edge"], 100 * ev["conv"]) + " | ring med th12 %.1f/%.1f th24 %.1f/%.1f" % (s[0]["ring_median_rel_peak_db"][12], s[1]["ring_median_rel_peak_db"][12], s[0]["ring_median_rel_peak_db"][24], s[1]["ring_median_rel_peak_db"][24])


ABL = [("B0 binary statics, 1 pair", {"binary_statics": True}),
       ("+I1 analog statics", {}),
       ("+I2 switched pairs", {"multi_pairs": True}),
       ("+I3 taper 0.5", {"multi_pairs": True, "taper": 0.5}),
       ("+I6 ring x1.5", {"multi_pairs": True, "taper": 0.5, "ring": 1.5}),
       ("+I7 sidelobe CD nu 100", {"multi_pairs": True, "taper": 0.5, "ring": 1.5, "nu": 100.0}),
       ("+I7 sidelobe CD nu 400", {"multi_pairs": True, "taper": 0.5, "ring": 1.5, "nu": 400.0}),
       ("I7 only (no I3/I6) nu 400", {"multi_pairs": True, "nu": 400.0})]

if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "ablation"
    if what == "lib":
        lib = load_lib(96)
        print(json.dumps({k: v for k, v in lib.items() if k != "levels"}, indent=1))
        print({k: round(float(v), 3) for k, v in lib["levels"].items()})
    elif what == "final":
        cfg = json.loads(sys.argv[2]) if len(sys.argv) > 2 else {"multi_pairs": True, "g": 3.0}
        cfg["keep"] = True
        ev = run_config(("final", cfg))
        print(summarize(ev), flush=True)
        (RES / f"ro3006_board_design_r{R_CELL}.json").write_text(json.dumps(ev, indent=1, default=float), encoding="utf-8", newline=chr(10))
    elif what == "ablation":
        with ProcessPoolExecutor(len(ABL)) as ex:
            res = list(ex.map(run_config, ABL))
        for ev in res:
            print(summarize(ev), flush=True)
        (RES / "ro3006_r222_ablation.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline="\n")
