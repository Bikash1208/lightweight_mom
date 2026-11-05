"""K2-default board designs (MODEL_ONLY): ONE diode group / one control line / no inverter. State 1 (all diodes OFF) = beam 1 (25, 0), state 2 (all ON) = beam 2 (25, 180), linear polarisation, 10.5 GHz, hexagon circumradius R (patch discs fully inside the outline).
Engine = sidequest_coded_k4 (per-cell least squares over option sets; HFSS-calibrated scalar cell: alpha 0.844 Np/m, beta/k0 1.0731). Static-level and family-B cell values: MoM library `k2_default_lib.json` (this script's normalisation: OFF/OFF = +1, ON/ON = -1 of the ls 2.8 cell).
  (i)  family A only: switched cell (+1 OFF, -1 ON) + analog statics (slit-length cells, strap).  Cells that would need ON -> OFF are compromised by the statics.
  (ii) adds family-B switched cells (through-slit, diode bridges islands: ON RAISES Y), several geometries (ls, pos) -> pairs (m_OFF, m_ON); still one line.
  old two-group design (K2_LP): GA/GB (two words), +-1 statics, 10.54 GHz, as built.
    python sidequest_k2_default_design.py        -> results/sidequest_hybrid/k2_default_designs.json"""
import json
import os
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(Path(__file__).resolve().parent))
os.environ.setdefault("K4_A0", "0.844")
os.environ.setdefault("K4_BETA_K0", "1.0731")
warnings.simplefilter("ignore")
import numpy as np

import sidequest_arch_study as st
import sidequest_coded_improve as I
import sidequest_coded_k4 as K

C0 = 299792458.0
LIB = ROOT / "results/sidequest_hybrid/k2_default_lib.json"
MARGIN = 3.2e-3  # disc radius 3.0 + 0.2 mm clearance to the board edge


# ---------------------------------------------------------------------------------------------- setup (frequency, outline)
def hexpoly(R):
    return np.array([[R * np.cos(np.radians(90 + 60 * i)), R * np.sin(np.radians(90 + 60 * i))] for i in range(6)])


def geom(R, margin):
    """Lattice cells (6.4 mm, integer index) whose disc (radius `margin`) lies inside the hexagon of circumradius R, rho >= 6.5 mm."""
    poly = hexpoly(R)
    n = int(np.ceil(R / st.PITCH)) + 3
    k = np.arange(-n, n + 1)
    xg, yg = np.meshgrid(st.PITCH * k, st.PITCH * k)
    x, y = xg.ravel(), yg.ravel()
    ok = np.ones(x.shape, bool)
    for a, b in zip(poly, np.roll(poly, -1, axis=0), strict=True):
        cr = (b[0] - a[0]) * (y - a[1]) - (b[1] - a[1]) * (x - a[0])
        ok &= cr / np.hypot(*(b - a)) >= margin - 1e-9
    ok &= np.hypot(x, y) >= st.RHO_KEEP
    x, y = x[ok], y[ok]
    psi = np.arctan2(y, x)
    return x, y, np.hypot(x, y), psi, st.redge(poly, None, psi)


def setup(f_hz, R, margin):
    st.K0 = 2 * np.pi * f_hz / C0
    K.BETA = K.BETA_K0 * st.K0
    K.LAM0 = C0 / f_hz
    area = 1.5 * np.sqrt(3) * R**2
    K.DMAX = 10 * np.log10(4 * np.pi * area / K.LAM0**2)
    I.hexgeom = lambda: geom(R, margin)


# ---------------------------------------------------------------------------------------------- library -> option levels
def library():
    rows = json.loads(LIB.read_text(encoding="utf-8"))["rows"]
    Y = {(r["family"], r["ls"], r["pos"], r["mode"]): complex(*r["Y_mS"][0]) * 0.5 + complex(*r["Y_mS"][3]) * 0.5 for r in rows}  # mean of Yxx, Yyy (mS)
    off, on = Y[("A", 2.8, 0.8, "open")], Y[("A", 2.8, 0.8, "ON")]
    c0, h0 = 0.5 * (off.imag + on.imag), 0.5 * (off.imag - on.imag)
    m = lambda y: (y.imag - c0) / h0  # noqa: E731
    stat = {}
    for (fam, ls, pos, mode), y in Y.items():
        if fam == "A" and mode in ("none", "open", "strap"):
            stat[(ls, "strap" if mode == "strap" else "open")] = (m(y), y.real)
    pairs = {}
    for (fam, ls, pos, mode), y in Y.items():
        if fam == "B" and mode == "open":
            yon = Y[(fam, ls, pos, "ON")]
            pairs[(ls, pos)] = (m(y), m(yon), yon.real / on.real)
    return stat, pairs, {"Y_OFF_mS": [off.real, off.imag], "Y_ON_mS": [on.real, on.imag], "c0_mS": c0, "h0_mS": h0}


def options(statics, switched, words):
    """statics: list of levels; switched: list of (m_off, m_on, g_on_rel); words: list of ON-patterns over the 2 states ((0,1) = ON in state 2). Scalar cell: m_x = m_y."""
    mx, nd, ld, wid = [], [], [], []
    z = np.zeros(2)
    for a in statics:
        mx.append(np.full(2, a)); nd.append(0); ld.append((z, z)); wid.append(-1)
    for wi, w in enumerate(words):
        on = np.array(w, float)
        for moff, mon, gon in switched:
            v = np.array([mon if b else moff for b in w])
            mx.append(v); nd.append(2); ld.append((on * gon, on * gon)); wid.append(wi)
    mx = np.array(mx)
    wid = np.array(wid)
    return mx, mx.copy(), np.array(nd, float), np.array(ld), wid, wid.copy(), [tuple(w) for w in words]


def design_null(bd, amp, opt, eps, lam, mu, g=1.0, sweeps=8, seed=0):
    """K.design plus an explicit wrong-beam penalty (K = 2): F = sum_k sum_i |tau_k,i - C_k,i m_k,i|^2 / s + lam*nd + mu * sum_k |W_kq|^2 / (N s), W_kq = sum_i C_q,i m_k,i (amplitude of the OTHER beam q in state k).
    Greedy coordinate descent over the cells (exact change of F for each option), start = K.design. Returns the same tuple as K.design."""
    mx, my, nd, ld, wxi, wyi, _W = opt
    _, c = bd.unit(amp)
    n, nopt = bd.n, mx.shape[0]
    C = [c["x"][q] + c["y"][q] for q in range(2)]  # scalar cell (mx = my): projection of a unit element on beam q
    s = np.mean(np.abs(g * amp) ** 2)
    own = np.zeros((n, nopt))
    for k in range(2):
        tau = g * amp * np.exp(-1j * eps[k])
        pred = c["x"][k][:, None] * mx[None, :, k] + c["y"][k][:, None] * my[None, :, k]
        own += np.abs(tau[:, None] - pred) ** 2
    own = own / s + lam * nd[None, :]
    cur = np.argmin(own, axis=1)
    W = [np.sum(C[1 - k] * mx[cur, k]) for k in range(2)]
    rng = np.random.default_rng(seed)
    pen = mu / (n * s)
    for _ in range(sweeps):
        changed = 0
        for i in rng.permutation(n):
            Wi = [W[k] - C[1 - k][i] * mx[cur[i], k] for k in range(2)]
            f = own[i].copy()
            for k in range(2):
                f += pen * np.abs(Wi[k] + C[1 - k][i] * mx[:, k]) ** 2
            o = int(np.argmin(f))
            if o != cur[i]:
                changed += 1
                for k in range(2):
                    W[k] = Wi[k] + C[1 - k][i] * mx[o, k]
                cur[i] = o
        if changed == 0:
            break
    best = cur
    wx = np.array([bd.rhat[0] ** 2, bd.rhat[1] ** 2])
    m = np.stack([np.stack([mx[best, k], my[best, k]]) for k in range(2)])
    loss = np.stack([ld[best, 0, k] * wx[0] + ld[best, 1, k] * wx[1] for k in range(2)])
    return m, nd[best], loss, (wxi[best], wyi[best])


# ---------------------------------------------------------------------------------------------- design + evaluation with isolation
def run_design(bd, opt, eps, lam, g, iters=3, mu=0.0):
    amp = np.sqrt(2 * K.A0 / bd.rho) * np.exp(-0.5 * (K.A0 + K.AA0) * (bd.rho - st.RHO_KEEP))
    for _ in range(iters):
        m, nd, loss, ids = design_null(bd, amp, opt, eps, lam, mu, g) if mu > 0 else K.design(bd, amp, opt, eps, lam, 0.0, g)
        amp, conv, edge = K.amp_profile(bd, m, loss)
    j, _ = bd.unit(amp)
    a_id = np.sqrt(np.mean(K.vec2(bd, m), axis=0)) * amp
    states = []
    for k in range(bd.K):
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
        gain = mi["D_dbi"] + 10 * np.log10(K.LAUNCH * conv) + 10 * np.log10(inten[i0] / inten_i[i0])
        scale = 10 ** (gain / 10) / inten[i0]  # dBi per intensity unit (linear): gain(d) = scale * inten(d)
        lv = {}
        for q in range(bd.K):  # level of every beam direction in this state (max within 4 deg)
            c = np.cos(st.TG) * np.cos(K.THETA) + np.sin(st.TG) * np.sin(K.THETA) * np.cos(st.PG - bd.phis[q])
            near = np.degrees(np.arccos(np.clip(c, -1, 1))) <= 4.0
            lv[q] = float(10 * np.log10(scale * inten[near].max() + 1e-30))
        peak = 10 * np.log10(scale * inten.max())
        states.append({"gain_dbi": float(gain), "peak_gain_dbi": float(peak), "peak_theta_deg": mb["peak_theta_deg"], "peak_phi_deg": mb["peak_phi_deg"], "sll_db": mb["SLL_db"], "xpd_db": mb["XPD_db"],
                       "level_at_beam_dbi": lv, "ap_eff": 10 ** ((gain - K.DMAX) / 10)})
    return {"states": states, "m": m, "nd": nd, "ids": ids, "conv": float(conv), "edge": float(edge), "amp": amp, "lam": lam, "g": g, "mu": mu}


def isolation(states):
    """beam-1 level in state 2 and beam-2 level in state 1, in dBi and relative to the state's own beam (dB)."""
    s1, s2 = states
    return {"beam2_in_state1_dbi": s1["level_at_beam_dbi"][1], "beam2_in_state1_rel_db": s1["level_at_beam_dbi"][1] - s1["level_at_beam_dbi"][0],
            "beam1_in_state2_dbi": s2["level_at_beam_dbi"][0], "beam1_in_state2_rel_db": s2["level_at_beam_dbi"][0] - s2["level_at_beam_dbi"][1]}


def summarise(r, tag, bd):
    s = r["states"]
    iso = isolation(s)
    nsw = int(np.sum(r["nd"] > 0))
    return {"design": tag, "diodes": int(r["nd"].sum()), "switched_cells": nsw, "cells": int(bd.n), "static_pct_cells": round(100 * (1 - nsw / bd.n), 1),
            "gain_dbi": [round(x["gain_dbi"], 2) for x in s], "peak_dir": [[round(x["peak_theta_deg"]), round(x["peak_phi_deg"])] for x in s], "sll_db": [round(x["sll_db"], 1) for x in s],
            "xpd_db": [round(x["xpd_db"], 1) for x in s], "ap_eff_pct": [round(100 * x["ap_eff"], 1) for x in s], "isolation": {k: round(v, 1) for k, v in iso.items()},
            "sw_to_beam_pct": round(100 * r["conv"], 1), "edge_pct": round(100 * r["edge"], 1), "lam": r["lam"], "g": r["g"], "mu": r["mu"]}


def choose(rows):
    """best = highest min-state gain among designs with ISOLATION (worse of the two) reported; ties irrelevant."""
    return max(rows, key=lambda t: min(t[0]["gain_dbi"]))


def scan(R, margin, f_hz, name, statics, switched, words_list, lam, g, pols=("LP-H", "LP-V"), epss=range(0, 360, 45), mu=0.0):
    setup(f_hz, R, margin)
    out = []
    for pol in pols:
        bd = K.KBoard([0, 180], [pol, pol])
        opt = options(statics, switched, words_list)
        for e in epss:
            r = run_design(bd, opt, [0.0, np.radians(e)], lam, g, mu=mu)
            sm = summarise(r, name, bd)
            sm.update({"pol": pol, "eps_B_deg": int(e)})
            out.append((sm, r, bd, opt))
            print(name, pol, e, sm["gain_dbi"], sm["sll_db"], sm["isolation"], sm["diodes"], flush=True)
    return out


if __name__ == "__main__":
    LAM, G, MU = 0.28, 2.0, 30.0  # LAM/G: same operating point as the old design (not swept); MU: wrong-beam suppression weight (mu = 0 gave only -4.6 dB isolation, see the trade-off rows)
    stat, pairs, norm = library()
    levels = sorted({round(v[0], 4) for v in stat.values()})
    print("statics", {k: round(v[0], 3) for k, v in stat.items()})
    print("B pairs", {k: tuple(round(x, 3) for x in v) for k, v in pairs.items()}, norm)
    A_sw = [(1.0, -1.0, 1.0)]
    ideal = [(a, min(a + 2.0, 1.0), 1.0) for a in (-1.7, -1.4, -1.0, -0.6, -0.2)]  # hypothetical family-B bound: ON raises m by +2 (mirror of the A swing), capped at +1, OFF anywhere in the static range
    B_lib = [p for k, p in pairs.items() if k[0] == 3.1]  # through-slit cells (the MoM'd B cells whose ON raises Y)
    res = {"norm": norm, "mu": MU, "lam": LAM, "g": G, "static_levels": {f"{k[0]}_{k[1]}": round(v[0], 4) for k, v in stat.items()}, "B_pairs": {f"ls{k[0]}_pos{k[1]}": [round(x, 4) for x in v] for k, v in pairs.items()}, "rows": {}}
    outs = {}
    for name, sw in (("i", A_sw), ("ii", A_sw + B_lib), ("ii_ideal", A_sw + ideal)):
        outs[name] = scan(0.150, MARGIN, 10.5e9, name, levels, sw, [(0, 1)], LAM, G, mu=MU)
        res["rows"][name] = {"all": [t[0] for t in outs[name]], "best": choose(outs[name])[0]}
    outs["i_R120"] = scan(0.120, MARGIN, 10.5e9, "i_R120", levels, A_sw, [(0, 1)], LAM, G, mu=MU)  # smaller board, same concept (cost fallback)
    res["rows"]["i_R120"] = {"all": [t[0] for t in outs["i_R120"]], "best": choose(outs["i_R120"])[0]}
    # trade-off of design (i) with the suppression weight (best pol / phase of the mu = MU scan)
    b = res["rows"]["i"]["best"]
    setup(10.5e9, 0.150, MARGIN)
    bd = K.KBoard([0, 180], [b["pol"], b["pol"]])
    res["mu_tradeoff_i"] = []
    for mu in (0.0, 10.0, 30.0, 100.0, 300.0):
        r = run_design(bd, options(levels, A_sw, [(0, 1)]), [0.0, np.radians(b["eps_B_deg"])], LAM, G, mu=mu)
        sm = summarise(r, "i", bd)
        sm.update({"pol": b["pol"], "eps_B_deg": b["eps_B_deg"]})
        res["mu_tradeoff_i"].append(sm)
        print("mu", mu, sm["gain_dbi"], sm["sll_db"], sm["isolation"], sm["diodes"], flush=True)
    # old two-group design as built (K2_LP): +-1 statics, words (1,0),(0,1), 10.54 GHz, cells by centres inside 150, LP-V, eps_B 180, no suppression term
    setup(10.54e9, 0.150, -1e-9)
    bd = K.KBoard([0, 180], ["LP-V", "LP-V"])
    r = run_design(bd, options([-1.0, 1.0], A_sw, [(1, 0), (0, 1)]), [0.0, np.pi], 0.28, 2.0)
    old = summarise(r, "two-group (as built, K2_LP)", bd)
    old.update({"pol": "LP-V", "eps_B_deg": 180, "freq_ghz": 10.54, "note": "statics +-1 (bit row), words GA/GB, outline by centres, mu 0"})
    res["rows"]["two_group"] = {"best": old}
    print("OLD", old)
    # old two-group topology (two words) at the new settings: 10.5 GHz, analog statics, discs inside 150, same suppression weight
    outs["tg_new"] = scan(0.150, MARGIN, 10.5e9, "two-group @ new settings", levels, A_sw, [(1, 0), (0, 1)], LAM, G, pols=("LP-V",), epss=(0, 180), mu=MU)
    res["rows"]["two_group_new_settings"] = {"all": [t[0] for t in outs["tg_new"]], "best": choose(outs["tg_new"])[0]}
    Path(ROOT / "results/sidequest_hybrid/k2_default_designs.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline="\n")
