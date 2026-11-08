"""SUB6 board model with ISLAND switched cells (MODEL_ONLY, PENDING_HFSS_CELL_CHECK).  Derived from sub6_board_model2.py (same engine / method: Phase 26 march_rays + cell_design
unmodified, real per-state alpha / beta along the rays, joint 90-degree slot rotations, taper, ring exclusion, placements a / b).  Changes:
* STATIC = NOCUT rows (rows_RT5880.jsonl, as before).  SWITCHED = island rows of results/sub6_island_cell (input copy `input_island_rows.json`, commit 6d53878), only rows
  WITHOUT the resonance flag, and (my filter) alpha_OFF >= 0.9 Np/m and sqrt(alpha_ON/alpha_OFF) <= 0.5 so that they are near-silent in ON.
* Architecture: static cells radiate beam A in BOTH states, island cells radiate beam B in OFF and little in ON.  General (round-1) weights with
  rho_d = mean signed amplitude ratio J_ON / J_OFF of the switched table, c = -rho_d / (1 - rho_d):
      static   J = (u_A + c u_B) / (1 - f)          switched J_OFF = u_B / ((1 - rho_d) f),   J_ON = rho_i J_OFF
  so beam B has coherent sum 1 in OFF and (ideally) 0 in ON.  rho_i = sign(Im Delta_ON / Im Delta_OFF) sqrt(alpha_ON / alpha_OFF) of the chosen row (the magnitude from alpha,
  the sign from the reactance polarity: my modelling choice), with the row's own beta in each state.
* No A_perp (island cells are x-mirror symmetric: |A_perp| < 1.5e-10 mS).  Launch fraction 0.656 (standard SMA pin of sub6-launcher).
    python sub6_board_model4.py run [--workers N]
"""

import argparse
import csv
import json
import math
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sub6_board_model as M
from src.hologram.tensor_leaky_cp import SupercellTable, cell_design, march_rays
from src.pattern.vector_aperture_pattern import copol_component

B = M.B
K0 = M.K0
RES = M.RES.parent / "sub6_board_island"
N_SW = (110, 120, 135, 180, 240, 300)
SIDES = (0.18, 0.2, 0.22)
LAUNCH = 0.656
ISL_JSON = M.RES.parent / "sub6_board_island" / "input_island_rows.json"
TAPERS = (0.5, 0.8, 1.0)
LS_ROWS = ["RT5880_p11.6_g0.73_q0.8_W2.9", "RT5880_p11.6_g1_q0.9_W2.9", "RT5880_p11.6_g1_q0.9_W2.4", "RT5880_p11.6_g0.8_q0.8_W3.6", "RT5880_p11.6_g0.73_q0.8_W3.6"]
XPD_DB, FLOOR_A_DB, FLOOR_B_DB, MU_B, MAX_FLIPS, SUPP_DB = 20.0, 0.5, 1.0, 1.0, 200, 17.0


def cplx(v):
    return complex(v[0], v[1]) * 1e-3


def tables():
    T = M.tables()
    J = json.load(open(ISL_JSON))
    rows = [r for r in J["rows"] if not r["resonance_flag"]]
    sel = []
    for r in rows:
        o, n = r["states"]["OFF_BAR50"], r["states"]["ON_BAR50"]
        ratio = math.sqrt(n["alpha_np_per_m"] / o["alpha_np_per_m"])
        if o["alpha_np_per_m"] >= 0.9 and ratio <= 0.5:
            sel.append(r)
    z = lambda v: complex(v[0], v[1])  # noqa: E731
    st = lambda r, s_, k: r["states"][s_][k]  # noqa: E731
    sgn = [1.0 if (z(st(r, "ON_BAR50", "Delta_mS")).imag / z(st(r, "OFF_BAR50", "Delta_mS")).imag) > 0 else -1.0 for r in sel]
    T["sw"] = dict(rows=sel, id=[r["row"]["id"] for r in sel], aoff=np.array([st(r, "OFF_BAR50", "alpha_np_per_m") for r in sel]), aon=np.array([st(r, "ON_BAR50", "alpha_np_per_m") for r in sel]),
                   boff=np.array([st(r, "OFF_BAR50", "beta_over_k0") for r in sel]), bon=np.array([st(r, "ON_BAR50", "beta_over_k0") for r in sel]),
                   Doff=np.array([z(st(r, "OFF_BAR50", "Delta_mS")) * 1e-3 for r in sel]), Don=np.array([z(st(r, "ON_BAR50", "Delta_mS")) * 1e-3 for r in sel]),
                   Apoff=np.zeros(len(sel), complex), Apon=np.zeros(len(sel), complex), Aqoff=np.zeros(len(sel), complex), Aqon=np.zeros(len(sel), complex))
    T["sw"]["ratio"] = np.array(sgn) * np.sqrt(T["sw"]["aon"] / T["sw"]["aoff"])
    T["n_island_total"], T["n_island_unflagged"] = len(J["rows"]), len(rows)
    return T


def tables_ls(T, c_pf):
    """Switched table restricted to the 5 ls_tol rows; OFF alpha / beta / Delta at the diode C of the case (A_perp nominal)."""
    import csv as _csv
    sw = T["sw"]
    idx = [sw["id"].index(i) for i in LS_ROWS if i in sw["id"]]
    out = dict(rows=[sw["rows"][i] for i in idx], id=[sw["id"][i] for i in idx])
    for k in ("aon", "bon", "Don", "Apon", "Aqon", "Apoff", "Aqoff"):
        out[k] = sw[k][idx]
    aoff, boff, doff = [], [], []
    allr = list(_csv.DictReader(open(RES / "input_ls_tol.csv")))
    for i in idx:
        rid = sw["id"][i]
        r = [x for x in allr if x["row"] == rid and x["case"] == f"OFF_C{c_pf:g}"][0]
        aoff.append(float(r["alpha_Np_per_m"]))
        boff.append(float(r["beta_over_k0"]))
        doff.append(complex(float(r["Delta_re_mS"]), float(r["Delta_im_mS"])) * 1e-3)
    out["aoff"], out["boff"], out["Doff"] = np.array(aoff), np.array(boff), np.array(doff)
    return out


def region_masks(S, lk):
    ring = np.maximum(np.abs(S.ij[:, 0]), np.abs(S.ij[:, 1]))
    outer = np.zeros(S.n, bool)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di or dj:
                outer |= lk.t[S.ij[:, 0] + di + lk.n, S.ij[:, 1] + dj + lk.n] < 0
    return ring, outer


def sites_uniform_allowed(S, n_b, allowed):
    g1, g2 = 0.7548776662466927, 0.5698402909980532
    sc = (S.ij[:, 0] * g1 + S.ij[:, 1] * g2) % 1.0
    sc = np.where(allowed, sc, 9.0)
    mask = np.zeros(S.n, bool)
    mask[np.argsort(sc)[:n_b]] = True
    return mask


def optimise_sites_allowed(S, dA, dB, f, rho_d, n_b, allowed, iters=80):
    """Copy of sub6_board_req.optimise_sites (unmodified algorithm) with the switched sites restricted to `allowed`."""
    c = -rho_d / (1.0 - rho_d)
    ax, ay = B.unit_fields(dA)
    bx, by = B.unit_fields(dB)
    amp = dA["amp"]
    aA = S.cell_field(B.BEAM_A, amp * ax / (1 - f), amp * ay / (1 - f))
    aBst = S.cell_field(B.BEAM_B, amp * c * bx / (1 - f), amp * c * by / (1 - f))
    aBsw = S.cell_field(B.BEAM_B, amp * bx / ((1 - rho_d) * f), amp * by / ((1 - rho_d) * f))
    mask = sites_uniform_allowed(S, n_b, allowed)

    def fields(mk):
        return aA[~mk].sum(), aBst[~mk].sum() + aBsw[mk].sum(), aBst[~mk].sum() + rho_d * aBsw[mk].sum()

    def obj(mk):
        FA, FB, R = fields(mk)
        return math.log(abs(FA)) + math.log(abs(FB)) - math.log(abs(R) ** 2 + 1e-3 * abs(FB) ** 2)

    cur = obj(mask)
    k = max(1, int(0.2 * n_b))
    for _ in range(iters):
        FA, FB, R = fields(mask)
        eps = 1e-3 * abs(FB) ** 2
        dFA, dFB, dR = -aA, aBsw - aBst, rho_d * aBsw - aBst
        score = (np.real(dFA * np.conj(FA)) / abs(FA) ** 2 + np.real(dFB * np.conj(FB)) / abs(FB) ** 2 - 2 * np.real(dR * np.conj(R)) / (abs(R) ** 2 + eps))
        st, sw = np.nonzero(~mask & allowed)[0], np.nonzero(mask)[0]
        done = False
        while k >= 1:
            new = mask.copy()
            new[st[np.argsort(-score[st])[:k]]] = True
            new[sw[np.argsort(score[sw])[:k]]] = False
            v = obj(new)
            if v > cur + 1e-12:
                mask, cur, done = new, v, True
                break
            k //= 2
        if not done:
            break
        k = max(1, 2 * k)
    return mask


def ar_beam(S, et, ep, beam, cone=12.0):
    """Worst AR (power formula of sub6_board_req.metrics) over the directions within the -3 dB RHCP contour of the beam (cone <= 12 deg around the design direction)."""
    rh = np.abs(copol_component(et, ep, "rhcp")) ** 2
    lh = np.abs(copol_component(et, ep, "lhcp")) ** 2
    tg, pg = S.tg, S.pg
    cd = np.cos(np.radians(tg)) * math.cos(math.radians(beam[0])) + np.sin(np.radians(tg)) * math.sin(math.radians(beam[0])) * np.cos(np.radians(pg - beam[1]))
    near = np.degrees(np.arccos(np.clip(cd, -1, 1))) <= cone
    pk = np.max(np.where(near, rh, 0.0))
    sel = near & (rh >= 0.5 * pk)
    a, b = np.sqrt(rh[sel]), np.sqrt(lh[sel])
    return float(np.max(20 * np.log10((a + b) / np.maximum(np.abs(a - b), 1e-12))))


def pol_fields(S, beam, jx, jy):
    """Per-cell complex (RHCP, LHCP) contribution at one direction."""
    import math as m
    th, ph = m.radians(beam[0]), m.radians(beam[1])
    kx, ky = K0 * m.sin(th) * m.cos(ph), K0 * m.sin(th) * m.sin(ph)
    d = np.exp(1j * (kx * S.x + ky * S.y)) * M.PITCH**2
    g = S.green.spectral_dyad(np.array([kx]), np.array([ky]))[0] * m.cos(th)
    ex = g[0, 0] * d * jx + g[0, 1] * d * jy
    ey = g[1, 0] * d * jx + g[1, 1] * d * jy
    et = ex * m.cos(th) * m.cos(ph) + ey * m.cos(th) * m.sin(ph)
    ep = -ex * m.sin(ph) + ey * m.cos(ph)
    return copol_component(et, ep, "rhcp"), copol_component(et, ep, "lhcp")


def design(side, n_sw, place, taper, excl, T, sw_rows=None):
    """Everything that is fixed by the build: sites, rows, slot axes (OFF phase), no flips yet."""
    S = B.get_side(side)
    lk = M.Lookup(S)
    sw, st = (sw_rows if sw_rows is not None else T["sw"]), T["st"]
    tab = M.static_table(T)
    f = n_sw / S.n
    b_fill = float(np.mean(st["b"]))
    kt = K0 * math.sin(math.radians(B.BEAM_A[0]))
    rho_d = float(np.mean(sw["ratio"]))  # mean signed amplitude ratio J_ON / J_OFF of the (filtered) switched table
    c_w = -rho_d / (1.0 - rho_d)
    p_static = (1.0 + c_w**2 + 2.0 * c_w * np.cos(2.0 * kt * S.x)) / (1 - f) ** 2  # |u_A + c u_B|^2 / (1-f)^2
    p_sw = 1.0 / ((1.0 - rho_d) ** 2 * f**2)
    ring, outer = region_masks(S, lk)
    allowed = (ring > 2) & ~outer if excl else np.ones(S.n, bool)
    if place == "a":
        mask0 = sites_uniform_allowed(S, n_sw, allowed)
    else:
        mi, dAi, dBi = B.march_for(side, 2.0)
        mask0 = optimise_sites_allowed(S, dAi, dBi, f, rho_d, n_sw, allowed)

    def base(table):
        m0 = march_rays(table, B.hexagon(side), (0.0, 0.0), K0, math.radians(B.BEAM_A[0]), math.radians(B.BEAM_A[1]), B.ETA_DESIGN, taper_c=taper, rho_min_m=B.RHO_CLEAR, beta0_over_k0=b_fill)
        d0 = cell_design(m0, table, S.x, S.y, (0.0, 0.0), K0, math.radians(B.BEAM_A[0]), math.radians(B.BEAM_A[1]), hand=+1)
        live = d0["amp"] > 1e-12
        mask = mask0 & live
        kap = float(np.sum(np.where(mask, 0.0, d0["amp"] ** 2 * p_static)) / max(np.sum(d0["amp"] ** 2), 1e-30))
        return m0, d0, live, mask, kap

    m0, dA0, live, mask, kappa = base(tab)
    dem = np.where(mask, p_sw, p_static) / kappa
    q90 = float(np.quantile(dem[live], 0.9))
    cap1 = min(float(st["a"].max()), float(sw["aoff"].max())) / max(q90, 1.0)
    keep = np.nonzero(st["a"] <= cap1)[0]
    if keep.size < 3:
        keep = np.argsort(st["a"])[:3]
    tab1 = SupercellTable(tuple(tab.slots[i] for i in keep), tab.k_grid, tab.alpha[keep], tab.beta[keep], tab.diameter_mm, tab.pitch_mm)
    m0, dA0, live, mask, kappa = base(tab1)
    alpha0 = dA0["alpha"]
    sc2 = 1.0 / kappa
    a_req_st = alpha0 * sc2 * p_static
    a_req_sw = alpha0 * sc2 * p_sw
    ist = M.nearest_sqrt(st["a"], a_req_st)
    isw = M.nearest_sqrt(sw["aoff"], a_req_sw)
    D = dict(S=S, lk=lk, m0=m0, mask=mask, live=live, ist=ist, isw=isw, a_req_st=a_req_st, a_req_sw=a_req_sw, alpha0=alpha0, b_fill=b_fill, kt=kt, f=f, q90=q90, cap1=cap1,
             base_rows=int(keep.size), rho_d=rho_d, c_w=c_w, n_excluded_allowed=int(np.sum(~allowed)), T=T, swt=sw, taper=taper, excl=excl, place=place, n_sw=n_sw, side=side)
    return D


def site_arrays(D, over=None):
    """alpha / beta / Delta / A per site and state from the chosen rows (`over`: dict of replacement switched-table arrays for the OFF state)."""
    T, sw, st, mask, ist, isw = D["T"], D["swt"], D["T"]["st"], D["mask"], D["ist"], D["isw"]
    aoff_t = (over or sw)["aoff"]
    boff_t = (over or sw)["boff"]
    doff_t = (over or sw)["Doff"]
    out = dict(
        alpha_off=np.where(mask, aoff_t[isw], st["a"][ist]), alpha_on=np.where(mask, sw["aon"][isw], st["a"][ist]),
        beta_off=np.where(mask, boff_t[isw], st["b"][ist]), beta_on=np.where(mask, sw["bon"][isw], st["b"][ist]),
        D_off=np.where(mask, doff_t[isw], 1.0 + 0j), D_on=np.where(mask, sw["Don"][isw], 1.0 + 0j),
        Ap_off=np.where(mask, sw["Apoff"][isw], 0j), Ap_on=np.where(mask, sw["Apon"][isw], 0j),
        Aq_off=np.where(mask, sw["Aqoff"][isw], 0j), Aq_on=np.where(mask, sw["Aqon"][isw], 0j))
    out["sgn_on"] = np.where(mask, np.sign(sw["ratio"][isw]), 1.0)
    live = D["live"]
    out["alpha_off"], out["alpha_on"] = np.where(live, out["alpha_off"], 0.0), np.where(live, out["alpha_on"], 0.0)
    return out


def physics(D, sa):
    """Real marches (per-state phase and weight) and the base directions/axes (designed with the OFF phase)."""
    S, lk, m0 = D["S"], D["lk"], D["m0"]
    smp = M.sampler(m0, S)
    b_fill = D["b_fill"]
    w_off, ph_off, eta_off = M.real_march(m0, lk, sa["alpha_off"], sa["beta_off"], b_fill)
    w_on, ph_on, eta_on = M.real_march(m0, lk, sa["alpha_on"], sa["beta_on"], b_fill)
    return dict(Woff=smp(w_off), Won=smp(w_on), Poff=smp(ph_off), Pon=smp(ph_on), eta_off=eta_off, eta_on=eta_on)


def directions(D, Poff):
    S, kt, mask = D["S"], D["kt"], D["mask"]
    chiA, chiB = Poff - kt * S.x, Poff + kt * S.x
    uA = np.stack([np.cos(chiA), np.sin(chiA)], 1)
    uB = np.stack([np.cos(chiB), np.sin(chiB)], 1)
    ssum = uA + D["c_w"] * uB
    nrm = np.linalg.norm(ssum, axis=1, keepdims=True)
    shat = np.where(nrm > 1e-9, ssum / np.maximum(nrm, 1e-12), uA)
    return np.where(mask[:, None], uB, shat)


def currents(D, sa, ph, direction, flip, state, with_a=True):
    """Site currents of one state (flip: array of 0/1, 1 = slot rotated by 90 deg)."""
    S, mask = D["S"], D["mask"]
    rho = np.hypot(S.x, S.y)
    phiR = np.arctan2(S.y, S.x)
    alpha, W, P = (sa["alpha_off"], ph["Woff"], ph["Poff"]) if state == "off" else (sa["alpha_on"], ph["Won"], ph["Pon"])
    a = np.sqrt(np.maximum(2.0 * alpha * np.maximum(W, 0.0) / np.maximum(rho, 1e-9), 0.0))
    sgn = (1.0 - 2.0 * flip) * (sa["sgn_on"] if state == "on" else 1.0)
    chi = np.arctan2(direction[:, 1], direction[:, 0])
    axis = 0.5 * (chi + phiR) + flip * (math.pi / 2)  # a = (chi + phi_R) / 2, rotated by 90 deg when flipped
    th = phiR - axis
    Dl = sa["D_" + state]
    Aeff = sa["Ap_" + state] * np.sin(th) ** 2 + sa["Aq_" + state] * np.cos(th) ** 2
    r = np.where(mask & with_a, Aeff / (0.5 * Dl), 0j)
    eA = np.stack([np.sin(phiR), -np.cos(phiR)], 1)
    jx = a * np.exp(-1j * P) * (sgn * direction[:, 0] + r * eA[:, 0])
    jy = a * np.exp(-1j * P) * (sgn * direction[:, 1] + r * eA[:, 1])
    return jx, jy, a, r


def optimise_flips(D, sa, ph, direction):
    S = D["S"]
    live = D["live"]
    X = {}
    for state in ("off", "on"):
        for fl in (0, 1):
            jx, jy, _, _ = currents(D, sa, ph, direction, np.full(S.n, fl, float), state)
            for bn, beam in (("A", B.BEAM_A), ("B", B.BEAM_B)):
                R, L = pol_fields(S, beam, jx, jy)
                X[(state, fl, bn, "R")], X[(state, fl, bn, "L")] = R * live, L * live
    flips = np.zeros(S.n, int)
    keys = [(s, b, p) for s in ("off", "on") for b in ("A", "B") for p in ("R", "L")]
    sums = {k: X[(k[0], 0, k[1], k[2])].sum() for k in keys}
    n0 = {(s, b): abs(sums[(s, b, "R")]) ** 2 for s in ("off", "on") for b in ("A", "B")}
    xr, flA, flB = 10 ** (-XPD_DB / 10), 10 ** (-FLOOR_A_DB / 10), 10 ** (-FLOOR_B_DB / 10)
    cap_bn = max(10 ** (-SUPP_DB / 10) * n0[("off", "B")], n0[("on", "B")])

    def cost(v):
        nA = {s: n0[(s, "A")] for s in ("off", "on")}
        nBo = n0[("off", "B")]
        ra = {s: np.abs(v[(s, "A", "R")]) ** 2 for s in ("off", "on")}
        la = {s: np.abs(v[(s, "A", "L")]) ** 2 for s in ("off", "on")}
        rbo, rbn = np.abs(v[("off", "B", "R")]) ** 2, np.abs(v[("on", "B", "R")]) ** 2
        lbo = np.abs(v[("off", "B", "L")]) ** 2
        c = -(ra["off"] / nA["off"] + ra["on"] / nA["on"]) - rbo / nBo + MU_B * rbn / nBo
        for s in ("off", "on"):
            c = c + 50.0 * (np.maximum(la[s] - xr * ra[s], 0.0) / nA[s] + np.maximum(flA * nA[s] - ra[s], 0.0) / nA[s])
        c = c + 50.0 * (np.maximum(lbo - xr * rbo, 0.0) / nBo + np.maximum(flB * nBo - rbo, 0.0) / nBo)
        c = c + 50.0 * np.maximum(rbn - cap_bn, 0.0) / nBo  # beam-B residual in ON stays at/below the unflipped level (floor 17 dB below B OFF)
        return c

    cur = float(cost(sums))
    start = cur
    cand = np.nonzero(live)[0]
    for it in range(600):
        v = {}
        for k in keys:
            d = np.where(flips[cand] == 0, X[(k[0], 1, k[1], k[2])][cand] - X[(k[0], 0, k[1], k[2])][cand], X[(k[0], 0, k[1], k[2])][cand] - X[(k[0], 1, k[1], k[2])][cand])
            v[k] = sums[k] + d
        cs = cost(v)
        if int(flips.sum()) >= MAX_FLIPS:
            cs = np.where(flips[cand] == 1, cs, np.inf)
        j = int(np.argmin(cs))
        if cs[j] >= cur - 1e-9:
            break
        i = cand[j]
        for k in keys:
            sums[k] = v[k][j]
        flips[i] ^= 1
        cur = float(cs[j])
    return flips, dict(flip_iters=it, flip_count=int(flips.sum()), flip_static=int(np.sum(flips & ~D["mask"])), flip_switched=int(np.sum(flips & D["mask"])), cost_start=start, cost_end=cur)


def evaluate(D, launch, over=None, flips=None, do_flips=True):
    S = D["S"]
    sa = site_arrays(D, over)
    ph = physics(D, sa)
    direction = directions(D, ph["Poff"])
    live, mask = D["live"], D["mask"]
    z = np.zeros(S.n)
    res = {}

    def met(flip, state, with_a=True, phase_off=False, far=False):
        jx, jy, _, _ = currents(D, sa, dict(ph, Pon=ph["Poff"]) if phase_off else ph, direction, flip, state, with_a)
        ff = S.far(jx, jy)
        r = B.metrics(S, *ff, ph["eta_" + state], launch)
        return (r, ff) if far else r

    # baselines: round-1 style (no A term), with A term unflipped, with A term + joint flips
    r_off0, r_on0 = met(z, "off", False), met(z, "on", False)
    r_off1, r_on1 = met(z, "off"), met(z, "on")
    info = {}
    if flips is None:
        flips, info = optimise_flips(D, sa, ph, direction) if do_flips else (np.zeros(S.n, int), {})
    fl = flips.astype(float)
    (r_off, ff_off), (r_on, ff_on) = met(fl, "off", far=True), met(fl, "on", far=True)
    ar3 = dict(off_AR_A_3dB=ar_beam(S, *ff_off, B.BEAM_A), off_AR_B_3dB=ar_beam(S, *ff_off, B.BEAM_B), on_AR_A_3dB=ar_beam(S, *ff_on, B.BEAM_A))
    r_cf = met(fl, "on", True, True)
    stat = ~mask & live
    dphi = np.degrees(ph["Pon"] - ph["Poff"])
    _, a0, _ = None, None, None
    jx0, jy0, a_off, r_off_ratio = currents(D, sa, ph, direction, fl, "off")
    _, _, a_on, r_on_ratio = currents(D, sa, ph, direction, fl, "on")
    wgt = (a_off**2) * stat
    sw_i = mask & live
    row = dict(side_mm=D["side"] * 1e3, N_cells=S.n, N_sw=D["n_sw"], f_sw=D["f"], placement=D["place"], taper=D["taper"], excl=int(D["excl"]),
               off_GA=r_off["gA"], off_GB=r_off["gB"], off_balance=abs(r_off["gA"] - r_off["gB"]), off_AR_A=r_off["arA"], off_AR_B=r_off["arB"],
               on_GA=r_on["gA"], on_GB=r_on["gB"], on_suppression=r_on["gA"] - r_on["gB"], on_AR_A=r_on["arA"], on_AR_B=r_on["arB"],
               sll_off=r_off["sll"], sll_on=r_on["sll"], eta_off=ph["eta_off"], eta_on=ph["eta_on"],
               noA_on_AR_A=r_on0["arA"], noA_sll_on=r_on0["sll"], noA_off_AR_A=r_off0["arA"], noA_on_supp=r_on0["gA"] - r_on0["gB"],
               preflip_on_AR_A=r_on1["arA"], preflip_sll_on=r_on1["sll"], preflip_off_AR_A=r_off1["arA"], preflip_off_AR_B=r_off1["arB"], preflip_on_supp=r_on1["gA"] - r_on1["gB"],
               preflip_off_bal=abs(r_off1["gA"] - r_off1["gB"]),
               on_AR_A_no_beta_change=r_cf["arA"], sll_on_no_beta_change=r_cf["sll"], A_loss_from_beta_db=r_cf["gA"] - r_on["gA"],
               dphi_rms_deg=float(np.sqrt(np.sum(wgt * dphi**2) / max(np.sum(wgt), 1e-30))), dphi_max_deg=float(np.max(np.abs(dphi[stat]))),
               sw_beta_on_minus_off=float(np.mean((sa["beta_on"] - sa["beta_off"])[sw_i])), sw_alpha_off_mean=float(np.mean(sa["alpha_off"][sw_i])), sw_alpha_on_mean=float(np.mean(sa["alpha_on"][sw_i])),
               rA_off_abs_mean=float(np.mean(np.abs(r_off_ratio)[sw_i])), rA_on_abs_mean=float(np.mean(np.abs(r_on_ratio)[sw_i])),
               **ar3, sw_alpha_saturated=int(np.sum(D["a_req_sw"][sw_i] > D["swt"]["aoff"].max())), st_alpha_saturated=int(np.sum(D["a_req_st"][stat] > D["T"]["st"]["a"].max())),
               q90_demand=D["q90"], rho_d=D["rho_d"], c_w=D["c_w"], ratio_sel_mean=float(np.mean(sa["alpha_on"][sw_i] ** 0.5 / np.maximum(sa["alpha_off"][sw_i] ** 0.5, 1e-9))),  cap1_Np_m=D["cap1"], alpha0_mean=float(np.mean(D["alpha0"][live])), n_not_allowed=D["n_excluded_allowed"], **info)
    ar = max(row["off_AR_A"], row["off_AR_B"], row["on_AR_A"])
    row["worst_margin_db"], row["gate"] = B.gate(row["off_balance"], ar, row["on_suppression"], r_off["sll"], r_on["sll"])
    row["worst_margin_3dB_beam_db"], row["gate_3dB_beam_AR"] = B.gate(row["off_balance"], max(row["off_AR_A_3dB"], row["off_AR_B_3dB"], row["on_AR_A_3dB"]), row["on_suppression"], r_off["sll"], r_on["sll"])
    return row, dict(flips=flips, sa=sa, ph=ph, direction=direction)


def layout_rows(D, ex):
    S, sa, flips, direction = D["S"], ex["sa"], ex["flips"], ex["direction"]
    sw, st, mask, live = D["swt"], D["T"]["st"], D["mask"], D["live"]
    phiR = np.arctan2(S.y, S.x)
    axis = np.mod(0.5 * (np.arctan2(direction[:, 1], direction[:, 0]) + phiR) + flips * math.pi / 2, math.pi)
    out = []
    for i in range(S.n):
        if live[i]:
            out.append(dict(site=i, x_mm=S.x[i] * 1e3, y_mm=S.y[i] * 1e3, cls="SWITCHED" if mask[i] else "STATIC", row_id=(sw["id"][D["isw"][i]] if mask[i] else st["id"][D["ist"][i]]),
                            axis_deg=math.degrees(axis[i]), rotated90=int(flips[i]), alpha_off=sa["alpha_off"][i], alpha_on=sa["alpha_on"][i], beta_off=sa["beta_off"][i], beta_on=sa["beta_on"][i]))
    return out


def task(a):
    side, n_sw, place, taper, excl, launch = a
    T = tables()
    D = design(side, n_sw, place, taper, excl, T)
    row, ex = evaluate(D, launch)
    if place in ("a", "b") and taper == 0.8 and excl:
        lp = RES / f"layout_{int(round(side * 1e3))}mm_N{n_sw}_{place}_t{taper:g}.csv"
        lr = layout_rows(D, ex)
        with open(lp, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(lr[0].keys()))
            w.writeheader()
            for r in lr:
                w.writerow({k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in r.items()})
    return row


def run(workers):
    RES.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    lf, lfd = LAUNCH, dict(note="standard SMA pin of sub6-launcher (5475c42), thin-wire estimate")
    tasks = [(s, n, pl, t, ex, lf) for s in SIDES for n in N_SW for pl in ("a", "b") for t in TAPERS for ex in (0, 1)]
    T = tables()
    print(f"{len(tasks)} tasks, {workers} workers, launch fraction {lf:.4f}; static rows {len(T['st']['a'])}, switched rows {len(T['sw']['aon'])}", flush=True)
    with Pool(min(workers, len(tasks))) as pool:
        rows = []
        for r in pool.imap_unordered(task, tasks, chunksize=1):
            rows.append(r)
            print(f"  {r['side_mm']:.0f} N{r['N_sw']} {r['placement']} t{r['taper']} ex{r['excl']} {r['gate']} {r['worst_margin_db']:.2f} flips {r.get('flip_count')}", flush=True)
    rows.sort(key=lambda r: (-r["side_mm"], r["N_sw"], r["placement"], r["taper"], r["excl"]))
    B.write_csv(RES / "points_island.csv", rows)
    prov = B.provenance(dict(stage="run_island", workers=min(workers, len(tasks)), n_tasks=len(tasks), seconds=time.time() - t0, launch_fraction=lf, launch_detail=lfd,
                             static_rows=len(T["st"]["a"]), switched_rows=len(T["sw"]["aon"]), pitch_m=M.PITCH, xpd_db=XPD_DB, floor_a_db=FLOOR_A_DB, floor_b_db=FLOOR_B_DB, mu_b=MU_B, max_flips=MAX_FLIPS, supp_db=SUPP_DB))
    (RES / "provenance_run_island.json").write_text(json.dumps(prov, indent=1, default=str), encoding="utf-8", newline="\n")
    print(f"done {time.time() - t0:.0f} s", flush=True)


def tol_task(a):
    best, c_pf, launch = a
    T = tables()
    T5 = tables_ls(T, c_pf)
    T0 = tables_ls(T, 0.15)
    D = design(best["side"], best["n_sw"], best["place"], best["taper"], best["excl"], T, sw_rows=T0)
    # nominal flips for C=0.15 and the same flips for the tolerance cases
    row0, ex0 = evaluate(D, launch, over=T0)
    row, ex = evaluate(D, launch, over=T5, flips=ex0["flips"])
    row["C_pF"] = c_pf
    return row


def tol(workers):
    pts = list(csv.DictReader(open(RES / "points_island.csv")))
    b = max(pts, key=lambda r: float(r["worst_margin_3dB_beam_db"]))
    best = dict(side=float(b["side_mm"]) / 1e3, n_sw=int(float(b["N_sw"])), place=b["placement"], taper=float(b["taper"]), excl=int(float(b["excl"])))
    lf, _ = B.launch_fraction()
    with Pool(3) as pool:
        rows = pool.map(tol_task, [(best, c, lf) for c in (0.12, 0.15, 0.20)])
    B.write_csv(RES / "tolerance_C.csv", rows)
    (RES / "tolerance_best_point.json").write_text(json.dumps(dict(best=best, point=b), indent=1), encoding="utf-8", newline="\n")
    for r in rows:
        print(r["C_pF"], r["gate"], round(r["worst_margin_db"], 2), "supp", round(r["on_suppression"], 1), "balance", round(r["off_balance"], 2), "AR off", round(r["off_AR_A"], 2), round(r["off_AR_B"], 2), "AR on", round(r["on_AR_A"], 2))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["run", "tol"])
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args()
    run(a.workers) if a.stage == "run" else tol(a.workers)
