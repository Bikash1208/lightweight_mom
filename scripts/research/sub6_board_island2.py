"""SUB6 board-island follow-up (MODEL_ONLY, PENDING_HFSS_CELL_CHECK).  Builds on sub6_board_model4.py (unchanged).

1. ON power budget and the beta-shift diagnosis: per state, fractions of the accepted TM0 power: radiated (model eta) split into the beam-A cone, the beam-B cone (<= 15 deg
   RHCP+LHCP total power) and the rest (sidelobes, cross-pol); unconverted (reaches the aperture edge, 1 - eta); and an ESTIMATE of the sheet-resistive absorption
   alpha_abs = Re(Y_mean) / (d Im Y_total / d beta) of every cell (Y_mean = mean of Y_par, Y_perp of the row in that state, 'mean_mS'), integrated along the rays with the
   leakage-only ray power (first order, NOT subtracted from the gains of the model).  G_A,ON vs G_A,OFF decomposed into the directivity difference and 10 log(eta_on/eta_off);
   G_A with the switched-cell beta shift off (ON amplitudes with the OFF phase) and on.
2. Beam-A design phase: the slot axes of the static cells are designed with P_d = P_OFF + w (P_ON - P_OFF) (w = 0: round 'island'; w = 0.5: averaged beta; w* = the w in
   {0, 0.25, 0.5, 0.75, 1} that maximises min(G_A,OFF, G_A,ON) after the single-direction rotation optimisation).  Switched cells keep the OFF phase (beam B is an OFF-state beam).
3. AR: rotation optimisation (90-degree flips of ALL sites, round-2 method) with the cross-pol (LHCP/RHCP <= -20 dB) required at SEVEN directions per beam (centre + 6 on a
   5-degree ring) instead of the single beam direction, for A in both states and B in OFF.
Variants: V0 = round 'island' (w = 0, single-direction flips); V1 = w = 0.5; V2 = w*; V3 = w* + multi-direction flips.
    python sub6_board_island2.py diag | run [--workers N]
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
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sub6_board_model4 as R4
from src.dispersion import sheet_model as SM
from src.frequency import angular_frequency

B, M = R4.B, R4.M
K0 = R4.K0
RES = R4.RES / "followup"
LAUNCH = R4.LAUNCH
W_GRID = (0.0, 0.25, 0.5, 0.75, 1.0)
RING_DEG, N_RING = 5.0, 6
CONE_DEG = 15.0
BEST4 = [(0.2, 135, "a", 1.0, 0), (0.2, 240, "b", 0.5, 1), (0.22, 135, "b", 0.5, 1), (0.22, 240, "b", 0.8, 0)]  # best per (side, N) of the island sweep (-3 dB margin)
POINTS = [(s, n) for s in (0.2, 0.22) for n in (120, 135, 240)]


# ------------------------------------------------------------------ directions with a design phase
def directions_w(D, ph, w):
    S, kt, mask = D["S"], D["kt"], D["mask"]
    pd = ph["Poff"] + w * (ph["Pon"] - ph["Poff"])
    pdx = np.where(mask, ph["Poff"], pd)
    chiA, chiB = pdx - kt * S.x, pdx + kt * S.x
    uA = np.stack([np.cos(chiA), np.sin(chiA)], 1)
    uB = np.stack([np.cos(chiB), np.sin(chiB)], 1)
    ssum = uA + D["c_w"] * uB
    nrm = np.linalg.norm(ssum, axis=1, keepdims=True)
    shat = np.where(nrm > 1e-9, ssum / np.maximum(nrm, 1e-12), uA)
    return np.where(mask[:, None], uB, shat)


# ------------------------------------------------------------------ budget helpers
def cone_fraction(S, et, ep, beam, cone=CONE_DEG):
    u = np.abs(et) ** 2 + np.abs(ep) ** 2
    w = np.sin(np.radians(S.tg))
    tg, pg = S.tg, S.pg
    cd = np.cos(np.radians(tg)) * math.cos(math.radians(beam[0])) + np.sin(np.radians(tg)) * math.sin(math.radians(beam[0])) * np.cos(np.radians(pg - beam[1]))
    near = np.degrees(np.arccos(np.clip(cd, -1, 1))) <= cone
    return float(np.sum((u * w)[near]) / np.sum(u * w))


def alpha_abs_site(sub, row_state_mean_mS, beta_over_k0):
    """Sheet-resistive attenuation estimate Re(Y_mean) / (d Im Y_tot / d beta) at the root, Np/m."""
    s = SM
    omega = angular_frequency(M.F)
    k0 = K0
    eps_c = sub.eps_complex

    def ytot(b):
        return (_y_air(b, omega, k0) + _y_slab(b, omega, eps_c, k0, sub.thickness)).imag

    b = beta_over_k0 * k0
    d = (ytot(b * (1 + 1e-4)) - ytot(b * (1 - 1e-4))) / (2e-4 * b)
    re = complex(row_state_mean_mS[0], row_state_mean_mS[1]).real * 1e-3
    return abs(re / d) if d != 0 else 0.0


_y_air, _y_slab = SM._y_air, SM._y_slab


def ray_absorbed(m, lk, a_leak, a_abs, b_site, b_fill):
    """Mean over rays of int 2 alpha_abs W drho, with the leakage-only ray power W (the model's march)."""
    psi, rho, drho = m["psi"], m["rho"], m["drho"]
    X, Y = rho[None, :] * np.cos(psi)[:, None], rho[None, :] * np.sin(psi)[:, None]
    idx = lk.at(X, Y)
    inside = rho[None, :] <= m["r_edge"][:, None]
    ok = (idx >= 0) & inside
    a = np.where(ok, a_leak[np.maximum(idx, 0)], 0.0)
    ab = np.where(ok, a_abs[np.maximum(idx, 0)], 0.0)
    ca = np.cumsum(a, axis=1) - a
    w = np.exp(-2.0 * ca * drho)
    return float(np.mean(np.sum(2.0 * ab * w, axis=1) * drho))


def budget(D, launch, w=0.0, flips=None):
    """ON/OFF power budget and gain decomposition at one design."""
    S = D["S"]
    sa = R4.site_arrays(D)
    ph = R4.physics(D, sa)
    direction = directions_w(D, ph, w)
    live, mask = D["live"], D["mask"]
    if flips is None:
        flips, _ = R4.optimise_flips(D, sa, ph, direction)
    fl = flips.astype(float)
    sub = M.B.SUB
    T, sw, st = D["T"], D["swt"], D["T"]["st"]
    rows_st = T["st"]["rows"]
    # per-site sheet absorption (Np/m)
    ab_off, ab_on = np.zeros(S.n), np.zeros(S.n)
    cache = {}
    for i in np.nonzero(live)[0]:
        if mask[i]:
            r = sw["rows"][D["isw"][i]]
            key = ("sw", int(D["isw"][i]))
            if key not in cache:
                cache[key] = (alpha_abs_site(sub, r["states"]["OFF_BAR50"]["mean_mS"], r["states"]["OFF_BAR50"]["beta_over_k0"]),
                              alpha_abs_site(sub, r["states"]["ON_BAR50"]["mean_mS"], r["states"]["ON_BAR50"]["beta_over_k0"]))
            ab_off[i], ab_on[i] = cache[key]
        else:
            r = rows_st[D["ist"][i]]
            key = ("st", int(D["ist"][i]))
            if key not in cache:
                v = alpha_abs_site(sub, r["states"]["NOCUT"]["mean_mS"], r["states"]["NOCUT"]["beta_over_k0"])
                cache[key] = (v, v)
            ab_off[i], ab_on[i] = cache[key]
    out = {}
    for state in ("off", "on"):
        jx, jy, _, _ = R4.currents(D, sa, ph, direction, fl, state)
        et, ep = S.far(jx, jy)
        r = B.metrics(S, et, ep, ph["eta_" + state], launch)
        eta = ph["eta_" + state]
        fa, fb = cone_fraction(S, et, ep, B.BEAM_A), cone_fraction(S, et, ep, B.BEAM_B)
        al = sa["alpha_" + state]
        bt = sa["beta_" + state]
        ab = ab_off if state == "off" else ab_on
        absorbed = ray_absorbed(D["m0"], D["lk"], al, ab, bt, D["b_fill"])
        absorbed_sw = ray_absorbed(D["m0"], D["lk"], al, np.where(mask, ab, 0.0), bt, D["b_fill"])
        # directivity of the A beam (referred to radiated power): G_A = D_A + 10 log10(launch eta)
        dA = r["gA"] - 10 * math.log10(launch * eta)
        out[state] = dict(G_A=r["gA"], G_B=r["gB"], eta=eta, D_A=dA, frac_A_cone=fa * eta, frac_B_cone=fb * eta, frac_other_radiated=(1 - fa - fb) * eta, unconverted=1 - eta,
                          absorbed_est=absorbed, absorbed_switched_est=absorbed_sw, SLL=r["sll"], suppression_or_balance=r["gA"] - r["gB"])
    jx, jy, _, _ = R4.currents(D, sa, dict(ph, Pon=ph["Poff"]), direction, fl, "on")
    rcf = B.metrics(S, *S.far(jx, jy), ph["eta_on"], launch)
    out["G_A_ON_no_beta_shift"] = rcf["gA"]
    out["G_A_ON_with_beta_shift"] = out["on"]["G_A"]
    out["dphi_rms_deg"] = float(np.degrees(np.sqrt(np.sum(((ph["Pon"] - ph["Poff"]) ** 2) * (~mask & live)) / max(np.sum(~mask & live), 1))))
    out["sw_beta_shift_mean"] = float(np.mean((sa["beta_on"] - sa["beta_off"])[mask & live]))
    out["delta_G_ON_minus_OFF"] = out["on"]["G_A"] - out["off"]["G_A"]
    out["delta_D"] = out["on"]["D_A"] - out["off"]["D_A"]
    out["delta_eta_db"] = 10 * math.log10(out["on"]["eta"] / out["off"]["eta"])
    return out


def diag():
    RES.mkdir(parents=True, exist_ok=True)
    T = R4.tables()
    res = []
    for a in BEST4:
        D = R4.design(*a[:5], T)
        b = budget(D, LAUNCH, 0.0)
        res.append(dict(point=dict(side=a[0], N=a[1], pl=a[2], taper=a[3], excl=a[4]), **b))
        o, n = b["off"], b["on"]
        print(a[:2], f"G_A OFF {o['G_A']:.2f} ON {n['G_A']:.2f} d {b['delta_G_ON_minus_OFF']:.2f} = dD {b['delta_D']:.2f} + deta {b['delta_eta_db']:.2f} | eta {o['eta']:.3f}/{n['eta']:.3f} | A cone {o['frac_A_cone']:.3f}/{n['frac_A_cone']:.3f} B cone {o['frac_B_cone']:.3f}/{n['frac_B_cone']:.4f} other {o['frac_other_radiated']:.3f}/{n['frac_other_radiated']:.3f} "
              f"| absorbed est {o['absorbed_est']:.4f}/{n['absorbed_est']:.4f} (switched {o['absorbed_switched_est']:.4f}/{n['absorbed_switched_est']:.4f}) | ON no beta shift {b['G_A_ON_no_beta_shift']:.2f} with {b['G_A_ON_with_beta_shift']:.2f} | dphi rms {b['dphi_rms_deg']:.1f} sw dbeta {b['sw_beta_shift_mean']:.3f}", flush=True)
    json.dump(res, open(RES / "budget_best4.json", "w"), indent=1, default=float)
    return res


# ------------------------------------------------------------------ multi-direction rotation optimiser
def ring_points(beam, rdeg=RING_DEG, n=N_RING):
    th, ph = math.radians(beam[0]), math.radians(beam[1])
    n0 = np.array([math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph), math.cos(th)])
    e1 = np.array([math.cos(th) * math.cos(ph), math.cos(th) * math.sin(ph), -math.sin(th)])
    e2 = np.array([-math.sin(ph), math.cos(ph), 0.0])
    pts = [tuple(beam)]
    for k in range(n):
        al = 2 * math.pi * k / n
        v = n0 * math.cos(math.radians(rdeg)) + math.sin(math.radians(rdeg)) * (math.cos(al) * e1 + math.sin(al) * e2)
        pts.append((math.degrees(math.acos(np.clip(v[2], -1, 1))), math.degrees(math.atan2(v[1], v[0])) % 360.0))
    return pts


def optimise_flips_multi(D, sa, ph, direction, ring_weight=1.0):
    S, live = D["S"], D["live"]
    P = {"A": ring_points(B.BEAM_A), "B": ring_points(B.BEAM_B)}
    X = {}
    for state in ("off", "on"):
        for fl in (0, 1):
            jx, jy, _, _ = R4.currents(D, sa, ph, direction, np.full(S.n, fl, float), state)
            for bn in ("A", "B"):
                for ip, beam in enumerate(P[bn]):
                    R, L = R4.pol_fields(S, beam, jx, jy)
                    X[(state, fl, bn, ip, "R")], X[(state, fl, bn, ip, "L")] = R * live, L * live
    keys = [(s_, b, ip, p) for s_ in ("off", "on") for b in ("A", "B") for ip in range(len(P["A"])) for p in ("R", "L")]
    flips = np.zeros(S.n, int)
    sums = {k: X[(k[0], 0, k[1], k[2], k[3])].sum() for k in keys}
    n0 = {(s_, b): abs(sums[(s_, b, 0, "R")]) ** 2 for s_ in ("off", "on") for b in ("A", "B")}
    xr, flA, flB = 10 ** (-R4.XPD_DB / 10), 10 ** (-R4.FLOOR_A_DB / 10), 10 ** (-R4.FLOOR_B_DB / 10)
    cap_bn = max(10 ** (-R4.SUPP_DB / 10) * n0[("off", "B")], n0[("on", "B")])
    nring = len(P["A"]) - 1

    def cost(v):
        a2 = lambda s_, b, ip, p: np.abs(v[(s_, b, ip, p)]) ** 2  # noqa: E731
        nA, nBo = {s_: n0[(s_, "A")] for s_ in ("off", "on")}, n0[("off", "B")]
        c = -(a2("off", "A", 0, "R") / nA["off"] + a2("on", "A", 0, "R") / nA["on"]) - a2("off", "B", 0, "R") / nBo + R4.MU_B * a2("on", "B", 0, "R") / nBo
        for s_ in ("off", "on"):
            c = c + 50.0 * (np.maximum(a2(s_, "A", 0, "L") - xr * a2(s_, "A", 0, "R"), 0.0) / nA[s_] + np.maximum(flA * nA[s_] - a2(s_, "A", 0, "R"), 0.0) / nA[s_])
            ring = sum(np.maximum(a2(s_, "A", ip, "L") - xr * a2(s_, "A", ip, "R"), 0.0) for ip in range(1, nring + 1)) / (nring * nA[s_])
            c = c + 50.0 * ring_weight * ring
        c = c + 50.0 * (np.maximum(a2("off", "B", 0, "L") - xr * a2("off", "B", 0, "R"), 0.0) / nBo + np.maximum(flB * nBo - a2("off", "B", 0, "R"), 0.0) / nBo)
        ringb = sum(np.maximum(a2("off", "B", ip, "L") - xr * a2("off", "B", ip, "R"), 0.0) for ip in range(1, nring + 1)) / (nring * nBo)
        c = c + 50.0 * ring_weight * ringb
        c = c + 50.0 * np.maximum(a2("on", "B", 0, "R") - cap_bn, 0.0) / nBo
        return c

    cur = float(cost(sums))
    cand = np.nonzero(live)[0]
    for it in range(600):
        v = {}
        for k in keys:
            d = np.where(flips[cand] == 0, X[(k[0], 1, k[1], k[2], k[3])][cand] - X[(k[0], 0, k[1], k[2], k[3])][cand], X[(k[0], 0, k[1], k[2], k[3])][cand] - X[(k[0], 1, k[1], k[2], k[3])][cand])
            v[k] = sums[k] + d
        cs = cost(v)
        if int(flips.sum()) >= R4.MAX_FLIPS:
            cs = np.where(flips[cand] == 1, cs, np.inf)
        j = int(np.argmin(cs))
        if cs[j] >= cur - 1e-9:
            break
        i = cand[j]
        for k in keys:
            sums[k] = v[k][j]
        flips[i] ^= 1
        cur = float(cs[j])
    return flips, dict(flip_iters=it, flip_count=int(flips.sum()))


def metrics_state(D, sa, ph, direction, flips, launch):
    S = D["S"]
    fl = flips.astype(float)
    out = {}
    for state in ("off", "on"):
        jx, jy, _, _ = R4.currents(D, sa, ph, direction, fl, state)
        ff = S.far(jx, jy)
        r = B.metrics(S, *ff, ph["eta_" + state], launch)
        r["ar3A"] = R4.ar_beam(S, *ff, B.BEAM_A)
        r["ar3B"] = R4.ar_beam(S, *ff, B.BEAM_B)
        out[state] = r
    ar3 = max(out["off"]["ar3A"], out["off"]["ar3B"], out["on"]["ar3A"])
    bal = abs(out["off"]["gA"] - out["off"]["gB"])
    sup = out["on"]["gA"] - out["on"]["gB"]
    margin, gate = B.gate(bal, ar3, sup, out["off"]["sll"], out["on"]["sll"])
    return dict(off_GA=out["off"]["gA"], off_GB=out["off"]["gB"], bal=bal, off_ARA3=out["off"]["ar3A"], off_ARB3=out["off"]["ar3B"], on_GA=out["on"]["gA"], supp=sup, on_ARA3=out["on"]["ar3A"],
                sll_off=out["off"]["sll"], sll_on=out["on"]["sll"], eta_off=ph["eta_off"], eta_on=ph["eta_on"], margin=margin, gate=gate, minGA=min(out["off"]["gA"], out["on"]["gA"]))


def point_task(a):
    side, n, pl, taper, excl, launch = a
    T = R4.tables()
    D = R4.design(side, n, pl, taper, excl, T)
    sa = R4.site_arrays(D)
    ph = R4.physics(D, sa)
    res = {}
    dirs = {w: directions_w(D, ph, w) for w in W_GRID}
    single = {}
    for w in W_GRID:
        fl, info = R4.optimise_flips(D, sa, ph, dirs[w])
        single[w] = (fl, metrics_state(D, sa, ph, dirs[w], fl, launch), info)
    res["V0"] = (0.0,) + single[0.0]
    res["V1"] = (0.5,) + single[0.5]
    wstar = max(W_GRID, key=lambda w: single[w][1]["minGA"])
    res["V2"] = (wstar,) + single[wstar]
    best = (wstar,) + single[wstar]  # V3 never worse than V2 by margin
    best = (best[0], best[1], best[2], dict(best[3], ring_weight=0.0))
    for w in sorted({0.0, wstar}):
        for rwt in (0.1, 0.3, 1.0):
            fl3, info3 = optimise_flips_multi(D, sa, ph, dirs[w], ring_weight=rwt)
            m3 = metrics_state(D, sa, ph, dirs[w], fl3, launch)
            info3 = dict(info3, ring_weight=rwt)
            if best is None or m3["margin"] > best[2]["margin"]:
                best = (w, fl3, m3, info3)
    res["V3"] = best
    rows = []
    for v, (w, fl, m, info) in res.items():
        rows.append(dict(variant=v, side_mm=side * 1e3, N_sw=n, placement=pl, taper=taper, excl=excl, w=w, rotated=int(fl.sum()), ring_weight=info.get("ring_weight", 0.0), **m))
    return rows


def run(workers):
    RES.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    tasks = [(s, n, pl, t, ex, LAUNCH) for (s, n) in POINTS for pl in ("a", "b") for t in (0.5, 0.8, 1.0) for ex in (0, 1)]
    allrows = []
    with Pool(min(workers, len(tasks))) as pool:
        for rows in pool.imap_unordered(point_task, tasks, chunksize=1):
            allrows += rows
    allrows.sort(key=lambda r: (-r["side_mm"], r["N_sw"], r["placement"], r["taper"], r["excl"], r["variant"]))
    B.write_csv(RES / "followup_points.csv", allrows)
    print(f"done {time.time() - t0:.0f}s, {len(allrows)} rows", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["diag", "run"])
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args()
    diag() if a.stage == "diag" else run(a.workers)
