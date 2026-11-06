"""Far-field-matrix coordinate-descent design for candidate B on SWHB.  MODEL_ONLY.
J = sum_s [ -beam power + mu_sl * sidelobe power (outside r_ex of every wanted/wrong beam) + mu_w * wrong-beam power ] / Pref + lam * diodes,
exact change of the far field per option (Born model, SW field of the current assignment, iterated)."""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_swhb as S
import sidequest_swhb_design as DS

F0 = DS.F0


def dirgrid(dth=3.0, dph=6.0):
    th = np.arange(1.0, 89.0, dth)
    ph = np.arange(0.0, 360.0, dph)
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    return TH.ravel(), PH.ravel()


def amat(env, xy, p, E, ysh, thd, phd):
    """Per-cell far-field (E_theta, E_phi) per unit susceptance on the directions (deg)."""
    th, ph = np.radians(thd), np.radians(phd)
    kx = env.k0 * np.sin(th) * np.cos(ph)
    ky = env.k0 * np.sin(th) * np.sin(ph)
    kr = np.sin(th)
    ytm = env.yup_tm(kr) + env.ydn_tm(kr) + ysh[0]
    yte = env.yup_te(kr) + env.ydn_te(kr) + ysh[1]
    rho = np.hypot(xy[:, 0], xy[:, 1])
    rh = xy / np.maximum(rho, 1e-9)[:, None]
    jtm = rh[:, 0:1] * np.cos(ph)[None, :] + rh[:, 1:2] * np.sin(ph)[None, :]
    jte = -rh[:, 0:1] * np.sin(ph)[None, :] + rh[:, 1:2] * np.cos(ph)[None, :]
    sinc = np.sinc(kx * p / (2 * math.pi)) * np.sinc(ky * p / (2 * math.pi))
    pref = (E[:, None] * np.exp(1j * (np.outer(xy[:, 0], kx) + np.outer(xy[:, 1], ky)))) * p * p * sinc[None, :]
    Ath = 1j * env.k0 / (2 * math.pi) * (-jtm / ytm[None, :]) * pref * 1j
    Aph = 1j * env.k0 * np.cos(th)[None, :] / (2 * math.pi) * (-jte / yte[None, :]) * pref * 1j
    return Ath, Aph


def design2(lib, R, variant, full=False, max_diodes=514, beams=((25.0, 0.0), (25.0, 180.0)), outer=3, sweeps=5, mu_sl=1.0, mu_w=1.0, w2=1.0, rho_keep=6.5e-3, seed=0, verbose=False, init=None, sq=False):
    env, p = lib.env, lib.p
    xy = DS.cells(R, p, rho_keep)
    n = len(xy)
    nopt = len(lib.opts)
    bs = (lib.y0.imag, lib.y1.imag)
    allowed = (lib.nd > 0) if full else np.ones(nopt, bool)
    rng = np.random.default_rng(seed)
    cur = rng.choice(np.where(allowed)[0], n) if init is None else init.copy()
    thd, phd = dirgrid()
    nd_dir = len(thd)
    ct, st_, cp, sp = np.cos(np.radians(thd)), np.sin(np.radians(thd)), np.cos(np.radians(phd)), np.sin(np.radians(phd))

    def ang(t1, p1):
        return np.degrees(np.arccos(np.clip(st_ * math.sin(math.radians(t1)) * np.cos(np.radians(phd - p1)) + ct * math.cos(math.radians(t1)), -1, 1)))

    hp = 0.886 * (S.C0 / F0) / (2 * R * math.cos(math.radians(beams[0][0]))) * 57.3
    rex = 1.3 * hp + 3.0
    near = [ang(*bm) for bm in beams]
    bidx = [int(np.argmin(near[k])) for k in range(2)]
    wanted = {0: [0], 1: [1]} if variant == "swap" else {0: [0], 1: [0, 1]}
    wrong = {0: [1], 1: [0]} if variant == "swap" else {0: [1], 1: []}
    SL = {}
    for s in range(2):
        m = np.ones(nd_dir, bool)
        for b in set(wanted[s]) | set(wrong[s]):
            m &= near[b] > rex
        SL[s] = m

    def ypol(Ft, Fp, k):
        return Ft[..., k] * ct[k] * sp[k] + Fp[..., k] * cp[k]

    for it in range(outer):
        ybr = [lib.y0[cur], lib.y1[cur]]
        A = []
        for s in range(2):
            r = S.run_state(env, xy, ybr[s], lib.ycom_sw, lib.ycom_rad(), p, 1.0, rho_keep, n_iter=3, want_field=False)
            if not r["ok"]:
                return {"ok": False, "why": f"no TM root state {s} (outer {it})"}
            E = DS.sw_field(env, xy, p, r["beta_k0"] * env.k0, r["Ecell2"])
            A.append(amat(env, xy, p, E, (lib.y_rad_tm + ybr[s].mean(), lib.y_rad_te + ybr[s].mean()), thd, phd))
        PREF = [float(np.sum(np.abs(A[s][0][:, bidx[0]]) ** 2 + np.abs(A[s][1][:, bidx[0]]) ** 2) * np.abs(bs[s]).max() ** 2) / 4 + 1e-30 for s in range(2)]

        def run_cd(lam, sw, cu):
            cu = cu.copy()
            for _ in range(sw):
                F = []
                for s in range(2):
                    bb = bs[s][cu]
                    F.append([A[s][0].T @ (bb - bb.mean()), A[s][1].T @ (bb - bb.mean())])
                changed = 0
                for i in rng.permutation(n):
                    cost = lam * lib.nd.astype(float)
                    for s in range(2):
                        db = (bs[s] - bs[s][cu[i]])[:, None]
                        Ft = F[s][0][None, :] + A[s][0][i][None, :] * db
                        Fp = F[s][1][None, :] + A[s][1][i][None, :] * db
                        u = np.abs(Ft) ** 2 + np.abs(Fp) ** 2
                        bp = 0.0
                        for b in wanted[s]:
                            bp = bp + (w2 if (variant == "add" and b == 1) else 1.0) * np.abs(ypol(Ft, Fp, bidx[b])) ** 2
                        sl = (((u[:, SL[s]] / PREF[s]) ** 2).sum(axis=1) * PREF[s]) if sq else u[:, SL[s]].sum(axis=1)
                        wr = 0.0
                        for b in wrong[s]:
                            wr = wr + np.abs(ypol(Ft, Fp, bidx[b])) ** 2
                        cost = cost + (-bp + mu_sl * sl / nd_dir * 8 + mu_w * wr) / PREF[s]
                    cost = np.where(allowed, cost, np.inf)
                    o = int(np.argmin(cost))
                    if o != cu[i]:
                        for s in range(2):
                            d = bs[s][o] - bs[s][cu[i]]
                            F[s][0] = F[s][0] + A[s][0][i] * d
                            F[s][1] = F[s][1] + A[s][1][i] * d
                        cu[i] = o
                        changed += 1
                if not changed:
                    break
            return cu

        cu = run_cd(0.0, sweeps, cur)
        if lib.nd[cu].sum() > max_diodes and not full:
            lo, hi = 0.0, 2.0
            best = None
            for _ in range(10):
                mid = 0.5 * (lo + hi)
                c2 = run_cd(mid, 2, cur)
                if lib.nd[c2].sum() > max_diodes:
                    lo = mid
                else:
                    hi = mid
                    if best is None or lib.nd[c2].sum() > lib.nd[best].sum():
                        best = c2
            cu = run_cd(hi, sweeps, best if best is not None else cur)
            if lib.nd[cu].sum() > max_diodes:
                cu = best if best is not None else cu
        elif full and lib.nd[cu].sum() > max_diodes:
            return {"ok": False, "why": "full needs more diodes than allowed"}
        cur = cu
        if verbose:
            print("outer", it, "diodes", int(lib.nd[cur].sum()), "of", n, flush=True)
    return {"ok": True, "xy": xy, "cur": cur, "lib": lib, "nd": int(lib.nd[cur].sum()), "n": n, "variant": variant, "beams": beams}
