"""Board design with candidate B on SWHB.  MODEL_ONLY (SWHB far field NOT validated against the HFSS set: see SIDEQUEST_M6_VALIDATION_SET.md).
Per cell option o: susceptance b_o^s in state s (OFF/ON), diode count nd.  Targets: holographic superposition,  b_target,i = kappa t_i sum_b w_b Re[j c_i^b e^{-j eps_b}]/|c_i^b|,  cost = sum_s |c_i|^2 (b_o^s - b_target^s)^2 + lam nd + wrong-beam penalty.
c_i^{s,b}: far-field copol (y-polarised) per unit susceptance at the beam direction, including the local SW field of state s (iterated with the SW march)."""
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_cellB_core as K
import sidequest_cellB_model as M
import sidequest_swhb as S

RES = ROOT / "results/sidequest_hybrid"
F0 = 10.5e9
W0 = 2 * math.pi * F0
K0 = W0 / S.C0


class CellLib:
    """Candidate B cell: common patch-grid sheet (MoM) + per-cell via/diode branch (lumped)."""

    def __init__(self, patch_rows, er, tand, h, p, L_ser=1.6e-9, diode="MA4AGP907_measured", dR=1.0, dC=1.0, dL=1.0, C_ap=0.0, ls_levels=None, polB=True):
        self.env = S.Env(er, tand, h, F0)
        self.p, self.h = p, h
        self.patch = M.Patch(patch_rows)
        rad = sorted(patch_rows, key=lambda r: r["kr"])[0]
        self.y_rad_tm = complex(*rad["Yuu_S"])
        self.y_rad_te = complex(rad["Y_full_S"][1][1][0], rad["Y_full_S"][1][1][1])
        kw = dict(diode=diode, L_ser=L_ser, dR=dR, dC=dC, dL=dL, C_ap=C_ap)
        yon = M.branch("ON", h, **kw)
        yoff = M.branch("OFF", h, **kw)
        self.opts = []  # (name, y_state0, y_state1, nd)
        self.opts.append(("SW_A(ON lowers)", yoff, yon, 1))
        if polB:
            # concept cell: permanent inductive via (L_p) + diode in series with a MIM capacitor: ON raises Y by about the same amount
            Lp = (M.MU0 * h + 1.6e-9 + 3e-9)
            R, C = M.DIODE[diode]
            yb_off = 1.0 / (1j * W0 * Lp) + 1j * W0 * C * dC * 0.0
            dy = (yoff - yon)  # A-pair contrast (positive imaginary)
            yb_on = yb_off + abs(dy.imag) * 1j * 1.0 + (yon.real if False else 0.0)
            self.opts.append(("SW_B(ON raises)", yb_off, yb_on + 0.2e-3, 1))
        levels = ls_levels if ls_levels is not None else [0.0] + list(np.geomspace(1.0e-9, 24e-9, 14))
        for Ls in levels:
            self.opts.append((f"ST_L{Ls * 1e9:.1f}", 1.0 / (1j * W0 * (M.MU0 * h + Ls)), 1.0 / (1j * W0 * (M.MU0 * h + Ls)), 0))
        self.y0 = np.array([o[1] for o in self.opts])
        self.y1 = np.array([o[2] for o in self.opts])
        self.nd = np.array([o[3] for o in self.opts])

    def ycom_sw(self, kr):
        return self.patch.y_sw(kr)

    def ycom_rad(self):
        return (self.y_rad_tm, self.y_rad_te)


def cells(R, p, rho_keep=6.5e-3, hexagon=False):
    n = int(R / p) + 2
    k = np.arange(-n, n + 1) * p
    xg, yg = np.meshgrid(k, k)
    x, y = xg.ravel(), yg.ravel()
    rho = np.hypot(x, y)
    ok = (rho <= R - p / 2) & (rho >= rho_keep)
    return np.stack([x[ok], y[ok]], axis=1)


def coupling(env, xy, p, E, beam, ysh_tm, ysh_te, cop="y"):
    """c_i: copol (y) far-field per unit susceptance-j at the beam direction, for the SW field E (complex, n). beam=(theta, phi) deg."""
    th, ph = math.radians(beam[0]), math.radians(beam[1])
    kx, ky = env.k0 * math.sin(th) * math.cos(ph), env.k0 * math.sin(th) * math.sin(ph)
    kr = math.sin(th)
    ytm = env.yup_tm(kr) + env.ydn_tm(kr) + ysh_tm
    yte = env.yup_te(kr) + env.ydn_te(kr) + ysh_te
    rho = np.hypot(xy[:, 0], xy[:, 1])
    rh = xy / np.maximum(rho, 1e-9)[:, None]
    jtm = rh[:, 0] * math.cos(ph) + rh[:, 1] * math.sin(ph)
    jte = -rh[:, 0] * math.sin(ph) + rh[:, 1] * math.cos(ph)
    eth = 1j * env.k0 / (2 * math.pi) * (-jtm / ytm)
    eph = 1j * env.k0 * math.cos(th) / (2 * math.pi) * (-jte / yte)
    # y-polarised component: theta_hat.y = cos(th) sin(ph), phi_hat.y = cos(ph)
    cy = eth * math.cos(th) * math.sin(ph) + eph * math.cos(ph)
    arr = np.exp(1j * (kx * xy[:, 0] + ky * xy[:, 1]))
    sinc = np.sinc(kx * p / (2 * math.pi)) * np.sinc(ky * p / (2 * math.pi))
    return cy * E * arr * p * p * sinc


def sw_field(env, xy, p, beta, amp2, rho_keep=6.5e-3):
    from scipy.special import hankel2

    rho = np.hypot(xy[:, 0], xy[:, 1])
    a = hankel2(1, beta * np.maximum(rho, 1e-4))
    return np.sqrt(amp2) * a / np.abs(a)


def design(lib, R, variant, full=False, max_diodes=514, beams=((25.0, 0.0), (25.0, 180.0)), iters=3, sweeps=6, mu=60.0, w2=1.0, eps_b=0.0, rho_keep=6.5e-3, seed=0, verbose=False):
    env, p = lib.env, lib.p
    xy = cells(R, p, rho_keep)
    n = len(xy)
    nopt = len(lib.opts)
    rho = np.hypot(xy[:, 0], xy[:, 1])
    ysh_tm, ysh_te = lib.y_rad_tm, lib.y_rad_te
    b0, b1 = lib.y0.imag, lib.y1.imag  # option susceptances (S) per state
    allowed = np.ones(nopt, bool)
    if full:
        allowed = lib.nd > 0
    cur = np.array([int(np.argmin(np.abs(b0 + b1))) for _ in range(n)])
    taper = np.ones(n)
    # initial assignment: pseudo-random switched
    rng = np.random.default_rng(seed)
    cur = rng.choice(np.where(allowed)[0], n)
    wanted = {0: [0], 1: [1]} if variant == "swap" else {0: [0], 1: [0, 1]}
    wrong = [(0, 1)] if variant == "swap" else []
    if variant == "swap":
        wrong.append((1, 0))
    bmax = max(np.abs(b0).max(), np.abs(b1).max())
    res = None
    for it in range(iters):
        ybr = [lib.y0[cur], lib.y1[cur]]
        runs = []
        E = []
        for s in range(2):
            r = S.run_state(env, xy, ybr[s], lib.ycom_sw, lib.ycom_rad(), p, 1.0, rho_keep, n_iter=3, want_field=False)
            if not r["ok"]:
                return {"ok": False, "why": f"no TM root state {s}"}
            runs.append(r)
            E.append(sw_field(env, xy, p, r["beta_k0"] * env.k0, r["Ecell2"]))
        C = {}
        for s in range(2):
            ysh = (lib.y_rad_tm + ybr[s].mean(), lib.y_rad_te + ybr[s].mean())
            for b in range(2):
                C[(s, b)] = coupling(env, xy, p, E[s], beams[b], *ysh)
        # target susceptance per state
        tgt = np.zeros((2, n))
        for s, bl in wanted.items():
            for b in bl:
                wb = (w2 if b == 1 and variant == "add" else 1.0)
                ph = np.exp(-1j * (eps_b if b == 1 else 0.0))
                a = np.real(1j * C[(s, b)] * ph) / np.maximum(np.abs(C[(s, b)]), 1e-30)
                tgt[s] += wb * a
        tgt *= 0.8 * bmax / max(np.abs(tgt).max(), 1e-30) * taper[None, :]
        wgt = [np.abs(C[(s, b)]) ** 2 for s in range(2) for b in range(2)]
        cw = np.zeros((2, n))
        for s in range(2):
            cw[s] = sum(np.abs(C[(s, b)]) ** 2 for b in wanted[s]) + 1e-30
        cw /= cw.mean()
        own = cw[0][:, None] * (b0[None, :] - tgt[0][:, None]) ** 2 + cw[1][:, None] * (b1[None, :] - tgt[1][:, None]) ** 2
        own = own / (bmax ** 2)
        own = np.where(allowed[None, :], own, np.inf)
        # coordinate descent with the diode budget
        def run_cd(lam):
            cu = np.argmin(own + lam * lib.nd[None, :], axis=1)
            Wsum = {sb: np.sum(C[sb] * 1j * (b0 if sb[0] == 0 else b1)[cu]) for sb in wrong}
            norm = np.sum(np.abs(C[(0, 0)]) ** 2) * bmax ** 2 + 1e-30
            for _ in range(sweeps):
                changed = 0
                for i in rng.permutation(n):
                    cost = own[i] + lam * lib.nd
                    for sb in wrong:
                        bsel = b0 if sb[0] == 0 else b1
                        Wi = Wsum[sb] - C[sb][i] * 1j * bsel[cu[i]]
                        cost = cost + mu * np.abs(Wi + C[sb][i] * 1j * bsel) ** 2 / norm * n
                    o = int(np.argmin(cost))
                    if o != cu[i]:
                        for sb in wrong:
                            bsel = b0 if sb[0] == 0 else b1
                            Wsum[sb] += C[sb][i] * 1j * (bsel[o] - bsel[cu[i]])
                        cu[i] = o
                        changed += 1
                if not changed:
                    break
            return cu
        lam = 0.0
        cu = run_cd(0.0)
        if lib.nd[cu].sum() > max_diodes:
            lo, hi = 0.0, 5.0
            for _ in range(12):
                mid = 0.5 * (lo + hi)
                cu = run_cd(mid)
                if lib.nd[cu].sum() > max_diodes:
                    lo = mid
                else:
                    hi = mid
            cu = run_cd(hi)
        cur = cu
        # taper feedback: uniform radiated density
        ybr = [lib.y0[cur], lib.y1[cur]]
        r0 = S.run_state(env, xy, ybr[0], lib.ycom_sw, lib.ycom_rad(), p, 1.0, rho_keep, n_iter=2, want_field=False)
        if verbose:
            print("iter", it, "diodes", int(lib.nd[cur].sum()), "edge", r0["P_edge"], flush=True)
        # push the illumination toward uniform: scale taper with ring radiated density
        q = np.abs(ybr[0] - ybr[0].mean()) ** 2 * r0["Ecell2"]
        bins = np.minimum((rho / (R / 8)).astype(int), 7)
        dens = np.array([q[bins == k].sum() / max((bins == k).sum(), 1) for k in range(8)])
        tgtd = dens[dens > 0].mean()
        fac = np.where(dens > 0, (tgtd / np.maximum(dens, 1e-30)) ** 0.25, 1.0)
        taper = np.clip(taper * fac[bins], 0.3, 1.0)
    return {"ok": True, "xy": xy, "cur": cur, "lib": lib, "nd": int(lib.nd[cur].sum()), "n": n, "variant": variant, "beams": beams}


def evaluate(des, eta_l=1.0, grid=(1.0, 3.0), rho_keep=6.5e-3):
    lib, xy, cur = des["lib"], des["xy"], des["cur"]
    env, p = lib.env, lib.p
    out = []
    for s in range(2):
        y = (lib.y0 if s == 0 else lib.y1)[cur]
        r = S.run_state(env, xy, y, lib.ycom_sw, lib.ycom_rad(), p, eta_l, rho_keep, n_iter=6, grid=grid)
        out.append(r)
    return out
