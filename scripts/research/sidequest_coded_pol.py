"""K = 2 polarisation-switched coded hologram: state A RHCP at (25, 0), state B LHCP at (25, 180). Whole 150 mm hexagon, centre probe, 1-bit cells (MA4AGP907).
Route R1: ONE diode per cell on a fixed, ROTATED axis chi_i (needs the Phase 28 rotated-load fix). Route R2: TWO diodes per cell on the x and y axes (0/90 deg only).
Cell: split disc, diameter 6.0 mm, diode gap 0.4 mm, footprint 1.2 mm (sidequest_coded_pol_cell.py, MoM mesh 32; |dY| ~4.9 mS, y_perp static ~5.1j, near-isotropic mean sheet).
Model: unit-sign cell current J = (pi/2) A(rho) (u . rho_hat) u e^{-j beta rho}, A = sqrt(2 alpha0/rho) e^{-(alpha_b + alpha_a)(rho - rho_keep)}; |J|^2 ~ cos^2(delta), delta = chi - azimuth of the cell
(alpha_b = alpha0 <cos^2>, global mean over the chosen axes). Far field: E_theta = cos(theta) (J . rho_hat_beam), E_phi = J . phi_hat (same scalar element factor as sidequest_arch_study); co-pol = conj(Jones) . E.
Per cell, option (chi, s_A, s_B): coherent contribution s_k Re(e^{j eps_k} c_k(chi)); a penalty lambda on s_A != s_B trades diodes for gain. MODEL_ONLY: first-order leakage, scalar element factor, static texture neglected."""
import json
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np

import sidequest_arch_study as st
import sidequest_coded_cell as cc
import sidequest_coded_model as M
import sidequest_gate_cell as gc

THETA = np.radians(25.0)
PHIS = np.radians([0.0, 180.0])
POLS = ("RHCP", "LHCP")
# disc 6.0 mm, mesh 32, mean over chi = 0 / 22.5 / 45 (coded_pol_cell.json)
Y_U_OFF = 2.22e-3j
Y_U_ON = (0.61 + 7.15j) * 1e-3
Y_V = 5.07e-3j
DELTA = Y_U_ON - Y_U_OFF
Y_ISO = 0.5 * (0.5 * (Y_U_ON + Y_U_OFF) + Y_V)
BETA_K0, A0, AABS0 = cc.alpha_pair(Y_ISO - DELTA / 2, Y_ISO + DELTA / 2)
BETA = BETA_K0 * st.K0
LAUNCH = M.LAUNCH
CHIS = np.radians(np.arange(0, 180, 22.5))


class PolBoard:
    def __init__(self, chis, fc):
        self.poly, self.x, self.y, _ = st.aperture("hex")
        keep = np.hypot(self.x, self.y) >= st.RHO_KEEP
        self.x, self.y = self.x[keep], self.y[keep]
        self.n = self.x.size
        self.rho = np.hypot(self.x, self.y)
        self.psi = np.arctan2(self.y, self.x)
        self.chis = np.asarray(chis)
        self.fc = fc
        self.amp = np.sqrt(2 * A0 / self.rho) * np.exp(-(A0 * fc + AABS0 * fc) * (self.rho - st.RHO_KEEP))
        # unit-sign contributions c[k, a, i] (co-pol Jones projection of the cell current, with the aperture phase to beam k)
        self.c = np.zeros((2, len(chis), self.n), complex)
        self.jx = np.zeros((len(chis), self.n), complex)
        self.jy = np.zeros((len(chis), self.n), complex)
        for a, ch in enumerate(chis):
            d = np.cos(ch - self.psi)
            j = (np.pi / 2) * self.amp * d * np.exp(-1j * BETA * self.rho)
            self.jx[a], self.jy[a] = j * np.cos(ch), j * np.sin(ch)
        for k in range(2):
            t = st.JONES[POLS[k]]
            ph = PHIS[k]
            arr = np.exp(1j * st.K0 * np.sin(THETA) * (self.x * np.cos(ph) + self.y * np.sin(ph)))
            for a in range(len(chis)):
                et = np.cos(THETA) * (np.cos(ph) * self.jx[a] + np.sin(ph) * self.jy[a])
                ep = -np.sin(ph) * self.jx[a] + np.cos(ph) * self.jy[a]
                self.c[k, a] = (np.conj(t[0]) * et + np.conj(t[1]) * ep) * arr

    def mean_cos2(self, chi_idx):
        return float(np.mean(np.cos(self.chis[chi_idx] - self.psi) ** 2))


def _elem(ra, rb, lam, scale):
    """Best (sA, sB) per cell for contributions ra, rb (n,): maximise sA ra + sB rb - lam*scale*[sA != sB]. Returns (score, sA, sB)."""
    cand = [(abs(ra) + abs(rb), np.sign(ra + (ra == 0)), np.sign(rb + (rb == 0)))]
    opts = []
    for sa, sb in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
        opts.append((sa * ra + sb * rb - (lam * scale if sa != sb else 0.0), sa, sb))
    sc = np.array([o[0] for o in opts])
    o = np.argmax(sc, axis=0)
    sa = np.array([1, 1, -1, -1])[o]
    sb = np.array([1, -1, 1, -1])[o]
    return np.max(sc, axis=0), sa, sb


def design(B, route, eps, lam):
    """R1: per cell choose one axis chi (and signs). R2: x and y elements decided independently. Returns a dict with per-cell element choices."""
    r = np.array([np.real(np.exp(1j * eps[k]) * B.c[k]) for k in range(2)])  # (2, nchi, n)
    scale = np.abs(r).max(axis=(0, 1))
    if route == "R2":
        out = [_elem(r[0, e], r[1, e], lam, scale) for e in range(2)]
        return {"idx": np.array([[0] * B.n, [1] * B.n]), "sa": np.array([o[1] for o in out]), "sb": np.array([o[2] for o in out])}
    sc, sa, sb = zip(*[_elem(r[0, a], r[1, a], lam, scale) for a in range(len(B.chis))])
    sc = np.array(sc)
    a = np.argmax(sc, axis=0)
    ii = np.arange(B.n)
    return {"idx": a[None, :], "sa": np.array(sa)[a, ii][None, :], "sb": np.array(sb)[a, ii][None, :]}


def realise(B, d):
    cur = []
    for sk in ("sa", "sb"):
        jx = np.zeros(B.n, complex)
        jy = np.zeros(B.n, complex)
        for e in range(d["idx"].shape[0]):
            jx += d[sk][e] * B.jx[d["idx"][e], np.arange(B.n)]
            jy += d[sk][e] * B.jy[d["idx"][e], np.arange(B.n)]
        cur.append((jx, jy))
    return cur


def amplitude(B, cur):
    return [abs(np.sum(np.exp(1j * st.K0 * np.sin(THETA) * (B.x * np.cos(PHIS[k]) + B.y * np.sin(PHIS[k]))) * (np.conj(st.JONES[POLS[k]][0]) * np.cos(THETA) * (np.cos(PHIS[k]) * cur[k][0] + np.sin(PHIS[k]) * cur[k][1]) + np.conj(st.JONES[POLS[k]][1]) * (-np.sin(PHIS[k]) * cur[k][0] + np.cos(PHIS[k]) * cur[k][1])))) for k in range(2)]


def currents(B, a, sa, sb):
    ii = np.arange(B.n)
    return [(s * B.jx[a, ii], s * B.jy[a, ii]) for s in (sa, sb)]


def summarize(d):
    sa, sb = d["sa"].ravel(), d["sb"].ravel()
    bits = np.array([sa, sb])
    static = sa == sb
    nonstatic = ~static
    words = {tuple(b) if b[0] > 0 else tuple(-v for v in b) for b in bits[:, nonstatic].T}
    return {"static_frac": float(static.mean()), "diodes": int(nonstatic.sum()), "lines_with_inverters": len(words), "lines_no_inverters": len({tuple(b) for b in bits[:, nonstatic].T})}


def pattern_metrics(B, cur, fc_ideal, label=""):
    res = []
    for k in range(2):
        et, ep = st.far_field(B.x, B.y, *cur[k])
        t = st.JONES[POLS[k]]
        j = st.jones_to_current(t, THETA, PHIS[k])
        amp_i = np.sqrt(2 * A0 * fc_ideal / B.rho) * np.exp(-(A0 + AABS0) * fc_ideal * (B.rho - st.RHO_KEEP))
        idl = amp_i * np.exp(-1j * st.K0 * np.sin(THETA) * (B.x * np.cos(PHIS[k]) + B.y * np.sin(PHIS[k])))
        ei, epi = st.far_field(B.x, B.y, j[0] * idl, j[1] * idl)
        mi = st.metrics(ei, epi, THETA, PHIS[k], POLS[k])
        mb = st.metrics(et, ep, THETA, PHIS[k], POLS[k])
        i0 = mi["_i"]
        inten_i = (np.abs(ei) ** 2 + np.abs(epi) ** 2)[i0]
        inten_b = (np.abs(et) ** 2 + np.abs(ep) ** 2)[i0]
        eta = LAUNCH * float(np.sum(2 * A0 * fc_ideal * np.exp(-2 * (A0 + AABS0) * fc_ideal * (B.rho - st.RHO_KEEP)) / B.rho * st.PITCH**2 / (2 * np.pi)))
        g_ideal = mi["D_dbi"] + 10 * np.log10(eta)
        gain = g_ideal + 10 * np.log10(inten_b / inten_i)
        res.append({"state": "AB"[k], "pol": POLS[k], "phi_deg": float(np.degrees(PHIS[k])), "peak_theta_deg": mb["peak_theta_deg"], "peak_phi_deg": mb["peak_phi_deg"],
                    "gain_dbi": float(gain), "G_ideal_dbi": float(g_ideal), "SLL_db": mb["SLL_db"], "AR_db": mb["AR_db"], "XPD_db": mb["XPD_db"], "converted_beam_power": eta})
    return res


def converge(route, eps, lam, iters=3):
    chis = CHIS if route == "R1" else np.radians([0.0, 90.0])
    fc = 0.5 if route == "R1" else 1.0
    for _ in range(iters if route == "R1" else 1):
        B = PolBoard(chis, fc)
        d = design(B, route, eps, lam)
        if route == "R1":
            fc = float(np.mean(np.cos(B.chis[d["idx"][0]] - B.psi) ** 2))
    return B, d, fc


def run(route):
    out = {"route": route}
    grid = np.radians(np.arange(0, 360, 22.5))
    B0, d0, fc0 = converge(route, [0.0, 0.0], 0.0)
    best_amp = np.array(amplitude(B0, realise(B0, d0)))
    out["lambda0_amp_rel"] = 1.0
    rows = []
    for lam in (0.0, 0.1, 0.2, 0.3, 0.5, 0.8, 1.2):
        for eb in grid:
            B, d, fc = converge(route, [0.0, eb], lam, iters=1 if route == "R1" else 1)
            amp = np.array(amplitude(B, realise(B, d)))
            sm = summarize(d)
            rows.append({"lam": lam, "eps_B_deg": float(np.degrees(eb)), "loss_db": float(-20 * np.log10(amp.min() / best_amp.min())), **sm, "fc": fc})
    out["scan"] = rows
    return out, (B0, d0, fc0), rows


if False:
    for route in ("R1", "R2"):
        o, base, rows = run(route)
        for lam in sorted({r["lam"] for r in rows}):
            rr = [r for r in rows if r["lam"] == lam]
            r = max(rr, key=lambda r: (r["static_frac"] if r["loss_db"] < 1.0 else -1))
            print(route, lam, "static %.3f diodes %d lines %d/%d loss %.2f eps_B %.0f fc %.2f" % (r["static_frac"], r["diodes"], r["lines_with_inverters"], r["lines_no_inverters"], r["loss_db"], r["eps_B_deg"], r["fc"]), flush=True)


def full(route, lam, eps_b_deg, label):
    B, d, fc = converge(route, [0.0, np.radians(eps_b_deg)], lam)
    cur = realise(B, d)
    res = pattern_metrics(B, cur, fc)
    sm = summarize(d)
    return {"label": label, "route": route, "lam": lam, "eps_B_deg": eps_b_deg, "fc": fc, **sm, "states": res, "axis_hist_deg": (np.bincount(d["idx"][0], minlength=len(B.chis)).tolist() if route == "R1" else None)}


def main():
    out = {"cell": {"beta_k0": BETA_K0, "alpha0": A0, "alpha_abs0": AABS0, "delta_mS": abs(DELTA) * 1e3, "y_v_mS": 5.07, "note": "split disc 6.0 mm, gap 0.4, footprint 1.2, MoM mesh 32, mean of chi 0/22.5/45"}, "routes": {}}
    for route in ("R1", "R2"):
        o, base, rows = run(route)
        sel = []
        for lam in (0.0, 0.3, 0.5, 0.8, 1.2):
            rr = [r for r in rows if r["lam"] == lam and r["loss_db"] < 1.0]
            if not rr:
                continue
            r = max(rr, key=lambda r: r["static_frac"])
            f = full(route, lam, r["eps_B_deg"], f"lam{lam}")
            sel.append(f)
            g = [x["gain_dbi"] for x in f["states"]]
            print(route, lam, "eps_B %.1f static %.3f diodes %d lines %d/%d gain %.1f/%.1f SLL %.1f/%.1f AR %.1f/%.1f XPD %.1f/%.1f" % (r["eps_B_deg"], f["static_frac"], f["diodes"], f["lines_with_inverters"], f["lines_no_inverters"], g[0], g[1], f["states"][0]["SLL_db"], f["states"][1]["SLL_db"], f["states"][0]["AR_db"], f["states"][1]["AR_db"], f["states"][0]["XPD_db"], f["states"][1]["XPD_db"]), flush=True)
        out["routes"][route] = {"scan": o["scan"], "designs": sel}
    Path("results/sidequest_hybrid/coded_pol_results.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8", newline="\n")
if __name__ == "__main__":
    main()
