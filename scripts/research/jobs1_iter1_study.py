"""ITERATION 1 design study (model, relative trends; MODEL_ONLY / NOT_HFSS_VALIDATED).  User 2026-11-10: diode = BAR50-02V, N = 40 / 50 diodes (one diode = one switched element), circular partitioned board, OFF = A + B, ON = A.
    SUB6_MESH=96 python jobs1_iter1_study.py <out.json> [r_on_amplitude]
Combines (i) gathered elements spanning k = 1 / 2 / 3 / 4 / 6 lattice cells per diode (k cells along the element's own slot axis, ONE common slot axis for the element: each member cell's current direction is the reflection of its own
ray about the common axis, chi_i = 2 a_c - phi_i; member amplitude = the single-cell law a = sqrt(2 alpha_cell W / rho); the phase quality is simulated, not assumed), (ii) placement of the element centres (inner / outer / sector_axis (+-x, the beam
axis) / sector_perp (+-y) / field_weighted (equal quantiles of the A-design power, golden-angle scattered) / coherence (best centres by simulated cluster coherence x sqrt(W/rho))), (iii) beam-A strength reduced by capping the static-cell
alpha (the cap that keeps B >= A - 3 dB with the largest A; else the best balance), (iv) board radius 110 / 140 / 180 mm.
Same engine as jobs1_diode_count.py (march_rays / cell_design, layered-medium far field).  alpha_cell = the MoM alpha_OFF of the island cell (argv / default 2.16 Np/m, mesh 96, C 0.15 pF); ON: the switched cells keep r_on of the OFF amplitude
(default sqrt(alpha_ON/alpha_OFF) of the same MoM row) AND deplete less.  Assumptions that need the HFSS cell/9-cell check: per-cell strength and ON suppression of a gathered element equal those of the single cell.
Level = 20 log10 |E_RHCP| relative to the all-static A of the same radius with its best cap; gains G = D x eta_rad (D from the pattern integral theta <= 90, eta_rad = 0.656 x converted power fraction).
"""
import json
import math
import os
import sys
from pathlib import Path

os.environ.setdefault("SUB6_MESH", "96")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1]))
import numpy as np

import sub6_board_model as M
import sub6_board_model4 as R4
from src.hologram.tensor_leaky_cp import SupercellTable, cell_design, march_rays
from src.pattern.vector_aperture_pattern import _vector_far_field, copol_component

B = M.B
K0, PITCH = B.K0, B.PITCH
BEAM_A, BEAM_B = B.BEAM_A, B.BEAM_B
R_EXCL = 0.025
ALPHA_CELL = 2.01
R_ON = math.sqrt(0.226 / 2.01)
CAPS = [0.5, 0.7, 1.0, 1.4, 2.0, 3.0, 4.0]
LAUNCH = 0.656


def circle(r, n=360):
    a = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.stack([r * np.cos(a), r * np.sin(a)], 1)


def sites(R):
    n = int(R / PITCH) + 1
    i, j = np.meshgrid(np.arange(-n, n + 1), np.arange(-n, n + 1), indexing="ij")
    x, y = i.ravel() * PITCH, j.ravel() * PITCH
    ok = (np.hypot(x, y) <= R - PITCH / 2) & (np.hypot(x, y) >= R_EXCL)
    return x[ok], y[ok]


def beam_field(green, x, y, beam, jx, jy):
    th, ph = math.radians(beam[0]), math.radians(beam[1])
    kx, ky = K0 * math.sin(th) * math.cos(ph), K0 * math.sin(th) * math.sin(ph)
    d = np.exp(1j * (kx * x + ky * y)) * PITCH**2
    g = green.spectral_dyad(np.array([kx]), np.array([ky]))[0] * math.cos(th)
    ex = g[0, 0] * d * jx + g[0, 1] * d * jy
    ey = g[1, 0] * d * jx + g[1, 1] * d * jy
    et = ex * math.cos(th) * math.cos(ph) + ey * math.cos(th) * math.sin(ph)
    ep = -ex * math.sin(ph) + ey * math.cos(ph)
    return copol_component(et, ep, "rhcp")


def far_grid(green, x, y, jx, jy, th, ph):
    """E_theta, E_phi on the (th, ph) grid (radians) from the layered-medium dyad (same expressions as diag_hfss_vs_model.model_field)"""
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    kx, ky = (K0 * np.sin(TH) * np.cos(PH)).ravel(), (K0 * np.sin(TH) * np.sin(PH)).ravel()
    g = green.spectral_dyad(kx, ky) * np.cos(TH).ravel()[:, None, None]
    Jx = np.zeros(kx.size, complex); Jy = np.zeros(kx.size, complex)
    for i0 in range(0, kx.size, 2000):
        sl = slice(i0, i0 + 2000)
        d = np.exp(1j * (np.outer(kx[sl], x) + np.outer(ky[sl], y))) * PITCH**2
        Jx[sl], Jy[sl] = d @ jx, d @ jy
    ex = g[:, 0, 0] * Jx + g[:, 0, 1] * Jy
    ey = g[:, 1, 0] * Jx + g[:, 1, 1] * Jy
    c, s_ = np.cos(PH).ravel(), np.sin(PH).ravel()
    ct = np.cos(TH).ravel()
    return (ex * ct * c + ey * ct * s_).reshape(TH.shape), (-ex * s_ + ey * c).reshape(TH.shape)


def pattern_D(green, x, y, jx, jy, beams):
    th = np.radians(np.arange(0, 91, 2.0)); ph = np.radians(np.arange(0, 360, 2.0))
    et, ep = far_grid(green, x, y, jx, jy, th, ph)
    u = np.abs(et) ** 2 + np.abs(ep) ** 2
    w = np.sin(th) * np.where((np.arange(th.size) == 0) | (np.arange(th.size) == th.size - 1), 0.5, 1.0)
    P = float(np.sum(u * w[:, None]) * math.radians(2.0) * math.radians(2.0))
    out = []
    for bm in beams:
        e1, e2 = far_grid(green, x, y, jx, jy, np.array([math.radians(bm[0])]), np.array([math.radians(bm[1])]))
        f = copol_component(e1, e2, "rhcp")[0, 0]
        out.append(10 * math.log10(4 * math.pi * abs(f) ** 2 / P))
    return out, P


class Board:
    def __init__(self, R, green, T, tab, b_fill):
        self.R = R
        self.green = green
        self.x, self.y = sites(R)
        self.rho, self.psi = np.hypot(self.x, self.y), np.arctan2(self.y, self.x)
        th0, ph0 = math.radians(BEAM_A[0]), math.radians(BEAM_A[1])
        thb, phb = math.radians(BEAM_B[0]), math.radians(BEAM_B[1])
        self.d = {}
        for cap in CAPS:
            keep = np.nonzero(T["st"]["a"] <= cap)[0]
            if keep.size < 3:
                keep = np.argsort(T["st"]["a"])[:3]
            t1 = SupercellTable(tuple(tab.slots[i] for i in keep), tab.k_grid, tab.alpha[keep], tab.beta[keep], tab.diameter_mm, tab.pitch_mm)
            m = march_rays(t1, circle(R), (0.0, 0.0), K0, th0, ph0, B.ETA_DESIGN, taper_c=0.5, rho_min_m=B.RHO_CLEAR, beta0_over_k0=b_fill)
            dA = cell_design(m, t1, self.x, self.y, (0.0, 0.0), K0, th0, ph0, hand=+1)
            dB = cell_design(m, t1, self.x, self.y, (0.0, 0.0), K0, thb, phb, hand=+1)
            live = dA["amp"] > 1e-12
            ip = np.clip(np.rint((self.psi - m["psi"][0]) / (m["psi"][1] - m["psi"][0])).astype(int), 0, m["psi"].size - 1)
            ir = np.clip(np.rint(self.rho / m["drho"]).astype(int), 0, m["rho"].size - 1)
            self.d[cap] = dict(live=live, W=m["w"][ip, ir], phi_ref=dB["phi_ref"], chiB=np.arctan2(dB["jy"], dB["jx"].real * 0 + dB["jx"]).real if False else None,
                               ampA=dA["amp"], jxA=dA["jx"] * live, jyA=dA["jy"] * live, cellA=beam_field(green, self.x, self.y, BEAM_A, dA["jx"] * live, dA["jy"] * live),
                               eta=float(m["eta_conv"]), axisB=dB["slot_angle_rad"])
            kdotr = K0 * math.sin(thb) * (self.x * math.cos(phb) + self.y * math.sin(phb))
            self.d[cap]["chiB"] = self.d[cap]["phi_ref"] - kdotr
        ref = max(CAPS, key=lambda c: abs(self.d[c]["cellA"].sum()))
        self.ref_cap = ref
        self.FA_ref = abs(self.d[ref]["cellA"].sum())

    # ------------------------------------------------------------------ elements
    def elements(self, centres, k, cap):
        d = self.d[cap]
        used = np.zeros(self.x.size, bool)
        els = []
        for c in centres:
            if used[c] or not d["live"][c]:
                continue
            a = 0.5 * (d["chiB"][c] + self.psi[c])
            if k == 1:
                idx = np.array([c])
            else:
                u = np.array([math.cos(a), math.sin(a)])
                dx, dy = self.x - self.x[c], self.y - self.y[c]
                t = dx * u[0] + dy * u[1]
                perp = np.abs(-dx * u[1] + dy * u[0])
                cand = np.nonzero((perp <= 0.75 * PITCH) & (np.abs(t) <= (k / 2 + 0.3) * PITCH) & d["live"] & ~used)[0]
                if cand.size < k or c not in cand:
                    continue
                # k contiguous: the k candidates closest to the target window centred so that c is inside
                order = cand[np.argsort(t[cand])]
                pos = int(np.where(order == c)[0][0])
                lo = max(0, min(pos - (k - 1) // 2, order.size - k))
                idx = order[lo:lo + k]
                if idx.size < k or np.any(np.diff(t[idx]) > 1.6 * PITCH):
                    continue
            used[idx] = True
            els.append((c, idx, a))
        return els

    def currents_unit(self, els, cap):
        """unit-amplitude (a = 1) B currents of the member cells with the common axis of each element"""
        d = self.d[cap]
        ux = np.zeros(self.x.size, complex); uy = np.zeros(self.x.size, complex)
        for c, idx, a in els:
            chi = 2 * a - self.psi[idx]
            ph = np.exp(-1j * d["phi_ref"][idx])
            ux[idx], uy[idx] = ph * np.cos(chi), ph * np.sin(chi)
        return ux, uy

    def centres(self, place, n, cap, k):
        d = self.d[cap]
        rho, psi, N = self.rho, self.psi, self.x.size
        live = d["live"]
        if place == "inner":
            return list(np.argsort(rho)[: n * 6])
        if place == "outer":
            return list(np.argsort(-rho)[: n * 6])
        if place in ("sector_axis", "sector_perp"):
            tgt = np.array([0.0, 1.0]) * math.pi if place == "sector_axis" else np.array([0.5, 1.5]) * math.pi
            dd = np.min(np.abs(((psi[:, None] - tgt[None, :]) + math.pi) % (2 * math.pi) - math.pi), axis=1)
            return list(np.argsort(dd + 1e-3 * rho)[: n * 8])
        if place == "field_weighted":
            order = np.argsort(rho)
            cum = np.cumsum(d["ampA"][order] ** 2); cum /= cum[-1]
            out = []
            for q in range(n * 3):
                qq = (q + 0.5) / (n * 3)
                pos = int(np.searchsorted(cum, qq))
                target = (q * 2.399963229728653) % (2 * math.pi)
                cand = [order[min(max(pos + o, 0), N - 1)] for o in range(-8, 9)]
                out.append(min(cand, key=lambda c: abs(((psi[c] % (2 * math.pi)) - target + math.pi) % (2 * math.pi) - math.pi)))
            return out
        if place == "coherence":
            # simulated coherence of a k-cluster at every live site x sqrt(W / rho)
            sc = np.zeros(N)
            for c in np.nonzero(live & (rho > 0.03))[0]:
                els = self.elements([c], k, cap)
                if not els:
                    continue
                ux, uy = self.currents_unit(els, cap)
                f = beam_field(self.green, self.x[els[0][1]], self.y[els[0][1]], BEAM_B, ux[els[0][1]], uy[els[0][1]])
                sc[c] = abs(f.sum()) / max(np.abs(f).sum(), 1e-30) * math.sqrt(d["W"][c] / rho[c])
            return list(np.argsort(-sc)[: n * 8])
        raise ValueError(place)

    # ------------------------------------------------------------------ evaluation
    def evaluate(self, place, n, k, alpha, r_on, cap, need_fields=False):
        d = self.d[cap]
        cen = self.centres(place, n, cap, k)
        els = self.elements(cen, k, cap)[:n]
        if len(els) < n:
            return None
        sw = np.concatenate([e[1] for e in els])
        mask = np.zeros(self.x.size, bool); mask[sw] = True
        ux, uy = self.currents_unit(els, cap)
        res = {}
        for state, al, amp_ratio in (("off", alpha, 1.0), ("on", alpha * r_on**2, r_on)):
            dxm = (np.cos(self.psi)[:, None] * self.y[sw][None, :] - np.sin(self.psi)[:, None] * self.x[sw][None, :])
            ahead = (self.rho[sw][None, :] < self.rho[:, None] - 1e-9) & (np.abs(dxm) <= PITCH / 2)
            dep = np.exp(-2.0 * al * PITCH * ahead.sum(axis=1))
            stat = np.where(mask, 0.0, np.sqrt(dep))
            a = np.sqrt(2.0 * alpha * np.maximum(d["W"], 0.0) * dep / np.maximum(self.rho, 1e-9)) * amp_ratio * mask   # B amplitude (OFF: alpha; ON: amplitude x r_on, depletion with alpha_on)
            jxB, jyB = a * ux, a * uy
            jxA, jyA = d["jxA"] * stat, d["jyA"] * stat
            FA = abs(beam_field(self.green, self.x, self.y, BEAM_A, jxA, jyA).sum())
            FB = abs(beam_field(self.green, self.x, self.y, BEAM_B, jxB, jyB).sum())
            # B-direction leakage of A currents and A-direction of B (negligible) is included in the sums below via the full currents
            Fa_full = abs(beam_field(self.green, self.x, self.y, BEAM_A, jxA + jxB, jyA + jyB).sum())
            Fb_full = abs(beam_field(self.green, self.x, self.y, BEAM_B, jxA + jxB, jyA + jyB).sum())
            p_conv = float((np.sum(d["ampA"] ** 2 * stat**2 * (~mask)) + np.sum(a**2)) * PITCH**2 / (2 * math.pi))
            res[state] = dict(A=20 * math.log10(Fa_full / self.FA_ref), B=20 * math.log10(max(Fb_full, 1e-30) / self.FA_ref), p_conv=p_conv, jx=jxA + jxB, jy=jyA + jyB)
        out = dict(place=place, N=n, k=k, cap=cap, alpha_cell=alpha, A_off=res["off"]["A"], B_off=res["off"]["B"], A_on=res["on"]["A"], B_on=res["on"]["B"],
                   balance_A_minus_B=res["off"]["A"] - res["off"]["B"], suppression_dB=res["off"]["B"] - res["on"]["B"], dA_off_to_on=res["on"]["A"] - res["off"]["A"],
                   p_conv_off=res["off"]["p_conv"], p_conv_on=res["on"]["p_conv"], n_switched_cells=int(mask.sum()), vias=2 * n, diodes=n)
        if need_fields:
            out["_cur"] = {s: (res[s]["jx"], res[s]["jy"]) for s in res}
        return out


def main():
    outp = sys.argv[1]
    r_on = float(sys.argv[2]) if len(sys.argv) > 2 else R_ON
    T = R4.tables(); tab = M.static_table(T)
    b_fill = float(np.mean(T["st"]["b"]))
    green = B.get_side(0.2).green
    radii = [float(v) for v in os.environ.get("JOB_RADII", "0.110,0.140,0.180").split(",")]
    ks = [int(v) for v in os.environ.get("JOB_KS", "1,2,3,4,6").split(",")]
    ns = [int(v) for v in os.environ.get("JOB_NS", "40,50").split(",")]
    places = os.environ.get("JOB_PLACES", "inner,outer,sector_axis,sector_perp,field_weighted,coherence").split(",")
    alphas = [float(v) for v in os.environ.get("JOB_ALPHAS", "2.16").split(",")]
    rows = []
    for R in radii:
        bd = Board(R, green, T, tab, b_fill)
        print("R", R, "sites", bd.x.size, "ref cap", bd.ref_cap, flush=True)
        for n in ns:
            for k in ks:
                for place in places:
                    for alpha in alphas:
                        best = None
                        for cap in CAPS:
                            r = bd.evaluate(place, n, k, alpha, r_on, cap)
                            if r is None:
                                continue
                            ok = abs(r["balance_A_minus_B"]) <= 3.0
                            key = (ok, min(r["A_off"], r["B_off"]) if ok else -abs(r["balance_A_minus_B"]))
                            if best is None or key > best[0]:
                                best = (key, r)
                        if best is None:
                            continue
                        r = best[1]; r["R_mm"] = R * 1e3; r["balanced_within_3dB"] = bool(best[0][0])
                        rows.append(r)
                print(f"R {R} N {n} k {k}: " + " ".join(f"{r['place'][:6]}:{r['A_off']:.1f}/{r['B_off']:.1f}{'*' if r['balanced_within_3dB'] else ''}" for r in rows if r["R_mm"] == R * 1e3 and r["N"] == n and r["k"] == k and r["alpha_cell"] == alphas[0]), flush=True)
    Path(outp).write_text(json.dumps(dict(label="MODEL_ONLY / NOT_HFSS_VALIDATED iteration 1 study", r_on_amplitude=r_on, alpha_cell=alphas, rows=rows), indent=1, default=float), encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
