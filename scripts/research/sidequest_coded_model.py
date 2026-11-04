"""Grouped coded hologram, model level. Whole 150 mm hexagon = 1-bit coded hologram (one MA4AGP907 per switched cell), single centre probe, K fixed states.
Cell pair (split patch 6.0 x 4.8 mm, diode gap 0.4 mm, footprint 1.2 mm, current along its axis; measured MA4AGP907): Y_OFF = 2.90j, Y_ON = 1.15 + 10.41j mS (mesh 32).
Sheet along the ray: m = cos^2(angle between ray and diode axis) scales the modulation depth, hence alpha_beam ~ m^2 and the ohmic alpha_abs ~ m (g_eff = g m). First-order alpha from the repo
`leakage_alpha_first_order` (theta 25 deg, zeta from the beam/ray azimuths). Currents J_i = (pi/2) a_i s_i exp(-j beta rho_i) pol, s_i = +-1 the diode state, a_i = sqrt(2 alpha/rho) exp(-(alpha+alpha_abs)(rho - rho_keep)).
Gain = G_ideal x |peak amplitude|^2 ratio (G_ideal = D_ideal x launch (0.627) x converted beam power). Polarisation: scalar TM cell => LP (theta_hat), ideal by construction. MODEL_ONLY."""
import json
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np

import sidequest_arch_study as st
import sidequest_coded_cell as cc
import sidequest_gate_cell as gc
from src.dispersion.leakage_model import leakage_alpha_first_order

Y_OFF = 2.90e-3j
Y_ON = (1.15 + 10.41j) * 1e-3
LAUNCH = 0.627
RHO_KEEP = st.RHO_KEEP
THETA = np.radians(25.0)

BETA_K0, A_BEAM0, A_ABS0 = cc.alpha_pair(Y_OFF, Y_ON)  # ray along the axis, zeta = 0
BETA = BETA_K0 * st.K0
_ZETA = np.radians(np.linspace(0, 90, 10))
_A_ZETA = np.array([cc.alpha_pair(Y_OFF, Y_ON, 25.0, z)[1] for z in np.degrees(_ZETA)])


def alpha_beam(zeta):
    z = np.abs(np.arctan2(np.sin(zeta), np.cos(zeta)))
    z = np.where(z > np.pi / 2, np.pi - z, z)
    return np.interp(z, _ZETA, _A_ZETA)


class Board:
    def __init__(self, beams_phi_deg, axis_mode="x", pitch=st.PITCH):
        self.poly, self.x, self.y, _ = st.aperture("hex")
        keep = np.hypot(self.x, self.y) >= RHO_KEEP
        self.x, self.y = self.x[keep], self.y[keep]
        self.n = self.x.size
        self.rho = np.hypot(self.x, self.y)
        self.psi = np.arctan2(self.y, self.x)
        if axis_mode == "x":
            self.axis = np.zeros(self.n)
        else:  # sectors: x-cells where |cos| >= |sin|, y-cells elsewhere
            self.axis = np.where(np.abs(np.cos(self.psi)) >= np.abs(np.sin(self.psi)), 0.0, np.pi / 2)
        self.phis = np.radians(beams_phi_deg)
        self.K = len(self.phis)
        self.m = np.cos(self.psi - self.axis) ** 2  # modulation depth factor along the ray
        self.a = np.zeros((self.K, self.n))
        self.w = np.zeros((self.K, self.n))
        self.cosphi = np.zeros((self.K, self.n))
        for k, ph in enumerate(self.phis):
            zeta = ph - self.psi
            ab = alpha_beam(zeta) * self.m**2 / 1.0
            aa = A_ABS0 * self.m
            tm = np.cos(zeta) ** 2  # co-polar (theta_hat) share of the radiated harmonic; the rest (TE part, E_phi) is cross-polar and not counted
            self.a[k] = np.sqrt(2 * ab / self.rho) * np.exp(-(ab + aa) * (self.rho - RHO_KEEP)) * np.sqrt(tm)
            self.w[k] = 2 * ab * np.exp(-2 * (ab + aa) * (self.rho - RHO_KEEP)) / self.rho * st.PITCH**2 / (2 * np.pi) * tm

    def phase(self, k, eps):
        return BETA * self.rho - st.K0 * np.sin(THETA) * (self.x * np.cos(self.phis[k]) + self.y * np.sin(self.phis[k])) + eps

    def ideal_bits(self, eps):
        bits = np.ones((self.K, self.n))
        for k in range(self.K):
            c = np.cos(self.phase(k, eps[k]))
            bits[k] = np.where(c >= 0, 1.0, -1.0)
        return bits

    def contrib(self, eps):
        return np.array([self.a[k] * (np.pi / 2) * np.abs(np.cos(self.phase(k, eps[k]))) for k in range(self.K)])

    def currents(self, bits, eps, pols=None, ideal=False):
        outs = []
        for k in range(self.K):
            th, ph = THETA, self.phis[k]
            j = st.jones_to_current(st.JONES["LP-H"], th, ph)
            if ideal:
                cur = self.a[k] * np.exp(-1j * st.K0 * np.sin(th) * (self.x * np.cos(ph) + self.y * np.sin(ph)))
            else:
                cur = (np.pi / 2) * self.a[k] * bits[k] * np.exp(-1j * BETA * self.rho + 1j * eps[k] * 0)
            outs.append((j[0] * cur, j[1] * cur))
        return outs


def classes(words):
    """Complement classes of +-1 words (rows): canonical word has first entry +1."""
    seen = {}
    for w in map(tuple, words):
        c = w if w[0] > 0 else tuple(-v for v in w)
        seen.setdefault(c, []).append(w)
    return seen


def greedy(board, eps, budget_db, frac_static=1.0, merge=True):
    """Static relaxation with budget `frac_static * budget` then group merging with the rest. Returns bits (K, n), per-state relative amplitude loss."""
    K, n = board.K, board.n
    s0 = board.ideal_bits(eps)
    c = board.contrib(eps)  # (K, n) contribution of a cell to the peak amplitude
    a0 = c.sum(axis=1)
    rmax = 1 - 10 ** (-budget_db / 20)
    bits = s0.copy()
    r = np.zeros(K)
    # --- phase 1: make cells static (all +1 or all -1), cheapest first
    cost_p = np.array([(2 * c[k] * (s0[k] < 0)) / a0[k] for k in range(K)])  # loss if the cell is forced to +1
    cost_m = np.array([(2 * c[k] * (s0[k] > 0)) / a0[k] for k in range(K)])  # forced to -1
    use_p = cost_p.sum(axis=0) <= cost_m.sum(axis=0)
    cost = np.where(use_p[None, :], cost_p, cost_m)
    order = np.argsort(cost.sum(axis=0), kind="stable")
    cum = np.cumsum(cost[:, order], axis=1)
    ok = np.all(cum <= frac_static * rmax, axis=0)
    nst = int(np.max(np.flatnonzero(ok)) + 1) if ok.any() else 0
    sel = order[:nst]
    bits[:, sel] = np.where(use_p[sel], 1.0, -1.0)[None, :]
    r = cost[:, sel].sum(axis=1)
    # --- phase 2: merge groups (complement classes)
    if merge:
        statics = [tuple([1.0] * K), tuple([-1.0] * K)]
        while True:
            nonstatic = np.flatnonzero(~np.all(bits == bits[0], axis=0))
            if nonstatic.size == 0:
                break
            cls = classes(bits[:, nonstatic].T)
            best = None
            for ckey, ws in cls.items():
                remaining = [w for ck, ww in cls.items() if ck != ckey for w in ww] + statics
                remaining = list({tuple(w) for w in remaining} | {tuple(-v for v in w) for w in remaining if w not in statics})
                members = [i for i in nonstatic if (tuple(bits[:, i]) in ws)]
                add = np.zeros(K)
                new = {}
                for i in members:
                    # added loss of moving cell i from its current bits to word w (relative to ideal s0)
                    best_c, best_w = None, None
                    for w in remaining:
                        wv = np.array(w)
                        loss = 2 * c[:, i] * (wv != s0[:, i]) / a0
                        cur = 2 * c[:, i] * (bits[:, i] != s0[:, i]) / a0
                        d = loss - cur
                        sc = d.sum()
                        if best_c is None or sc < best_c:
                            best_c, best_w, best_d = sc, wv, d
                    add += best_d
                    new[i] = best_w
                if np.all(r + add <= rmax + 1e-12):
                    score = float(np.max(add))
                    if best is None or score < best[0]:
                        best = (score, add, new)
            if best is None:
                break
            _, add, new = best
            for i, w in new.items():
                bits[:, i] = w
            r = r + add
    return bits, s0, r, a0


def summarize(board, bits):
    static = np.all(bits == bits[0], axis=0)
    nonstatic = np.flatnonzero(~static)
    words = classes(bits[:, nonstatic].T) if nonstatic.size else {}
    n_words = len({tuple(w) for w in bits[:, nonstatic].T}) if nonstatic.size else 0
    return {"static_frac": float(static.mean()), "diodes": int(nonstatic.size), "lines_with_inverters": len(words), "lines_no_inverters": n_words}


def evaluate(board, bits, eps, ideal_cache=None):
    res = []
    cur = board.currents(bits, eps)
    idl = ideal_cache or board.currents(bits, eps, ideal=True)
    for k in range(board.K):
        et, ep = st.far_field(board.x, board.y, *cur[k])
        ei, epi = st.far_field(board.x, board.y, *idl[k])
        th, ph = THETA, board.phis[k]
        mi = st.metrics(ei, epi, th, ph, "LP-H")
        mb = st.metrics(et, ep, th, ph, "LP-H")
        i0 = mi["_i"]
        inten_i = (np.abs(ei) ** 2 + np.abs(epi) ** 2)[i0]
        inten_b = (np.abs(et) ** 2 + np.abs(ep) ** 2)[i0]
        eta_beam = LAUNCH * float(board.w[k].sum())
        g_ideal = mi["D_dbi"] + 10 * np.log10(eta_beam)
        gain = g_ideal + 10 * np.log10(inten_b / inten_i)
        res.append({"state": k + 1, "phi_deg": float(np.degrees(ph)), "peak_theta_deg": mb["peak_theta_deg"], "peak_phi_deg": mb["peak_phi_deg"], "gain_dbi": float(gain), "G_ideal_dbi": float(g_ideal),
                    "D_ideal_dbi": mi["D_dbi"], "SLL_db": mb["SLL_db"], "converted_beam_power": eta_beam})
    return res, idl
