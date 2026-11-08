"""SUB6_SPLITDISC_1BIT stage 1b: board requirement map with IDEAL cells.   MODEL_ONLY / IDEAL_CELLS / NOT_HFSS_VALIDATED.

Board-model path (no frozen code is edited, nothing under src/ is touched):
  * Phase 26 tensor engine `src/hologram/tensor_leaky_cp`: `march_rays` and `cell_design` are called UNMODIFIED with a SYNTHETIC
    `SupercellTable` (ideal cells: alpha available 0..cap Np/m, constant beta/k0, no K dependence).  The module-level 10.5 GHz
    constants of `tensor_cp_design` (F, PITCH, SUB, R_HEX, hexagon(), lattice_cells()) are NOT used; frequency, pitch, hexagon, lattice and
    substrate are passed as arguments / rebuilt here.
  * Far field: the `_vector_far_field` algebra (aperture DFT x LayeredGreenFunction dyad) written as a matrix product so that all sweep
    points share one phase matrix (checked against `_vector_far_field` by `selftest`).
  * Why not `sidequest_cp_hybrid.py` build/pattern_maps/analyse/joint_opt2: they are hard-wired to the 10.5 GHz board (`cells.csv`,
    `k2_default_lib.json`, GB library, module constants), so idealised cells would need edits of that script.  Its `placement()` idea
    (projection of each site onto the beam field, reference-phase iteration) is re-implemented here in `optimise_sites`.

Model (generalised architecture, Amendment 2026-11-07).  Every site is a tensor cell on one ray march (alpha0, W, phi_ref).  B sites (N_B,
the diodes) are switched, the other sites are static.  rho = Delta_ON/Delta_OFF of a switched cell (real nominal value, may be negative).
Currents are GROUP-NORMALISED: the static group alone, with weight 1, builds a complete beam, so each group carries 1/(1-f) resp. 1/f
times the full-aperture density, f = N_B/N:
    static site   J = (e_A + c e_B)/(1-f),         c = -rho/(1-rho)
    switched site J_OFF = e_B/((1-rho) f),         J_ON = rho_actual J_OFF
so the beam-B coherent sum is  c + 1/(1-rho) = 1 in OFF and  c + rho/(1-rho) = 0 in ON; beam A comes from the static sites only.
(e_A, e_B: unit-amplitude Phase 26 tensor-cell currents for beams A and B.)  All J are scaled by one number so that the mean radiated power
per area equals the march target (mean alpha = alpha0), amplitude ~ sqrt(alpha), and every site is clipped at its alpha <= alpha_max limit
(the march table is capped at alpha_max / (90th percentile of the site power factor)).  rho_actual = rho x scale x exp(j dphi) for the
tolerance cases.  ON: the power the switched sites stop radiating stays on the ray, W_on/W_off = exp(2 Lambda (1 - P_on/P_off)),
Lambda = int alpha0 drho, P = sum |J|^2; the mean beta of the switched sites shifts by delta_beta, so every site downstream picks up the
phase k0 delta_beta f rho (global f; local density NOT used).  eta (conversion) = march eta x clip loss, so absolute gains are model
numbers only.  Not modelled: coupling between neighbouring cells of different group, the cell's own phase change on switching beyond rho,
scattering at group boundaries.  Placement (b) is chosen with the march-0 amplitudes (cap alpha_max); the final march is then run for it.
"""
import argparse
import csv
import json
import math
import os
import platform
import subprocess
import sys
import time
from multiprocessing import Pool
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_sll_def as SLD
from src.frequency import free_space_wavenumber
from src.green_function import LayeredGreenFunction
from src.hologram.tensor_leaky_cp import SupercellTable, cell_design, march_rays
from src.materials import Substrate
from src.pattern.vector_aperture_pattern import _vector_far_field, copol_component

RES = ROOT / "results/sub6_board_req"
F = 5.8e9
K0 = free_space_wavenumber(F)
PITCH = 10e-3
RHO_CLEAR = 10e-3
H_SUB = 3.175e-3
SUB = Substrate("rt5880", 2.2, 0.0009, H_SUB)
BETA = 1.06  # mean beta/k0 of the ideal cells (BETA_NB of tensor_cp_design; not computed for this stack)
BEAM_A, BEAM_B = (25.0, 0.0), (25.0, 180.0)
ETA_DESIGN, TAPER = 0.9, 0.5
THETA = np.arange(0.0, 90.0 + 1e-9, 1.0)
PHI = np.arange(0.0, 360.0, 2.0)
OFFS = {"A_db": -1.8, "B_db": +1.5, "sll_db": +1.0}  # 10.5 GHz hybrid HFSS offsets, SIDEQUEST_CP_HYBRID.md (beam 1, beam 2, SLL)

N_B_LIST = (70, 100, 150, 200, 300, 400, 500, 600)
RHO_LIST = (0.0, 0.3, -1.0, -1.4, -2.0, -3.5)
DB_LIST = (0.0, 0.01, 0.03)
ALPHA_LIST = (1.0, 2.0, 4.0)
TOL = (("mag0.9", 0.9, 0.0), ("mag1.1", 1.1, 0.0), ("ph+5", 1.0, 5.0), ("ph+10", 1.0, 10.0))


def launch_fraction():
    """Pin launch fraction at 5.8 GHz: the Phase 26 10.5 GHz pin (top 1.219 mm above ground in h = 1.575 mm, radius 0.305 mm; 0.592)
    scaled geometrically by the h ratio, evaluated with src/dispersion/pin_launch.py."""
    from src.dispersion.pin_launch import Launcher, PinLaunch

    s = H_SUB / 1.575e-3
    pl = PinLaunch(F, 2.2, H_SUB, Launcher(1.219e-3 * s - H_SUB, 0.305e-3 * s, 0.0), tan_delta=0.0, n_seg=24)
    sp = pl.split()
    return float(sp["launch_fraction"]), {k: float(v) for k, v in sp.items()}


def hexagon(side):
    return np.array([[side * np.cos(np.radians(90 + 60 * i)), side * np.sin(np.radians(90 + 60 * i))] for i in range(6)])


def lattice(side):
    poly = hexagon(side)
    n = int(np.ceil(side / PITCH)) + 3
    k = np.arange(-n, n + 1)
    xg, yg = np.meshgrid(PITCH * k, PITCH * k)
    x, y = xg.ravel(), yg.ravel()
    ok = np.ones(x.shape, bool)
    for a, b in zip(poly, np.roll(poly, -1, axis=0), strict=True):
        ok &= (b[0] - a[0]) * (y - a[1]) - (b[1] - a[1]) * (x - a[0]) >= -1e-9
    ok &= np.hypot(x, y) >= RHO_CLEAR
    ij = np.stack([np.round(x / PITCH), np.round(y / PITCH)], 1).astype(int)
    return x[ok], y[ok], ij[ok]


def ideal_table(cap, n_slots=48):
    al = np.linspace(0.02, max(cap, 0.03), n_slots)
    kg = np.linspace(0.0, 3.0, 31)
    return SupercellTable(tuple((float(i), 1.0) for i in range(n_slots)), kg, np.repeat(al[:, None], kg.size, 1),
                          np.full((n_slots, kg.size), BETA), 6.0, PITCH * 1e3)


class Side:
    """Per-hexagon-side cache: cells, phase matrix, Green dyad."""

    def __init__(self, side):
        self.side = side
        self.x, self.y, self.ij = lattice(side)
        self.n = self.x.size
        tg, pg = np.meshgrid(THETA, PHI, indexing="ij")
        self.tg, self.pg = tg, pg
        th, ph = np.radians(tg.ravel()), np.radians(pg.ravel())
        self.th, self.ph = th, ph
        kx, ky = K0 * np.sin(th) * np.cos(ph), K0 * np.sin(th) * np.sin(ph)
        self.D = np.exp(1j * (np.outer(kx, self.x) + np.outer(ky, self.y))).astype(np.complex64) * np.float32(PITCH**2)
        self.green = LayeredGreenFunction(SUB, F, backing="pec")
        self.G = self.green.spectral_dyad(kx, ky) * np.cos(th)[:, None, None]

    def far(self, jx, jy):
        a = self.D @ jx.astype(np.complex64)
        b = self.D @ jy.astype(np.complex64)
        ex = self.G[:, 0, 0] * a + self.G[:, 0, 1] * b
        ey = self.G[:, 1, 0] * a + self.G[:, 1, 1] * b
        et = ex * np.cos(self.th) * np.cos(self.ph) + ey * np.cos(self.th) * np.sin(self.ph)
        ep = -ex * np.sin(self.ph) + ey * np.cos(self.ph)
        s = self.tg.shape
        return et.reshape(s), ep.reshape(s)

    def cell_field(self, beam, jx, jy):
        """Per-cell RHCP complex contribution at one direction."""
        th, ph = math.radians(beam[0]), math.radians(beam[1])
        kx, ky = K0 * math.sin(th) * math.cos(ph), K0 * math.sin(th) * math.sin(ph)
        d = np.exp(1j * (kx * self.x + ky * self.y)) * PITCH**2
        g = self.green.spectral_dyad(np.array([kx]), np.array([ky]))[0] * math.cos(th)
        ex = g[0, 0] * d * jx + g[0, 1] * d * jy
        ey = g[1, 0] * d * jx + g[1, 1] * d * jy
        et = ex * math.cos(th) * math.cos(ph) + ey * math.cos(th) * math.sin(ph)
        ep = -ex * math.sin(ph) + ey * math.cos(ph)
        return copol_component(et, ep, "rhcp")


_CACHE = {}


def get_side(side):
    if side not in _CACHE:
        _CACHE[side] = Side(side)
    return _CACHE[side]


def unit_fields(d):
    amp = np.maximum(d["amp"], 1e-30)
    ok = d["amp"] > 1e-12
    return np.where(ok, d["jx"] / amp, 0.0), np.where(ok, d["jy"] / amp, 0.0)


def site_units(dA, dB, mask, f, rho_d):
    """Group-normalised unit currents (jx, jy) of every site for the OFF state (before power normalisation)."""
    c = -rho_d / (1.0 - rho_d)
    ax, ay = unit_fields(dA)
    bx, by = unit_fields(dB)
    sx, sy = (ax + c * bx) / (1 - f), (ay + c * by) / (1 - f)
    wx, wy = bx / ((1 - rho_d) * f), by / ((1 - rho_d) * f)
    return np.where(mask, wx, sx), np.where(mask, wy, sy)


def march_for(side, cap):
    tab = ideal_table(cap)
    S = get_side(side)
    m = march_rays(tab, hexagon(side), (0.0, 0.0), K0, math.radians(BEAM_A[0]), math.radians(BEAM_A[1]), ETA_DESIGN, taper_c=TAPER,
                   rho_min_m=RHO_CLEAR, beta0_over_k0=BETA)
    dA = cell_design(m, tab, S.x, S.y, (0.0, 0.0), K0, math.radians(BEAM_A[0]), math.radians(BEAM_A[1]), hand=+1)
    dB = cell_design(m, tab, S.x, S.y, (0.0, 0.0), K0, math.radians(BEAM_B[0]), math.radians(BEAM_B[1]), hand=+1)
    return m, dA, dB


def cell_index(S, m):
    phi = np.arctan2(S.y, S.x)
    ip = np.argmin(np.abs(((phi[:, None] - m["psi"][None, :] + np.pi) % (2 * np.pi)) - np.pi), axis=1)
    ir = np.minimum(np.round(np.hypot(S.x, S.y) / m["drho"]).astype(int), m["rho"].size - 1)
    return ip, ir


def sites_uniform(S, n_b):
    """Placement (a): uniform interleave, low-discrepancy (R2 sequence): aperiodic, no grating lobes from a regular sublattice."""
    g1, g2 = 0.7548776662466927, 0.5698402909980532
    sc = (S.ij[:, 0] * g1 + S.ij[:, 1] * g2) % 1.0
    mask = np.zeros(S.n, bool)
    mask[np.argsort(sc)[:n_b]] = True
    return mask


def optimise_sites(S, dA, dB, f, rho_d, n_b, iters=80):
    """Placement (b): choose the N_B switched sites by projection / pair swaps (idea of sidequest_cp_hybrid.placement), maximising
    L = log|F_A| + log|F_B,OFF| - log(|R_ON|^2 + 1e-3 |F_B|^2), R_ON = beam-B residual in ON for the nominal rho (coherent RHCP sums
    of per-site fields from the march-0 amplitudes).  Monotone; number of swaps halved on failure."""
    c = -rho_d / (1.0 - rho_d)
    ax, ay = unit_fields(dA)
    bx, by = unit_fields(dB)
    amp = dA["amp"]
    aA = S.cell_field(BEAM_A, amp * ax / (1 - f), amp * ay / (1 - f))
    aBst = S.cell_field(BEAM_B, amp * c * bx / (1 - f), amp * c * by / (1 - f))
    aBsw = S.cell_field(BEAM_B, amp * bx / ((1 - rho_d) * f), amp * by / ((1 - rho_d) * f))
    mask = sites_uniform(S, n_b)

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
        dFA, dFB, dR = -aA, aBsw - aBst, rho_d * aBsw - aBst  # static -> switched
        score = (np.real(dFA * np.conj(FA)) / abs(FA) ** 2 + np.real(dFB * np.conj(FB)) / abs(FB) ** 2
                 - 2 * np.real(dR * np.conj(R)) / (abs(R) ** 2 + eps))
        st, sw = np.nonzero(~mask)[0], np.nonzero(mask)[0]
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


def build(side, alpha_max, n_b, rho_d, place):
    """March-0 (cap alpha_max) -> switched-site mask -> power factor -> final march with cap alpha_max/q90 -> clipped OFF currents."""
    S = get_side(side)
    f = n_b / S.n
    m0, dA0, dB0 = march_for(side, alpha_max)
    mask = sites_uniform(S, n_b) if place == "a" else optimise_sites(S, dA0, dB0, f, rho_d, n_b)
    ux, uy = site_units(dA0, dB0, mask, f, rho_d)
    p_unit = np.abs(ux) ** 2 + np.abs(uy) ** 2
    amp0 = dA0["amp"]
    live = amp0 > 1e-12
    kappa = float(np.sum(amp0**2 * p_unit) / max(np.sum(amp0**2), 1e-30))
    sc = 1.0 / math.sqrt(kappa)
    q90 = float(np.quantile((sc**2 * p_unit)[live], 0.9))
    cap1 = alpha_max / max(q90, 1.0)
    m, dA, dB = march_for(side, cap1)
    ux, uy = site_units(dA, dB, mask, f, rho_d)
    jx, jy = dA["amp"] * ux * sc, dA["amp"] * uy * sc
    amp_max = np.sqrt(2.0 * alpha_max * np.maximum(dA["w"], 0.0) / np.maximum(dA["rho"], 1e-9))
    jmag = np.sqrt(np.abs(jx) ** 2 + np.abs(jy) ** 2)
    clip = np.minimum(1.0, amp_max / np.maximum(jmag, 1e-30))
    r_clip = float(np.sum((jmag * clip) ** 2) / max(np.sum(jmag**2), 1e-30))
    ip, ir = cell_index(S, m)
    lam = np.cumsum(m["alpha"], axis=1) * m["drho"]
    re = np.minimum(np.floor(m["r_edge"] / m["drho"]).astype(int), m["rho"].size - 1)
    return dict(f=f, rho_d=rho_d, mask=mask, jx=jx * clip, jy=jy * clip, r_clip=r_clip, kappa=kappa, q90=q90, cap=cap1,
                lam_cell=lam[ip, ir], lam_edge=lam[np.arange(m["psi"].size), re], w_edge=m["w_edge"], eta_march=float(m["eta_conv"]),
                rho_cell=np.hypot(S.x, S.y), alpha_mean=float(np.mean(m["alpha"][m["alpha"] > 0])), n_clipped=int(np.sum(clip < 0.999)))


def state_currents(d, state_on, rho_act, dbeta):
    if not state_on:
        return d["jx"], d["jy"], d["eta_march"] * d["r_clip"]
    sw = np.where(d["mask"], rho_act, 1.0)
    jx, jy = d["jx"] * sw, d["jy"] * sw
    pr = float(np.sum(np.abs(jx) ** 2 + np.abs(jy) ** 2) / max(np.sum(np.abs(d["jx"]) ** 2 + np.abs(d["jy"]) ** 2), 1e-30))
    gain = np.exp(d["lam_cell"] * (1.0 - pr))
    phase = np.exp(-1j * K0 * dbeta * d["f"] * d["rho_cell"])
    we = d["w_edge"] * np.exp(2.0 * d["lam_edge"] * (1.0 - pr))
    eta = float(np.mean(1.0 - np.minimum(we, 1.0))) * d["r_clip"]
    return jx * gain * phase, jy * gain * phase, eta


def metrics(S, et, ep, eta, launch):
    tg, pg = S.tg, S.pg
    u = np.abs(et) ** 2 + np.abs(ep) ** 2
    rh = np.abs(copol_component(et, ep, "rhcp")) ** 2
    lh = np.abs(copol_component(et, ep, "lhcp")) ** 2
    dth, dph = np.radians(THETA[1] - THETA[0]), np.radians(PHI[1] - PHI[0])
    prad = float(np.sum(u * np.sin(np.radians(tg))) * dth * dph)
    g0 = 4 * math.pi / prad * launch * eta
    gt, gr, gl = g0 * u, g0 * rh, g0 * lh

    def peak(b):
        cd = np.cos(np.radians(tg)) * math.cos(math.radians(b[0])) + np.sin(np.radians(tg)) * math.sin(math.radians(b[0])) * np.cos(np.radians(pg - b[1]))
        near = np.degrees(np.arccos(np.clip(cd, -1, 1))) <= 6.0
        i = np.unravel_index(np.argmax(np.where(near, gr, -1)), gr.shape)
        ar = 20 * math.log10((math.sqrt(gr[i]) + math.sqrt(gl[i])) / max(abs(math.sqrt(gr[i]) - math.sqrt(gl[i])), 1e-12))
        return 10 * math.log10(gr[i] + 1e-30), ar

    gA, arA = peak(BEAM_A)
    gB, arB = peak(BEAM_B)
    ins = np.zeros(gt.shape, bool)
    u0 = 0.0
    for bm in (BEAM_A, BEAM_B):
        try:
            mk, uu = SLD.analyze(THETA, PHI, gt, peak=bm, return_mask=True)
            ins |= mk
            u0 = max(u0, uu)
        except Exception:
            pass
    outs = np.where(ins | (tg > 90.0), -np.inf, gt)
    sll = float(10 * math.log10(max(outs.max(), 1e-300) / max(u0, g0 * rh.max())))
    return dict(gA=gA, gB=gB, arA=arA, arB=arB, sll=sll, eta=eta)


def gate(bal, ar, sup, s_off, s_on):
    mg = dict(BALANCE=3.0 - bal, AR=3.0 - ar, SUPPRESSION=sup - 12.0, SLL=-13.0 - max(s_off, s_on))
    miss = [k for k, v in mg.items() if v < 0]
    return min(mg.values()), ("MEETS_ALL" if not miss else "MISSES_" + "+".join(miss))


def evaluate(side, alpha_max, n_b, place, rho_d, dbeta, launch, d):
    S = get_side(side)
    r = {}
    for st in (False, True):
        jx, jy, eta = state_currents(d, st, rho_d, dbeta)
        r[st] = metrics(S, *S.far(jx, jy), eta, launch)
    off, on = r[False], r[True]
    row = dict(side_mm=side * 1e3, alpha_max=alpha_max, N_B=n_b, N_cells=S.n, placement=place, rho_B=rho_d, abs_1_minus_rho=abs(1 - rho_d),
               c_static=-rho_d / (1 - rho_d), dbeta_over_k0=dbeta, launch=launch, alpha_mean_Np_m=d["alpha_mean"], cap_Np_m=d["cap"],
               sites_clipped=d["n_clipped"], eta_off=off["eta"], eta_on=on["eta"],
               off_GA=off["gA"], off_GB=off["gB"], off_balance=abs(off["gA"] - off["gB"]), off_AR_A=off["arA"], off_AR_B=off["arB"],
               on_GA=on["gA"], on_GB=on["gB"], on_suppression=on["gA"] - on["gB"], on_AR_A=on["arA"], sll_off=off["sll"], sll_on=on["sll"],
               beamA_change_on_minus_off=on["gA"] - off["gA"])
    ar = max(row["off_AR_A"], row["off_AR_B"], row["on_AR_A"])
    row["worst_margin_db"], row["gate"] = gate(row["off_balance"], ar, row["on_suppression"], off["sll"], on["sll"])
    # OFFSET_FROM_10GHz_HYBRID (not validated at 5.8 GHz): beam 1 -1.8 dB, beam 2 +1.5 dB, SLL +1.0 dB
    row["adj_off_balance"] = abs((off["gA"] + OFFS["A_db"]) - (off["gB"] + OFFS["B_db"]))
    row["adj_on_suppression"] = (on["gA"] + OFFS["A_db"]) - (on["gB"] + OFFS["B_db"])
    row["adj_sll_off"], row["adj_sll_on"] = off["sll"] + OFFS["sll_db"], on["sll"] + OFFS["sll_db"]
    row["adj_worst_margin_db"], row["adj_gate"] = gate(row["adj_off_balance"], ar, row["adj_on_suppression"], row["adj_sll_off"], row["adj_sll_on"])
    # tolerance of rho (ON only; OFF unchanged): |rho| x0.9 / x1.1, arg(rho) +5 / +10 deg
    for name, mag, deg in TOL:
        if rho_d == 0.0:
            row[f"sup_{name}"], row[f"gate_{name}"] = float("nan"), "NOT_APPLICABLE_RHO0"
            continue
        ra = rho_d * mag * np.exp(1j * math.radians(deg))
        jx, jy, eta = state_currents(d, True, ra, dbeta)
        mt = metrics(S, *S.far(jx, jy), eta, launch)
        sup = mt["gA"] - mt["gB"]
        _, g = gate(row["off_balance"], max(row["off_AR_A"], row["off_AR_B"], mt["arA"]), sup, off["sll"], mt["sll"])
        row[f"sup_{name}"], row[f"gate_{name}"] = sup, g
    return row


def task(a):
    side, alpha_max, n_b, place, rho_d, db_list, launch = a
    d = build(side, alpha_max, n_b, rho_d, place)
    return [evaluate(side, alpha_max, n_b, place, rho_d, db, launch, d) for db in db_list]


def provenance(extra):
    def g(*a):
        try:
            return subprocess.check_output(["git", *a], cwd=ROOT, text=True).strip()
        except Exception:
            return "unknown"
    import contextlib
    import io
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            np.show_config()
        blas = buf.getvalue()
    except Exception:
        blas = "unknown"
    return dict(machine=platform.node(), os=platform.platform(), git_commit=g("rev-parse", "HEAD"),
                git_dirty=bool(g("status", "--porcelain", "--untracked-files=no")), python=sys.version, numpy=np.__version__,
                blas_show_config=blas[-1500:], threads={k: os.environ.get(k) for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")},
                board_model_path="Phase 26 tensor engine src/hologram/tensor_leaky_cp (march_rays, cell_design, UNMODIFIED) with a synthetic ideal SupercellTable; "
                                 "far field = _vector_far_field algebra as a matrix product (scripts/research/sub6_board_req.py); no frozen code edited",
                labels="MODEL_ONLY IDEAL_CELLS NOT_HFSS_VALIDATED", **extra)


def write_csv(path, rows):
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})


def selftest():
    S = get_side(0.1)
    rng = np.random.default_rng(1)
    jx = rng.normal(size=S.n) + 1j * rng.normal(size=S.n)
    jy = rng.normal(size=S.n) + 1j * rng.normal(size=S.n)
    et, ep = S.far(jx, jy)
    worst = 0.0
    for a, b in zip([3, 25, 60, 88], [10, 90, 170, 40]):
        e1, e2 = _vector_far_field(np.array([THETA[a]]), np.array([PHI[b]]), S.x, S.y, jx, jy, F, S.green, PITCH, PITCH)
        worst = max(worst, abs(e1[0] - et[a, b]) / abs(e1[0]), abs(e2[0] - ep[a, b]) / abs(e2[0]))
    print("selftest far field vs _vector_far_field: max rel err %.2e (complex64 phase matrix)" % worst)
    cf = S.cell_field((THETA[25], PHI[90]), jx, jy).sum()
    ref = copol_component(et[25, 90], ep[25, 90], "rhcp")
    print("selftest cell_field sum vs grid: rel err %.2e" % (abs(cf - ref) / abs(ref)))
    n = S.n
    print("cells (side 100 mm):", n, " (side 200 mm):", get_side(0.2).n, " (side 150 mm):", get_side(0.15).n)
    # group-normalisation identity: B coherent sum = 1 in OFF and 0 in ON for rho = -1.4 (sampling fluctuation only)
    side = 0.2
    d = build(side, 2.0, 300, -1.4, "a")
    S2 = get_side(side)
    print("selftest build: kappa %.3f q90 %.3f cap %.3f eta_march %.3f r_clip %.3f clipped %d" % (d["kappa"], d["q90"], d["cap"], d["eta_march"], d["r_clip"], d["n_clipped"]))
    return worst


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["selftest", "coarse", "fine", "hex150"])
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args()
    RES.mkdir(parents=True, exist_ok=True)
    if a.stage == "selftest":
        selftest()
        sys.exit(0)
    t0 = time.time()
    lf, lfd = launch_fraction()
    side = 0.2
    if a.stage == "coarse":
        tasks = [(side, al, nb, "a", rho, DB_LIST, lf) for al in ALPHA_LIST for nb in N_B_LIST for rho in RHO_LIST]
    elif a.stage == "fine":
        # placement (b) only where (a) is within 3 dB of all targets (nominal worst margin >= -3 dB), for each (alpha, N_B, rho) all dbeta
        coarse = list(csv.DictReader(open(RES / "points_coarse.csv")))
        keys = sorted({(float(r["alpha_max"]), int(r["N_B"]), float(r["rho_B"])) for r in coarse if float(r["worst_margin_db"]) >= -3.0})
        tasks = [(side, al, nb, "b", rho, DB_LIST, lf) for al, nb, rho in keys]
    else:
        best = json.load(open(RES / "best_point_200mm.json"))
        nb = int(round(best["N_B"] * get_side(0.15).n / get_side(0.2).n))
        tasks = [(0.15, best["alpha_max"], nb, best["placement"], best["rho_B"], (best["dbeta_over_k0"],), lf)]
    print(f"{a.stage}: {len(tasks)} tasks, {a.workers} workers, launch fraction {lf:.4f}", flush=True)
    rows = []
    with Pool(min(a.workers, max(1, len(tasks)))) as pool:
        for rr in pool.imap_unordered(task, tasks):
            rows.extend(rr)
    dt = time.time() - t0
    rows.sort(key=lambda r: (r["side_mm"], r["alpha_max"], r["placement"], r["N_B"], r["rho_B"], r["dbeta_over_k0"]))
    write_csv(RES / f"points_{a.stage}.csv", rows)
    prov = provenance(dict(stage=a.stage, workers=a.workers, n_tasks=len(tasks), n_points=len(rows), seconds=dt, launch_fraction=lf, launch_detail=lfd,
                           start=time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(t0)), end=time.strftime("%Y-%m-%d %H:%M:%S"),
                           fixed=dict(f_hz=F, pitch_m=PITCH, substrate="RT5880 er 2.2 tan 0.0009 h 3.175 mm PEC backed", beta_over_k0=BETA,
                                      beams=[BEAM_A, BEAM_B], eta_design=ETA_DESIGN, taper=TAPER)))
    (RES / f"provenance_{a.stage}.json").write_text(json.dumps(prov, indent=1, default=str), encoding="utf-8", newline="\n")
    print(f"{a.stage}: {len(rows)} points in {dt:.1f} s", flush=True)
