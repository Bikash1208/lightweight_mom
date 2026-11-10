"""JOB D (model, relative trends; MODEL_ONLY / NOT_HFSS_VALIDATED): diode-count study on a CIRCULAR partitioned board (beam A static, beam B switched).
    SUB6_MESH=96 python jobs1_diode_count.py <out.json>
Engine: Phase 26 march_rays / cell_design (unmodified) on a circular aperture (R 180 mm ~ same area as the 200 mm hexagon, square lattice 11.6 mm, ring exclusion 25 mm), static NOCUT table of the board
(sub6_board_model4.tables), RHCP far field at beam A (25, 0) and B (25, 180) from the same layered-medium dyad as the board model.
Architecture (simplest, as asked): the N_sw switched cells radiate beam B in OFF and are quiet in ON (current ratio r_on); the other cells radiate beam A in both states
(designed amplitude x 1/(1 - f_p), f_p = share of the A-design power sitting in the switched cells, i.e. the static cells are strengthened to compensate; reported as 'A_boost_dB').
Switched-cell strength = alpha_sw [Np/m] (same for all switched cells): amplitude a = sqrt(2 alpha_sw W / rho) with the ray power W of the design march, attenuated by the switched cells ahead of it on the same ray
(exp(-2 alpha_sw L), L = pitch per switched cell within half a pitch of the ray).  Realised-gain ratio B/A (OFF) = |E_B|^2 / |E_A|^2 (same accepted power).
For each placement and N_sw: alpha_need = smallest alpha_sw with B >= A - 3 dB (None if not reached up to ALPHA_MAX), the B-A level at the island-cell strength 2.16 Np/m, and the ON requirement r_on <= 10^(-12/20) (power ratio alpha_on/alpha_off <= 0.063).
Placements: uniform (the board's low-discrepancy 'a'), inner (nearest the feed), outer (furthest), sector_B (wedge nearest the phi = 180 axis), sectors4 (4 wedges at 0/90/180/270), field_weighted (equal-quantile of the A-design power, golden-angle scattered).
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
from src.hologram.tensor_leaky_cp import cell_design, march_rays
from src.pattern.vector_aperture_pattern import copol_component

B = M.B
K0 = B.K0
PITCH = B.PITCH
R_BOARD, R_EXCL = 0.180, 0.025
BEAM_A, BEAM_B = B.BEAM_A, B.BEAM_B
ALPHA_MAX = 30.0
ALPHA_ISL = 2.16
NS = tuple(int(v) for v in os.environ.get("JOB_NS", "40,60,80,100,150").split(","))
PLACES = tuple(os.environ.get("JOB_PLACES", "uniform,inner,outer,sector_B,sectors4,field_weighted").split(","))
R_ON_MAX = 10 ** (-12 / 20)


def circle(r, n=360):
    a = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.stack([r * np.cos(a), r * np.sin(a)], 1)


def sites():
    n = int(R_BOARD / PITCH) + 1
    i, j = np.meshgrid(np.arange(-n, n + 1), np.arange(-n, n + 1), indexing="ij")
    x, y = i.ravel() * PITCH, j.ravel() * PITCH
    ok = (np.hypot(x, y) <= R_BOARD - PITCH / 2) & (np.hypot(x, y) >= R_EXCL)
    return x[ok], y[ok]


def beam_cell_field(green, x, y, beam, jx, jy):
    th, ph = math.radians(beam[0]), math.radians(beam[1])
    kx, ky = K0 * math.sin(th) * math.cos(ph), K0 * math.sin(th) * math.sin(ph)
    d = np.exp(1j * (kx * x + ky * y)) * PITCH**2
    g = green.spectral_dyad(np.array([kx]), np.array([ky]))[0] * math.cos(th)
    ex = g[0, 0] * d * jx + g[0, 1] * d * jy
    ey = g[1, 0] * d * jx + g[1, 1] * d * jy
    et = ex * math.cos(th) * math.cos(ph) + ey * math.cos(th) * math.sin(ph)
    ep = -ex * math.sin(ph) + ey * math.cos(ph)
    return copol_component(et, ep, "rhcp")


def select(place, n, x, y, rho, psi, powA):
    N = x.size
    if place == "uniform":
        g1, g2 = 0.7548776662466927, 0.5698402909980532
        sc = (np.round(x / PITCH) * g1 + np.round(y / PITCH) * g2) % 1.0
        idx = np.argsort(sc)[:n]
    elif place == "inner":
        idx = np.argsort(rho)[:n]
    elif place == "outer":
        idx = np.argsort(-rho)[:n]
    elif place == "sector_B":
        d = np.abs(((psi - math.pi) + math.pi) % (2 * math.pi) - math.pi)
        idx = np.argsort(d)[:n]
    elif place == "sectors4":
        d = np.min(np.abs(((psi[:, None] - np.array([0, 0.5, 1, 1.5])[None, :] * math.pi) + math.pi) % (2 * math.pi) - math.pi), axis=1)
        idx = np.argsort(d)[:n]
    elif place == "field_weighted":
        order = np.argsort(rho)
        cum = np.cumsum(powA[order]); cum /= cum[-1]
        taken, idx = set(), []
        for k in range(n):
            q = (k + 0.5) / n
            pos = int(np.searchsorted(cum, q))
            target = (k * 2.399963229728653) % (2 * math.pi)
            cand = [order[min(max(pos + o, 0), N - 1)] for o in range(-6, 7)]
            cand = [c for c in cand if c not in taken]
            if not cand:
                cand = [c for c in range(N) if c not in taken]
            best = min(cand, key=lambda c: abs(((psi[c] % (2 * math.pi)) - target + math.pi) % (2 * math.pi) - math.pi))
            taken.add(best); idx.append(best)
        idx = np.array(idx)
    else:
        raise ValueError(place)
    mask = np.zeros(N, bool)
    mask[idx] = True
    return mask


def main():
    outp = sys.argv[1]
    T = R4.tables()
    tab = M.static_table(T)
    b_fill = float(np.mean(T["st"]["b"]))
    x, y = sites()
    rho, psi = np.hypot(x, y), np.arctan2(y, x)
    poly = circle(R_BOARD)
    from src.hologram.tensor_leaky_cp import SupercellTable
    th0, ph0 = math.radians(BEAM_A[0]), math.radians(BEAM_A[1])
    green = B.get_side(0.2).green
    caps = [float(v) for v in os.environ.get("JOB_CAPS", "0.5,0.7,1.0,1.4,2.0,3.0,4.0").split(",")]
    designs = {}
    for cap in caps:
        keep = np.nonzero(T["st"]["a"] <= cap)[0]
        if keep.size < 3:
            keep = np.argsort(T["st"]["a"])[:3]
        t1 = SupercellTable(tuple(tab.slots[i] for i in keep), tab.k_grid, tab.alpha[keep], tab.beta[keep], tab.diameter_mm, tab.pitch_mm)
        m = march_rays(t1, poly := circle(R_BOARD), (0.0, 0.0), K0, th0, ph0, B.ETA_DESIGN, taper_c=0.5, rho_min_m=B.RHO_CLEAR, beta0_over_k0=b_fill)
        dA = cell_design(m, t1, x, y, (0.0, 0.0), K0, th0, ph0, hand=+1)
        dB = cell_design(m, t1, x, y, (0.0, 0.0), K0, math.radians(BEAM_B[0]), math.radians(BEAM_B[1]), hand=+1)
        live = dA["amp"] > 1e-12
        ip = np.clip(np.rint((psi - m["psi"][0]) / (m["psi"][1] - m["psi"][0])).astype(int), 0, m["psi"].size - 1)
        ir = np.clip(np.rint(rho / m["drho"]).astype(int), 0, m["rho"].size - 1)
        ub = np.where(live, 1.0 / np.maximum(dB["amp"], 1e-30), 0.0)
        designs[cap] = dict(live=live, W=m["w"][ip, ir], cellA=beam_cell_field(green, x, y, BEAM_A, dA["jx"] * live, dA["jy"] * live),
                            cellB_unit=beam_cell_field(green, x, y, BEAM_B, dB["jx"] * ub, dB["jy"] * ub), ampA=dA["amp"], eta=float(m["eta_conv"]), a_mean=float(np.mean(dA["alpha"][live])))
        d_ = designs[cap]
        d_["FA_ref"] = abs(d_["cellA"].sum())
        print("cap", cap, "slots", keep.size, "eta_conv %.3f" % m["eta_conv"], "A_ref level dB %.2f" % (20 * math.log10(d_["FA_ref"])), flush=True)
    ref_cap = max(caps, key=lambda c: designs[c]["FA_ref"])
    res = dict(label="MODEL_ONLY / NOT_HFSS_VALIDATED diode-count study", n_sites=int(x.size), R_board_mm=R_BOARD * 1e3, caps=caps, ref_all_static_best_cap=ref_cap,
               ref_all_static_A_level_dB=float(20 * math.log10(designs[ref_cap]["FA_ref"])), r_on_max_for_12dB=R_ON_MAX, power_ratio_on_off_max=R_ON_MAX**2,
               static_cap_levels={str(c): dict(A_ref_dB=float(20 * math.log10(designs[c]["FA_ref"])), eta_conv=designs[c]["eta"]) for c in caps}, rows=[])
    FAref = designs[ref_cap]["FA_ref"]
    alphas = [float(v) for v in os.environ.get("JOB_ALPHAS", "1,2.16,3,4,6,8,12,20").split(",")]
    for place in PLACES:
        for n in NS:
            mask0 = select(place, n, x, y, rho, psi, designs[ref_cap]["ampA"] ** 2)
            for alpha in alphas:
                best = None
                for cap in caps:
                    d = designs[cap]
                    mask = mask0 & d["live"]
                    sw = np.nonzero(mask)[0]
                    if sw.size == 0:
                        continue
                    # switched cells ahead of every site on its ray (perpendicular distance <= PITCH/2): depletion of W for both switched and downstream static cells
                    dxm = (np.cos(psi)[:, None] * y[sw][None, :] - np.sin(psi)[:, None] * x[sw][None, :])
                    ahead = (rho[sw][None, :] < rho[:, None] - 1e-9) & (np.abs(dxm) <= PITCH / 2)
                    dep = np.exp(-2.0 * alpha * PITCH * ahead.sum(axis=1))
                    FA = (d["cellA"] * np.where(mask, 0.0, np.sqrt(dep))).sum()
                    a = np.sqrt(2.0 * alpha * np.maximum(d["W"][sw], 0.0) * dep[sw] / np.maximum(rho[sw], 1e-9))
                    FB = abs((d["cellB_unit"][sw] * a).sum())
                    lA, lB = 20 * math.log10(abs(FA) / FAref), 20 * math.log10(max(FB, 1e-30) / FAref)
                    row = dict(cap=cap, A_dB=lA, B_dB=lB, B_minus_A=lB - lA)
                    if best is None or (row["B_minus_A"] >= -3.0 and (best["B_minus_A"] < -3.0 or row["A_dB"] > best["A_dB"])) or (best["B_minus_A"] < -3.0 and row["B_minus_A"] > best["B_minus_A"]):
                        best = row
                res["rows"].append(dict(placement=place, N_sw=int(mask0.sum()), alpha_sw=alpha, **best, ok_B_within_3dB=bool(best["B_minus_A"] >= -3.0)))
            tail = [r for r in res["rows"] if r["placement"] == place and r["N_sw"] == int(mask0.sum())]
            print(place, n, " ".join("a%.2f:%.1f/%.1f%s" % (r["alpha_sw"], r["A_dB"], r["B_minus_A"], "*" if r["ok_B_within_3dB"] else "") for r in tail), flush=True)
    Path(outp).write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
