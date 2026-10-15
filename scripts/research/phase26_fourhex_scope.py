"""Phase 26: scope of the 4-hexagon tensor-CP board (layout of `four_port_hex150_gap0`) -- model only, no export.

    python scripts/research/phase26_fourhex_scope.py <out.json>

(1) per-port design/pattern for beams phi = 45/135/225/315 (theta 30), same cells, pin and board; gain anchored as the
    aperture-size study (offset to the 150 mm HFSS run).
(2) cell lattices: per-port lattices (each centred on its feed, as `four_port_hex150_gap0`) vs one shared lattice:
    minimum centre distance of cells of different hexagons (D 6.0 circles overlap below 6.0 mm).
(3) unconverted-wave crossing: the surface wave of the driven port continues across the shared edge into the
    neighbours; their cells (slot angles of THEIR design) are excited by it. Current of a neighbour cell
    `A exp(-j phi_in) (cos chi, sin chi)`, chi = 2 a_cell - psi_in (the same reflection of the wave E about the slot
    axis), `A = sqrt(2 alpha W / rho)` with alpha the cell's own design alpha (ASSUMPTION: the counter/oblique wave sees
    the cell's designed alpha), W decayed along the straight ray through the neighbour. Far field of the board =
    driven hexagon + the pickup.
(4) port-to-port coupling estimate: capture width `l lambda_sw / (2 pi)` of a launcher with TM0 share `l`, power per
    radian of the driven wave decayed to the neighbour's feed; the same estimate applied to the 1-bit
    `four_port_hex150_gap0` board (HFSS `.s4p`) gives the bias carried to the tensor board.
Labels: MODEL_ONLY, NOT_HFSS_VALIDATED, SOLVER_SUPERCELL.
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import phase26_tensor_cp_hexagon as cp  # noqa: E402
from src.hologram.tensor_leaky_cp import SupercellTable, cell_design, march_rays  # noqa: E402

FEEDS = {1: (64.952, 112.5), 2: (-194.856, 112.5), 3: (-64.952, -112.5), 4: (194.856, -112.5)}
PHI = {1: 45.0, 2: 135.0, 3: 225.0, 4: 315.0}
P_MM = 6.4
LAUNCH, OFFSET_DB = 0.62, -0.6826632782724573
BETA_NB = 1.06


def hex_poly_mm(c):
    return np.array([[c[0] + 150 * np.cos(np.radians(90 + 60 * i)), c[1] + 150 * np.sin(np.radians(90 + 60 * i))] for i in range(6)])


def inside(poly, x, y, tol=1e-9):
    ok = np.ones(np.shape(x), dtype=bool)
    for a, b in zip(poly, np.roll(poly, -1, axis=0), strict=True):
        ok &= (b[0] - a[0]) * (y - a[1]) - (b[1] - a[1]) * (x - a[0]) >= -tol
    return ok


def lattice_cells_mm(centre, shared, keepout=6.5):
    n = 70  # +/- 448 mm: covers the whole 650 x 525 mm board for the shared lattice
    k = np.arange(-n, n + 1)
    if shared:
        xg, yg = np.meshgrid(P_MM * k, P_MM * k)
    else:
        xg, yg = np.meshgrid(centre[0] + P_MM * k, centre[1] + P_MM * k)
    x, y = xg.ravel(), yg.ravel()
    m = inside(hex_poly_mm(centre), x, y) & (np.hypot(x - centre[0], y - centre[1]) >= keepout)
    return x[m], y[m]


def main(out):
    table = SupercellTable.from_files(cp.TABLE_GLOB)
    green = cp.LayeredGreenFunction(cp.SUB, cp.F, backing="pec")
    theta, phi_g = np.arange(0.0, 90.0 + 1e-9, 1.0), np.arange(0.0, 360.0, 2.0)
    poly0 = cp.hexagon()
    res = {"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "SOLVER_SUPERCELL", "GAMMA_0_COLLINEAR_ONLY"], "ports": {}}

    # ---- (2) lattices -------------------------------------------------------------------------
    per = {j: lattice_cells_mm(FEEDS[j], shared=False) for j in FEEDS}
    sh = {j: lattice_cells_mm(FEEDS[j], shared=True) for j in FEEDS}
    lat = {}
    for tag, cells in (("per_port_lattices", per), ("shared_lattice", sh)):
        dmins, n_close, n_over = [], 0, 0
        for a in FEEDS:
            ta = cKDTree(np.c_[cells[a][0], cells[a][1]])
            for b in FEEDS:
                if b <= a:
                    continue
                d, _ = ta.query(np.c_[cells[b][0], cells[b][1]])
                dmins.append(float(d.min()))
                n_close += int((d < 6.4 - 1e-6).sum())
                n_over += int((d < 6.0).sum())
        lat[tag] = {"min_center_distance_between_hexagons_mm": float(min(dmins)),
                    "cells_of_one_hexagon_closer_than_6.4mm_to_another_hexagon": n_close, "cells_overlapping_D6_circles": n_over,
                    "n_cells": int(sum(len(c[0]) for c in cells.values()))}
    # shared lattice: also drop cells within 3.2 mm of any feed (pin through a cell)
    res["lattice"] = lat
    print("lattice", json.dumps(lat), flush=True)

    # ---- (1) per-port designs on the shared lattice ------------------------------------------------
    allx, ally, owner, alpha_cell, slot_ang = [], [], [], [], []
    designs = {}
    for j, f in FEEDS.items():
        xs, ys = sh[j]
        keep = np.hypot(xs - f[0], ys - f[1]) >= 3.2
        xs, ys = xs[keep], ys[keep]
        th0, ph0 = np.radians(30.0), np.radians(PHI[j])
        m = march_rays(table, poly0, (0.0, 0.0), cp.K0, th0, ph0, 0.9, taper_c=0.5, n_psi=720)
        xr, yr = (xs - f[0]) * 1e-3, (ys - f[1]) * 1e-3
        d = cell_design(m, table, xr, yr, (0.0, 0.0), cp.K0, th0, ph0, hand=+1)
        designs[j] = (m, d, xr, yr)
        et, ep = cp.pattern(d, xr, yr, green, theta, phi_g)
        mt = cp.metrics(et, ep, theta, phi_g, beam=(30.0, PHI[j]))
        g_un = mt["directivity_ap_rhcp_dbi"] + 10 * np.log10(LAUNCH * m["eta_conv"])
        res["ports"][str(j)] = {
            "beam_phi_deg": PHI[j], "n_cells": int(xs.size), "eta_conv": m["eta_conv"], "d_ap_dbi": mt["directivity_ap_rhcp_dbi"],
            "gain_anchored_dbi": float(g_un + OFFSET_DB), "xpd_db": mt["xpd_at_peak_db"], "ar_peak_db": mt["ar_at_peak_db"],
            "ar_beam_max_db": mt["ar_max_over_minus3db_beam_db"], "sll_db": mt["rhcp_sll_db"],
            "peak_theta_phi": [mt["peak_rhcp_theta_deg"], mt["peak_rhcp_phi_deg"]], "rhcp_fraction": mt["rhcp_fraction"],
        }
        print("port", j, {k: (round(v, 2) if isinstance(v, float) else v) for k, v in res["ports"][str(j)].items()}, flush=True)
        allx.append(xs); ally.append(ys); owner.append(np.full(xs.size, j))
        alpha_cell.append(d["alpha"]); slot_ang.append(d["slot_angle_rad"])
    allx, ally = np.concatenate(allx), np.concatenate(ally)
    owner, alpha_cell, slot_ang = np.concatenate(owner), np.concatenate(alpha_cell), np.concatenate(slot_ang)
    tree = cKDTree(np.c_[allx, ally])

    # ---- (3) crossing wave -----------------------------------------------------------------------
    cross = {}
    for j, f in FEEDS.items():
        m, d, xr, yr = designs[j]
        th0, ph0 = np.radians(30.0), np.radians(PHI[j])
        nb = owner != j
        xn, yn = allx[nb] - f[0], ally[nb] - f[1]
        rho = np.hypot(xn, yn) * 1e-3
        psi = np.arctan2(yn, xn)
        r_edge_ray = np.interp(psi, m["psi"], m["r_edge"], period=2 * np.pi)
        w_edge = np.interp(psi, m["psi"], m["w_edge"], period=2 * np.pi)
        iph = np.minimum(np.round(r_edge_ray / m["drho"]).astype(int), m["phi_ref"].shape[1] - 1)
        ip = np.argmin(np.abs(((psi[:, None] - m["psi"][None, :] + np.pi) % (2 * np.pi)) - np.pi), axis=1)
        phi_edge = m["phi_ref"][ip, iph]
        # decay beyond the edge through the neighbours (alpha of the nearest cell, 0 over bare strips)
        w = np.empty(rho.size)
        for i in range(rho.size):
            n_s = max(int((rho[i] - r_edge_ray[i]) / 1e-3), 1)
            s = np.linspace(r_edge_ray[i], rho[i], n_s + 1)[1:]
            pts = np.c_[f[0] + 1e3 * s * np.cos(psi[i]), f[1] + 1e3 * s * np.sin(psi[i])]
            dist, idx = tree.query(pts)
            a = np.where(dist < 4.5, alpha_cell[idx], 0.0)
            w[i] = w_edge[i] * np.exp(-2.0 * np.sum(a) * (s[1] - s[0] if s.size > 1 else 1e-3))
        phi_in = phi_edge + BETA_NB * cp.K0 * (rho - r_edge_ray)
        a_nb = alpha_cell[nb]
        amp = np.sqrt(np.maximum(2.0 * a_nb * w / np.maximum(rho, 1e-9), 0.0))
        chi = 2.0 * slot_ang[nb] - psi
        jj = amp * np.exp(-1j * phi_in)
        cur = {"jx": jj * np.cos(chi), "jy": jj * np.sin(chi)}
        # patterns in the global frame (positions in m)
        xg, yg = allx[nb] * 1e-3, ally[nb] * 1e-3
        et_n, ep_n = cp.pattern(cur, xg, yg, green, theta, phi_g)
        # the driven hexagon in the same frame (its cells + design currents, phase origin = its feed: shift)
        xm, ym = designs[j][2] + f[0] * 1e-3, designs[j][3] + f[1] * 1e-3
        # phase reference: the neighbour currents carry the propagation phase from feed j, so use positions relative to feed j for both
        et_m, ep_m = cp.pattern(d, xm - f[0] * 1e-3, ym - f[1] * 1e-3, green, theta, phi_g)
        et_n, ep_n = cp.pattern(cur, xg - f[0] * 1e-3, yg - f[1] * 1e-3, green, theta, phi_g)
        w_t = np.sin(np.radians(theta))[:, None] * np.radians(1.0) * np.radians(2.0)
        pw = lambda a, b: float(np.sum((np.abs(a) ** 2 + np.abs(b) ** 2) * w_t))  # noqa: E731
        p_main, p_nb = pw(et_m, ep_m), pw(et_n, ep_n)
        from src.pattern.vector_aperture_pattern import copol_component

        rh_n = np.abs(copol_component(et_n, ep_n, "rhcp")) ** 2
        lh_n = np.abs(copol_component(et_n, ep_n, "lhcp")) ** 2
        rh_m = np.abs(copol_component(et_m, ep_m, "rhcp")) ** 2
        i_n = np.unravel_index(np.argmax(lh_n), lh_n.shape)
        i_r = np.unravel_index(np.argmax(rh_n), rh_n.shape)
        tot_t, tot_p = et_m + et_n, ep_m + ep_n
        mt_all = cp.metrics(tot_t, tot_p, theta, phi_g, beam=(30.0, PHI[j]))
        mt_one = cp.metrics(et_m, ep_m, theta, phi_g, beam=(30.0, PHI[j]))
        eta = res["ports"][str(j)]["eta_conv"]
        cross[str(j)] = {
            "radiated_power_ratio_pickup_over_main": p_nb / p_main,
            "pickup_share_of_accepted": float(p_nb / p_main * LAUNCH * eta),
            "pickup_peak_lhcp_theta_phi": [float(theta[i_n[0]]), float(phi_g[i_n[1]])],
            "pickup_peak_lhcp_rel_main_peak_db": float(10 * np.log10(lh_n.max() / rh_m.max())),
            "pickup_peak_rhcp_theta_phi": [float(theta[i_r[0]]), float(phi_g[i_r[1]])],
            "pickup_peak_rhcp_rel_main_peak_db": float(10 * np.log10(rh_n.max() / rh_m.max())),
            "pickup_lhcp_fraction": float(np.sum(lh_n * w_t) / max(np.sum((rh_n + lh_n) * w_t), 1e-30)),
            "with_pickup": {k: mt_all[k] for k in ("directivity_ap_rhcp_dbi", "xpd_at_peak_db", "ar_at_peak_db", "ar_max_over_minus3db_beam_db", "rhcp_sll_db", "within_8deg_fraction")},
            "main_only": {k: mt_one[k] for k in ("directivity_ap_rhcp_dbi", "xpd_at_peak_db", "ar_at_peak_db", "ar_max_over_minus3db_beam_db", "rhcp_sll_db", "within_8deg_fraction")},
            "neighbour_cells_excited": int(nb.sum()),
        }
        print("cross", j, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in cross[str(j)].items() if not isinstance(v, dict)}, flush=True)
    res["crossing"] = cross

    # ---- (4) coupling estimate ----------------------------------------------------------------------
    lam_sw = 299792458.0 / cp.F / BETA_NB
    w_eff = LAUNCH * lam_sw / (2 * np.pi)
    est = {}
    for a_ in FEEDS:
        for b_ in FEEDS:
            if b_ <= a_:
                continue
            f, g = FEEDS[a_], FEEDS[b_]
            dist = float(np.hypot(g[0] - f[0], g[1] - f[1])) * 1e-3
            if dist > 0.3:
                est[f"{a_}-{b_}"] = {"distance_mm": dist * 1e3, "adjacent": False}
                continue
            m, d, xr, yr = designs[a_]
            psi = np.arctan2(g[1] - f[1], g[0] - f[0])
            r_e = float(np.interp(psi, m["psi"], m["r_edge"], period=2 * np.pi))
            w_e = float(np.interp(psi, m["psi"], m["w_edge"], period=2 * np.pi))
            # through the neighbour's cells to its feed
            s = np.linspace(r_e, dist, 200)[1:]
            pts = np.c_[f[0] + 1e3 * s * np.cos(psi), f[1] + 1e3 * s * np.sin(psi)]
            dd, idx = tree.query(pts)
            al = np.where(dd < 4.5, alpha_cell[idx], 0.0)
            w_f = w_e * np.exp(-2.0 * np.sum(al) * (s[1] - s[0]))
            flux = LAUNCH / (2 * np.pi) * w_f / dist  # accepted-power-normalised W/m
            raw = flux * w_eff
            est[f"{a_}-{b_}"] = {"distance_mm": dist * 1e3, "adjacent": True, "w_edge_unconverted": w_e,
                                 "w_at_neighbour_feed": float(w_f), "raw_estimate_db": float(10 * np.log10(raw))}
    # the same estimate on the 1-bit board (alpha 1.03 Np/m everywhere, launch 0.6), vs HFSS .s4p at 10.54 GHz
    from src.phase9.hfss_reference import read_touchstone

    net = read_touchstone(ROOT / "hfss_results/hex150_gap0/HEX_HFSSDesign1.s4p")
    fgh = net.frequency_hz / 1e9
    i0 = int(np.argmin(np.abs(fgh - 10.54)))
    hf = {f"{a_}-{b_}": float(20 * np.log10(abs(net.component(a_, b_)[i0]))) for a_ in range(1, 5) for b_ in range(a_ + 1, 5)}
    a_1b, l_1b = 1.03, 0.6
    dist = 0.2598
    w_f1 = np.exp(-2 * a_1b * (0.1299 + 0.1299))
    raw_1b = l_1b / (2 * np.pi) * w_f1 / dist * (l_1b * lam_sw / (2 * np.pi))
    bias = float(10 * np.log10(raw_1b)) - float(np.mean([v for v in hf.values() if np.isfinite(v)]))
    res["coupling"] = {"estimates": est, "capture_width_mm": float(w_eff * 1e3),
                       "one_bit_board": {"hfss_s_ij_db_at_10p54": hf, "raw_estimate_db_same_method": float(10 * np.log10(raw_1b)),
                                         "bias_db_raw_minus_mean_hfss": bias},
                       "note": "calibrated = raw - bias (the 1-bit board's HFSS coupling is far below the raw estimate)"}
    for k, v in est.items():
        if v.get("adjacent"):
            v["calibrated_estimate_db"] = v["raw_estimate_db"] - bias
    print("coupling", json.dumps(res["coupling"], indent=0, default=float)[:1500], flush=True)
    Path(out).write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main(sys.argv[1])
