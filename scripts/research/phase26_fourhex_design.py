"""Phase 26: designs and model predictions of the two 4-hexagon tensor-CP boards (shared 6.4 mm lattice).

    python scripts/research/phase26_fourhex_design.py <out_dir> [--quick] [--scale-s 1.0,1.7]

Board A (static): beams (theta 30; phi 45/135/225/315) on ports 1..4, each hexagon with its own phase reference.
Board B (Butler-ready): all four hexagons aimed at (30, 225), ONE global phase reference (the board centre), so that
equal-phase drive of the four ports is a single coherent aperture and the Butler phase steps steer it.

Per port the model gives: main pattern (cell currents of the design, grounded-slab Green function) and the
`pickup` pattern of the unconverted surface wave crossing the shared edges into the neighbours' cells (cell current
`A exp(-j phi_in) (cos chi, sin chi)`, chi = 2 a_cell - psi_in, alpha of the cell's own design: ASSUMPTION), with a
sensitivity scale `s` on the unconverted edge power (s = 1: the design march; s = 1.7: upper reading of the HFSS
|E| line of the single hexagon, `hfss_results/hex150_tensor_cp_single/efield_line_magE_phi225_1mm.csv`).
Gain anchoring: E_j is normalised to the fraction `LAUNCH eta_j (1 + r_j)` of the accepted power radiated by the cells,
+ the HFSS offset of the single hexagon (-0.68 dB). Labels: MODEL_ONLY, NOT_HFSS_VALIDATED, SOLVER_SUPERCELL.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import phase26_fourhex_scope as sc  # noqa: E402
import phase26_tensor_cp_hexagon as cp  # noqa: E402
from src.hologram.tensor_leaky_cp import SupercellTable, cell_design, march_rays  # noqa: E402
from src.pattern.vector_aperture_pattern import copol_component  # noqa: E402

FEEDS, PHI_A, LAUNCH, OFFSET_DB, BETA_NB = sc.FEEDS, sc.PHI, sc.LAUNCH, sc.OFFSET_DB, sc.BETA_NB
PHI_B = {1: 225.0, 2: 225.0, 3: 225.0, 4: 225.0}
KEEP_FEED_MM = 3.2


def hex_of_point(x, y):
    for j, c in FEEDS.items():
        if sc.inside(sc.hex_poly_mm(c), np.asarray(x), np.asarray(y), tol=1e-6):
            return j
    return 0


def design_board(board, table):
    """cells (global mm) + per-port march/design, `board` 'A' or 'B'."""
    phis = PHI_A if board == "A" else PHI_B
    poly0 = cp.hexagon()
    out = {}
    for j, f in FEEDS.items():
        xs, ys = sc.lattice_cells_mm(f, shared=True)
        keep = np.hypot(xs - f[0], ys - f[1]) >= KEEP_FEED_MM
        xs, ys = xs[keep], ys[keep]
        th0, ph0 = np.radians(30.0), np.radians(phis[j])
        m = march_rays(table, poly0, (0.0, 0.0), cp.K0, th0, ph0, 0.9, taper_c=0.5, n_psi=720)
        xr, yr = (xs - f[0]) * 1e-3, (ys - f[1]) * 1e-3
        d = cell_design(m, table, xr, yr, (0.0, 0.0), cp.K0, th0, ph0, hand=+1)
        if board == "B":  # global phase reference: the board centre instead of this feed
            kt = cp.K0 * np.sin(th0)
            kdotr = kt * (xs * 1e-3 * np.cos(ph0) + ys * 1e-3 * np.sin(ph0))
            chi = d["phi_ref"] - kdotr
            axis = 0.5 * (chi + d["phi"])
            jj = d["amp"] * np.exp(-1j * d["phi_ref"])
            d = dict(d)
            d["slot_angle_rad"] = np.mod(axis, np.pi)
            d["jx"], d["jy"] = jj * np.cos(chi), jj * np.sin(chi)
        out[j] = {"x_mm": xs, "y_mm": ys, "march": m, "design": d, "beam_phi": phis[j]}
    return out


def crossing_edges(board_design):
    """Unconverted power (fraction of the launched surface-wave power) leaving through each shared edge / free edge."""
    res = {}
    for j, f in FEEDS.items():
        m = board_design[j]["march"]
        psi, r_e, w_e = m["psi"], m["r_edge"], m["w_edge"]
        dpsi = 2 * np.pi / psi.size
        tgt = np.array([hex_of_point(f[0] + (r * 1e3 + 1.0) * np.cos(p), f[1] + (r * 1e3 + 1.0) * np.sin(p)) for p, r in zip(psi, r_e, strict=True)])
        row = {}
        for k in [0] + [q for q in FEEDS if q != j]:
            s = tgt == k
            row["free_edge" if k == 0 else f"to_hex{k}"] = float(np.sum(w_e[s]) * dpsi / (2 * np.pi))
        row["total_unconverted"] = float(np.sum(w_e) * dpsi / (2 * np.pi))
        res[j] = row
    return res


def far_fields(bd, scale_s, green, theta, phi_g):
    """Main + pickup far fields (global frame, positions relative to the BOARD CENTRE)."""
    ax = np.concatenate([bd[j]["x_mm"] for j in FEEDS])
    ay = np.concatenate([bd[j]["y_mm"] for j in FEEDS])
    owner = np.concatenate([np.full(bd[j]["x_mm"].size, j) for j in FEEDS])
    alpha = np.concatenate([bd[j]["design"]["alpha"] for j in FEEDS])
    ang = np.concatenate([bd[j]["design"]["slot_angle_rad"] for j in FEEDS])
    tree = cKDTree(np.c_[ax, ay])
    out = {}
    for j, f in FEEDS.items():
        d, m = bd[j]["design"], bd[j]["march"]
        xg, yg = bd[j]["x_mm"] * 1e-3, bd[j]["y_mm"] * 1e-3
        em = cp.pattern(d, xg, yg, green, theta, phi_g)
        nb = owner != j
        xn, yn = ax[nb] - f[0], ay[nb] - f[1]
        rho = np.hypot(xn, yn) * 1e-3
        psi = np.arctan2(yn, xn)
        r_edge = np.interp(psi, m["psi"], m["r_edge"], period=2 * np.pi)
        w_edge = np.interp(psi, m["psi"], m["w_edge"], period=2 * np.pi) * scale_s
        iph = np.minimum(np.round(r_edge / m["drho"]).astype(int), m["phi_ref"].shape[1] - 1)
        ip = np.argmin(np.abs(((psi[:, None] - m["psi"][None, :] + np.pi) % (2 * np.pi)) - np.pi), axis=1)
        phi_in = m["phi_ref"][ip, iph] + BETA_NB * cp.K0 * (rho - r_edge)
        w = np.empty(rho.size)
        for i in range(rho.size):
            n_s = max(int((rho[i] - r_edge[i]) / 1e-3), 1)
            s = np.linspace(r_edge[i], rho[i], n_s + 1)[1:]
            pts = np.c_[f[0] + 1e3 * s * np.cos(psi[i]), f[1] + 1e3 * s * np.sin(psi[i])]
            dist, idx = tree.query(pts)
            a = np.where(dist < 4.5, alpha[idx], 0.0)
            w[i] = w_edge[i] * np.exp(-2.0 * np.sum(a) * (s[1] - s[0] if s.size > 1 else 1e-3))
        amp = np.sqrt(np.maximum(2.0 * alpha[nb] * w / np.maximum(rho, 1e-9), 0.0))
        chi = 2.0 * ang[nb] - psi
        jj = amp * np.exp(-1j * phi_in)
        cur = {"jx": jj * np.cos(chi), "jy": jj * np.sin(chi)}
        ep = cp.pattern(cur, ax[nb] * 1e-3, ay[nb] * 1e-3, green, theta, phi_g)
        out[j] = {"main": em, "pickup": ep, "eta_conv": m["eta_conv"], "n_nb": int(nb.sum())}
    return out


def intensity(et, ep):
    return np.abs(et) ** 2 + np.abs(ep) ** 2


def prad(et, ep, theta, phi):
    tg, _ = np.meshgrid(theta, phi, indexing="ij")
    w = np.sin(np.radians(tg)) * np.radians(theta[1] - theta[0]) * np.radians(phi[1] - phi[0])
    return float(np.sum(intensity(et, ep) * w))


def normalised(ff, theta, phi):
    """E_j scaled so that its radiated power is LAUNCH eta (1 + r) of 1 W accepted (r = pickup/main power)."""
    em, pk = ff["main"], ff["pickup"]
    p_main = prad(*em, theta, phi)
    et, ep = em[0] + pk[0], em[1] + pk[1]
    r = prad(*pk, theta, phi) / p_main
    rho = LAUNCH * ff["eta_conv"] * (1.0 + r)
    sc_ = np.sqrt(rho / prad(et, ep, theta, phi))
    return et * sc_, ep * sc_, r


def pol_metrics(et, ep, theta, phi, beam):
    g = 4 * np.pi * intensity(et, ep) * 10 ** (OFFSET_DB / 10)
    rh = np.abs(copol_component(et, ep, "rhcp")) ** 2 * 4 * np.pi * 10 ** (OFFSET_DB / 10)
    lh = np.abs(copol_component(et, ep, "lhcp")) ** 2 * 4 * np.pi * 10 ** (OFFSET_DB / 10)
    tg, pg = np.meshgrid(theta, phi, indexing="ij")
    t0, p0 = np.radians(beam[0]), np.radians(beam[1])
    cosd = np.cos(np.radians(tg)) * np.cos(t0) + np.sin(np.radians(tg)) * np.sin(t0) * np.cos(np.radians(pg) - p0)
    near = np.degrees(np.arccos(np.clip(cosd, -1, 1))) <= 12.0
    i = np.unravel_index(np.argmax(np.where(near, rh, -1)), rh.shape)
    # sub-grid peak by a parabola through the three samples in theta and phi (periodic in phi)
    def para(vm, v0, vp):
        den = vm - 2 * v0 + vp
        return 0.0 if abs(den) < 1e-30 else 0.5 * (vm - vp) / den
    it, ip = i
    dt = para(rh[max(it - 1, 0), ip], rh[it, ip], rh[min(it + 1, rh.shape[0] - 1), ip])
    dp = para(rh[it, (ip - 1) % rh.shape[1]], rh[it, ip], rh[it, (ip + 1) % rh.shape[1]])
    th_pk = float(theta[it] + dt * (theta[1] - theta[0]))
    ph_pk = float(phi[ip] + dp * (phi[1] - phi[0]))
    cop, cr = np.sqrt(rh[i]), np.sqrt(lh[i])
    ar = 20 * np.log10((cop + cr) / max(abs(cop - cr), 1e-12))
    dpk = np.degrees(np.arccos(np.clip(np.cos(np.radians(tg)) * np.cos(np.radians(tg[i])) + np.sin(np.radians(tg)) * np.sin(np.radians(tg[i]))
                                       * np.cos(np.radians(pg) - np.radians(pg[i])), -1, 1)))
    # local maxima of the RHCP pattern outside a 9 deg cone about the peak (lobes), top 4
    from scipy.ndimage import maximum_filter

    mx = (rh == maximum_filter(rh, size=(5, 9), mode="wrap")) & (dpk > 9.0) & (rh > 10 ** (-30 / 10) * rh[i])
    lobes = sorted(((float(rh[k] / rh[i]), float(tg[k]), float(pg[k]), float(dpk[k])) for k in zip(*np.nonzero(mx), strict=True)), reverse=True)[:4]
    beam_pts = (rh >= 0.5 * rh[i]) & (dpk <= 10.0)
    ar_map = 20 * np.log10((np.sqrt(rh) + np.sqrt(lh)) / np.maximum(np.abs(np.sqrt(rh) - np.sqrt(lh)), 1e-12))
    # strongest LHCP lobe (off-beam) relative to the RHCP peak
    lm = (lh == maximum_filter(lh, size=(5, 9), mode="wrap")) & (dpk > 9.0)
    lk = max(zip(*np.nonzero(lm), strict=True), key=lambda k: lh[k])
    return {
        "peak_theta_deg": th_pk, "peak_phi_deg": ph_pk, "rhcp_gain_dbi": float(10 * np.log10(rh[i])), "lhcp_at_peak_dbi": float(10 * np.log10(max(lh[i], 1e-30))),
        "xpd_db": float(10 * np.log10(rh[i] / max(lh[i], 1e-30))), "ar_peak_db": float(ar), "ar_worst_in_minus3db_beam_db": float(ar_map[beam_pts].max()),
        "lobes_rel_db_theta_phi_deg_from_peak": [[float(10 * np.log10(r)), t, p, dd] for r, t, p, dd in lobes],
        "strongest_lhcp_lobe": {"rel_to_rhcp_peak_db": float(10 * np.log10(lh[lk] / rh[i])), "theta_deg": float(tg[lk]), "phi_deg": float(pg[lk]), "deg_from_peak": float(dpk[lk])},
        "total_gain_peak_dbi": float(10 * np.log10(g.max())),
    }


def main(out_dir, quick=False, scales=(1.0, 1.7)):
    t0 = time.time()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    table = SupercellTable.from_files(cp.TABLE_GLOB)
    green = cp.LayeredGreenFunction(cp.SUB, cp.F, backing="pec")
    theta = np.arange(0.0, 90.0 + 1e-9, 3.0 if quick else 1.0)
    phi = np.arange(0.0, 360.0, 6.0 if quick else 1.0)
    res = {"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "SOLVER_SUPERCELL", "GAMMA_0_COLLINEAR_ONLY"], "grid": {"theta_step": float(theta[1] - theta[0]), "phi_step": float(phi[1] - phi[0])}}
    fields = {}
    boards = {}
    for board in ("A", "B"):
        bd = design_board(board, table)
        boards[board] = bd
        n = sum(bd[j]["x_mm"].size for j in FEEDS)
        res[f"board_{board}"] = {"n_cells": int(n), "crossing_unconverted_power_fraction_of_launched_wave": {str(j): v for j, v in crossing_edges(bd).items()}}
        cells = []
        for j in FEEDS:
            d = bd[j]["design"]
            for i in range(bd[j]["x_mm"].size):
                cells.append({"x_mm": float(bd[j]["x_mm"][i]), "y_mm": float(bd[j]["y_mm"][i]), "slot_L_mm": float(d["slot_L_mm"][i]),
                              "slot_W_mm": float(d["slot_W_mm"][i]), "slot_angle_deg": float(np.degrees(d["slot_angle_rad"][i])),
                              "alpha_np_per_m": float(d["alpha"][i]), "port": int(j)})
        (out_dir / f"fourhex_{board}_cells.json").write_text(json.dumps({"board": board, "cells": cells}))
        for s in scales:
            ff = far_fields(bd, s, green, theta, phi)
            fields[(board, s)] = {j: normalised(ff[j], theta, phi) for j in FEEDS}
            fields[(board, s, "raw")] = ff
            print(f"board {board} s={s}: fields done ({time.time() - t0:.0f}s)", flush=True)
    # B vs A: the slot-angle offset per hexagon (the same lattice, same cells; the phase reference only rotates the slots)
    off = {}
    for j in FEEDS:
        da = np.degrees(boards["B"][j]["design"]["slot_angle_rad"] - boards["A"][j]["design"]["slot_angle_rad"])
        da = (da + 90.0) % 180.0 - 90.0  # slots are undirected: modulo 180 deg
        mean = float(np.degrees(0.5 * np.angle(np.mean(np.exp(2j * np.radians(da))))))
        off[str(j)] = {"mean_deg": mean, "max_abs_deviation_from_mean_deg": float(np.max(np.abs(((da - mean) + 90.0) % 180.0 - 90.0)))}
        # board A's beam phi differs from B's for ports 1, 2, 4: the offset is only a pure rotation for port 3 (phi 225 in both)
    res["board_B_minus_A_slot_angle_offset_deg"] = {"by_port": off, "note": "only port 3 (phi 225 in both boards) is a pure uniform rotation; "
                                                      "for the others the beam azimuth differs between A and B"}
    kt = cp.K0 * np.sin(np.radians(30.0))
    res["board_B_slot_rotation_vs_single_hexagon_design_deg"] = {
        str(j): float((np.degrees(-0.5 * kt * (FEEDS[j][0] * 1e-3 * np.cos(np.radians(225.0)) + FEEDS[j][1] * 1e-3 * np.sin(np.radians(225.0)))) + 90.0) % 180.0 - 90.0)
        for j in FEEDS}
    # ---- per-port-alone metrics ---------------------------------------------------------------
    for board, phis in (("A", PHI_A), ("B", PHI_B)):
        res[f"board_{board}"]["ports_alone"] = {}
        for s in scales:
            row = {}
            for j in FEEDS:
                et, ep, r = fields[(board, s)][j]
                mt = pol_metrics(et, ep, theta, phi, (30.0, phis[j]))
                mt["pickup_over_main_radiated_power"] = r
                # main-only (no pickup), same normalisation by the main pattern
                em = fields[(board, s, "raw")][j]["main"]
                sc_ = np.sqrt(LAUNCH * fields[(board, s, "raw")][j]["eta_conv"] / prad(*em, theta, phi))
                mo = pol_metrics(em[0] * sc_, em[1] * sc_, theta, phi, (30.0, phis[j]))
                mt["main_only"] = {k: mo[k] for k in ("rhcp_gain_dbi", "xpd_db", "ar_peak_db", "ar_worst_in_minus3db_beam_db")}
                pk = fields[(board, s, "raw")][j]["pickup"]
                lh_p = np.abs(copol_component(pk[0], pk[1], "lhcp")) ** 2
                rh_p = np.abs(copol_component(pk[0], pk[1], "rhcp")) ** 2
                rh_m = np.abs(copol_component(em[0], em[1], "rhcp")) ** 2
                tg_, pg_ = np.meshgrid(theta, phi, indexing="ij")
                k_l, k_r = np.unravel_index(np.argmax(lh_p), lh_p.shape), np.unravel_index(np.argmax(rh_p), rh_p.shape)
                t0_, p0_ = np.radians(30.0), np.radians(phis[j])
                dist_ = lambda k: float(np.degrees(np.arccos(np.clip(np.cos(np.radians(tg_[k])) * np.cos(t0_) + np.sin(np.radians(tg_[k])) * np.sin(t0_) * np.cos(np.radians(pg_[k]) - p0_), -1, 1))))  # noqa: E731
                mt["pickup_only"] = {
                    "lhcp_peak_rel_main_rhcp_peak_db": float(10 * np.log10(lh_p.max() / rh_m.max())), "lhcp_peak_theta_phi_deg": [float(tg_[k_l]), float(pg_[k_l])],
                    "lhcp_peak_deg_from_beam": dist_(k_l), "rhcp_peak_rel_main_rhcp_peak_db": float(10 * np.log10(rh_p.max() / rh_m.max())),
                    "rhcp_peak_theta_phi_deg": [float(tg_[k_r]), float(pg_[k_r])], "rhcp_peak_deg_from_beam": dist_(k_r),
                    "lhcp_fraction_of_pickup": float(np.sum(lh_p) / np.sum(lh_p + rh_p))}
                row[str(j)] = mt
            res[f"board_{board}"]["ports_alone"][f"s={s}"] = row
    # ---- Butler inputs for board B -----------------------------------------------------------------
    order = [2, 3, 1, 4]  # Butler outputs 1..4 = ports by increasing x
    steps = {"1L": -45.0, "2R": 135.0, "2L": -135.0, "1R": 45.0}
    res["butler"] = {"output_order_ports": order, "x_mm_of_outputs": [FEEDS[p][0] for p in order], "steps_deg": steps, "by_scale": {}}
    for s in scales:
        rows = {"equal_phase": None}
        cases = {"equal_phase": 0.0, **{k: v for k, v in steps.items()}}
        for name, step in cases.items():
            et = np.zeros_like(fields[("B", s)][1][0])
            ep = np.zeros_like(et)
            phases = {}
            for n, p in enumerate(order):
                ph = np.exp(1j * np.radians(n * step))
                phases[p] = float(((n * step + 180.0) % 360.0) - 180.0)
                a_t, a_p, _ = fields[("B", s)][p]
                et = et + 0.5 * ph * a_t
                ep = ep + 0.5 * ph * a_p
            # nominal beam from the array theory (u shift = -step/(k d)), d = 129.904 mm
            d_m = 0.129904
            du = -np.radians(step) / (cp.K0 * d_m)
            u0, v0 = np.sin(np.radians(30)) * np.cos(np.radians(225)), np.sin(np.radians(30)) * np.sin(np.radians(225))
            u1 = u0 + du
            rr = np.hypot(u1, v0)
            nominal = (float(np.degrees(np.arcsin(rr))), float(np.degrees(np.arctan2(v0, u1)) % 360.0))
            mt = pol_metrics(et, ep, theta, phi, nominal)
            mt["nominal_theta_phi_deg_array_theory"] = nominal
            mt["port_phases_deg"] = phases
            mt["accepted_power_fraction_radiated_by_cells"] = prad(et, ep, theta, phi)
            rows[name] = mt
        res["butler"]["by_scale"][f"s={s}"] = rows
    res["seconds"] = time.time() - t0
    (out_dir / "fourhex_design.json").write_text(json.dumps(res, indent=1, default=float))
    print("done", res["seconds"])
    np.savez_compressed(out_dir / "_fields_scratch.npz", theta=theta, phi=phi,
                        **{f"{b}_{s}_{j}_{c}": fields[(b, s)][j][i] for (b, s) in [(k[0], k[1]) for k in fields if len(k) == 2] for j in FEEDS for i, c in enumerate(("et", "ep"))})


if __name__ == "__main__":
    a = sys.argv[1:]
    sc_arg = tuple(float(v) for v in a[a.index("--scale-s") + 1].split(",")) if "--scale-s" in a else (1.0, 1.7)
    main(a[0], quick="--quick" in a, scales=sc_arg)
