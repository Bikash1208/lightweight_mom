"""Tensor-CP (rotating-slot) monopole-fed hologram design: the Phase 26 flow as an importable module.

Single 150 mm hexagon (`design_single_hexagon`), and the 4-hexagon boards on one shared 6.4 mm lattice (`design_fourhex`:
board A = static beams phi 45/135/225/315, board B = combined beam with a global phase reference). Moved from
`scripts/research/phase26_*.py` without a change of numerics (the scripts are thin wrappers). Labels of every result:
MODEL_ONLY, NOT_HFSS_VALIDATED (except the single hexagon's HFSS check, see `docs/validation/README.md` section 7),
SOLVER_SUPERCELL. Gain is anchored to the single-hexagon HFSS run (`OFFSET_DB`).
"""

import functools
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from ..frequency import free_space_wavenumber
from ..green_function import LayeredGreenFunction
from ..materials import Substrate
from ..pattern.vector_aperture_pattern import _vector_far_field, circular_components, copol_component
from .leaky_forward import pattern_budget
from .tensor_leaky_cp import SupercellTable, cell_design, march_rays

F = 10.54e9
K0 = free_space_wavenumber(F)
PITCH = 6.4e-3
SUB = Substrate("rt5880", 2.2, 0.0009, 1.575e-3)
R_HEX = 0.150
THETA0, PHI0 = np.radians(30.0), np.radians(225.0)
LAUNCH_BARE, LAUNCH_LOADED = 0.592, (0.66, 0.88)
RHO_CLEAR = 6.5e-3
TABLE_GLOB = str(Path(__file__).resolve().parents[2] / "results/phase26/supercell_table_*.json")


def hexagon():
    return np.array(
        [[R_HEX * np.cos(np.radians(90 + 60 * i)), R_HEX * np.sin(np.radians(90 + 60 * i))] for i in range(6)]
    )


def lattice_cells(poly):
    n = int(np.ceil(R_HEX / PITCH)) + 3  # lattice wide enough for any hexagon size
    k = np.arange(-n, n + 1)
    xg, yg = np.meshgrid(PITCH * k, PITCH * k)
    x, y = xg.ravel(), yg.ravel()
    ok = np.ones(x.shape, dtype=bool)
    for a, b in zip(poly, np.roll(poly, -1, axis=0), strict=True):
        cross = (b[0] - a[0]) * (y - a[1]) - (b[1] - a[1]) * (x - a[0])
        ok &= cross >= -1e-9
    ok &= np.hypot(x, y) >= RHO_CLEAR
    return x[ok], y[ok]


def pattern(design, x, y, green, theta_deg, phi_deg):
    tg, pg = np.meshgrid(theta_deg, phi_deg, indexing="ij")
    et = np.zeros(tg.size, dtype=complex)
    ep = np.zeros(tg.size, dtype=complex)
    th, ph = tg.ravel(), pg.ravel()
    for lo in range(0, th.size, 4000):
        sl = slice(lo, lo + 4000)
        et[sl], ep[sl] = _vector_far_field(th[sl], ph[sl], x, y, design["jx"], design["jy"], F, green, PITCH, PITCH)
    return et.reshape(tg.shape), ep.reshape(tg.shape)


def metrics(et, ep, theta, phi, label_launch=LAUNCH_BARE, eta_conv=1.0, beam=(30.0, 225.0)):
    b = pattern_budget(theta, phi, et, ep, beam[0], beam[1], (8.0, 20.0))
    er, el = circular_components(et, ep)
    rh = np.abs(copol_component(et, ep, "rhcp")) ** 2
    lh = np.abs(copol_component(et, ep, "lhcp")) ** 2
    tg, pg = np.meshgrid(theta, phi, indexing="ij")
    t0, p0 = np.radians(beam[0]), np.radians(beam[1])
    cosd = np.cos(np.radians(tg)) * np.cos(t0) + np.sin(np.radians(tg)) * np.sin(t0) * np.cos(np.radians(pg) - p0)
    dist = np.degrees(np.arccos(np.clip(cosd, -1, 1)))
    near = dist <= 6.0
    i = np.unravel_index(np.argmax(np.where(near, rh, -1)), rh.shape)
    peak = rh[i]
    # -3 dB beam region around the peak (within 10 deg of the beam direction)
    pk_dir = (np.radians(tg[i]), np.radians(pg[i]))
    cd = np.cos(np.radians(tg)) * np.cos(pk_dir[0]) + np.sin(np.radians(tg)) * np.sin(pk_dir[0]) * np.cos(
        np.radians(pg) - pk_dir[1]
    )
    dpk = np.degrees(np.arccos(np.clip(cd, -1, 1)))
    beam = (rh >= 0.5 * peak) & (dpk <= 10.0)
    cop, cross = np.sqrt(rh), np.sqrt(lh)
    ar_lin = (cop + cross) / np.maximum(np.abs(cop - cross), 1e-12)
    ar_db = 20 * np.log10(ar_lin)
    out = dict(b)
    out["ar_at_peak_db"] = float(ar_db[i])
    out["ar_max_over_minus3db_beam_db"] = float(ar_db[beam].max())
    out["minus3db_beam_points"] = int(beam.sum())
    # SLL: RHCP outside a 15 deg cone about the peak
    side = rh[dpk > 15.0]
    out["rhcp_sll_db"] = float(10 * np.log10(side.max() / peak))
    out["peak_rhcp_theta_deg"], out["peak_rhcp_phi_deg"] = float(tg[i]), float(pg[i])
    return out


FEEDS = {1: (64.952, 112.5), 2: (-194.856, 112.5), 3: (-64.952, -112.5), 4: (194.856, -112.5)}
PHI = {1: 45.0, 2: 135.0, 3: 225.0, 4: 315.0}
P_MM = 6.4
LAUNCH, OFFSET_DB = 0.62, -0.6826632782724573
BETA_NB = 1.06


def hex_poly_mm(c):
    return np.array(
        [[c[0] + 150 * np.cos(np.radians(90 + 60 * i)), c[1] + 150 * np.sin(np.radians(90 + 60 * i))] for i in range(6)]
    )


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


PHI_A = PHI
PHI_B = {1: 225.0, 2: 225.0, 3: 225.0, 4: 225.0}
KEEP_FEED_MM = 3.2


def hex_of_point(x, y):
    for j, c in FEEDS.items():
        if inside(hex_poly_mm(c), np.asarray(x), np.asarray(y), tol=1e-6):
            return j
    return 0


def design_board(board, table, taper_c=0.5, eta_design=0.9, density_factory=None):
    """cells (global mm) + per-port march/design, `board` 'A' or 'B'."""
    phis = PHI_A if board == "A" else PHI_B
    poly0 = hexagon()
    out = {}
    for j, f in FEEDS.items():
        xs, ys = lattice_cells_mm(f, shared=True)
        keep = np.hypot(xs - f[0], ys - f[1]) >= KEEP_FEED_MM
        xs, ys = xs[keep], ys[keep]
        th0, ph0 = np.radians(30.0), np.radians(phis[j])
        m = march_rays(
            table,
            poly0,
            (0.0, 0.0),
            K0,
            th0,
            ph0,
            eta_design,
            taper_c=taper_c,
            n_psi=720,
            density_fn=None if density_factory is None else density_factory(j),
        )
        xr, yr = (xs - f[0]) * 1e-3, (ys - f[1]) * 1e-3
        d = cell_design(m, table, xr, yr, (0.0, 0.0), K0, th0, ph0, hand=+1)
        if board == "B":  # global phase reference: the board centre instead of this feed
            kt = K0 * np.sin(th0)
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
        tgt = np.array(
            [
                hex_of_point(f[0] + (r * 1e3 + 1.0) * np.cos(p), f[1] + (r * 1e3 + 1.0) * np.sin(p))
                for p, r in zip(psi, r_e, strict=True)
            ]
        )
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
        em = pattern(d, xg, yg, green, theta, phi_g)
        nb = owner != j
        xn, yn = ax[nb] - f[0], ay[nb] - f[1]
        rho = np.hypot(xn, yn) * 1e-3
        psi = np.arctan2(yn, xn)
        r_edge = np.interp(psi, m["psi"], m["r_edge"], period=2 * np.pi)
        w_edge = np.interp(psi, m["psi"], m["w_edge"], period=2 * np.pi) * scale_s
        iph = np.minimum(np.round(r_edge / m["drho"]).astype(int), m["phi_ref"].shape[1] - 1)
        ip = np.argmin(np.abs(((psi[:, None] - m["psi"][None, :] + np.pi) % (2 * np.pi)) - np.pi), axis=1)
        phi_in = m["phi_ref"][ip, iph] + BETA_NB * K0 * (rho - r_edge)
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
        ep = pattern(cur, ax[nb] * 1e-3, ay[nb] * 1e-3, green, theta, phi_g)
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
    dpk = np.degrees(
        np.arccos(
            np.clip(
                np.cos(np.radians(tg)) * np.cos(np.radians(tg[i]))
                + np.sin(np.radians(tg)) * np.sin(np.radians(tg[i])) * np.cos(np.radians(pg) - np.radians(pg[i])),
                -1,
                1,
            )
        )
    )
    # local maxima of the RHCP pattern outside a 9 deg cone about the peak (lobes), top 4
    from scipy.ndimage import maximum_filter

    mx = (rh == maximum_filter(rh, size=(5, 9), mode="wrap")) & (dpk > 9.0) & (rh > 10 ** (-30 / 10) * rh[i])
    lobes = sorted(
        ((float(rh[k] / rh[i]), float(tg[k]), float(pg[k]), float(dpk[k])) for k in zip(*np.nonzero(mx), strict=True)),
        reverse=True,
    )[:4]
    beam_pts = (rh >= 0.5 * rh[i]) & (dpk <= 10.0)
    ar_map = 20 * np.log10((np.sqrt(rh) + np.sqrt(lh)) / np.maximum(np.abs(np.sqrt(rh) - np.sqrt(lh)), 1e-12))
    # strongest LHCP lobe (off-beam) relative to the RHCP peak
    lm = (lh == maximum_filter(lh, size=(5, 9), mode="wrap")) & (dpk > 9.0)
    lk = max(zip(*np.nonzero(lm), strict=True), key=lambda k: lh[k])
    return {
        "peak_theta_deg": th_pk,
        "peak_phi_deg": ph_pk,
        "rhcp_gain_dbi": float(10 * np.log10(rh[i])),
        "lhcp_at_peak_dbi": float(10 * np.log10(max(lh[i], 1e-30))),
        "xpd_db": float(10 * np.log10(rh[i] / max(lh[i], 1e-30))),
        "ar_peak_db": float(ar),
        "ar_worst_in_minus3db_beam_db": float(ar_map[beam_pts].max()),
        "lobes_rel_db_theta_phi_deg_from_peak": [[float(10 * np.log10(r)), t, p, dd] for r, t, p, dd in lobes],
        "strongest_lhcp_lobe": {
            "rel_to_rhcp_peak_db": float(10 * np.log10(lh[lk] / rh[i])),
            "theta_deg": float(tg[lk]),
            "phi_deg": float(pg[lk]),
            "deg_from_peak": float(dpk[lk]),
        },
        "total_gain_peak_dbi": float(10 * np.log10(g.max())),
    }


def design_single_hexagon(quick=False, table_glob=None, log=print):
    """The Phase 26 C single-hexagon design scan. Returns `(res, cells)`; `res` is the content of `tensor_cp_hexagon.json`
    (without `seconds`-independent file I/O), `cells` that of `tensor_cp_hexagon_cells.json`."""
    t0 = time.time()
    table = SupercellTable.from_files(table_glob or TABLE_GLOB)
    poly = hexagon()
    x, y = lattice_cells(poly)
    green = LayeredGreenFunction(SUB, F, backing="pec")
    theta = np.arange(0.0, 90.0 + 1e-9, 2.0 if quick else 1.0)
    phi = np.arange(0.0, 360.0, 3.0 if quick else 2.0)
    res = {
        "labels": ["SOLVER_SUPERCELL", "MODEL_ONLY", "NOT_HFSS_VALIDATED", "GAMMA_0_COLLINEAR_ONLY"],
        "n_cells": int(x.size),
        "hexagon_circumradius_mm": R_HEX * 1e3,
        "pitch_mm": PITCH * 1e3,
        "slots_mm": table.slots,
        "k_grid_over_k0": table.k_grid.tolist(),
        "launch_fraction_bare": LAUNCH_BARE,
        "launch_fraction_loaded_range": LAUNCH_LOADED,
        "scan": [],
    }
    best = None
    for hand in (+1,):
        for eta_design in (0.5, 0.7, 0.9):
            for taper in (0.0, 0.5, 0.8):
                m = march_rays(
                    table, poly, (0.0, 0.0), K0, THETA0, PHI0, eta_design, taper_c=taper, n_psi=360 if quick else 720
                )
                d = cell_design(m, table, x, y, (0.0, 0.0), K0, THETA0, PHI0, hand=hand)
                et, ep = pattern(d, x, y, green, theta, phi)
                mt = metrics(et, ep, theta, phi)
                eta = m["eta_conv"]
                rec = {
                    "hand": hand,
                    "eta_design": eta_design,
                    "taper_c": taper,
                    "eta_conv": eta,
                    "eta_conv_by_ray_min_max": [float(m["eta_conv_by_ray"].min()), float(m["eta_conv_by_ray"].max())],
                    "metrics_ap": mt,
                    "gain_rhcp_dbi_launch_bare": mt["directivity_ap_rhcp_dbi"] + 10 * np.log10(LAUNCH_BARE * eta),
                    "gain_rhcp_dbi_launch_loaded": [
                        mt["directivity_ap_rhcp_dbi"] + 10 * np.log10(l * eta) for l in LAUNCH_LOADED
                    ],
                }
                res["scan"].append(rec)
                log(
                    f"hand {hand:+d} eta_d {eta_design} taper {taper}: eta_conv {eta:.3f}  D_ap {mt['directivity_ap_rhcp_dbi']:.2f}  "
                    f"G {rec['gain_rhcp_dbi_launch_bare']:.2f}  XPD {mt['xpd_at_peak_db']:.1f}  RHCPfrac {mt['rhcp_fraction']:.3f}  "
                    f"SLL {mt['rhcp_sll_db']:.1f}  peak ({mt['peak_rhcp_theta_deg']:.0f},{mt['peak_rhcp_phi_deg']:.0f})"
                )
                key = (mt["rhcp_fraction"] > 0.5, mt["rhcp_sll_db"] <= -13.0, rec["gain_rhcp_dbi_launch_bare"])
                if best is None or key > best[0]:
                    best = (key, rec, d, m)
    _, rec, d, m = best
    res["selected"] = {k: v for k, v in rec.items()}
    # ---- budget of the ACCEPTED power (assumptions stated, nothing invented) ---------------------
    mt, eta = rec["metrics_ap"], rec["eta_conv"]
    bud = {}
    for tag, lf in (
        ("launch_bare_0.592", LAUNCH_BARE),
        ("launch_loaded_0.66", LAUNCH_LOADED[0]),
        ("launch_loaded_0.88", LAUNCH_LOADED[1]),
    ):
        bud[tag] = {
            "pin_direct_space_wave_assumed_outside_beam": 1 - lf,
            "unconverted_reaches_edge": lf * (1 - eta),
            "aperture_radiated": lf * eta,
            "main_beam_8deg": lf * eta * mt["within_8deg_fraction"],
            "main_beam_8deg_rhcp": lf * eta * mt["within_8deg_rhcp_fraction"],
            "lhcp_total_in_aperture_pattern": lf * eta * mt["lhcp_fraction"],
            "grazing_plus_back_upper_bound_unmodelled": (1 - lf) + lf * (1 - eta),
        }
    res["accepted_power_budget"] = bud
    # ---- element-factor check: ideal polarization pre-compensation (NOT realised by the cell family) ----
    kx0 = K0 * np.sin(THETA0) * np.cos(PHI0)
    ky0 = K0 * np.sin(THETA0) * np.sin(PHI0)
    g = np.asarray(green.spectral_dyad(np.array(kx0), np.array(ky0)))
    mmat = np.linalg.inv(g)
    jc = np.stack([d["jx"], d["jy"]], axis=-1) @ mmat.T
    dc = dict(d)
    dc["jx"], dc["jy"] = jc[:, 0], jc[:, 1]
    etc, epc = pattern(dc, x, y, green, theta, phi)
    mc = metrics(etc, epc, theta, phi)
    res["ideal_polarization_compensation"] = {
        "label": "UPPER_BOUND_NOT_REALISED_BY_THE_CELL_FAMILY",
        "metrics_ap": mc,
        "note": "J' = g(k_t,beam)^-1 J per cell: removes the grounded slab's TM/TE element-factor ellipticity at the beam "
        "direction; needs an independently controlled second modulation (not designed here).",
    }
    log(
        f"uncompensated XPD {mt['xpd_at_peak_db']:.1f} AR {mt['ar_at_peak_db']:.2f} | ideal-compensated XPD {mc['xpd_at_peak_db']:.1f} "
        f"AR {mc['ar_at_peak_db']:.2f} dB RHCPfrac {mc['rhcp_fraction']:.3f}"
    )
    sel = {
        "x_mm": (x * 1e3).tolist(),
        "y_mm": (y * 1e3).tolist(),
        "slot_L_mm": d["slot_L_mm"].tolist(),
        "slot_W_mm": d["slot_W_mm"].tolist(),
        "slot_angle_deg": np.degrees(d["slot_angle_rad"]).tolist(),
        "alpha_np_per_m": d["alpha"].tolist(),
        "beta_over_k0": d["beta_over_k0"].tolist(),
        "amp": d["amp"].tolist(),
    }
    res["seconds"] = time.time() - t0
    log("selected:", {k: rec[k] for k in ("hand", "eta_design", "taper_c", "eta_conv", "gain_rhcp_dbi_launch_bare")})
    return res, sel


def _deg_from_beam(tg, pg, t0, p0, k):
    c = np.cos(np.radians(tg[k])) * np.cos(t0) + np.sin(np.radians(tg[k])) * np.sin(t0) * np.cos(np.radians(pg[k]) - p0)
    return float(np.degrees(np.arccos(np.clip(c, -1, 1))))


def design_fourhex(quick=False, scales=(1.0, 1.7), table_glob=None, log=print):
    """Boards A and B on the shared lattice. Returns `(res, cells, fields)`: `fourhex_design.json` content, `{board: cell list}`, per-port fields."""
    t0 = time.time()
    table = SupercellTable.from_files(table_glob or TABLE_GLOB)
    green = LayeredGreenFunction(SUB, F, backing="pec")
    theta = np.arange(0.0, 90.0 + 1e-9, 3.0 if quick else 1.0)
    phi = np.arange(0.0, 360.0, 6.0 if quick else 1.0)
    res = {
        "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "SOLVER_SUPERCELL", "GAMMA_0_COLLINEAR_ONLY"],
        "grid": {"theta_step": float(theta[1] - theta[0]), "phi_step": float(phi[1] - phi[0])},
    }
    fields = {}
    boards = {}
    cells_by_board = {}
    for board in ("A", "B"):
        bd = design_board(board, table)
        boards[board] = bd
        n = sum(bd[j]["x_mm"].size for j in FEEDS)
        res[f"board_{board}"] = {
            "n_cells": int(n),
            "crossing_unconverted_power_fraction_of_launched_wave": {str(j): v for j, v in crossing_edges(bd).items()},
        }
        cells = []
        for j in FEEDS:
            d = bd[j]["design"]
            for i in range(bd[j]["x_mm"].size):
                cells.append(
                    {
                        "x_mm": float(bd[j]["x_mm"][i]),
                        "y_mm": float(bd[j]["y_mm"][i]),
                        "slot_L_mm": float(d["slot_L_mm"][i]),
                        "slot_W_mm": float(d["slot_W_mm"][i]),
                        "slot_angle_deg": float(np.degrees(d["slot_angle_rad"][i])),
                        "alpha_np_per_m": float(d["alpha"][i]),
                        "port": int(j),
                    }
                )
        cells_by_board[board] = cells
        for s in scales:
            ff = far_fields(bd, s, green, theta, phi)
            fields[(board, s)] = {j: normalised(ff[j], theta, phi) for j in FEEDS}
            fields[(board, s, "raw")] = ff
            log(f"board {board} s={s}: fields done ({time.time() - t0:.0f}s)")
    # B vs A: the slot-angle offset per hexagon (the same lattice, same cells; the phase reference only rotates the slots)
    off = {}
    for j in FEEDS:
        da = np.degrees(boards["B"][j]["design"]["slot_angle_rad"] - boards["A"][j]["design"]["slot_angle_rad"])
        da = (da + 90.0) % 180.0 - 90.0  # slots are undirected: modulo 180 deg
        mean = float(np.degrees(0.5 * np.angle(np.mean(np.exp(2j * np.radians(da))))))
        off[str(j)] = {
            "mean_deg": mean,
            "max_abs_deviation_from_mean_deg": float(np.max(np.abs(((da - mean) + 90.0) % 180.0 - 90.0))),
        }
        # board A's beam phi differs from B's for ports 1, 2, 4: the offset is only a pure rotation for port 3 (phi 225 in both)
    res["board_B_minus_A_slot_angle_offset_deg"] = {
        "by_port": off,
        "note": "only port 3 (phi 225 in both boards) is a pure uniform rotation; "
        "for the others the beam azimuth differs between A and B",
    }
    kt = K0 * np.sin(np.radians(30.0))
    res["board_B_slot_rotation_vs_single_hexagon_design_deg"] = {
        str(j): float(
            (
                np.degrees(
                    -0.5
                    * kt
                    * (FEEDS[j][0] * 1e-3 * np.cos(np.radians(225.0)) + FEEDS[j][1] * 1e-3 * np.sin(np.radians(225.0)))
                )
                + 90.0
            )
            % 180.0
            - 90.0
        )
        for j in FEEDS
    }
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
                mt["main_only"] = {
                    k: mo[k] for k in ("rhcp_gain_dbi", "xpd_db", "ar_peak_db", "ar_worst_in_minus3db_beam_db")
                }
                pk = fields[(board, s, "raw")][j]["pickup"]
                lh_p = np.abs(copol_component(pk[0], pk[1], "lhcp")) ** 2
                rh_p = np.abs(copol_component(pk[0], pk[1], "rhcp")) ** 2
                rh_m = np.abs(copol_component(em[0], em[1], "rhcp")) ** 2
                tg_, pg_ = np.meshgrid(theta, phi, indexing="ij")
                k_l, k_r = np.unravel_index(np.argmax(lh_p), lh_p.shape), np.unravel_index(np.argmax(rh_p), rh_p.shape)
                t0_, p0_ = np.radians(30.0), np.radians(phis[j])
                dist_ = functools.partial(_deg_from_beam, tg_, pg_, t0_, p0_)
                mt["pickup_only"] = {
                    "lhcp_peak_rel_main_rhcp_peak_db": float(10 * np.log10(lh_p.max() / rh_m.max())),
                    "lhcp_peak_theta_phi_deg": [float(tg_[k_l]), float(pg_[k_l])],
                    "lhcp_peak_deg_from_beam": dist_(k_l),
                    "rhcp_peak_rel_main_rhcp_peak_db": float(10 * np.log10(rh_p.max() / rh_m.max())),
                    "rhcp_peak_theta_phi_deg": [float(tg_[k_r]), float(pg_[k_r])],
                    "rhcp_peak_deg_from_beam": dist_(k_r),
                    "lhcp_fraction_of_pickup": float(np.sum(lh_p) / np.sum(lh_p + rh_p)),
                }
                row[str(j)] = mt
            res[f"board_{board}"]["ports_alone"][f"s={s}"] = row
    # ---- Butler inputs for board B -----------------------------------------------------------------
    order = [2, 3, 1, 4]  # Butler outputs 1..4 = ports by increasing x
    steps = {"1L": -45.0, "2R": 135.0, "2L": -135.0, "1R": 45.0}
    res["butler"] = {
        "output_order_ports": order,
        "x_mm_of_outputs": [FEEDS[p][0] for p in order],
        "steps_deg": steps,
        "by_scale": {},
    }
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
            du = -np.radians(step) / (K0 * d_m)
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
    return res, cells_by_board, fields


def single_hexagon_cells(table, eta_design=0.9, taper_c=0.5, beam_deg=(30.0, 225.0), hand=+1, n_psi=720):
    """The selected single-hexagon design only (no pattern scan): cell list `[{x_mm, y_mm, slot_L_mm, slot_W_mm, slot_angle_deg, alpha_np_per_m}]`
    identical to `design_single_hexagon`'s `selected` cells for the same parameters."""
    poly = hexagon()
    x, y = lattice_cells(poly)
    th0, ph0 = np.radians(beam_deg[0]), np.radians(beam_deg[1])
    m = march_rays(table, poly, (0.0, 0.0), K0, th0, ph0, eta_design, taper_c=taper_c, n_psi=n_psi)
    d = cell_design(m, table, x, y, (0.0, 0.0), K0, th0, ph0, hand=hand)
    return [
        {
            "x_mm": float(x[i] * 1e3),
            "y_mm": float(y[i] * 1e3),
            "slot_L_mm": float(d["slot_L_mm"][i]),
            "slot_W_mm": float(d["slot_W_mm"][i]),
            "slot_angle_deg": float(np.degrees(d["slot_angle_rad"][i])),
            "alpha_np_per_m": float(d["alpha"][i]),
        }
        for i in range(x.size)
    ]


def single_hex_predictions(res):
    """`predictions.json` content of the single-hexagon export from the design-scan result `res` (`tensor_cp_hexagon.json`)."""
    sel = res["selected"]
    pred = {
        "labels": res["labels"],
        "state": "PREDICTION_NOT_HFSS_VALIDATED",
        "beam_target_theta_phi_deg": [30.0, 225.0],
        "frequency_ghz": 10.54,
        "selected_design": {
            k: sel[k] for k in ("hand", "eta_design", "taper_c", "eta_conv", "eta_conv_by_ray_min_max")
        },
        "aperture_pattern": sel["metrics_ap"],
        "realized_gain_rhcp_dbi": {
            "launch_bare_0.592": sel["gain_rhcp_dbi_launch_bare"],
            "launch_loaded_0.66_to_0.88": sel["gain_rhcp_dbi_launch_loaded"],
        },
        "accepted_power_budget": res["accepted_power_budget"],
        "ideal_polarization_compensation_upper_bound": res["ideal_polarization_compensation"]["metrics_ap"],
        "frequencies_10p33_10p75": "NOT_AVAILABLE_FROM_CURRENT_DATA (supercell table is at 10.54 GHz only)",
    }
    return pred


def board_a_predictions(design):
    a = design["board_A"]
    pa = a["ports_alone"]
    cond = {}
    for s_key in pa:
        cond[s_key] = {}
        for j, v in pa[s_key].items():
            cond[s_key][j] = {
                "beam_target_theta_phi_deg": [30.0, PHI[int(j)]],
                "model_peak_theta_phi_deg": [round(v["peak_theta_deg"], 1), round(v["peak_phi_deg"], 1)],
                "realized_gain_rhcp_dbi_anchored": round(v["rhcp_gain_dbi"], 2),
                "xpd_at_peak_db": round(v["xpd_db"], 1),
                "ar_at_peak_db": round(v["ar_peak_db"], 2),
                "ar_worst_in_minus3db_beam_db": round(v["ar_worst_in_minus3db_beam_db"], 2),
                "main_only_gain_dbi": round(v["main_only"]["rhcp_gain_dbi"], 2),
                "pickup_over_main_radiated_power": round(v["pickup_over_main_radiated_power"], 3),
                "check_item_neighbour_lhcp_lobe": {
                    "rel_to_main_rhcp_peak_db": round(v["pickup_only"]["lhcp_peak_rel_main_rhcp_peak_db"], 1),
                    "theta_phi_deg": v["pickup_only"]["lhcp_peak_theta_phi_deg"],
                    "deg_from_beam": round(v["pickup_only"]["lhcp_peak_deg_from_beam"], 1),
                },
            }
    return cond


def fourhex_a_predictions(design, scope, n_cells):
    """`predictions.json` content of the board A export (`design`: `fourhex_design.json`, `scope`: `fourhex_scope.json`)."""
    cond = board_a_predictions(design)
    pred = {
        "labels": design["labels"],
        "state": "PREDICTION_NOT_HFSS_VALIDATED",
        "board": "A: static beams phi 45/135/225/315, theta 30, RHCP",
        "frequency_ghz": 10.54,
        "n_cells": n_cells,
        "feed": {"pin_top_above_surface_mm": 6.0, "mono_h_from_ground_mm": 7.575, "r_pin_mm": 0.635},
        "per_port_alone_by_edge_power_scale": cond,
        "edges_unconverted_power_fraction_of_launched_wave": design["board_A"][
            "crossing_unconverted_power_fraction_of_launched_wave"
        ],
        "edge_note": "multiply by the launch fraction 0.62 for the fraction of the accepted power; E-line of the single hexagon: rim/inner power ratio "
        "0.7 +- factor 2 vs model 0.54 (non-discriminating): s=1.7 is the upper reading",
        "coupling": {
            "adjacent_pairs": {
                k: {
                    kk: vv for kk, vv in v.items() if kk in ("distance_mm", "raw_estimate_db", "calibrated_estimate_db")
                }
                for k, v in scope["coupling"]["estimates"].items()
                if v.get("adjacent")
            },
            "one_bit_board_hfss_range_db": [-73.5, -51.5],
            "state": "ESTIMATE_BIAS_FROM_ONE_BOARD",
        },
        "hfss_size": {
            "cells": n_cells,
            "ratio_vs_single_hexagon_cells": round(n_cells / 1430.0, 2),
            "note": "scale the user's single-hexagon tetrahedra/RAM/time by ~4x / ~4x / 5-8x; NOT measured",
        },
        "check_items": [
            "gain per port within +-1 dB of the table (23.0 dBi)",
            "neighbour LHCP lobe level/direction (table)",
            "adjacent-pair coupling vs estimate",
            "AR at peak <= 3 dB per port",
            "S11 per port vs the single-hexagon -15 dB at 10.54 GHz",
        ],
    }
    return pred
