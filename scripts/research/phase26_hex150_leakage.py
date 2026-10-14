"""Phase 26 A: re-predict the existing hex150 board's port 3 with the leakage model and set it against the
user's HFSS run (`hfss_results/hex150_gap0/`).

    python scripts/research/phase26_hex150_leakage.py <out.json> [--freq-ghz 10.5] [--quick]

Everything printed is MODEL_ONLY except the HFSS column, which is read from the committed antenna-parameter
table and realized-gain CSV and integrated here (2 deg grid, +-5 %).
"""

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.constants import ETA0  # noqa: E402
from src.dispersion.analytic_tm0 import tm0_beta  # noqa: E402
from src.dispersion.leakage_model import (  # noqa: E402
    bound_beta,
    floquet_complex_beta,
    leakage_alpha_first_order,
    radiated_harmonic_geometry,
    state_surface_reactances,
    surface_wave_power_per_h2,
    two_level_c1,
    two_level_fourier,
)
from src.dispersion.sheet_model import sheet_admittance_from_gamma  # noqa: E402
from src.frequency import free_space_wavenumber  # noqa: E402
from src.hologram.auto_states import auto_states  # noqa: E402
from src.hologram.butler_forward import board_from_cells, far_field  # noqa: E402
from src.hologram.leaky_forward import leakage_profile, pattern_budget, scaled_board  # noqa: E402
from src.materials import Substrate  # noqa: E402

EXPORT = ROOT / "hfss_exports" / "four_port_hex150_gap0"
HFSS = ROOT / "hfss_results" / "hex150_gap0"
EPS_R, TAND, H_M, PERIOD, SQUARE_MM = 2.2, 0.0009, 1.575e-3, 6.4e-3, 5.85
PORT, BEAM_THETA, BEAM_PHI = 3, 30.0, 225.0
HFSS_CONES = (8.0, 20.0)


# ----------------------------------------------------------------------------------- the board
def load_port(port=PORT):
    """Port `port` of the exported board: its owned cells (centres inside the hexagon, on the lattice the
    patches sit on) and the state of each (1 = a square patch is drawn there)."""
    layout = json.loads((EXPORT / "layout.json").read_text())
    outline = np.array([[float(a), float(b)] for a, b in layout["port_outlines"][port - 1]]) * 1e-3
    feed = next(f for f in layout["feeds"] if f["port"] == port)
    feed_xy = (float(feed["x_mm"]) * 1e-3, float(feed["y_mm"]) * 1e-3)
    patches = [p for p in layout["patches"] if p["port"] == port]
    centres = {
        (round((float(p["x_min_mm"]) + float(p["x_max_mm"])) / 2, 3), round((float(p["y_min_mm"]) + float(p["y_max_mm"])) / 2, 3))
        for p in patches
    }
    cx0, cy0 = next(iter(centres))
    n = 80
    ks = np.arange(-n, n + 1)
    xg, yg = np.meshgrid(cx0 * 1e-3 + PERIOD * ks, cy0 * 1e-3 + PERIOD * ks)
    inside = _inside(outline, xg, yg)
    x, y = xg[inside], yg[inside]
    keys = [(round(a * 1e3, 3), round(b * 1e3, 3)) for a, b in zip(x, y, strict=True)]
    missing = centres - set(keys)
    if missing:
        raise RuntimeError(f"{len(missing)} patches are not on the reconstructed lattice (first {list(missing)[:3]})")
    level = np.array([k in centres for k in keys], dtype=int)
    return {"x": x, "y": y, "level": level, "feed": feed_xy, "outline": outline}


def _inside(poly, x, y):
    n = len(poly)
    res = np.ones(x.shape, dtype=bool)
    for i in range(n):
        a, b = poly[i], poly[(i + 1) % n]
        cross = (b[0] - a[0]) * (y - a[1]) - (b[1] - a[1]) * (x - a[0])
        res &= cross <= 1e-9  # clockwise or counter-clockwise: fixed below
    if res.any():
        return res
    res = np.ones(x.shape, dtype=bool)
    for i in range(n):
        a, b = poly[i], poly[(i + 1) % n]
        cross = (b[0] - a[0]) * (y - a[1]) - (b[1] - a[1]) * (x - a[0])
        res &= cross >= -1e-9
    return res


# ----------------------------------------------------------------------------------- HFSS side
def hfss_numbers():
    rows = list(csv.reader(open(HFSS / "port3_realized_gain_RHCP_LHCP_10p5GHz.csv")))[1:]
    a = np.array([[float(v) for v in r] for r in rows])
    ph, th, gl, gr = a[:, 1], a[:, 2], a[:, 3], a[:, 4]
    u, v = np.unique(th), np.unique(ph)
    pl = np.zeros((v.size, u.size))
    pr = np.zeros_like(pl)
    iv, iu = np.searchsorted(v, ph), np.searchsorted(u, th)
    pl[iv, iu] = 10 ** (gl / 10)
    pr[iv, iu] = 10 ** (gr / 10)
    tg, pg = np.meshgrid(np.radians(u), np.radians(v))
    t0, p0 = np.radians(BEAM_THETA), np.radians(224.0)
    dist = np.degrees(np.arccos(np.clip(np.cos(tg) * np.cos(t0) + np.sin(tg) * np.sin(t0) * np.cos(pg - p0), -1, 1)))
    w = np.sin(tg) * np.radians(2) ** 2
    w[0, :] *= 0.5
    w[-1, :] *= 0.5
    tot = float(np.sum((pl + pr) * w) / (4 * np.pi))
    out = {
        "peak_realized_gain_rhcp_dbi": float(10 * np.log10(pr.max())),
        "peak_at_theta_phi_deg": [float(u[np.unravel_index(pr.argmax(), pr.shape)[1]]), float(v[np.unravel_index(pr.argmax(), pr.shape)[0]])],
        "integrated_total_over_accepted": tot,
        "peak_directivity_dbi_table": 10 * np.log10(99.40236),
        "total_efficiency_table": 0.977420,
        "peak_xpd_at_rhcp_peak_db": float(
            10 * np.log10(pr.max() / pl[np.unravel_index(pr.argmax(), pr.shape)])
        ),
        "lhcp_fraction": float(np.sum(pl * w) / (4 * np.pi) / tot),
        "above_70deg_fraction": float(np.sum((pl + pr) * w * (np.degrees(tg) > 70)) / (4 * np.pi) / tot),
        "back_fraction": float(np.sum((pl + pr) * w * (np.degrees(tg) > 90)) / (4 * np.pi) / tot),
    }
    for c in HFSS_CONES:
        m = dist <= c
        out[f"within_{c:g}deg_fraction"] = float(np.sum((pl + pr) * w * m) / (4 * np.pi) / tot)
        out[f"within_{c:g}deg_rhcp_fraction"] = float(np.sum(pr * w * m) / (4 * np.pi) / tot)
    return out


# ----------------------------------------------------------------------------------- the model
def mom_context_ratios():
    """kappa = (in-context contrast)/(isolated contrast) from the Phase 25 supercell MoM: the specular order
    of the 1x2 stripe and 2x2 checkerboard supercells of the design pair (50 % square) against the
    isolated-cell admittance average."""
    d = json.loads((ROOT / "results" / "phase25" / "block_locality.json").read_text())
    sub = Substrate("user", EPS_R, TAND, H_M)
    f = 10.54e9
    out = {}
    for mesh, m in d["meshes"].items():
        g0 = m["isolated_empty"][0] * np.exp(1j * np.radians(m["isolated_empty"][1]))
        g1 = m["isolated_square"][0] * np.exp(1j * np.radians(m["isolated_square"][1]))
        y0, y1 = (sheet_admittance_from_gamma(g, sub, f) for g in (g0, g1))
        for case in ("1x2_stripes", "2x2_checkerboard"):
            mag, ph = m["cases"][case]["orders"]["0,0"]["mom"]
            gm = mag * np.exp(1j * np.radians(ph))
            ym = sheet_admittance_from_gamma(gm, sub, f)
            kappa = float(np.real(2.0 * (ym - y0) / (y1 - y0)))
            out[f"{mesh}:{case}"] = kappa
    return out


def run(freq_hz, quick=False):
    t0 = time.time()
    sub = Substrate("user", EPS_R, TAND, H_M)
    k0 = free_space_wavenumber(freq_hz)
    port = load_port()
    states = auto_states(1, sub, freq_hz, PERIOD, fixed_L_mm=SQUARE_MM)
    g0, g1 = (st.gamma for st in states.states)
    y0 = sheet_admittance_from_gamma(g0, sub, freq_hz)
    y1 = sheet_admittance_from_gamma(g1, sub, freq_hz)
    duty = float(port["level"].mean())
    kappas = mom_context_ratios()
    ctx_path = ROOT / "results" / "phase26" / "context_contrast.json"
    if ctx_path.exists():
        for k, v in json.loads(ctx_path.read_text())["cases"].items():
            kappas[f"mesh24:{k}"] = v["kappa"]
    # the homogenisation is only meaningful for supercells well below a wavelength (period < ~0.7 lambda):
    valid = {k: v for k, v in kappas.items() if "period4" not in k}
    kappa_ctx = float(min(valid.values()))
    kappa_hi = float(max(valid.values()))
    res = {
        "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED"],
        "frequency_ghz": freq_hz / 1e9,
        "port": PORT,
        "cells": int(port["level"].size),
        "patch_fraction": duty,
        "state_phase_deg": [float(np.degrees(np.angle(g))) for g in (g0, g1)],
        "y_sheet_empty_state_S": [y0.real, y0.imag],
        "y_sheet_square_state_S": [y1.real, y1.imag],
        "context_contrast_ratios_phase25_mom": kappas,
        "context_contrast_ratio_range_used": [kappa_ctx, kappa_hi],
        "context_contrast_note": "period-4 (25.6 mm = 0.9 lambda) stripes excluded: not a sub-wavelength supercell",
    }
    beta_bare = tm0_beta(freq_hz, EPS_R, H_M)
    res["beta_bare_over_k0"] = beta_bare / k0

    th0, ph0 = np.radians(BEAM_THETA), np.radians(BEAM_PHI)
    kt_beam = k0 * np.sin(th0)
    kvec_beam = kt_beam * np.array([np.cos(ph0), np.sin(ph0)])
    cases = {}
    for name, kappa in (
        ("ISOLATED_CELL_CONTRAST", 1.0),
        ("MOM_CONTEXT_CONTRAST_LOW", kappa_ctx),
        ("MOM_CONTEXT_CONTRAST_HIGH", kappa_hi),
    ):
        # the contrast of the square state vs empty is scaled; the empty state stays the bare slab
        y1c = y0 + kappa * (y1 - y0)
        y_mean = (1.0 - duty) * y0 + duty * y1c
        beta = bound_beta(y_mean, sub, freq_hz)
        x_states = state_surface_reactances([y0, y1c], beta, sub, freq_hz)
        c1 = two_level_c1(y0, y1c, duty)
        psi_g = np.radians(np.arange(0.0, 360.0, 15.0))
        a_closed, a_hb, k_mods, gam = [], [], [], []
        for psi in psi_g:
            rhat = np.array([np.cos(psi), np.sin(psi)])
            kmod = beta * rhat - kvec_beam
            km = float(np.hypot(*kmod))
            gamma = float(np.arccos(np.clip(np.dot(rhat, kmod) / km, -1, 1)))
            zeta = float(ph0 - psi)
            a_closed.append(
                leakage_alpha_first_order(beta, y_mean, c1, th0, zeta, sub, freq_hz)
            )
            k_mods.append(km)
            gam.append(gamma)
            if not quick:
                yc = two_level_fourier(y0, y1c, 6, duty)
                b, _ = floquet_complex_beta(yc, 2 * np.pi / km, gamma, sub, freq_hz, beta, n_max=6)
                a_hb.append(float(-b.imag))
        a_closed = np.array(a_closed)
        profile = leakage_profile(
            port["x"], port["y"], port["feed"], port["outline"],
            lambda p, a=a_closed, g=psi_g: np.interp(np.mod(p, 2 * np.pi), np.append(g, 2 * np.pi), np.append(a, a[0])),
            PERIOD,
        )
        cases[name] = {
            "kappa": kappa,
            "beta_over_k0": beta / k0,
            "x_bar_over_eta": float(np.real(np.mean(x_states)) / ETA0),
            "x_states_over_eta": [float(v / ETA0) for v in x_states],
            "surface_wave_w_ohm_m": surface_wave_power_per_h2(beta, EPS_R, H_M, freq_hz),
            "alpha_closed_np_per_m": {"min": float(a_closed.min()), "mean": float(a_closed.mean()), "max": float(a_closed.max())},
            "alpha_floquet_total_np_per_m": None
            if not a_hb
            else {"min": float(min(a_hb)), "mean": float(np.mean(a_hb)), "max": float(max(a_hb))},
            "eta_conv_closed": profile.eta_conv,
            "edge_fraction": profile.edge_fraction,
            "k_mod_over_k0": [float(min(k_mods) / k0), float(max(k_mods) / k0)],
            "alpha_closed_by_psi_deg": {f"{int(round(np.degrees(p)))}": float(a) for p, a in zip(psi_g, a_closed, strict=True)},
            "alpha_floquet_by_psi_deg": {f"{int(round(np.degrees(p)))}": float(a) for p, a in zip(psi_g, a_hb, strict=True)}
            if a_hb
            else None,
            "eta_conv_floquet_total": float(np.mean(1.0 - np.exp(-2.0 * np.interp(np.mod(profile.psi_grid_rad, 2 * np.pi), np.append(psi_g, 2 * np.pi), np.append(a_hb, a_hb[0])) * profile.rho_edge_grid_m)))
            if a_hb
            else None,
            "_profile": profile,
            "_alpha_psi": a_closed,
        }
    # ---- forward patterns (aperture only) ---------------------------------------------------
    gamma_cell = np.where(port["level"] == 1, g1, g0)
    gamma_cell = gamma_cell / np.abs(gamma_cell)  # the phase-only convention of every earlier pattern
    board = board_from_cells(
        port["x"], port["y"], gamma_cell, [port["feed"]], PERIOD, freq_hz, sub, beta_bare, "rhcp", "phase_only"
    )
    theta = np.arange(0.0, 90.0 + 1e-9, 1.0 if not quick else 2.0)
    phi = np.arange(0.0, 360.0, 1.0 if not quick else 2.0)
    tg, pg = np.meshgrid(theta, phi, indexing="ij")

    def pattern(b):
        et, ep = far_field(b, [1.0], tg.ravel(), pg.ravel())
        return pattern_budget(theta, phi, et.reshape(tg.shape), ep.reshape(tg.shape), BEAM_THETA, 224.0, HFSS_CONES)

    res["baseline_phase_only"] = pattern(board)
    for name, c in cases.items():
        # the illumination taper changes the aperture directivity by < 0.1 dB between these cases
        c["pattern_ap"] = pattern(scaled_board(board, c["_profile"].amplitude))
    res["cases"] = {}
    for name, c in cases.items():
        d_ap = c["pattern_ap"]["directivity_ap_rhcp_dbi"]
        eta = c["eta_conv_closed"]
        gain_ub = d_ap + 10 * np.log10(max(eta, 1e-12))
        entry = {k: v for k, v in c.items() if not k.startswith("_")}
        entry["gain_upper_bound_dbi"] = float(gain_ub)
        entry["note_gain"] = "directivity_ap + 10 log10(eta_conv); launch_fraction = 1 (upper bound)"
        res["cases"][name] = entry
    # ---- what alpha would the HFSS gain need? ----------------------------------------------
    hf = hfss_numbers()
    res["hfss"] = hf
    c = cases["ISOLATED_CELL_CONTRAST"]
    d_ap = c["pattern_ap"]["directivity_ap_rhcp_dbi"]
    need = 10 ** ((hf["peak_realized_gain_rhcp_dbi"] - d_ap) / 10.0)  # eta_launch * eta_conv
    pg_edge = c["_profile"].rho_edge_grid_m
    a_mean_needed = -np.log(max(1.0 - need, 1e-12)) / (2.0 * float(np.mean(pg_edge)))
    res["implied_by_hfss"] = {
        "eta_launch_times_eta_conv": need,
        "alpha_needed_if_launch_1_np_per_m": a_mean_needed,
        "alpha_ratio_to_isolated_closed": a_mean_needed / float(c["_alpha_psi"].mean()),
    }
    d_iso = c["pattern_ap"]["directivity_ap_rhcp_dbi"]
    sweep = {}
    for kap in (0.8, 0.86, 0.9, 0.97, 1.0, 1.04, 1.1, 1.2):
        y1c = y0 + kap * (y1 - y0)
        y_mean = (1.0 - duty) * y0 + duty * y1c
        beta = bound_beta(y_mean, sub, freq_hz)
        c1 = two_level_c1(y0, y1c, duty)
        psi_g = np.radians(np.arange(0.0, 360.0, 15.0))
        a = np.array([leakage_alpha_first_order(beta, y_mean, c1, th0, float(ph0 - p), sub, freq_hz) for p in psi_g])
        prof = leakage_profile(
            port["x"], port["y"], port["feed"], port["outline"],
            lambda p, a=a, g=psi_g: np.interp(np.mod(p, 2 * np.pi), np.append(g, 2 * np.pi), np.append(a, a[0])),
            PERIOD,
        )
        sweep[f"{kap:g}"] = {
            "alpha_mean_np_per_m": float(a.mean()),
            "eta_conv": prof.eta_conv,
            "gain_upper_bound_dbi": float(d_iso + 10 * np.log10(prof.eta_conv)),
        }
    res["gain_upper_bound_vs_contrast_ratio"] = sweep
    res["gain_vs_launch_fraction_at_isolated_contrast"] = {
        f"{lf:g}": float(res["cases"]["ISOLATED_CELL_CONTRAST"]["gain_upper_bound_dbi"] + 10 * np.log10(lf))
        for lf in (1.0, 0.75, 0.5, 0.35)
    }
    res["seconds"] = time.time() - t0
    return res


def main():
    args = sys.argv[1:]
    freq = 10.5e9
    if "--freq-ghz" in args:
        freq = float(args[args.index("--freq-ghz") + 1]) * 1e9
    out = Path(args[0])
    res = run(freq, quick="--quick" in args)
    out.write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
