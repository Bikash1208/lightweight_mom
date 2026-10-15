"""Phase 26 C: pure-RHCP rotating-slot tensor hologram on ONE 150 mm hexagon (RT/duroid 5880 1.575 mm).

    python scripts/research/phase26_tensor_cp_hexagon.py <out.json> [--quick]

Design from the SOLVER_SUPERCELL table (`supercell_table_*.json`, nearest grid, mesh-stable slots only: bridges
>= 0.43 mm), mesh pitch 6.4 mm, D 6.0 mm circle, hexagon circumradius 150 mm (vertex at +/-y, the orientation of
port 3 of `four_port_hex150_gap0`), feed at the centre, beam theta 30 deg / phi 225 deg, 10.54 GHz.
A scan over the design conversion `eta_design` and amplitude taper picks the design of the highest predicted
gain whose RHCP SLL is acceptable; every scanned design is reported. Everything is MODEL_ONLY.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.frequency import free_space_wavenumber  # noqa: E402
from src.green_function import LayeredGreenFunction  # noqa: E402
from src.hologram.leaky_forward import pattern_budget  # noqa: E402
from src.hologram.tensor_leaky_cp import SupercellTable, cell_design, march_rays  # noqa: E402
from src.materials import Substrate  # noqa: E402
from src.pattern.vector_aperture_pattern import _vector_far_field, circular_components, copol_component  # noqa: E402

F = 10.54e9
K0 = free_space_wavenumber(F)
PITCH = 6.4e-3
SUB = Substrate("rt5880", 2.2, 0.0009, 1.575e-3)
R_HEX = 0.150
THETA0, PHI0 = np.radians(30.0), np.radians(225.0)
LAUNCH_BARE, LAUNCH_LOADED = 0.592, (0.66, 0.88)
RHO_CLEAR = 6.5e-3
TABLE_GLOB = str(ROOT / "results/phase26/supercell_table_*.json")


def hexagon():
    return np.array([[R_HEX * np.cos(np.radians(90 + 60 * i)), R_HEX * np.sin(np.radians(90 + 60 * i))] for i in range(6)])


def lattice_cells(poly):
    n = 30
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


def metrics(et, ep, theta, phi, label_launch=LAUNCH_BARE, eta_conv=1.0):
    b = pattern_budget(theta, phi, et, ep, 30.0, 225.0, (8.0, 20.0))
    er, el = circular_components(et, ep)
    rh = np.abs(copol_component(et, ep, "rhcp")) ** 2
    lh = np.abs(copol_component(et, ep, "lhcp")) ** 2
    tg, pg = np.meshgrid(theta, phi, indexing="ij")
    t0, p0 = np.radians(30.0), np.radians(225.0)
    cosd = np.cos(np.radians(tg)) * np.cos(t0) + np.sin(np.radians(tg)) * np.sin(t0) * np.cos(np.radians(pg) - p0)
    dist = np.degrees(np.arccos(np.clip(cosd, -1, 1)))
    near = dist <= 6.0
    i = np.unravel_index(np.argmax(np.where(near, rh, -1)), rh.shape)
    peak = rh[i]
    # -3 dB beam region around the peak (within 10 deg of the beam direction)
    pk_dir = (np.radians(tg[i]), np.radians(pg[i]))
    cd = np.cos(np.radians(tg)) * np.cos(pk_dir[0]) + np.sin(np.radians(tg)) * np.sin(pk_dir[0]) * np.cos(np.radians(pg) - pk_dir[1])
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


def main(out_path, quick=False):
    t0 = time.time()
    table = SupercellTable.from_files(TABLE_GLOB)
    poly = hexagon()
    x, y = lattice_cells(poly)
    green = LayeredGreenFunction(SUB, F, backing="pec")
    theta = np.arange(0.0, 90.0 + 1e-9, 2.0 if quick else 1.0)
    phi = np.arange(0.0, 360.0, 3.0 if quick else 2.0)
    res = {"labels": ["SOLVER_SUPERCELL", "MODEL_ONLY", "NOT_HFSS_VALIDATED", "GAMMA_0_COLLINEAR_ONLY"],
           "n_cells": int(x.size), "hexagon_circumradius_mm": R_HEX * 1e3, "pitch_mm": PITCH * 1e3, "slots_mm": table.slots,
           "k_grid_over_k0": table.k_grid.tolist(), "launch_fraction_bare": LAUNCH_BARE, "launch_fraction_loaded_range": LAUNCH_LOADED,
           "scan": []}
    best = None
    for hand in (+1,):
        for eta_design in (0.5, 0.7, 0.9):
            for taper in (0.0, 0.5, 0.8):
                m = march_rays(table, poly, (0.0, 0.0), K0, THETA0, PHI0, eta_design, taper_c=taper, n_psi=360 if quick else 720)
                d = cell_design(m, table, x, y, (0.0, 0.0), K0, THETA0, PHI0, hand=hand)
                et, ep = pattern(d, x, y, green, theta, phi)
                mt = metrics(et, ep, theta, phi)
                eta = m["eta_conv"]
                rec = {"hand": hand, "eta_design": eta_design, "taper_c": taper, "eta_conv": eta,
                       "eta_conv_by_ray_min_max": [float(m["eta_conv_by_ray"].min()), float(m["eta_conv_by_ray"].max())],
                       "metrics_ap": mt,
                       "gain_rhcp_dbi_launch_bare": mt["directivity_ap_rhcp_dbi"] + 10 * np.log10(LAUNCH_BARE * eta),
                       "gain_rhcp_dbi_launch_loaded": [mt["directivity_ap_rhcp_dbi"] + 10 * np.log10(l * eta) for l in LAUNCH_LOADED]}
                res["scan"].append(rec)
                print(f"hand {hand:+d} eta_d {eta_design} taper {taper}: eta_conv {eta:.3f}  D_ap {mt['directivity_ap_rhcp_dbi']:.2f}  "
                      f"G {rec['gain_rhcp_dbi_launch_bare']:.2f}  XPD {mt['xpd_at_peak_db']:.1f}  RHCPfrac {mt['rhcp_fraction']:.3f}  "
                      f"SLL {mt['rhcp_sll_db']:.1f}  peak ({mt['peak_rhcp_theta_deg']:.0f},{mt['peak_rhcp_phi_deg']:.0f})", flush=True)
                key = (mt["rhcp_fraction"] > 0.5, mt["rhcp_sll_db"] <= -13.0, rec["gain_rhcp_dbi_launch_bare"])
                if best is None or key > best[0]:
                    best = (key, rec, d, m)
        Path(out_path).write_text(json.dumps(res, indent=1, default=float))
    _, rec, d, m = best
    res["selected"] = {k: v for k, v in rec.items()}
    # ---- budget of the ACCEPTED power (assumptions stated, nothing invented) ---------------------
    mt, eta = rec["metrics_ap"], rec["eta_conv"]
    bud = {}
    for tag, lf in (("launch_bare_0.592", LAUNCH_BARE), ("launch_loaded_0.66", LAUNCH_LOADED[0]), ("launch_loaded_0.88", LAUNCH_LOADED[1])):
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
        "label": "UPPER_BOUND_NOT_REALISED_BY_THE_CELL_FAMILY", "metrics_ap": mc,
        "note": "J' = g(k_t,beam)^-1 J per cell: removes the grounded slab's TM/TE element-factor ellipticity at the beam "
                "direction; needs an independently controlled second modulation (not designed here)."}
    print("uncompensated XPD %.1f AR %.2f | ideal-compensated XPD %.1f AR %.2f dB RHCPfrac %.3f"
          % (mt["xpd_at_peak_db"], mt["ar_at_peak_db"], mc["xpd_at_peak_db"], mc["ar_at_peak_db"], mc["rhcp_fraction"]), flush=True)
    sel = {
        "x_mm": (x * 1e3).tolist(), "y_mm": (y * 1e3).tolist(), "slot_L_mm": d["slot_L_mm"].tolist(), "slot_W_mm": d["slot_W_mm"].tolist(),
        "slot_angle_deg": np.degrees(d["slot_angle_rad"]).tolist(), "alpha_np_per_m": d["alpha"].tolist(),
        "beta_over_k0": d["beta_over_k0"].tolist(), "amp": d["amp"].tolist(),
    }
    Path(out_path).with_name(Path(out_path).stem + "_cells.json").write_text(json.dumps(sel))
    res["seconds"] = time.time() - t0
    Path(out_path).write_text(json.dumps(res, indent=1, default=float))
    print("selected:", {k: rec[k] for k in ("hand", "eta_design", "taper_c", "eta_conv", "gain_rhcp_dbi_launch_bare")})


if __name__ == "__main__":
    main(sys.argv[1], "--quick" in sys.argv)
