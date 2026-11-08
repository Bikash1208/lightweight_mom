"""SUB6 board model with HFSS-calibrated switched (island) cells (MODEL_ONLY, PENDING more HFSS data).  Job sub6-board-island-hfsscal.

Base: sub6_board_island2.py / sub6_board_model4.py (V3 board synthesis, unchanged and imported).  Changes (module-level patches, applied in every worker):
* the switched table is restricted to the ONE row that has HFSS data: ISL_W2.9_q0.8_Li0.5_c0.6 (mode `base` = uncalibrated MoM numbers with this single row; mode `cal` = calibrated);
* CALIBRATION per state and direction d (along the slot / across it) from `hfss_data.json` (ISL_vias2 = vias + bias stubs, FR4 0.8): the MoM dispersion f_M(phase) of the island cell (pred_table, 7 frequencies,
  mesh 64) is corrected by a fit of the HFSS points,  f_cal(phi) = f_M(phi) + a + b (phi - phi_ref)   (ONE point: shift a only; >= 2 points with different phases: shift + slope b);
  beta_d/k0 at 5.8 GHz = phi_cal(5.8 GHz) / (p k0);  Q mean of the HFSS Q values;  v_g from the slope of f_cal;  extra loss  alpha_x = omega / (2 Q v_g) = f (d beta / d f) / (2 Q)  [Np/m]
  (absorbed in the ray march, NOT radiated: the MoM leakage alpha of the row is kept and is MoM-only);
* direction dependence on the board: for every switched site the slot axis a follows from the hologram direction (before the 90-degree rotation flips) and theta' = phi_ray - a; beta(theta') =
  beta_along cos^2 theta' + beta_across sin^2 theta', alpha_x likewise; iterated 3 times (axes depend on the phase, the phase on beta).  The 90-degree flips of the V0-V3 optimisation are NOT fed back into beta
  (a flipped slot would swap along / across): flagged, the number of rotated switched cells is small;
* state ON: not calibrated yet. `placeholder_on = True`: ON uses the MoM beta and alpha of the row and alpha_x = 0 (labelled in every output).
Radiated power per ray = sum 2 alpha_leak W drho; eta = mean of that; the extra loss only lowers W (so the gain is referred to the accepted power as before).
    python sub6_board_hfsscal.py run [--workers N]
"""
import argparse
import csv
import json
import math
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from scipy.optimize import brentq

import sub6_board_island2 as X
import sub6_board_model4 as R4

B, M = R4.B, R4.M
K0 = R4.K0
RES = M.RES.parent / "sub6_board_island_hfsscal"
HFSS = RES / "hfss_data.json"
MOM_TABLE = RES / "input_pred_table_mom_nodata.json"
ROW_ID = "ISL_W2.9_q0.8_Li0.5_c0.6"
P_MM = 11.6
LAUNCH = R4.LAUNCH
MODE = "base"  # set per task: "base" | "cal"
CAL = None  # loaded lazily
_orig = dict(tables=R4.tables, site_arrays=R4.site_arrays, physics=R4.physics)


# ------------------------------------------------------------------ calibration from the HFSS points
def mom_curve(state, d):
    t = json.load(open(MOM_TABLE))
    rows = t[f"ISL|{state}|{d}"]
    f = np.array([r["f_GHz"] for r in rows])
    ph = np.array([r["phase_deg"] for r in rows])
    o = np.argsort(f)
    return f[o], ph[o]


def calibrate(state="OFF_BAR50"):
    """-> {direction: dict(beta58, alpha_x58, a, b, n_points, Q, slope_deg_per_GHz, ...)} or None for a state without HFSS data."""
    data = json.load(open(HFSS))
    out = {}
    for d, key in (("x", "along"), ("y", "across")):
        pts = [p for p in data["points"] if p["state"] == state and p["direction"] == key and p.get("f_GHz") is not None]
        if not pts:
            return None
        f, ph = mom_curve(state, d)
        fM = lambda phi: float(np.interp(phi, ph, f))  # noqa: E731
        phis = np.array([p["phase_deg"] for p in pts], float)
        df = np.array([p["f_GHz"] - fM(p["phase_deg"]) for p in pts])
        ref = float(phis.mean())
        if len(pts) >= 2 and np.ptp(phis) > 1.0:
            A = np.stack([np.ones_like(phis), phis - ref], 1)
            a, b = np.linalg.lstsq(A, df, rcond=None)[0]
        else:
            a, b = float(df.mean()), 0.0
        fcal = lambda phi: fM(phi) + a + b * (phi - ref)  # noqa: E731
        lo, hi = float(ph.min()), float(ph.max())
        phi58 = brentq(lambda q: fcal(q) - 5.8, lo, hi)
        h = 0.05
        phi_p, phi_m = brentq(lambda q: fcal(q) - (5.8 + h), lo, hi), brentq(lambda q: fcal(q) - (5.8 - h), lo, hi)
        dphi_df = (phi_p - phi_m) / (2 * h)  # deg / GHz
        k0 = 2 * math.pi * 5.8e9 / 299792458.0
        beta58 = math.radians(phi58) / (P_MM * 1e-3) / k0
        dbeta_df = math.radians(dphi_df) / (P_MM * 1e-3)  # rad/m per GHz
        Q = float(np.mean([p["Q"] for p in pts]))
        out[d] = dict(direction=key, n_points=len(pts), shift_a_GHz=float(a), slope_b_GHz_per_deg=float(b), phi_ref_deg=ref, beta_over_k0_5p8=float(beta58), beta_over_k0_MoM_5p8=float(np.interp(5.8, f, [math.radians(p) / (P_MM * 1e-3) / k0 for p in ph])),
                      dphi_df_deg_per_GHz=float(dphi_df), Q=Q, alpha_x_Np_per_m=float(5.8 * dbeta_df / (2 * Q)), points=pts, f_MoM_at_phase_GHz=[fM(p["phase_deg"]) for p in pts])
    return out


def get_cal():
    global CAL
    if CAL is None:
        CAL = {"OFF": calibrate("OFF_BAR50"), "ON": calibrate("ON_BAR50")}
    return CAL


# ------------------------------------------------------------------ patched pieces
def tables_single():
    T = _orig["tables"]()
    sw = T["sw"]
    i = sw["id"].index(ROW_ID)
    new = {}
    for k, v in sw.items():
        if isinstance(v, list):
            new[k] = [v[i]]
        elif isinstance(v, np.ndarray):
            new[k] = v[i:i + 1]
        else:
            new[k] = v
    T["sw"] = new
    return T


def site_arrays_cal(D, over=None):
    sa = _orig["site_arrays"](D, over)
    z = np.zeros_like(sa["alpha_off"])
    sa["loss_off"], sa["loss_on"] = z.copy(), z.copy()
    return sa


def real_march_loss(m, lk, a_leak, a_loss, b_site, b_fill):
    """as sub6_board_model.real_march (alpha decays the ray, beta gives the phase) plus an absorbed extra loss; returns W, phi_ref, radiated fraction, absorbed fraction."""
    psi, rho, drho = m["psi"], m["rho"], m["drho"]
    Xg, Yg = rho[None, :] * np.cos(psi)[:, None], rho[None, :] * np.sin(psi)[:, None]
    idx = lk.at(Xg, Yg)
    inside = rho[None, :] <= m["r_edge"][:, None]
    ok = (idx >= 0) & inside
    a = np.where(ok, a_leak[np.maximum(idx, 0)], 0.0)
    ls = np.where(ok, a_loss[np.maximum(idx, 0)], 0.0)
    b = np.where(idx >= 0, b_site[np.maximum(idx, 0)], b_fill)
    ca, cb = np.cumsum(a + ls, axis=1) - (a + ls), np.cumsum(b, axis=1) - b
    w = np.exp(-2.0 * ca * drho)
    phi = cb * K0 * drho
    rad = float(np.mean(np.sum(2.0 * a * w, axis=1) * drho))
    absd = float(np.mean(np.sum(2.0 * ls * w, axis=1) * drho))
    return w, phi, rad, absd


def physics_run(D, sa):
    S, lk, m0 = D["S"], D["lk"], D["m0"]
    smp = M.sampler(m0, S)
    w_off, ph_off, rad_off, abs_off = real_march_loss(m0, lk, sa["alpha_off"], sa["loss_off"], sa["beta_off"], D["b_fill"])
    w_on, ph_on, rad_on, abs_on = real_march_loss(m0, lk, sa["alpha_on"], sa["loss_on"], sa["beta_on"], D["b_fill"])
    return dict(Woff=smp(w_off), Won=smp(w_on), Poff=smp(ph_off), Pon=smp(ph_on), eta_off=rad_off, eta_on=rad_on, abs_off=abs_off, abs_on=abs_on)


def physics_cal(D, sa):
    """MODE base: the MoM beta / alpha, loss 0 (same march function).  MODE cal: iterate beta(theta') and alpha_x(theta') of the switched sites for the OFF state (ON = MoM placeholder)."""
    if MODE == "cal":
        cal = get_cal()["OFF"]
        mask = D["mask"] & D["live"]
        S = D["S"]
        phiR = np.arctan2(S.y, S.x)
        b_al, b_ac = cal["x"]["beta_over_k0_5p8"], cal["y"]["beta_over_k0_5p8"]
        l_al, l_ac = cal["x"]["alpha_x_Np_per_m"], cal["y"]["alpha_x_Np_per_m"]
        sa["beta_off"] = np.where(mask, 0.5 * (b_al + b_ac), sa["beta_off"])
        for it in range(3):
            ph = physics_run(D, sa)
            direction = X.directions_w(D, ph, 0.0)
            chi = np.arctan2(direction[:, 1], direction[:, 0])
            axis = np.mod(0.5 * (chi + phiR), math.pi)
            c2 = np.cos(phiR - axis) ** 2
            sa["beta_off"] = np.where(mask, b_al * c2 + b_ac * (1 - c2), sa["beta_off"])
            sa["loss_off"] = np.where(mask, l_al * c2 + l_ac * (1 - c2), 0.0)
        sa["theta_stats"] = dict(cos2_mean=float(c2[mask].mean()), beta_off_mean=float(sa["beta_off"][mask].mean()), beta_off_min=float(sa["beta_off"][mask].min()), beta_off_max=float(sa["beta_off"][mask].max()),
                                 alpha_x_mean=float(sa["loss_off"][mask].mean()))
    return physics_run(D, sa)


R4.tables = tables_single
R4.site_arrays = site_arrays_cal
R4.physics = physics_cal


# ------------------------------------------------------------------ tasks
def task(a):
    global MODE
    mode, side, n, pl, taper, excl, launch = a
    MODE = mode
    rows = X.point_task((side, n, pl, taper, excl, launch))
    # extras of the design (same for all variants): switched-cell beta / loss statistics and the radiated / absorbed fractions
    T = R4.tables()
    D = R4.design(side, n, pl, taper, excl, T)
    sa = R4.site_arrays(D)
    ph = R4.physics(D, sa)
    mm = D["mask"] & D["live"]
    ex = dict(mode=mode, placeholder_on=(mode == "cal"), sw_beta_off_mean=float(sa["beta_off"][mm].mean()), sw_beta_on_mean=float(sa["beta_on"][mm].mean()), sw_alpha_leak_off=float(sa["alpha_off"][mm].mean()), sw_alpha_x_off_mean=float(sa["loss_off"][mm].mean()),
              absorbed_frac_off=ph["abs_off"], absorbed_frac_on=ph["abs_on"], N_live=int(mm.sum()))
    for r in rows:
        r.update(ex)
    return rows


def run(workers):
    RES.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    cal = get_cal()
    tasks = [(mode, 0.2, n, pl, t, ex, LAUNCH) for mode in ("base", "cal") for n in (120, 135, 240) for pl in ("a", "b") for t in (0.5, 0.8, 1.0) for ex in (0, 1)]
    allrows = []
    with Pool(min(workers, len(tasks))) as pool:
        for rows in pool.imap_unordered(task, tasks, chunksize=1):
            allrows += rows
    allrows.sort(key=lambda r: (r["mode"], r["N_sw"], r["placement"], r["taper"], r["excl"], r["variant"]))
    B.write_csv(RES / "points_hfsscal.csv", allrows)
    (RES / "calibration.json").write_text(json.dumps(dict(OFF=cal["OFF"], ON=cal["ON"], on_status="MoM placeholder (no HFSS ON data yet)" if cal["ON"] is None else "calibrated"), indent=1, default=float), encoding="utf-8", newline="\n")
    (RES / "provenance_run.json").write_text(json.dumps(B.provenance(dict(stage="hfsscal", workers=min(workers, len(tasks)), n_tasks=len(tasks), seconds=time.time() - t0, launch_fraction=LAUNCH, row=ROW_ID)), indent=1, default=str), encoding="utf-8", newline="\n")
    print(f"done {time.time() - t0:.0f}s, {len(allrows)} rows", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["run", "cal"])
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args()
    if a.stage == "cal":
        print(json.dumps(get_cal(), indent=1, default=float))
    else:
        run(a.workers)
