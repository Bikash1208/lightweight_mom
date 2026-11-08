"""SUB6 FINAL BOARD, steps 1-2: HFSS calibration of the switched island cell (OFF and ON) and the V3 board synthesis, calibrated vs uncalibrated (MODEL_ONLY until board HFSS).
Cell: island cell ISL_W2.9_q0.8_Li0.5_c0.6 with PLAIN vias (island -1.78,-0.53; disc -3.63,-3.63), no chokes, bottom traces on FR4 0.8.  Static cells = the NOCUT library.
Extends sub6_board_hfsscal.py (unchanged):  ON is now calibrated too.  Method per state and direction d (along the slot = x of the MoM, across = y):
  * the MoM dispersion f_M(phase) of the island cell (pred_table, mesh 64; for ON across the UPPER branch only, f >= 5.7 GHz, the lower branch is the island resonance), linearly extrapolated beyond the last point;
  * HFSS points (top surface wave only; low-Q via/island resonances and the isotropic beta/k0 1.013 FR4 underside mode are not used): f_cal(phi) = f_M(phi) + a   (one point) or  + b (phi - phi_ref) (two points);
    ON: ONE point per direction (95 deg): the 80 deg ON points are ambiguous in HFSS, so the MoM ON branch is used there, anchored to the 95 deg points (constant shift a);
  * beta_d/k0 at 5.8 GHz from f_cal;  alpha_x = f (d beta / d f) / (2 Q)  with Q = mean of the recorded Q of that state and direction (OFF 95 deg Q not recorded: the 85 deg Q is used);
  * on the board beta(theta') = beta_along cos^2 theta' + beta_across sin^2 theta' (theta' = ray angle - slot axis, axes before the 90-degree flips), alpha_x likewise, for OFF and ON; MoM leakage alpha kept (MoM-only).
    python sub6_final_cal.py cal | run [--workers N]
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
RES = M.RES.parent / "sub6_final_board"
HFSS = RES / "hfss_data.json"
MOM_TABLE = RES / "input_pred_table_mom_nodata.json"
ROW_ID = "ISL_W2.9_q0.8_Li0.5_c0.6"
P_MM = 11.6
LAUNCH = R4.LAUNCH
MODE = "base"
CAL = None
MESH = os.environ.get("SUB6_MESH", "64")  # "96": the switched-cell row ROW_ID takes its OFF/ON states (alpha, beta, Delta, J) from the mesh-96 table (AEM-LAB1, compute/sub6-island-meshcheck d789d9a)
OUTD = RES / "mesh96" if MESH == "96" else RES
if MESH == "96":
    OUTD.mkdir(parents=True, exist_ok=True)
    _j64 = json.load(open(R4.ISL_JSON))
    _j96 = json.load(open(OUTD / "input_island_rows_n96_d789d9a.json"))
    _new = {r["row"]["id"]: r for r in _j96["rows"]}[ROW_ID]
    for _r in _j64["rows"]:
        if _r["row"]["id"] == ROW_ID:
            _r["states"] = _new["states"]
            _r["row"]["mesh96_real_dims"] = _new["row"]["real"]
    (OUTD / "input_island_rows_merged_m96.json").write_text(json.dumps(_j64), encoding="utf-8")
    R4.ISL_JSON = OUTD / "input_island_rows_merged_m96.json"
_orig = dict(tables=R4.tables, site_arrays=R4.site_arrays, physics=R4.physics)


def mom_curve(state, d):
    """MoM (f, phase) of the branch that contains 5.8 GHz (ON across: the upper branch)."""
    rows = json.load(open(MOM_TABLE))[f"ISL|{state}|{d}"]
    f = np.array([r["f_GHz"] for r in rows])
    ph = np.array([r["phase_deg"] for r in rows])
    o = np.argsort(f)
    f, ph = f[o], ph[o]
    cuts = [0] + [i for i in range(1, len(f)) if ph[i] < ph[i - 1] - 5.0] + [len(f)]
    for i0, i1 in zip(cuts[:-1], cuts[1:]):
        if f[i0] - 1e-9 <= 5.8 <= f[i1 - 1] + 1e-9:
            return f[i0:i1], ph[i0:i1]
    return f, ph


def interp_ext(x, xp, fp):
    """np.interp with linear extrapolation from the two end points."""
    if x < xp[0]:
        return float(fp[0] + (x - xp[0]) * (fp[1] - fp[0]) / (xp[1] - xp[0]))
    if x > xp[-1]:
        return float(fp[-1] + (x - xp[-1]) * (fp[-1] - fp[-2]) / (xp[-1] - xp[-2]))
    return float(np.interp(x, xp, fp))


def calibrate(state):
    data = json.load(open(HFSS))
    out = {}
    for d, key in (("x", "along"), ("y", "across")):
        pts = [p for p in data["points"] if p["state"] == state and p["direction"] == key]
        if not pts:
            return None
        f, ph = mom_curve(state, d)
        fM = lambda phi: interp_ext(phi, ph, f)  # noqa: E731
        phi_of_f = lambda ff: interp_ext(ff, f, ph)  # noqa: E731
        phis = np.array([p["phase_deg"] for p in pts], float)
        df = np.array([p["f_GHz"] - fM(p["phase_deg"]) for p in pts])
        ref = float(phis.mean())
        if len(pts) >= 2 and np.ptp(phis) > 1.0:
            A = np.stack([np.ones_like(phis), phis - ref], 1)
            a, b = np.linalg.lstsq(A, df, rcond=None)[0]
        else:
            a, b = float(df.mean()), 0.0
        fcal = lambda phi: fM(phi) + a + b * (phi - ref)  # noqa: E731
        lo, hi = float(ph.min()) - 25.0, float(ph.max()) + 25.0
        phi58 = brentq(lambda q: fcal(q) - 5.8, lo, hi)
        h = 0.05
        phi_p, phi_m = brentq(lambda q: fcal(q) - (5.8 + h), lo, hi), brentq(lambda q: fcal(q) - (5.8 - h), lo, hi)
        dphi_df = (phi_p - phi_m) / (2 * h)
        k0 = 2 * math.pi * 5.8e9 / 299792458.0
        beta58 = math.radians(phi58) / (P_MM * 1e-3) / k0
        dbeta_df = math.radians(dphi_df) / (P_MM * 1e-3)
        qs = [p["Q"] for p in pts if p.get("Q")]
        if not qs:
            qs = [p["Q"] for p in data["points"] if p["state"] == state and p["direction"] == key and p.get("Q")]
        Q = float(np.mean(qs))
        mom58 = math.radians(float(phi_of_f(5.8))) / (P_MM * 1e-3) / k0
        light_phase_58 = math.degrees(k0 * P_MM * 1e-3)
        out[d] = dict(direction=key, n_points=len(pts), shift_a_GHz=float(a), slope_b_GHz_per_deg=float(b), phi_ref_deg=ref, beta_over_k0_5p8=float(beta58), beta_over_k0_MoM_5p8=float(mom58), phase58_deg=float(phi58), light_line_phase_5p8_deg=float(light_phase_58),
                      dphi_df_deg_per_GHz=float(dphi_df), Q=Q, alpha_x_Np_per_m=float(5.8 * dbeta_df / (2 * Q)), points=pts, f_MoM_at_phase_GHz=[float(fM(p["phase_deg"])) for p in pts],
                      extrapolated_mom=[bool(p["phase_deg"] > ph.max() or p["phase_deg"] < ph.min()) for p in pts], mom_branch_f_GHz=[float(f.min()), float(f.max())])
    return out


def get_cal():
    global CAL
    if CAL is None:
        CAL = {"OFF": calibrate("OFF_BAR50"), "ON": calibrate("ON_BAR50")}
    return CAL


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
    return w, phi, float(np.mean(np.sum(2.0 * a * w, axis=1) * drho)), float(np.mean(np.sum(2.0 * ls * w, axis=1) * drho))


def physics_run(D, sa):
    S, lk, m0 = D["S"], D["lk"], D["m0"]
    smp = M.sampler(m0, S)
    w_off, ph_off, rad_off, abs_off = real_march_loss(m0, lk, sa["alpha_off"], sa["loss_off"], sa["beta_off"], D["b_fill"])
    w_on, ph_on, rad_on, abs_on = real_march_loss(m0, lk, sa["alpha_on"], sa["loss_on"], sa["beta_on"], D["b_fill"])
    return dict(Woff=smp(w_off), Won=smp(w_on), Poff=smp(ph_off), Pon=smp(ph_on), eta_off=rad_off, eta_on=rad_on, abs_off=abs_off, abs_on=abs_on)


def physics_cal(D, sa):
    if MODE == "cal":
        cal = get_cal()
        mask = D["mask"] & D["live"]
        S = D["S"]
        phiR = np.arctan2(S.y, S.x)
        c0 = {s: (cal[s]["x"], cal[s]["y"]) for s in ("OFF", "ON")}
        for s in ("OFF", "ON"):
            sa["beta_" + s.lower()] = np.where(mask, 0.5 * (max(1.0, c0[s][0]["beta_over_k0_5p8"]) + max(1.0, c0[s][1]["beta_over_k0_5p8"])), sa["beta_" + s.lower()])
        for it in range(3):
            ph = physics_run(D, sa)
            direction = X.directions_w(D, ph, 0.0)
            chi = np.arctan2(direction[:, 1], direction[:, 0])
            axis = np.mod(0.5 * (chi + phiR), math.pi)
            c2 = np.cos(phiR - axis) ** 2
            for s in ("OFF", "ON"):
                al, ac = c0[s]
                # beta < 1 cannot be a guided wave: the calibrated ON-along value 0.990 (rigid shift crosses the light line) is floored at 1.000 (flagged)
                sa["beta_" + s.lower()] = np.where(mask, max(1.0, al["beta_over_k0_5p8"]) * c2 + max(1.0, ac["beta_over_k0_5p8"]) * (1 - c2), sa["beta_" + s.lower()])
                sa["loss_" + s.lower()] = np.where(mask, al["alpha_x_Np_per_m"] * c2 + ac["alpha_x_Np_per_m"] * (1 - c2), 0.0)
    return physics_run(D, sa)


R4.tables = tables_single
R4.site_arrays = site_arrays_cal
R4.physics = physics_cal


def task(a):
    global MODE
    mode, side, n, pl, taper, excl, launch = a
    MODE = mode
    rows = X.point_task((side, n, pl, taper, excl, launch))
    T = R4.tables()
    D = R4.design(side, n, pl, taper, excl, T)
    sa = R4.site_arrays(D)
    ph = R4.physics(D, sa)
    mm = D["mask"] & D["live"]
    ex = dict(mode=mode, N_live=int(mm.sum()), sw_beta_off_mean=float(sa["beta_off"][mm].mean()), sw_beta_on_mean=float(sa["beta_on"][mm].mean()), sw_beta_off_min=float(sa["beta_off"][mm].min()), sw_beta_on_min=float(sa["beta_on"][mm].min()),
              sw_alpha_leak_off=float(sa["alpha_off"][mm].mean()), sw_alpha_leak_on=float(sa["alpha_on"][mm].mean()), sw_alpha_x_off_mean=float(sa["loss_off"][mm].mean()), sw_alpha_x_on_mean=float(sa["loss_on"][mm].mean()),
              absorbed_frac_off=ph["abs_off"], absorbed_frac_on=ph["abs_on"])
    for r in rows:
        r.update(ex)
    return rows


def run(workers):
    OUTD.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    cal = get_cal()
    tasks = [(mode, 0.2, n, pl, t, ex, LAUNCH) for mode in ("base", "cal") for n in (240, 135) for pl in ("a", "b") for t in (0.5, 0.8, 1.0) for ex in (0, 1)]
    allrows = []
    with Pool(min(workers, len(tasks))) as pool:
        for rows in pool.imap_unordered(task, tasks, chunksize=1):
            allrows += rows
    allrows.sort(key=lambda r: (r["mode"], r["N_sw"], r["placement"], r["taper"], r["excl"], r["variant"]))
    B.write_csv(OUTD / "points_final.csv", allrows)
    (OUTD / "calibration.json").write_text(json.dumps(cal, indent=1, default=float), encoding="utf-8", newline="\n")
    (OUTD / "provenance_final.json").write_text(json.dumps(B.provenance(dict(stage="final_cal_synthesis", workers=min(workers, len(tasks)), n_tasks=len(tasks), seconds=time.time() - t0, launch_fraction=LAUNCH, row=ROW_ID)), indent=1, default=str), encoding="utf-8", newline="\n")
    print(f"done {time.time() - t0:.0f}s, {len(allrows)} rows", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["run", "cal"])
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args()
    if a.stage == "cal":
        c = get_cal()
        for s in ("OFF", "ON"):
            for d, v in c[s].items():
                print(s, v["direction"], {k: (round(x, 4) if isinstance(x, float) else x) for k, x in v.items() if k not in ("points",)})
    else:
        run(a.workers)
