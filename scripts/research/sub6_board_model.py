"""SUB6 board model with REAL cell tables (MODEL_ONLY, PENDING_HFSS_CELL_CHECK).  5.8 GHz, RT5880 3.175, p 11.6, centre probe, beams A (25,0) / B (25,180) RHCP.

No frozen code is edited.  `sub6_board_req.py` (copied unmodified from compute/sub6-board-requirements 7694971) supplies the lattice, the far-field matrix `Side`, the
metrics and the gate; its module constant PITCH is set to 11.6 mm here (all functions read the module global at call time).  Phase 26 `march_rays` / `cell_design`
are called unmodified.

Architecture (user, 2026-11-08): OFF = static cells radiate A+B; ON = static cells unchanged + switched cells radiate -B (near-silent in OFF).  Group-normalised:
    static site   J = a e^{-j phi} s_hat,  s_hat = (uA + uB)/|uA + uB|, target amplitude amp0 sc |uA + uB|/(1-f)      (one linear current: amplitude via the cell alpha, axis psi)
    switched site J_ON = -a e^{-j phi} uB, target amplitude amp0 sc/f ;  J_OFF = same direction, amplitude from its (tiny) OFF alpha
so the beam-B coherent sum is 1 in OFF and 0 in ON (this is the rho -> -infinity limit of the earlier ideal-cell architecture).  f = N_sw / N.
Real cells: every site is assigned one ROW of results/sub6_slotdisc/rows_RT5880.jsonl:
    static   : NOCUT rows (alpha_NOCUT, beta_NOCUT), nearest in sqrt(alpha) to the required alpha (alpha_req = alpha0 sc^2 |uA+uB|^2/(1-f)^2, alpha0 = base-march alpha)
    switched : notch rows (alpha_OFF, beta_OFF, alpha_ON, beta_ON), nearest in sqrt(alpha_ON) to alpha_req = alpha0 sc^2/f^2
Filters (mine, flagged): no resonance_flag row, static beta_NOCUT <= 1.2, switched alpha_OFF <= 0.1 Np/m, alpha_ON >= 0.3, beta_ON <= 1.2.
Per-state propagation ("real march", my bookkeeping outside the engine): along every ray (720 rays, 0.25 mm) the cell under the ray supplies its OWN alpha and beta for that
state, W = exp(-2 int alpha), phi = k0 int beta; amplitude a = sqrt(2 alpha W / rho).  The slot axes (chi) are designed with the OFF phase; the ON state keeps those axes but
the wave phase is the ON phase (switched cells beta 1.05 -> 1.11-1.14), which is the beam-A phase error reported.
Not modelled: coupling between neighbouring cells of different type (rows differ in gap and W), the A+B current being an ideal linear current of the chosen cell,
the switched cell's OFF polarisation (assumed the ON direction), diode loss beyond the BAR50 model.  Placement (b) is chosen with the ideal-cell objective of the earlier job
(rho = -1000) and then evaluated with the real cells.
    python sub6_board_model.py run [--workers N]  |  report
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
os.environ["MKL_NUM_THREADS"] = "1"
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sub6_board_req as B

B.PITCH = 11.6e-3
PITCH = B.PITCH
RES = ROOT / "results" / "sub6_board_model"
ROWS = ROOT / "results" / "sub6_slotdisc" / "rows_RT5880.jsonl"
from scipy.interpolate import RegularGridInterpolator

from src.hologram.tensor_leaky_cp import SupercellTable, cell_design, march_rays

K0, F = B.K0, B.F
N_SW = (150, 200, 250, 300)  # at side 200 mm; side 150 mm: scaled by the cell-count ratio
SIDES = (0.2, 0.15)
RHO_LIM = -1.0e3


def tables():
    rows = [json.loads(line) for line in ROWS.read_text().splitlines() if line.strip()]
    g = lambda r, s, k: r["states"][s][k]  # noqa: E731
    st = [r for r in rows if not r["resonance_flag"] and g(r, "NOCUT", "beta_over_k0") <= 1.2]
    sw = [r for r in rows if not r["resonance_flag"] and g(r, "OFF_BAR50", "alpha_np_per_m") <= 0.1 and g(r, "ON_BAR50", "alpha_np_per_m") >= 0.3 and g(r, "ON_BAR50", "beta_over_k0") <= 1.2]
    T = dict(
        st=dict(rows=st, id=[r["row"]["id"] for r in st], a=np.array([g(r, "NOCUT", "alpha_np_per_m") for r in st]), b=np.array([g(r, "NOCUT", "beta_over_k0") for r in st])),
        sw=dict(rows=sw, id=[r["row"]["id"] for r in sw], aoff=np.array([g(r, "OFF_BAR50", "alpha_np_per_m") for r in sw]), boff=np.array([g(r, "OFF_BAR50", "beta_over_k0") for r in sw]),
                aon=np.array([g(r, "ON_BAR50", "alpha_np_per_m") for r in sw]), bon=np.array([g(r, "ON_BAR50", "beta_over_k0") for r in sw])))
    T["n_rows_total"], T["n_flagged"] = len(rows), sum(bool(r["resonance_flag"]) for r in rows)
    return T


def static_table(T):
    st = T["st"]
    n = len(st["a"])
    return SupercellTable(tuple((r["row"]["L"], r["row"]["W"]) for r in st["rows"]), np.array([0.0, 3.0]), np.repeat(st["a"][:, None], 2, 1), np.repeat(st["b"][:, None], 2, 1), 6.0, PITCH * 1e3)


def nearest_sqrt(avail, req):
    return np.argmin(np.abs(np.sqrt(avail)[None, :] - np.sqrt(np.maximum(req, 0.0))[:, None]), axis=1)


class Lookup:
    def __init__(self, S):
        self.n = int(np.abs(S.ij).max()) + 2
        self.t = -np.ones((2 * self.n + 1, 2 * self.n + 1), int)
        self.t[S.ij[:, 0] + self.n, S.ij[:, 1] + self.n] = np.arange(S.n)

    def at(self, x, y):
        i = np.clip(np.rint(x / PITCH).astype(int), -self.n, self.n)
        j = np.clip(np.rint(y / PITCH).astype(int), -self.n, self.n)
        return self.t[i + self.n, j + self.n]


def real_march(m, lk, a_site, b_site, b_fill):
    """Per-state march along the engine's rays with the alpha / beta of the cell under each ray point.  Returns W, phi_ref (n_psi, n_rho) and the edge conversion."""
    psi, rho, drho = m["psi"], m["rho"], m["drho"]
    X, Y = rho[None, :] * np.cos(psi)[:, None], rho[None, :] * np.sin(psi)[:, None]
    idx = lk.at(X, Y)
    inside = rho[None, :] <= m["r_edge"][:, None]
    ok = (idx >= 0) & inside
    a = np.where(ok, a_site[np.maximum(idx, 0)], 0.0)
    b = np.where(idx >= 0, b_site[np.maximum(idx, 0)], b_fill)
    ca, cb = np.cumsum(a, axis=1) - a, np.cumsum(b, axis=1) - b
    w = np.exp(-2.0 * ca * drho)
    phi = cb * K0 * drho
    re = np.minimum(np.floor(m["r_edge"] / drho).astype(int), rho.size - 1)
    w_edge = w[np.arange(psi.size), re]
    return w, phi, float(np.mean(1.0 - w_edge))


def sampler(m, S):
    psi, rg = m["psi"], m["rho"]
    pp = np.concatenate([[psi[-1] - 2 * np.pi], psi, [psi[0] + 2 * np.pi]])
    phi_s, rho_s = np.arctan2(S.y, S.x), np.minimum(np.hypot(S.x, S.y), rg[-1])

    def f(field):
        return RegularGridInterpolator((pp, rg), np.vstack([field[-1:], field, field[:1]]), bounds_error=False, fill_value=None)((phi_s, rho_s))

    return f


def design_and_evaluate(side, n_sw, place, launch, T):
    S = B.get_side(side)
    lk = Lookup(S)
    tab = static_table(T)
    f = n_sw / S.n
    b_fill = float(np.mean(T["st"]["b"]))
    kt = K0 * math.sin(math.radians(B.BEAM_A[0]))
    p_static = 4.0 * np.cos(kt * S.x) ** 2 / (1 - f) ** 2  # |uA + uB|^2 / (1-f)^2, |uA + uB|^2 = 4 cos^2(kt x)
    if place == "a":
        mask0 = B.sites_uniform(S, n_sw)
    else:
        mi, dAi, dBi = B.march_for(side, 2.0)
        mask0 = B.optimise_sites(S, dAi, dBi, f, RHO_LIM, n_sw)

    def base(table):
        m0 = march_rays(table, B.hexagon(side), (0.0, 0.0), K0, math.radians(B.BEAM_A[0]), math.radians(B.BEAM_A[1]), B.ETA_DESIGN, taper_c=B.TAPER, rho_min_m=B.RHO_CLEAR, beta0_over_k0=b_fill)
        d0 = cell_design(m0, table, S.x, S.y, (0.0, 0.0), K0, math.radians(B.BEAM_A[0]), math.radians(B.BEAM_A[1]), hand=+1)
        live = d0["amp"] > 1e-12
        mask = mask0 & live
        kap = float(np.sum(np.where(mask, 0.0, d0["amp"] ** 2 * p_static)) / max(np.sum(d0["amp"] ** 2), 1e-30))
        return m0, d0, live, mask, kap

    # pass 0: full static table -> demand quantile q90 of the power factors (static sc^2 p, switched sc^2/f^2) -> pass 1: base march with the rows alpha <= cap1,
    # exactly the cap alpha_max / q90 of the earlier ideal-cell build, so the sites that need the most amplitude fit inside the available alpha
    m0, dA0, live, mask, kappa = base(tab)
    dem = np.where(mask, 1.0 / f**2, p_static) / kappa
    q90 = float(np.quantile(dem[live], 0.9))
    cap1 = min(float(T["st"]["a"].max()), float(T["sw"]["aon"].max())) / max(q90, 1.0)
    keep = np.nonzero(T["st"]["a"] <= cap1)[0]
    if keep.size < 3:
        keep = np.argsort(T["st"]["a"])[:3]
    tab1 = SupercellTable(tuple(tab.slots[i] for i in keep), tab.k_grid, tab.alpha[keep], tab.beta[keep], tab.diameter_mm, tab.pitch_mm)
    m0, dA0, live, mask, kappa = base(tab1)
    alpha0, amp0 = dA0["alpha"], dA0["amp"]
    sc2 = 1.0 / kappa
    a_req_st = alpha0 * sc2 * p_static
    a_req_sw = alpha0 * sc2 / f**2
    ist = nearest_sqrt(T["st"]["a"], a_req_st)
    isw = nearest_sqrt(T["sw"]["aon"], a_req_sw)
    sw, st = T["sw"], T["st"]
    alpha_off = np.where(mask, sw["aoff"][isw], st["a"][ist])
    alpha_on = np.where(mask, sw["aon"][isw], st["a"][ist])
    beta_off = np.where(mask, sw["boff"][isw], st["b"][ist])
    beta_on = np.where(mask, sw["bon"][isw], st["b"][ist])
    alpha_off, alpha_on = np.where(live, alpha_off, 0.0), np.where(live, alpha_on, 0.0)
    # real marches, per state
    smp = sampler(m0, S)
    w_off, ph_off, eta_off = real_march(m0, lk, alpha_off, beta_off, b_fill)
    w_on, ph_on, eta_on = real_march(m0, lk, alpha_on, beta_on, b_fill)
    Woff, Won, Poff, Pon = smp(w_off), smp(w_on), smp(ph_off), smp(ph_on)
    rho = np.hypot(S.x, S.y)
    dx, dy = S.x, S.y
    chiA = Poff - kt * dx
    chiB = Poff + kt * dx  # beam B: phi0 = 180 deg
    uA = np.stack([np.cos(chiA), np.sin(chiA)], 1)
    uB = np.stack([np.cos(chiB), np.sin(chiB)], 1)
    ssum = uA + uB
    nrm = np.linalg.norm(ssum, axis=1, keepdims=True)
    shat = np.where(nrm > 1e-9, ssum / np.maximum(nrm, 1e-12), uA)
    d_st, d_sw = shat, -uB
    direction = np.where(mask[:, None], d_sw, d_st)

    def currents(alpha, W, ph):
        a = np.sqrt(np.maximum(2.0 * alpha * np.maximum(W, 0.0) / np.maximum(rho, 1e-9), 0.0))
        j = a * np.exp(-1j * ph)
        return j * direction[:, 0], j * direction[:, 1], a

    jx0, jy0, a0 = currents(alpha_off, Woff, Poff)
    jx1, jy1, a1 = currents(alpha_on, Won, Pon)
    jx2, jy2, _ = currents(alpha_on, Won, Poff)  # counterfactual: ON amplitudes, OFF phase (no beta change of the switched cells)
    r_off = B.metrics(S, *S.far(jx0, jy0), eta_off, launch)
    r_on = B.metrics(S, *S.far(jx1, jy1), eta_on, launch)
    r_cf = B.metrics(S, *S.far(jx2, jy2), eta_on, launch)
    stat = ~mask & live
    dphi = np.degrees(Pon - Poff)
    wgt = a0**2 * stat
    row = dict(side_mm=side * 1e3, N_cells=S.n, N_sw=n_sw, f_sw=f, placement=place,
               off_GA=r_off["gA"], off_GB=r_off["gB"], off_balance=abs(r_off["gA"] - r_off["gB"]), off_AR_A=r_off["arA"], off_AR_B=r_off["arB"],
               on_GA=r_on["gA"], on_GB=r_on["gB"], on_suppression=r_on["gA"] - r_on["gB"], on_AR_A=r_on["arA"], on_AR_B=r_on["arB"],
               sll_off=r_off["sll"], sll_on=r_on["sll"], eta_off=eta_off, eta_on=eta_on,
               dphi_rms_deg=float(np.sqrt(np.sum(wgt * dphi**2) / max(np.sum(wgt), 1e-30))), dphi_max_deg=float(np.max(np.abs(dphi[stat]))), dphi_mean_deg=float(np.sum(wgt * dphi) / max(np.sum(wgt), 1e-30)),
               on_GA_no_beta_change=r_cf["gA"], A_loss_from_beta_db=r_cf["gA"] - r_on["gA"], on_supp_no_beta_change=r_cf["gA"] - r_cf["gB"], on_AR_A_no_beta_change=r_cf["arA"], sll_on_no_beta_change=r_cf["sll"],
               sw_alpha_req_over_avail_max=float(np.max(a_req_sw[mask]) / sw["aon"].max()), sw_sites_alpha_saturated=int(np.sum(a_req_sw[mask] > sw["aon"].max())),
               st_sites_alpha_saturated=int(np.sum(a_req_st[stat] > st["a"].max())),
               sw_amp_err_rms_pct=float(100 * np.sqrt(np.mean((np.sqrt(sw["aon"][isw][mask]) / np.sqrt(np.maximum(a_req_sw[mask], 1e-12)) - 1) ** 2))),
               q90_demand=q90, cap1_Np_m=cap1, base_rows=int(keep.size), sw_alpha_on_mean=float(np.mean(sw["aon"][isw][mask])), sw_beta_on_minus_off_mean=float(np.mean((sw["bon"] - sw["boff"])[isw][mask])), alpha0_mean=float(np.mean(alpha0[live])))
    ar = max(row["off_AR_A"], row["off_AR_B"], row["on_AR_A"])
    row["worst_margin_db"], row["gate"] = B.gate(row["off_balance"], ar, row["on_suppression"], r_off["sll"], r_on["sll"])
    axis = np.mod(0.5 * (np.arctan2(direction[:, 1], direction[:, 0]) + np.arctan2(S.y, S.x)), np.pi)
    layout = [dict(site=i, x_mm=S.x[i] * 1e3, y_mm=S.y[i] * 1e3, cls="SWITCHED" if mask[i] else "STATIC",
                   row_id=(sw["id"][isw[i]] if mask[i] else st["id"][ist[i]]), axis_deg=math.degrees(axis[i]),
                   alpha_off=alpha_off[i], alpha_on=alpha_on[i], beta_off=beta_off[i], beta_on=beta_on[i], alpha_required=(a_req_sw[i] if mask[i] else a_req_st[i]),
                   dphi_on_minus_off_deg=dphi[i]) for i in range(S.n) if live[i]]
    return row, layout


def task(a):
    side, n_sw, place, launch = a
    T = tables()
    row, layout = design_and_evaluate(side, n_sw, place, launch, T)
    lp = RES / f"layout_{int(side * 1e3)}mm_N{n_sw}_{place}.csv"
    with open(lp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(layout[0].keys()))
        w.writeheader()
        for r in layout:
            w.writerow({k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in r.items()})
    return row


def run(workers):
    RES.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    lf, lfd = B.launch_fraction()
    n200 = B.get_side(0.2).n
    tasks = []
    for side in SIDES:
        scale = B.get_side(side).n / n200
        for n in N_SW:
            for pl in ("a", "b"):
                tasks.append((side, int(round(n * scale)) if side != 0.2 else n, pl, lf))
    T = tables()
    print(f"{len(tasks)} tasks, {workers} workers, launch fraction {lf:.4f}; static rows {len(T['st']['a'])}, switched rows {len(T['sw']['aon'])} of {T['n_rows_total']} ({T['n_flagged']} resonance-flagged)", flush=True)
    with Pool(min(workers, len(tasks))) as pool:
        rows = []
        for r in pool.imap_unordered(task, tasks):
            rows.append(r)
            print(f"  {r['side_mm']:.0f} N{r['N_sw']} {r['placement']} {r['gate']} margin {r['worst_margin_db']:.2f}", flush=True)
    rows.sort(key=lambda r: (-r["side_mm"], r["placement"], r["N_sw"]))
    B.write_csv(RES / "points.csv", rows)
    prov = B.provenance(dict(stage="run", workers=min(workers, len(tasks)), n_tasks=len(tasks), seconds=time.time() - t0, launch_fraction=lf, launch_detail=lfd,
                             static_rows=len(T["st"]["a"]), switched_rows=len(T["sw"]["aon"]), rows_total=T["n_rows_total"], rows_resonance_flagged=T["n_flagged"], pitch_m=PITCH))
    (RES / "provenance_run.json").write_text(json.dumps(prov, indent=1, default=str), encoding="utf-8", newline="\n")
    print(f"done {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["run"])
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args()
    run(a.workers)
