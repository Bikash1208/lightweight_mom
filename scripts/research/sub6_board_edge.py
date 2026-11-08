"""SUB6 board model with EDGE-DIODE switched cells in straight series strings (MODEL_ONLY, PENDING_HFSS_CELL_CHECK).

Same engine / V3 method as sub6_board_island2.py (unchanged, imported): static cells (NOCUT rows) carry beam A in both states, switched cells carry beam B in OFF and are near-quiet in
ON, general weights with rho_d = mean signed ratio J_ON / J_OFF of the switched table, joint 90-degree rotations (V0 single-direction; V1 w = 0.5; V2 w*; V3 multi-direction, best of V2 and the
multi-direction runs by gate margin).  Changes:
* SWITCHED table = rows of results/sub6_edge_diode (input copy `input_edge_diode.json`, branch compute/sub6-edge-diode-cell 34c4a37) with the slot at 90 deg to the diode axis, real OFF / ON alpha, beta
  and Delta; my filter: alpha_OFF >= 0.9 Np/m and sqrt(alpha_ON / alpha_OFF) <= 0.5 -> 5 rows (q 0.7 W 2.4 / 2.9 / 3.6, q 0.8 W 2.4 / 2.9).  The edge-diode scan has NO resonance proxy (resonance_flip = null).
  NOTE: the table rows are periodic-lattice values with a diode to the +x neighbour in EVERY cell and a 1.00 mm gap (disc r 5.30); the board cells next to a string end or beside a string differ.
* PLACEMENT: the switched cells form straight lines of consecutive neighbours along x or y (each line = one series string, diode axis along the line, n_cells - 1 diodes in series, length 10-24 cells,
  no two lines of the same axis touching end to end, no overlaps).  Placement 'a' = random compatible lines (seeded); placement 'b' = greedy line selection + remove-and-replace local search with the
  ideal-cell projection objective of sub6_board_req.optimise_sites (exact additive per-site deltas, prefix sums over the line candidates), total cell count exactly N_sw.
  The slot axis of every cell follows the hologram; the table row (slot at 90 deg to the diode axis) is used for every switched cell whatever the actual angle between the slot and the line
  (flagged; the angle distribution and the amplitude factor of a theta-aware alpha are reported).
* String voltage = (n - 1) x 0.95 V at 10 mA (diode forward drop given by the user).
    python sub6_board_edge.py run [--workers N]
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

import sub6_board_island2 as X
import sub6_board_model4 as R4
from src.hologram.tensor_leaky_cp import SupercellTable, cell_design, march_rays

B, M = R4.B, R4.M
K0 = R4.K0
RES = M.RES.parent / "sub6_board_edge"
EDGE_JSON = RES / "input_edge_diode.json"
LMIN, LMAX = 10, 24
V_DIODE = 0.95
SIDES = (0.2, 0.22)
N_LIST = (135, 240)
LAUNCH = R4.LAUNCH


# ------------------------------------------------------------------ tables
def edge_tables():
    T = M.tables()
    J = json.load(open(EDGE_JSON))["results"]
    z = lambda v: complex(v[0], v[1]) * 1e-3  # noqa: E731
    sel = [r for r in J if r["theta"] == 90 and r["alpha_OFF"] >= 0.9 and math.sqrt(r["alpha_ON"] / r["alpha_OFF"]) <= 0.5]
    n = len(sel)
    sw = dict(rows=sel, id=[r["id"] for r in sel], aoff=np.array([r["alpha_OFF"] for r in sel]), aon=np.array([r["alpha_ON"] for r in sel]),
              boff=np.array([r["beta_OFF"] for r in sel]), bon=np.array([r["beta_ON"] for r in sel]),
              Doff=np.array([z(r["Delta_mS"]["OFF_BAR50"]) for r in sel]), Don=np.array([z(r["Delta_mS"]["ON_BAR50"]) for r in sel]),
              Apoff=np.zeros(n, complex), Apon=np.zeros(n, complex), Aqoff=np.zeros(n, complex), Aqon=np.zeros(n, complex))
    sgn = np.where(np.array([z(r["Delta_mS"]["ON_BAR50"]).imag / z(r["Delta_mS"]["OFF_BAR50"]).imag for r in sel]) > 0, 1.0, -1.0)
    sw["ratio"] = sgn * np.sqrt(sw["aon"] / sw["aoff"])
    T["sw"] = sw
    T["n_edge_rows_total"], T["n_edge_theta90"] = len(J), sum(r["theta"] == 90 for r in J)
    return T


# ------------------------------------------------------------------ line placement
class Lines:
    def __init__(self, S, lk, allowed):
        self.S, self.lk, self.n = S, lk, S.n
        self.runs, self.axis = [], []
        n = lk.n
        for ax in ("x", "y"):
            for fixed in range(-n, n + 1):
                cur = []
                for var in range(-n, n + 1):
                    i, j = (var, fixed) if ax == "x" else (fixed, var)
                    idx = lk.t[i + n, j + n]
                    if idx >= 0 and allowed[idx]:
                        cur.append(idx)
                    else:
                        if len(cur) >= LMIN:
                            self.runs.append(np.array(cur)), self.axis.append(ax)
                        cur = []
                if len(cur) >= LMIN:
                    self.runs.append(np.array(cur)), self.axis.append(ax)
        rid, st, ln = [], [], []
        for k, r in enumerate(self.runs):
            for s in range(0, len(r) - LMIN + 1):
                for L in range(LMIN, min(LMAX, len(r) - s) + 1):
                    rid.append(k), st.append(s), ln.append(L)
        self.crun, self.cst, self.cln = np.array(rid), np.array(st), np.array(ln)
        self.rstart = np.zeros(len(self.runs) + 1, int)
        np.add.at(self.rstart, np.array(rid) + 1, 1)
        self.rstart = np.cumsum(self.rstart)

    def deltas(self, dvec):
        """dvec (n_sites, 3) complex per-site deltas -> (n_cand, 3) line deltas via prefix sums."""
        out = np.zeros((len(self.crun), 3), complex)
        for k, r in enumerate(self.runs):
            P = np.vstack([np.zeros((1, 3), complex), np.cumsum(dvec[r], axis=0)])
            a, b = self.rstart[k], self.rstart[k + 1]
            out[a:b] = P[self.cst[a:b] + self.cln[a:b]] - P[self.cst[a:b]]
        return out

    def compat(self, used_any, used_ax):
        ok = np.zeros(len(self.crun), bool)
        for k, r in enumerate(self.runs):
            a, b = self.rstart[k], self.rstart[k + 1]
            ua = np.concatenate([[0], np.cumsum(used_any[r])])
            ux = used_ax[self.axis[k]][r]
            s, L = self.cst[a:b], self.cln[a:b]
            free = (ua[s + L] - ua[s]) == 0
            left = np.where(s > 0, ux[np.maximum(s - 1, 0)], False)
            right = np.where(s + L < len(r), ux[np.minimum(s + L, len(r) - 1)], False)
            ok[a:b] = free & ~left & ~right
        return ok

    def cells(self, c):
        return self.runs[self.crun[c]][self.cst[c]: self.cst[c] + self.cln[c]]


def obj_vec(sums):
    FA, FB, R = sums[..., 0], sums[..., 1], sums[..., 2]
    return np.log(np.abs(FA)) + np.log(np.abs(FB)) - np.log(np.abs(R) ** 2 + 1e-3 * np.abs(FB) ** 2)


def place_lines(S, lk, n_sw, place, allowed, rho_d, side, seed):
    LN = Lines(S, lk, allowed)
    used_any = np.zeros(S.n, bool)
    used_ax = {"x": np.zeros(S.n, bool), "y": np.zeros(S.n, bool)}
    rng = np.random.default_rng(seed)
    if place == "b":
        mi, dAi, dBi = B.march_for(side, 2.0)
        c = -rho_d / (1.0 - rho_d)
        f = n_sw / S.n
        ax_, ay_ = B.unit_fields(dAi)
        bx_, by_ = B.unit_fields(dBi)
        amp = dAi["amp"]
        aA = S.cell_field(B.BEAM_A, amp * ax_ / (1 - f), amp * ay_ / (1 - f))
        aBst = S.cell_field(B.BEAM_B, amp * c * bx_ / (1 - f), amp * c * by_ / (1 - f))
        aBsw = S.cell_field(B.BEAM_B, amp * bx_ / ((1 - rho_d) * f), amp * by_ / ((1 - rho_d) * f))
        base = np.array([aA.sum(), aBst.sum(), aBst.sum()])
        dvec = np.stack([-aA, aBsw - aBst, rho_d * aBsw - aBst], 1)
        cand_d = LN.deltas(dvec)
    chosen, total = [], 0
    sums = base.copy() if place == "b" else None

    def feasible_len(rem):
        return (LN.cln <= rem) & ((rem - LN.cln == 0) | (rem - LN.cln >= LMIN))

    def mark(c, val):
        cells = LN.cells(c)
        used_any[cells] = val
        used_ax[LN.axis[LN.crun[c]]][cells] = val

    while total < n_sw:
        rem = n_sw - total
        ok = LN.compat(used_any, used_ax) & feasible_len(rem)
        if not ok.any():
            break
        if place == "b":
            sc = np.where(ok, obj_vec(sums[None, :] + cand_d), -np.inf)
            c = int(np.argmax(sc))
        else:
            idx = np.nonzero(ok)[0]
            c = int(rng.choice(idx))
        chosen.append(c), mark(c, True)
        total += int(LN.cln[c])
        if place == "b":
            sums = sums + cand_d[c]
    if place == "b" and chosen:
        for _ in range(6):
            improved = False
            for pos in range(len(chosen)):
                c0 = chosen[pos]
                mark(c0, False)
                s_wo = sums - cand_d[c0]
                ok = LN.compat(used_any, used_ax) & (LN.cln == LN.cln[c0])
                sc = np.where(ok, obj_vec(s_wo[None, :] + cand_d), -np.inf)
                c1 = int(np.argmax(sc))
                if sc[c1] > obj_vec(s_wo + cand_d[c0]) + 1e-12:
                    chosen[pos] = c1
                    mark(c1, True)
                    sums = s_wo + cand_d[c1]
                    improved = True
                else:
                    mark(c0, True)
            if not improved:
                break
    mask = np.zeros(S.n, bool)
    lines = []
    for c in chosen:
        cells = LN.cells(c)
        mask[cells] = True
        lines.append(dict(axis=LN.axis[LN.crun[c]], cells=cells, L=int(LN.cln[c])))
    return mask, lines


# ------------------------------------------------------------------ design (copy of sub6_board_model4.design with the line placement)
def design_edge(side, n_sw, place, taper, excl, T):
    S = B.get_side(side)
    lk = M.Lookup(S)
    sw, st = T["sw"], T["st"]
    tab = M.static_table(T)
    f = n_sw / S.n
    b_fill = float(np.mean(st["b"]))
    kt = K0 * math.sin(math.radians(B.BEAM_A[0]))
    rho_d = float(np.mean(sw["ratio"]))
    c_w = -rho_d / (1.0 - rho_d)
    p_static = (1.0 + c_w**2 + 2.0 * c_w * np.cos(2.0 * kt * S.x)) / (1 - f) ** 2
    p_sw = 1.0 / ((1.0 - rho_d) ** 2 * f**2)
    ring, outer = R4.region_masks(S, lk)
    allowed = (ring > 2) & ~outer if excl else np.ones(S.n, bool)
    mask0, lines = place_lines(S, lk, n_sw, place, allowed, rho_d, side, seed=int(side * 1000) * 1000 + n_sw)

    def base(table):
        m0 = march_rays(table, B.hexagon(side), (0.0, 0.0), K0, math.radians(B.BEAM_A[0]), math.radians(B.BEAM_A[1]), B.ETA_DESIGN, taper_c=taper, rho_min_m=B.RHO_CLEAR, beta0_over_k0=b_fill)
        d0 = cell_design(m0, table, S.x, S.y, (0.0, 0.0), K0, math.radians(B.BEAM_A[0]), math.radians(B.BEAM_A[1]), hand=+1)
        live = d0["amp"] > 1e-12
        mask = mask0 & live
        kap = float(np.sum(np.where(mask, 0.0, d0["amp"] ** 2 * p_static)) / max(np.sum(d0["amp"] ** 2), 1e-30))
        return m0, d0, live, mask, kap

    m0, dA0, live, mask, kappa = base(tab)
    dem = np.where(mask, p_sw, p_static) / kappa
    q90 = float(np.quantile(dem[live], 0.9))
    cap1 = min(float(st["a"].max()), float(sw["aoff"].max())) / max(q90, 1.0)
    keep = np.nonzero(st["a"] <= cap1)[0]
    if keep.size < 3:
        keep = np.argsort(st["a"])[:3]
    tab1 = SupercellTable(tuple(tab.slots[i] for i in keep), tab.k_grid, tab.alpha[keep], tab.beta[keep], tab.diameter_mm, tab.pitch_mm)
    m0, dA0, live, mask, kappa = base(tab1)
    alpha0 = dA0["alpha"]
    sc2 = 1.0 / kappa
    a_req_st = alpha0 * sc2 * p_static
    a_req_sw = alpha0 * sc2 * p_sw
    ist = M.nearest_sqrt(st["a"], a_req_st)
    isw = M.nearest_sqrt(sw["aoff"], a_req_sw)
    return dict(S=S, lk=lk, m0=m0, mask=mask, live=live, ist=ist, isw=isw, a_req_st=a_req_st, a_req_sw=a_req_sw, alpha0=alpha0, b_fill=b_fill, kt=kt, f=f, q90=q90, cap1=cap1,
                base_rows=int(keep.size), rho_d=rho_d, c_w=c_w, n_excluded_allowed=int(np.sum(~allowed)), T=T, swt=sw, taper=taper, excl=excl, place=place, n_sw=n_sw, side=side, lines=lines)


R4.design = design_edge  # point_task of sub6_board_island2 calls R4.design at call time (also in spawned workers: this module is imported there)
R4.tables = edge_tables


# ------------------------------------------------------------------ extras
def line_stats(D):
    lines = [l for l in D["lines"] if D["mask"][l["cells"]].any()]
    L = np.array([l["L"] for l in lines])
    volts = (L - 1) * V_DIODE
    return dict(n_lines=len(lines), n_lines_x=sum(l["axis"] == "x" for l in lines), n_lines_y=sum(l["axis"] == "y" for l in lines), line_len_min=int(L.min()) if len(L) else 0, line_len_max=int(L.max()) if len(L) else 0,
                line_len_mean=float(L.mean()) if len(L) else 0.0, n_diodes_in_strings=int((L - 1).sum()), string_V_min=float(volts.min()) if len(L) else 0.0, string_V_max=float(volts.max()) if len(L) else 0.0,
                string_V_total_if_series=float(volts.sum()))


def theta_stats(D, ph, direction, flips, tol_deg=7.5):
    """Angle between the hologram slot axis of every switched cell and its string (0 = slot along the string, 90 = across = the tabulated row).  The edge-diode table has theta = 0 / 45 / 90 only:
    cells whose angle is NOT within tol_deg of one of them need an angle table that does not exist yet (NEED_ANGLE_TABLE); nothing is interpolated or guessed."""
    S = D["S"]
    phiR = np.arctan2(S.y, S.x)
    chi = np.arctan2(direction[:, 1], direction[:, 0])
    axis = np.mod(0.5 * (chi + phiR) + flips * math.pi / 2, math.pi)
    ang_line = np.zeros(S.n)
    for l in D["lines"]:
        ang_line[l["cells"]] = 0.0 if l["axis"] == "x" else math.pi / 2
    rel = np.degrees(np.abs(((axis - ang_line + math.pi / 2) % math.pi) - math.pi / 2))
    m = D["mask"] & D["live"]
    r = rel[m]
    near = np.min(np.abs(r[:, None] - np.array([0.0, 45.0, 90.0])[None, :]), axis=1) <= tol_deg
    near90 = np.abs(r - 90.0) <= tol_deg
    hist, _ = np.histogram(r, bins=[0, 15, 30, 45, 60, 75, 90.001])
    return dict(theta_rel_mean_deg=float(r.mean()), theta_frac_within_tol_of_0_45_90=float(near.mean()), theta_frac_within_tol_of_90=float(near90.mean()), theta_frac_need_angle_table=float(1.0 - near.mean()),
                theta_hist_15deg_bins="/".join(str(int(h)) for h in hist), theta_tol_deg=tol_deg, NEED_ANGLE_TABLE=bool(near.mean() < 1.0))


def beta_effect(D, sa, ph, direction, flips, launch):
    S = D["S"]
    fl = flips.astype(float)
    out = {}
    jx, jy, _, _ = R4.currents(D, sa, ph, direction, fl, "on")
    r_on = B.metrics(S, *S.far(jx, jy), ph["eta_on"], launch)
    jx, jy, _, _ = R4.currents(D, sa, dict(ph, Pon=ph["Poff"]), direction, fl, "on")
    r_cf = B.metrics(S, *S.far(jx, jy), ph["eta_on"], launch)
    stat = ~D["mask"] & D["live"]
    dphi = np.degrees(ph["Pon"] - ph["Poff"])
    a_off = np.sqrt(np.maximum(2.0 * sa["alpha_off"] * np.maximum(ph["Woff"], 0) / np.maximum(np.hypot(S.x, S.y), 1e-9), 0))
    w = a_off**2 * stat
    out["on_GA_beta_shift_on"] = r_on["gA"]
    out["on_GA_beta_shift_off"] = r_cf["gA"]
    out["beta_shift_cost_dB"] = r_cf["gA"] - r_on["gA"]
    out["static_dphi_rms_deg"] = float(np.sqrt(np.sum(w * dphi**2) / max(np.sum(w), 1e-30)))
    out["static_dphi_max_deg"] = float(np.max(np.abs(dphi[stat])))
    m = D["mask"] & D["live"]
    out["sw_dbeta_on_minus_off"] = float(np.mean((sa["beta_on"] - sa["beta_off"])[m]))
    return out


DESIGN = None  # sub6_board_edge2 sets this to its chord-aware design function


def edge_task(a):
    side, n, pl, taper, excl, launch = a
    T = edge_tables()
    D = (DESIGN or design_edge)(side, n, pl, taper, excl, T)
    sa = R4.site_arrays(D)
    ph = R4.physics(D, sa)
    dirs = {w: X.directions_w(D, ph, w) for w in X.W_GRID}
    single = {}
    for w in X.W_GRID:
        fl, info = R4.optimise_flips(D, sa, ph, dirs[w])
        single[w] = (fl, X.metrics_state(D, sa, ph, dirs[w], fl, launch), info)
    res = {"V0": (0.0,) + single[0.0], "V1": (0.5,) + single[0.5]}
    wstar = max(X.W_GRID, key=lambda w: single[w][1]["minGA"])
    res["V2"] = (wstar,) + single[wstar]
    best = (wstar,) + single[wstar]
    best = (best[0], best[1], best[2], dict(best[3], ring_weight=0.0))
    for w in sorted({0.0, wstar}):
        for rwt in (0.1, 0.3, 1.0):
            fl3, info3 = X.optimise_flips_multi(D, sa, ph, dirs[w], ring_weight=rwt)
            m3 = X.metrics_state(D, sa, ph, dirs[w], fl3, launch)
            info3 = dict(info3, ring_weight=rwt)
            if m3["margin"] > best[2]["margin"]:
                best = (w, fl3, m3, info3)
    res["V3"] = best
    ls = line_stats(D)
    rows = []
    for v, (w, fl, m, info) in res.items():
        row = dict(variant=v, side_mm=side * 1e3, N_sw=n, N_sw_live=int((D["mask"] & D["live"]).sum()), placement=pl, taper=taper, excl=excl, w=w, rotated=int(fl.sum()), ring_weight=info.get("ring_weight", 0.0), rho_d=D["rho_d"], **m, **ls)
        if v in ("V0", "V3"):
            row.update(beta_effect(D, sa, ph, dirs[w], fl, launch))
            row.update(theta_stats(D, ph, dirs[w], fl))
        rows.append(row)
    return rows


def run(workers):
    RES.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    tasks = [(s, n, pl, t, ex, LAUNCH) for s in SIDES for n in N_LIST for pl in ("a", "b") for t in (0.5, 0.8, 1.0) for ex in (0, 1)]
    T = edge_tables()
    print(f"{len(tasks)} tasks, {workers} workers; static rows {len(T['st']['a'])}, switched (edge theta 90, filtered) rows {len(T['sw']['aoff'])} of {T['n_edge_theta90']}", flush=True)
    allrows = []
    with Pool(min(workers, len(tasks))) as pool:
        for rows in pool.imap_unordered(edge_task, tasks, chunksize=1):
            allrows += rows
    allrows.sort(key=lambda r: (-r["side_mm"], r["N_sw"], r["placement"], r["taper"], r["excl"], r["variant"]))
    B.write_csv(RES / "points_edge.csv", allrows)
    prov = B.provenance(dict(stage="run_edge", workers=min(workers, len(tasks)), n_tasks=len(tasks), seconds=time.time() - t0, launch_fraction=LAUNCH, static_rows=len(T["st"]["a"]), switched_rows=len(T["sw"]["aoff"]),
                             switched_ids=T["sw"]["id"], line_len=[LMIN, LMAX], v_per_diode=V_DIODE))
    (RES / "provenance_run_edge.json").write_text(json.dumps(prov, indent=1, default=str), encoding="utf-8", newline="\n")
    print(f"done {time.time() - t0:.0f}s, {len(allrows)} rows", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["run", "test"])
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args()
    if a.stage == "run":
        run(a.workers)
    else:
        t = time.time()
        rows = edge_task((0.2, 135, "b", 0.5, 1, LAUNCH))
        for r in rows:
            print(r["variant"], r["gate"], round(r["margin"], 2), "GA", round(r["off_GA"], 2), round(r["on_GA"], 2), "AR3", round(r["off_ARA3"], 2), round(r["off_ARB3"], 2), round(r["on_ARA3"], 2), "supp", round(r["supp"], 1), "sll", round(r["sll_off"], 1), round(r["sll_on"], 1),
                  "lines", r["n_lines"], r["line_len_min"], r["line_len_max"], "N", r["N_sw_live"])
        print({k: round(v, 3) for k, v in rows[3].items() if isinstance(v, float) and ("beta" in k or "dphi" in k or "theta" in k)}, round(time.time() - t, 1))
