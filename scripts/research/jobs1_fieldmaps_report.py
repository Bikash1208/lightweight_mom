"""JOB 1 A report: maps (PNG), score table, 5 lowest-field DC-reachable via candidates per body, current via spots.  MODEL_ONLY / NOT_HFSS_VALIDATED.
    python jobs1_fieldmaps_report.py
Score of a point (per case = state x propagation direction): s = (|Et|_top/rms + |Ez|_top/rms + |Ez|_gnd/rms) / 3, each map normalised to its own rms over the cell; the point score is the MAX over the cases
(island: OFF/ON x along/across, static: along/across).  Candidates: on copper of the body inside BOTH the MoM mask and the STEP polygon, >= MARGIN from every copper edge (pad r 0.30 + 0.05),
>= 1.0 mm apart.  Island = 'island' body (DC-isolated, own via); disc = disc body without the slot (own via).
"""
import csv
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import diag_full_tensor_lib as L

FM = ROOT / "results" / "jobs1" / "fieldmaps"
MARGIN = 0.35
meta = json.load(open(FM / "meta.json"))
cells = list(csv.DictReader(open(ROOT / "results/sub6_final_board/mesh96/fab/cell_table.csv")))
sw = [c for c in cells if c["class"] == "SWITCHED"]
c0 = sw[0]
th0 = math.radians(float(c0["slot_axis_deg"]))


def local(c, key):
    x, y, th = float(c["x_mm"]), float(c["y_mm"]), math.radians(float(c["slot_axis_deg"]))
    dx, dy = float(c[key + "_via_x"]) - x, float(c[key + "_via_y"]) - y
    return (math.cos(th) * dx + math.sin(th) * dy, -math.sin(th) * dx + math.cos(th) * dy)


cur = {k: [local(c, k) for c in sw] for k in ("disc", "island")}
cur_spread = {k: (float(np.ptp([p[0] for p in v])), float(np.ptp([p[1] for p in v]))) for k, v in cur.items()}
cur_via = {k: cur[k][0] for k in cur}
R_DISC, LH, WH = float(c0["disc_r_mm"]), float(c0["slot_L_mm"]) / 2, float(c0["slot_W_mm"]) / 2


def load(tag, state, d):
    z = np.load(FM / f"maps_{tag}_{state}_{d}.npz")
    et = np.hypot(np.abs(z["Etx_top"]), np.abs(z["Ety_top"]))
    return z["xs_mm"], et, np.abs(z["Ez_top"]), np.abs(z["Ez_gnd"])


def norm(a):
    return a / math.sqrt(np.mean(a ** 2))


def cases(tag, states):
    out = {}
    for st in states:
        for d in "xy":
            xs, et, ezt, ezg = load(tag, st, d)
            out[(st, d)] = (norm(et), norm(ezt), norm(ezg))
    return xs, out


def score_map(cs):
    return np.max([(a + b + c) / 3 for a, b, c in cs.values()], axis=0)


# ---- masks on the map grid (x = first index)
def island_masks(xs):
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    dxm = L.P / 96
    n_l, n_w, n_i, cb = meta["island_cells"]
    gap = L.GAP_CELLS
    y_top, y_bot = (n_w - gap) * dxm, (-n_w + cb) * dxm
    mom_isl = (np.abs(X) <= n_i * dxm) & (Y >= y_bot) & (Y <= y_top)
    slot = (np.abs(X) <= n_l * dxm) & (np.abs(Y) <= n_w * dxm)
    mom_disc = (np.hypot(X, Y) <= 5.435) & ~slot
    step_isl = (np.abs(X) <= 2.174) & (Y >= -0.85) & (Y <= 0.45)
    step_disc = (np.hypot(X, Y) <= R_DISC) & ~((np.abs(X) <= LH) & (np.abs(Y) <= WH))
    return X, Y, mom_isl & step_isl, mom_disc & step_disc


def erode(mask, X, Y, r, finer):
    """points of `mask` at least r from the mask boundary (distance on the fine grid)"""
    from scipy.ndimage import distance_transform_edt
    d = distance_transform_edt(mask, sampling=float(X[1, 0] - X[0, 0]))
    return mask & (d >= r)


def pick(score, region, X, Y, n=5, sep=1.0):
    idx = np.argwhere(region)
    order = np.argsort(score[region])
    chosen = []
    for k in order:
        i, j = idx[k]
        p = (float(X[i, j]), float(Y[i, j]))
        if all(math.hypot(p[0] - q[0], p[1] - q[1]) >= sep for q, _ in chosen):
            chosen.append((p, float(score[i, j])))
        if len(chosen) == n:
            break
    return chosen


def at(score, xs, p):
    i = int(np.clip(np.rint(p[0] / (xs[1] - xs[0])) + len(xs) // 2, 0, len(xs) - 1))
    j = int(np.clip(np.rint(p[1] / (xs[1] - xs[0])) + len(xs) // 2, 0, len(xs) - 1))
    return float(score[i, j])


def draw(tag, states, cs, xs, sc, outline, marks, fn, title):
    nrow = len(states) * 2
    fig, ax = plt.subplots(nrow, 3, figsize=(12, 3.7 * nrow))
    r = 0
    names = ["|E_t| top (z=+0.1 mm)", "|E_z| top (z=+0.1 mm)", "|E_z| ground plane (z=-h)"]
    for st in states:
        for d in "xy":
            for k in range(3):
                a = ax[r, k]
                im = a.imshow(np.minimum(cs[(st, d)][k], 4.0).T, origin="lower", extent=[xs[0], xs[-1], xs[0], xs[-1]], cmap="magma", vmin=0, vmax=3)
                a.contour(xs, xs, outline.T.astype(float), levels=[0.5], colors="c", linewidths=0.6)
                for lab, p, col, mk in marks:
                    a.plot(p[0], p[1], mk, color=col, ms=7, mfc="none")
                a.set_title(f"{st} kt||{'x (along)' if d == 'x' else 'y (across)'}: {names[k]}", fontsize=8)
                a.set_xlim(-6, 6); a.set_ylim(-6, 6)
            r += 1
    fig.colorbar(im, ax=ax, shrink=0.6, label="normalised to cell rms (clipped at 3)")
    fig.suptitle(title + "  [cyan: copper outline; red o: current vias; green o + number: low-field candidates]", fontsize=9)
    fig.savefig(fn, dpi=110)
    plt.close(fig)


def main():
    res = dict(label="MODEL_ONLY / NOT_HFSS_VALIDATED field maps", margin_mm=MARGIN, current_via_local_uv_mm=cur_via, current_via_spread_over_cells_mm=cur_spread)
    xsi, csi = cases("island", ["OFF_BAR50", "ON_BAR50"])
    X, Y, isl, disc = island_masks(xsi)
    s_isl = score_map(csi)
    cand = {}
    for nm, region in (("island", isl), ("disc", disc)):
        reg = erode(region, X, Y, MARGIN, None)
        cand[nm] = pick(s_isl, reg, X, Y)
    marks = [("cur", cur_via["island"], "red", "o"), ("cur", cur_via["disc"], "red", "o")]
    for nm in cand:
        for k, (p, s) in enumerate(cand[nm]):
            marks.append((f"{nm}{k + 1}", p, "lime", "o"))
    outline = isl | disc
    for st in ("OFF_BAR50", "ON_BAR50"):
        draw("island", [st], {k: v for k, v in csi.items() if k[0] == st}, xsi, s_isl, outline, marks, FM / f"fig_island_{st}.png", f"Island cell {st}")
    res["island_cell"] = dict(score_current_island_via=at(s_isl, xsi, cur_via["island"]), score_current_disc_via=at(s_isl, xsi, cur_via["disc"]),
                              rms_normalised_score_median_island_body=float(np.median(s_isl[isl])), median_disc_body=float(np.median(s_isl[disc])),
                              candidates={nm: [dict(rank=k + 1, u_mm=p[0], v_mm=p[1], score=s) for k, (p, s) in enumerate(cand[nm])] for nm in cand},
                              per_case_at_current_vias={f"{st}|{d}": dict(island_via=[float(at(v, xsi, cur_via['island'])) for v in cs_], disc_via=[float(at(v, xsi, cur_via['disc'])) for v in cs_]) for (st, d), cs_ in csi.items()},
                              per_case_at_best={f"{st}|{d}": dict(island=[float(at(v, xsi, cand['island'][0][0])) for v in cs_], disc=[float(at(v, xsi, cand['disc'][0][0])) for v in cs_]) for (st, d), cs_ in csi.items()})
    # static cell
    xss, css = cases("static", ["NOCUT"])
    Xs, Ys = np.meshgrid(xss, xss, indexing="ij")
    rs = meta["static_row"]["r"]; Ls, Ws = meta["static_row"]["L"] / 2, meta["static_row"]["W"] / 2
    dmask = (np.hypot(Xs, Ys) <= rs) & ~((np.abs(Xs) <= Ls) & (np.abs(Ys) <= Ws))
    s_st = score_map(css)
    cst = pick(s_st, erode(dmask, Xs, Ys, MARGIN, None), Xs, Ys)
    marks_s = [("c", p, "lime", "o") for p, _ in cst]
    draw("static", ["NOCUT"], css, xss, s_st, dmask, marks_s, FM / "fig_static_NOCUT.png", f"Static NOCUT cell ({meta['static_row']['id']}); no via on the board today")
    res["static_cell"] = dict(row=meta["static_row"]["id"], median_score_disc=float(np.median(s_st[dmask])), candidates=[dict(rank=k + 1, u_mm=p[0], v_mm=p[1], score=s) for k, (p, s) in enumerate(cst)])
    # field statistics: where is the field strong (island vs disc body, along/across)
    res["notes"] = dict(maps_normalised="each map divided by its own rms over the 11.6 mm cell")
    Path(FM / "fieldmaps_summary.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline="\n")
    best = dict(island=cand["island"][0][0], disc=cand["disc"][0][0])
    Path(FM / "best_spots_uv.json").write_text(json.dumps(best), encoding="utf-8")
    print(json.dumps(res, indent=1, default=float)[:6000])


if __name__ == "__main__":
    main()
