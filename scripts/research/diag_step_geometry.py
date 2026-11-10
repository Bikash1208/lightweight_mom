"""DIAGNOSTIC item 1: measure the cells in BOARD_N240_RF_full_win.step itself and compare with the design (fd = sub6_final_design.final_design: site x/y, slot axis psi, class).  MODEL_ONLY.
    SUB6_MESH=96 python diag_step_geometry.py <step> <out.json>
Per cell: disc centre (outer circle of the top face z = 0.035), slot polygon (inner wire of the top face): long-axis angle (edge lengths grouped by direction) and area centroid;
switched cells: island top-face centroid / long axis, diode sheet D<string>_<pos> centre / direction.  All 758 cells are measured; mirror / rotation / swap hypotheses are tested on the residuals.
"""
import json
import math
import os
import sys
from pathlib import Path

os.environ.setdefault("SUB6_MESH", "96")
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import numpy as np

import sub6_final_design as FD
from sidequest_board_k2_lp_package import read_names


def poly_info(pts):
    """pts: (n,2) ordered polygon vertices -> (centroid, long-axis angle deg mod 180, area)"""
    p = np.asarray(pts)
    q = np.roll(p, -1, axis=0)
    cr = p[:, 0] * q[:, 1] - q[:, 0] * p[:, 1]
    a = 0.5 * cr.sum()
    cx = ((p[:, 0] + q[:, 0]) * cr).sum() / (6 * a)
    cy = ((p[:, 1] + q[:, 1]) * cr).sum() / (6 * a)
    ang = {}
    for u, v in zip(p, q):
        d = v - u
        ln = math.hypot(*d)
        if ln < 1e-9:
            continue
        th = round(math.degrees(math.atan2(d[1], d[0])) % 180.0, 1)
        ang[th] = ang.get(th, 0.0) + ln
    # merge the two directions that are within 0.2 deg (rounding of the slightly bent edges)
    keys = sorted(ang)
    best = max(keys, key=lambda k: sum(v for kk, v in ang.items() if abs(((kk - k + 90) % 180) - 90) < 0.3))
    num = sum(v * (kk - best if abs(kk - best) < 90 else kk - best + (180 if kk < best else -180)) for kk, v in ang.items() if abs(((kk - best + 90) % 180) - 90) < 0.3)
    den = sum(v for kk, v in ang.items() if abs(((kk - best + 90) % 180) - 90) < 0.3)
    return (float(cx), float(cy)), float((best + num / den) % 180.0), float(abs(a))


def wire_pts(w):
    pts = []
    for e in w.Edges():
        pts.append((e.startPoint().x, e.startPoint().y))
    # order edges by connectivity: cadquery returns them in wire order for a polygon
    return pts


def main():
    step, outp = sys.argv[1], sys.argv[2]
    names = read_names(step)
    metal = names["METAL_PEC"]
    lumps = metal.Solids()
    fd = FD.final_design("cal")
    D = fd["D"]
    S = D["S"]
    xs, ys = S.x * 1e3, S.y * 1e3
    mask = D["mask"] & D["live"]
    axis = np.mod(np.degrees(fd["axis"]), 180.0)
    rt = json.load(open(Path(__file__).resolve().parents[2] / "results" / "sub6_final_board" / "mesh96" / "routing.json")) if False else None
    discs, islands = [], []
    for l in lumps:
        tops = [f for f in l.Faces() if abs(f.Center().z - 0.035) < 1e-6 and f.normalAt().z > 0.5]
        if not tops:
            continue
        f = max(tops, key=lambda f: f.Area())
        circ = [e for e in f.outerWire().Edges() if e.geomType() == "CIRCLE"]
        if circ and f.innerWires():
            c = circ[0].arcCenter()
            r = circ[0].radius()
            if r < 4.0:   # via pad / barrel top, not a cell disc
                continue
            iw = f.innerWires()[0]
            vs = wire_pts(iw)
            cen, ang, area = poly_info(vs)
            discs.append(dict(cx=c.x, cy=c.y, r=r, slot_c=cen, slot_ang=ang, slot_area=area, n_slot_vertices=len(vs)))
        elif not circ and not f.innerWires():
            vs = wire_pts(f.outerWire())
            if len(vs) in (4, 6) and f.Area() < 8:
                cen, ang, area = poly_info(vs)
                islands.append(dict(c=cen, ang=ang, area=area, n=len(vs)))
    sheets = {}
    for n, s in names.items():
        if n.startswith("D") and "_" in n and n[1:2].isdigit():
            ctr = s.Center() if hasattr(s, "Center") else None
            vs = [(v.X, v.Y) for v in s.Vertices()]
            sheets[n] = np.mean(np.array(vs), axis=0), vs
    print("measured discs", len(discs), "islands", len(islands), "sheets", len(sheets), flush=True)
    # ---- match discs to design sites
    dc = np.array([[d["cx"], d["cy"]] for d in discs])
    res = []
    used = set()
    for i in range(S.n):
        j = int(np.argmin(np.hypot(dc[:, 0] - xs[i], dc[:, 1] - ys[i])))
        d = discs[j]
        dpos = math.hypot(d["cx"] - xs[i], d["cy"] - ys[i])
        dax = ((d["slot_ang"] - axis[i] + 90) % 180) - 90
        res.append(dict(site=i, cls="SWITCHED" if mask[i] else "STATIC", design_xy=(float(xs[i]), float(ys[i])), step_disc_xy=(round(d["cx"], 4), round(d["cy"], 4)), pos_err_mm=dpos,
                        design_axis_deg=float(axis[i]), step_axis_deg=round(d["slot_ang"], 3), axis_err_deg=float(dax), slot_centre_offset_mm=math.hypot(d["slot_c"][0] - d["cx"], d["slot_c"][1] - d["cy"]),
                        disc_r=d["r"], slot_area=d["slot_area"]))
        used.add(j)
    pe = np.array([r["pos_err_mm"] for r in res])
    ae = np.array([abs(r["axis_err_deg"]) for r in res])
    summ = dict(n_discs=len(discs), n_sites=int(S.n), n_unmatched_discs=len(discs) - len(used), pos_err_max_mm=float(pe.max()), pos_err_mean_mm=float(pe.mean()), axis_err_max_deg=float(ae.max()), axis_err_mean_deg=float(ae.mean()),
                slot_centre_offset_max_mm=float(max(r["slot_centre_offset_mm"] for r in res)), n_switched=int(mask.sum()), n_static=int((~mask).sum()))
    # ---- hypotheses on the axis residual and on the positions (all 758 cells)
    hyp = {}
    stepx = np.array([r["step_disc_xy"][0] for r in res])
    stepy = np.array([r["step_disc_xy"][1] for r in res])
    stepa = np.array([r["step_axis_deg"] for r in res])
    cand = {"identity": (stepx, stepy, stepa), "x-mirror": (-stepx, stepy, 180 - stepa), "y-mirror": (stepx, -stepy, -stepa), "rot180": (-stepx, -stepy, stepa), "x<->y swap": (stepy, stepx, 90 - stepa), "axis+90": (stepx, stepy, stepa + 90)}
    for k, (cx, cy, ca) in cand.items():
        hyp[k] = dict(pos_err_max_mm=float(np.max(np.hypot(cx - xs, cy - ys))), axis_err_max_deg=float(np.max(np.abs(((ca - axis + 90) % 180) - 90))), axis_err_mean_deg=float(np.mean(np.abs(((ca - axis + 90) % 180) - 90))))
    # ---- switched cells: island and diode sheet vs the design (table angle = slot axis, +180 where the string optimisation rotated the cell)
    rt = json.load(open(Path(__file__).resolve().parents[2] / "results" / "sub6_final_board" / "mesh96" / "routing.json"))
    cells = rt["cells"]
    ic = np.array([i["c"] for i in islands])
    sw_rows = []
    for i in np.nonzero(mask)[0]:
        c = cells[int(i)]
        thb = c.get("theta_board", c["theta"])
        exp_isl = np.array([xs[i], ys[i]]) + np.array([-math.sin(thb), math.cos(thb)]) * (-0.20)
        k = int(np.argmin(np.hypot(ic[:, 0] - exp_isl[0], ic[:, 1] - exp_isl[1])))
        sw_rows.append(dict(site=int(i), island_err_mm=float(math.hypot(*(ic[k] - exp_isl))), island_axis_err_deg=float(((islands[k]["ang"] - math.degrees(thb) + 90) % 180) - 90), board_rotation_deg=float(math.degrees(thb) % 360)))
    summ["island_pos_err_max_mm"] = float(max(r["island_err_mm"] for r in sw_rows))
    summ["island_axis_err_max_deg"] = float(max(abs(r["island_axis_err_deg"]) for r in sw_rows))
    # diode sheets: centre and direction (from the sheet corners) vs design (centre at cell-local (0, 0.95) rotated by the board angle)
    shr = []
    for i in np.nonzero(mask)[0]:
        c = cells[int(i)]
        thb = c.get("theta_board", c["theta"])
        exp_c = np.array([xs[i], ys[i]]) + np.array([-math.sin(thb), math.cos(thb)]) * 0.95
        best = min(sheets.items(), key=lambda kv: np.hypot(*(kv[1][0] - exp_c)))
        shr.append(dict(name=best[0], err_mm=float(np.hypot(*(best[1][0] - exp_c)))))
    summ["sheet_centre_err_max_mm"] = float(max(r["err_mm"] for r in shr))
    summ["sheets_matched_unique"] = len({r["name"] for r in shr})
    # ---- >= 30 cells spread over the board (12 switched + 24 static, by azimuth sectors)
    order = np.argsort(np.degrees(np.arctan2(ys, xs)))
    pick = []
    for cls, n in (("SWITCHED", 12), ("STATIC", 24)):
        ids = [i for i in order if res[i]["cls"] == cls]
        pick += [ids[int(k * len(ids) / n)] for k in range(n)]
    table = [res[i] for i in pick]
    out = dict(label="MODEL_ONLY diagnostic: STEP geometry vs design", step=str(step), summary=summ, hypotheses_on_all_cells=hyp, sample_cells=table, island_rows_head=sw_rows[:5])
    Path(outp).write_text(json.dumps(out, indent=1, default=str), encoding="utf-8", newline="\n")
    print(json.dumps(dict(summary=summ, hypotheses=hyp), indent=1, default=str))


if __name__ == "__main__":
    main()
