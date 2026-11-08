"""SUB6 FINAL BOARD step 3/4 export: routes the strings (sub6_final_layout.build), verifies clearances with an independent sampler, writes
DXF (layers), drill CSV (Excellon-like), cell table CSV, string table CSV/JSON.  PRELIMINARY_EXTRAPOLATED_CELL / MODEL_ONLY / NOT_HFSS_VALIDATED.
    python sub6_final_export.py
All coordinates rounded to 0.001 mm (sub6_step_exact grid).
"""
import csv
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ezdxf
import numpy as np
from scipy.spatial import cKDTree

import sub6_final_layout as L

OUT = L.OUT
FAB = OUT / "fab"
R3 = lambda v: round(float(v) + 0.0, 3)


def rect_pts(cx, cy, hx, hy, th, ox=0.0, oy=0.0):
    """rectangle (local centre (ox,oy), half sizes) in a local frame rotated by th about (cx,cy)."""
    pts = []
    for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
        v = L.rot(np.array([ox + sx * hx, oy + sy * hy]), th)
        pts.append((R3(cx + v[0]), R3(cy + v[1])))
    return pts


def verify_full(res):
    """Independent clearance check, per class (copper EDGE to copper edge, mm; requirement >= 0.15 for every class)."""
    nets = [n for n in res["nets"] if n["poly"] is not None]
    hw = L.W_TR / 2
    samples = [np.array(L.resample(n["poly"], 0.02)) for n in nets]
    pts = np.vstack(samples)
    lab = np.concatenate([np.full(len(s), k) for k, s in enumerate(samples)])
    tree = cKDTree(pts)
    pairs = tree.query_pairs(1.0, output_type="ndarray")
    pairs = pairs[lab[pairs[:, 0]] != lab[pairs[:, 1]]]
    dmin = float(np.linalg.norm(pts[pairs[:, 0]] - pts[pairs[:, 1]], axis=1).min()) if len(pairs) else 1.0
    vp = np.array([v for ci, (vi, vd) in res["vias"].items() for v in (vi, vd)])
    pad_c = np.array(res["pads"])
    vt, pt = cKDTree(vp), cKDTree(pad_c)
    d_via, d_pad = 9.0, 9.0
    for k, n in enumerate(nets):
        s = samples[k]
        ends = [n["p0"], n["p1"]]
        d, j = vt.query(s, k=3)
        for r in range(3):
            ok = np.array([min(np.linalg.norm(vp[jj] - e) for e in ends) > 1e-6 for jj in j[:, r]])
            if ok.any():
                d_via = min(d_via, float(d[ok, r].min()))
        # connector pads are 1.0 x 1.0 squares: distance from the sample to the square
        idx = pt.query_ball_point(s, 2.0)
        for i, lst in enumerate(idx):
            for jj in lst:
                if min(np.linalg.norm(pad_c[jj] - e) for e in ends) <= 1e-6:
                    continue
                q = np.abs(s[i] - pad_c[jj]) - 0.5
                d_pad = min(d_pad, float(np.hypot(max(q[0], 0.0), max(q[1], 0.0))) if max(q) > 0 else 0.0)
    out_poly = res["outline"]
    edge_min = 9.0
    for s in samples:
        for a, b in zip(out_poly, np.roll(out_poly, -1, axis=0)):
            ln = math.hypot(*(b - a))
            cr = ((b[0] - a[0]) * (s[:, 1] - a[1]) - (b[1] - a[1]) * (s[:, 0] - a[0])) / ln
            edge_min = min(edge_min, float(cr.min()))
    r_orig = min(float(np.hypot(s[:, 0], s[:, 1]).min()) for s in samples)
    cl = dict(
        trace_to_trace_gap_mm=dmin - L.W_TR,
        trace_edge_to_foreign_antipad_edge_mm=d_via - hw - L.ANTIPAD_R,
        trace_edge_to_foreign_via_pad_edge_mm=d_via - hw - L.PAD_R,
        trace_edge_to_foreign_connector_pad_edge_mm=d_pad - hw,
        trace_edge_to_board_outline_mm=edge_min - hw,
        trace_edge_to_probe_opening_edge_mm=r_orig - hw - 2.05,
    )
    return dict(min_clearances_mm={k: round(v, 4) for k, v in cl.items()}, requirement_mm=0.15, all_ge_0p15=bool(all(v >= 0.15 - 1e-9 for v in cl.values())), antipad_crossings=int(cl["trace_edge_to_foreign_antipad_edge_mm"] < 0),
                min_trace_distance_to_SMA_axis_mm=r_orig, n_samples=int(len(pts)))


def main():
    res = L.build()
    if res["failed"]:
        raise SystemExit(f"UNROUTED {res['failed']}")
    FAB.mkdir(parents=True, exist_ok=True)
    cells = res["cells"]
    ver = verify_full(res)
    ver['antipad_check_AEMLAB1_patch'] = {k: v for k, v in L.antipad_check(res).items() if k != 'worst_hits'}
    L.save_result(res)
    print(ver, flush=True)
    # ---- cell table
    with open(FAB / "cell_table.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["site", "x_mm", "y_mm", "class", "row_id", "disc_r_mm", "slot_L_mm", "slot_W_mm", "slot_axis_deg", "flipped90", "string", "pos_in_string", "island_via_x", "island_via_y", "disc_via_x", "disc_via_y"])
        sinfo = {}
        for st in res["strings"]:
            for k, ci in enumerate(st["cells"]):
                sinfo[ci] = (st["id"], k + 1)
        for i, c in enumerate(cells):
            th = c.get("theta_board", c["theta"])
            if c["cls"] == "SWITCHED":
                vi, vd = res["vias"][i]
                s_, k_ = sinfo[i]
                vv = [R3(vi[0]), R3(vi[1]), R3(vd[0]), R3(vd[1])]
            else:
                s_, k_, vv = "", "", ["", "", "", ""]
            w.writerow([c["site"], R3(c["x"]), R3(c["y"]), c["cls"], c["row_id"], c["r"], c["L"], c["W"], R3(math.degrees(th) % 360), c["flipped"], s_, k_] + vv)
    # ---- strings
    srows = []
    for st in res["strings"]:
        links = [n for n in res["nets"] if n["kind"] == "link" and n["string"] == st["id"]]
        leads = [n for n in res["nets"] if n["kind"] == "lead" and n["string"] == st["id"]]
        Ll = sum(L.poly_len(n["poly"]) for n in links)
        Ld = sum(L.poly_len(n["poly"]) for n in leads)
        srows.append(dict(string=st["id"], n_diodes=len(st["cells"]), link_trace_mm=R3(Ll), lead_trace_mm=R3(Ld), total_trace_mm=R3(Ll + Ld), V_forward_at_10mA=round(len(st["cells"]) * L.V_DIODE, 2), I_mA=L.I_DIODE * 1e3, pad_start=[n["pad"] for n in leads if n["pos"] == "start"][0] + 1, pad_end=[n["pad"] for n in leads if n["pos"] == "end"][0] + 1))
    with open(FAB / "string_table.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(srows[0]))
        w.writeheader()
        w.writerows(srows)
    # ---- drill
    holes = []
    for ci, (vi, vd) in res["vias"].items():
        holes += [("VIA_ISLAND", ci, vi), ("VIA_DISC", ci, vd)]
    holes.append(("SMA_PIN", -1, np.zeros(2)))
    with open(FAB / "drill.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["type", "site", "x_mm", "y_mm", "drill_dia_mm", "pad_dia_mm", "antipad_dia_mm"])
        for t, ci, p in holes:
            if t == "SMA_PIN":
                w.writerow([t, "", 0.0, 0.0, 1.27, "", ""])
            else:
                w.writerow([t, cells[ci]["site"], R3(p[0]), R3(p[1]), 0.3, 0.6, 1.0])
    with open(FAB / "drill.drl", "w") as f:
        f.write("M48\nMETRIC,TZ\nT1C0.300\nT2C1.270\n%\nT1\n")
        for t, ci, p in holes:
            if t != "SMA_PIN":
                f.write(f"X{R3(p[0]):.3f}Y{R3(p[1]):.3f}\n")
        f.write("T2\nX0.000Y0.000\nM30\n")
    # ---- DXF
    doc = ezdxf.new("R2010", setup=True)
    doc.units = 4
    msp = doc.modelspace()
    for name, col in (("TOP_DISC", 1), ("TOP_SLOT", 2), ("TOP_ISLAND", 3), ("DIODE_SHEET", 6), ("VIA_DRILL", 5), ("VIA_PAD", 4), ("ANTIPAD", 8), ("BOT_TRACE", 30), ("BOT_PAD", 40), ("OUTLINE", 7), ("SMA_HOLE", 7), ("CONNECTOR", 50)):
        doc.layers.add(name, color=col)
    for i, c in enumerate(cells):
        th = c.get("theta_board", c["theta"])
        cx, cy = R3(c["x"]), R3(c["y"])
        msp.add_circle((cx, cy), c["r"], dxfattribs={"layer": "TOP_DISC"})
        msp.add_lwpolyline(rect_pts(cx, cy, c["L"] / 2, c["W"] / 2, th), close=True, dxfattribs={"layer": "TOP_SLOT"})
        if c["cls"] == "SWITCHED":
            msp.add_lwpolyline(rect_pts(cx, cy, 2.174, 0.65, th, 0.0, -0.20), close=True, dxfattribs={"layer": "TOP_ISLAND"})
            msp.add_lwpolyline(rect_pts(cx, cy, 0.175, 0.5, th, 0.0, 0.95), close=True, dxfattribs={"layer": "DIODE_SHEET"})
            vi, vd = res["vias"][i]
            for p in (vi, vd):
                msp.add_circle((R3(p[0]), R3(p[1])), 0.15, dxfattribs={"layer": "VIA_DRILL"})
                msp.add_circle((R3(p[0]), R3(p[1])), 0.30, dxfattribs={"layer": "VIA_PAD"})
                msp.add_circle((R3(p[0]), R3(p[1])), 0.50, dxfattribs={"layer": "ANTIPAD"})
                msp.add_circle((R3(p[0]), R3(p[1])), 0.30, dxfattribs={"layer": "BOT_PAD"})
    for n in res["nets"]:
        pl = [(R3(a), R3(b)) for a, b in n["poly"]]
        msp.add_lwpolyline(pl, dxfattribs={"layer": "BOT_TRACE", "const_width": L.W_TR})
    for k, p in enumerate(res["pads"]):
        msp.add_lwpolyline(rect_pts(R3(p[0]), R3(p[1]), 0.5, 0.5, 0.0), close=True, dxfattribs={"layer": "CONNECTOR"})
        msp.add_text(f"P{k+1}", dxfattribs={"layer": "CONNECTOR", "height": 0.8, "insert": (R3(p[0] + 1.2), R3(p[1] - 0.4))})
    msp.add_lwpolyline([(R3(a), R3(b)) for a, b in res["outline"]], close=True, dxfattribs={"layer": "OUTLINE"})
    msp.add_circle((0.0, 0.0), 0.635, dxfattribs={"layer": "SMA_HOLE"})
    doc.saveas(FAB / "sub6_final_board.dxf")
    meta = dict(label="PRELIMINARY_EXTRAPOLATED_CELL / MODEL_ONLY / NOT_HFSS_VALIDATED", verification=ver, strings=srows, n_cells=len(cells), n_switched=sum(c["cls"] == "SWITCHED" for c in cells), n_vias=len(holes) - 1, outline_hex_circumradius_mm=R3(res["r_out"]), connector_x_mm=R3(res["conn"][0]), n_nets=len(res["nets"]),
                note="Cell local frame: slot along x, island at local y -0.85..0.45, 1.00 mm diode gap at local +y (y 0.45..1.45), vias island (-1.78,-0.53), disc (-3.63,-3.63); board rotation = slot axis (+180 deg where the string DP chose it). Diode anode on island pad, cathode on disc pad of the same cell.")
    (FAB / "layout_summary.json").write_text(json.dumps(meta, indent=1))
    # trace polylines for the later HFSS/STEP stage
    (FAB / "traces.json").write_text(json.dumps([dict(kind=n["kind"], string=n["string"], pos=n["pos"], poly=[(R3(a), R3(b)) for a, b in n["poly"]]) for n in res["nets"]]))
    print(json.dumps(srows[:3]), flush=True)
    print("EXPORT_DONE", flush=True)


if __name__ == "__main__":
    main()
