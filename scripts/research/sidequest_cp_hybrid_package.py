"""HFSS package: CP option (i) hybrid board.  Fixed RHCP beam 1 (25, 0) from Phase 26 rotating-slot tensor cells + LP beam 2 (25, 180) from 180 GB crossed-slit cells (2 diodes each, ONE group, S_w).  MODEL_ONLY.
Same hexagon, stack, probe and port as `hfss_exports/sidequest_board_k2_lp` (RT5880 1.575, 150 mm hexagon, lattice 6.4).  Tensor cells: disc r 3.0 with ONE slot L x W (centred, long axis at the cell angle), 35 um;
the slot corners are rounded to 1 um (exact grid), so L, W are met to +-0.001 mm.  GB cells: disc r 3.0 with a cross slit 0.4 x 5.6 and the two sheets GB_*_DX / GB_*_DY (0.4 x 1.0 at (0, +0.8) / (+0.8, 0), z 0.0175).
    python sidequest_cp_hybrid_package.py [design.npz]"""
import csv
import json
import math
import shutil
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import cadquery as cq
import numpy as np

import sidequest_board_k2_lp_package as K2P
import sidequest_cp_hybrid as H
from src.transmitarray.export import inexact_literals, snap_step_text, verify_snap

mm = K2P.mm
OUT = ROOT / "hfss_exports/sidequest_board_CP_LP_hybrid"
SYNC = Path(r"C:\hmsa_sync\data\sidequest\board_CP_LP_hybrid")
RES = ROOT / "results/sidequest_hybrid"
PITCH_UM, T_UM, R_UM = K2P.PITCH_UM, K2P.T_UM, K2P.R_DISC_UM


def slot_corners_um(cx_um, cy_um, L_mm, W_mm, ang_rad):
    u = np.array([math.cos(ang_rad), math.sin(ang_rad)])
    v = np.array([-u[1], u[0]])
    pts = []
    for sl, sw in ((1, 1), (-1, 1), (-1, -1), (1, -1)):
        p = np.array([cx_um, cy_um]) + 500.0 * L_mm * sl * u + 500.0 * W_mm * sw * v
        pts.append((int(round(p[0])), int(round(p[1]))))
    return pts


def make_cells(design_file):
    x, y, kind0, grp0 = H.load_cells()
    dz = np.load(design_file)
    gb, flip = dz["gb"], dz["flip"]
    table, m, d = H.tensor_design(x, y)
    ix = np.rint(x / 6.4e-3).astype(int)
    iy = np.rint(y / 6.4e-3).astype(int)
    cells = []
    for i in range(len(x)):
        if gb[i]:
            cells.append(dict(ix=int(ix[i]), iy=int(iy[i]), kind="GB", grp="GB", L=5.6, W=0.4, ang=0.0, pts=None))
        else:
            ang = float(d["slot_angle_rad"][i]) + (math.pi / 2 if flip[i] < 0 else 0.0)
            ang = ang % math.pi
            cells.append(dict(ix=int(ix[i]), iy=int(iy[i]), kind="T", grp="", L=float(d["slot_L_mm"][i]), W=float(d["slot_W_mm"][i]), ang=ang, flip=int(flip[i] < 0),
                              pts=slot_corners_um(ix[i] * PITCH_UM, iy[i] * PITCH_UM, float(d["slot_L_mm"][i]), float(d["slot_W_mm"][i]), ang)))
    return cells, (gb, flip, d, m)


def build(cells):
    t, r = mm(T_UM), mm(R_UM)
    t_open, _ = K2P.templates()
    H_, GT = mm(K2P.HSUB_UM), mm(K2P.GT_UM)
    metal = [K2P.hex_wp(-H_ - GT, -H_).cut(cq.Workplane("XY").workplane(offset=-H_ - GT).circle(mm(K2P.PORT_R_UM)).extrude(GT)).val(),
             cq.Workplane("XY").workplane(offset=-H_).circle(mm(K2P.PIN_R_UM)).extrude(H_ + mm(K2P.PIN_TOP_UM)).val()]
    sheets = []
    for c in cells:
        cx, cy = mm(c["ix"] * PITCH_UM), mm(c["iy"] * PITCH_UM)
        if c["kind"] == "GB":
            metal.append(t_open.moved(cq.Location(cq.Vector(cx, cy, 0.0))))
            tag = f"{c['ix']:+03d}_{c['iy']:+03d}".replace("+", "p").replace("-", "m")
            sheets.append((f"GB_{tag}_DX", K2P.sheet(cx, cy + mm(K2P.POS_UM), mm(K2P.SLOT_G_UM), mm(K2P.WD_UM), t / 2), "GB"))
            sheets.append((f"GB_{tag}_DY", K2P.sheet(cx + mm(K2P.POS_UM), cy, mm(K2P.WD_UM), mm(K2P.SLOT_G_UM), t / 2), "GB"))
        else:
            pts = [(mm(px), mm(py)) for px, py in c["pts"]]
            wp = cq.Workplane("XY").center(cx, cy).circle(r).polyline([(px - cx, py - cy) for px, py in pts]).close().extrude(t)
            metal.append(wp.val())
    sub = K2P.hex_wp(-H_, 0.0).cut(cq.Workplane("XY").workplane(offset=-H_).circle(mm(K2P.PIN_R_UM)).extrude(H_)).val()
    air = K2P.box(0, 0, 2 * mm(K2P.HEX_X_UM + K2P.PAD_XY_UM), 2 * mm(K2P.HEX_R_UM + K2P.PAD_XY_UM), -H_ - GT, mm(K2P.AIR_TOP_UM)).val()
    port = cq.Face.makeFromWires(cq.Wire.makeCircle(mm(K2P.PORT_R_UM), cq.Vector(0, 0, -H_), cq.Vector(0, 0, 1)), [cq.Wire.makeCircle(mm(K2P.PIN_R_UM), cq.Vector(0, 0, -H_), cq.Vector(0, 0, 1))])
    return {"METAL_PEC": cq.Compound.makeCompound(metal), "SUBSTRATE": sub, "AIR": air, "PORT": port}, sheets


def intended(cells):
    offs = sorted({o for o in [0, 200, 300, 1300, 2800, 3000]} | {-o for o in [0, 200, 300, 1300, 2800, 3000]})
    xs, ys = set(), set()
    for c in cells:
        for o in offs:
            xs.add(c["ix"] * PITCH_UM + o)
            ys.add(c["iy"] * PITCH_UM + o)
        if c["pts"]:
            for px, py in c["pts"]:
                xs.add(px)
                ys.add(py)
    glob = [0, K2P.PIN_R_UM, -K2P.PIN_R_UM, K2P.PORT_R_UM, -K2P.PORT_R_UM, K2P.HEX_X_UM, -K2P.HEX_X_UM, K2P.HEX_R_UM, -K2P.HEX_R_UM, K2P.HEX_Y_UM, -K2P.HEX_Y_UM, K2P.HEX_X_UM + K2P.PAD_XY_UM, -K2P.HEX_X_UM - K2P.PAD_XY_UM, K2P.HEX_R_UM + K2P.PAD_XY_UM, -K2P.HEX_R_UM - K2P.PAD_XY_UM]
    xy = sorted({mm(v) for v in xs | ys | set(glob)})
    zs = sorted({mm(v) for v in (0, T_UM, T_UM / 2.0, -K2P.HSUB_UM, -K2P.HSUB_UM - K2P.GT_UM, K2P.PIN_TOP_UM, K2P.AIR_TOP_UM)})
    return dict(xy=np.array(xy), z=np.array(zs), dir=np.array([-1.0, 0.0, 1.0]))


def check(path, cells, sheets_expected):
    parts = K2P.read_names(path)
    res, ok_all = {}, True

    def rec(name, ok, **kw):
        nonlocal ok_all
        ok_all &= bool(ok)
        res[name] = {"ok": bool(ok), **kw}

    names = set(parts)
    sh_names = {n for n in names if n.startswith("GB_")}
    exp_sh = {n for n, _f, _g in sheets_expected}
    rec("bodies_and_names", {"METAL_PEC", "SUBSTRATE", "AIR", "PORT"} <= names and sh_names == exp_sh and len(names) == 4 + len(exp_sh), n_objects=len(names), n_sheets=len(sh_names))
    solids = parts["METAL_PEC"].Solids()
    rec("metal_single_part_lumps", len(solids) == 2 + len(cells), lumps=len(solids), expected=2 + len(cells))
    t_open, _ = K2P.templates()
    vol_gb = t_open.Volume()
    got = {}
    for s in solids:
        bb = s.BoundingBox()
        if abs(bb.zmin) < 1e-6 and abs(bb.zmax - mm(T_UM)) < 1e-6:
            cx, cy = (bb.xmin + bb.xmax) / 2, (bb.ymin + bb.ymax) / 2
            got[(round(cx / 6.4), round(cy / 6.4))] = (s.Volume(), abs(bb.xlen - 6.0) < 1e-6 and abs(bb.ylen - 6.0) < 1e-6 and abs(cx - mm(round(cx / 6.4) * PITCH_UM)) < 1e-6 and abs(cy - mm(round(cy / 6.4) * PITCH_UM)) < 1e-6)
    bad = 0
    for c in cells:
        k = (c["ix"], c["iy"])
        if k not in got or not got[k][1]:
            bad += 1
            continue
        if c["kind"] == "GB":
            exp = vol_gb
        else:
            px = np.array([p[0] for p in c["pts"]]) / 1000.0
            py = np.array([p[1] for p in c["pts"]]) / 1000.0
            area = 0.5 * abs(np.dot(px, np.roll(py, -1)) - np.dot(py, np.roll(px, -1)))
            exp = (math.pi * 9.0 - area) * 0.035
        bad += abs(got[k][0] - exp) > 2e-6
    rec("patch_vs_design_map", bad == 0 and len(got) == len(cells), n_patches=len(got), bad=int(bad))
    # tensor slot sizes: the rounded polygon sides against L, W
    worst = 0.0
    for c in cells:
        if c["pts"]:
            p = np.array(c["pts"], float) / 1000.0
            sides = [np.hypot(*(p[(i + 1) % 4] - p[i])) for i in range(4)]
            worst = max(worst, abs(sides[0] - c["L"]), abs(sides[2] - c["L"]), abs(sides[1] - c["W"]), abs(sides[3] - c["W"]))
    rec("tensor_slot_sizes_vs_table", worst <= 0.0016, max_deviation_mm=round(worst, 4), note="corners rounded to 1 um; L x W of the Phase 26 table within +-0.0015 mm")
    sub = parts["SUBSTRATE"]
    gp = cq.Compound.makeCompound([s for s in solids if s.BoundingBox().zmin < -1e-6])
    pat = cq.Compound.makeCompound([s for s in solids if s.BoundingBox().zmin > -1e-6 and s.BoundingBox().zmax < 0.04])
    c1, c2 = sub.intersect(gp).Volume(), sub.intersect(pat).Volume()
    rec("no_intersections", c1 < 1e-9 and c2 < 1e-9, substrate_x_ground_probe=c1, substrate_x_patches=c2)
    area = 2 * mm(K2P.HEX_X_UM) * 2 * mm(K2P.HEX_Y_UM) + 2 * (0.5 * 2 * mm(K2P.HEX_X_UM) * mm(K2P.HEX_R_UM - K2P.HEX_Y_UM))
    rec("substrate", abs(sub.Volume() - (area - np.pi * 0.635 ** 2) * 1.575) < 1e-6 and abs(sub.BoundingBox().zmin + 1.575) < 1e-9, thickness=1.575)
    ab = parts["AIR"].BoundingBox()
    rec("air_box", len(parts["AIR"].Faces()) == 6 and abs(ab.zmin + 1.61) < 1e-9 and abs(ab.zmax - 20.0) < 1e-9, z=[ab.zmin, ab.zmax], faces=len(parts["AIR"].Faces()))
    pin = [s for s in solids if abs(s.BoundingBox().zmax - 6.0) < 1e-9]
    pb = pin[0].BoundingBox() if len(pin) == 1 else None
    port = parts["PORT"]
    pbx = port.BoundingBox()
    rec("probe_and_port", pb is not None and abs(pb.xlen - 1.27) < 1e-9 and abs(pb.zmin + 1.575) < 1e-9 and abs(pbx.xlen - 4.1) < 1e-9 and abs(port.Area() - np.pi * (2.05 ** 2 - 0.635 ** 2)) < 1e-6, pin_r=0.635, pin_top_above_substrate=6.0, port_r=[0.635, 2.05])
    byname = {n: s for n, s in parts.items() if n in exp_sh}
    bad_sh = 0
    for n, face, grp in sheets_expected:
        bb = byname[n].BoundingBox()
        a = n.split("_")
        ix, iy = int(a[1].replace("p", "+").replace("m", "-")), int(a[2].replace("p", "+").replace("m", "-"))
        cx, cy = mm(ix * PITCH_UM), mm(iy * PITCH_UM)
        ex = (cx - 0.2, cx + 0.2, cy + 0.3, cy + 1.3) if n.endswith("DX") else (cx + 0.3, cx + 1.3, cy - 0.2, cy + 0.2)
        bad_sh += not (abs(bb.zmin - 0.0175) < 1e-9 and abs(bb.xmin - ex[0]) < 1e-9 and abs(bb.xmax - ex[1]) < 1e-9 and abs(bb.ymin - ex[2]) < 1e-9 and abs(bb.ymax - ex[3]) < 1e-9 and abs(byname[n].Area() - 0.4) < 1e-9)
    rec("sheet_geometry", bad_sh == 0, n_checked=len(sheets_expected), size="0.4 x 1.0", z=0.0175, bad=int(bad_sh))
    return res, ok_all


def main(design_file):
    cells, (gb, flip, d, m) = make_cells(design_file)
    n_gb = sum(1 for c in cells if c["kind"] == "GB")
    print("cells", len(cells), "GB", n_gb, "diodes", 2 * n_gb, "rotated tensor slots", sum(c.get("flip", 0) for c in cells), flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    parts, sheets = build(cells)
    asm = cq.Assembly()
    cols = {"METAL_PEC": (0.9, 0.7, 0.2), "SUBSTRATE": (0.25, 0.55, 0.25), "AIR": (0.7, 0.8, 1.0), "PORT": (0.9, 0.1, 0.9)}
    for nm, shp in parts.items():
        asm.add(shp, name=nm, color=cq.Color(*cols[nm], 1))
    for nm, f, g in sheets:
        asm.add(f, name=nm, color=cq.Color(0.1, 0.1, 0.9, 1))
    path = OUT / "board_CP_LP_hybrid.step"
    asm.export(str(path), exportType="STEP")
    raw = path.read_text(encoding="utf-8")
    it = intended(cells)
    snapped, stats = snap_step_text(raw, it, description="sidequest CP+LP hybrid board, mm, 0.001 mm grid")
    path.write_text(snapped, encoding="utf-8", newline="\n")
    bad = inexact_literals(path.read_text(encoding="utf-8"), it)
    vs = verify_snap(raw, snapped)
    res, ok_all = check(path, cells, sheets)
    exact_ok = len(bad) == 0 and vs.get("collapsed_edges", 1) == 0 and not vs.get("distinct_merge_violations")
    res["EXACT_DIMS"] = {"ok": bool(exact_ok), "inexact_literals": len(bad), "n_snapped": stats["n_snapped"], "max_snap_mm": stats["max_snap_mm"], "grid_mm": 0.001, "positions": "integer index x 6.4; tensor slot corners rounded to 1 um"}
    ok_all &= exact_ok
    size_mb = path.stat().st_size / 1e6
    (OUT / "geometry_checks.json").write_text(json.dumps({"label": "MODEL_ONLY geometry checks on the re-imported STEP", "all_pass": bool(ok_all), "step_mb": round(size_mb, 1), "checks": res}, indent=1, default=float), encoding="utf-8", newline="\n")
    with open(OUT / "cells.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["ix", "iy", "x_mm", "y_mm", "kind", "group", "slot_L_mm", "slot_W_mm", "slot_angle_deg", "rotated_90_for_beam2"])
        for c in cells:
            w.writerow([c["ix"], c["iy"], f"{mm(c['ix'] * PITCH_UM):.3f}", f"{mm(c['iy'] * PITCH_UM):.3f}", "GB_switched_crossed_slit" if c["kind"] == "GB" else "tensor_rotating_slot", c["grp"], f"{c['L']:.3f}", f"{c['W']:.3f}",
                        "" if c["kind"] == "GB" else f"{math.degrees(c['ang']):.4f}", "" if c["kind"] == "GB" else c.get("flip", 0)])
    with open(OUT / "sheets.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["sheet", "group", "cx_mm", "cy_mm", "z_mm", "size_x_mm", "size_y_mm", "current_line_start_x", "current_line_start_y", "current_line_end_x", "current_line_end_y"])
        for n, f, g in sheets:
            bb = f.BoundingBox()
            cx, cy = (bb.xmin + bb.xmax) / 2, (bb.ymin + bb.ymax) / 2
            a, b = ((bb.xmin, cy), (bb.xmax, cy)) if n.endswith("DX") else ((cx, bb.ymin), (cx, bb.ymax))
            w.writerow([n, g, f"{cx:.3f}", f"{cy:.3f}", "0.0175", f"{bb.xlen:.3f}", f"{bb.ylen:.3f}", f"{a[0]:.3f}", f"{a[1]:.3f}", f"{b[0]:.3f}", f"{b[1]:.3f}"])
    with open(OUT / "states.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["state", "S_w_ohm_per_sq_on_all_GB_sheets", "diodes", "expected"])
        w.writerow(["GB_OFF (default)", "1e6", "OFF", "beam 1 RHCP (25, 0) only"])
        w.writerow(["GB_ON", "6.1", "ON", "beam 1 RHCP (25, 0) + beam 2 LP-V (25, 180)"])
    return cells, ok_all, res, size_mb, (gb, flip, d, m)


if __name__ == "__main__":
    df = Path(sys.argv[1]) if len(sys.argv) > 1 else RES / "cp_hybrid_design_n180_lb2.0.npz"
    cells, ok_all, res, size_mb, _ = main(df)
    print(json.dumps({k: v["ok"] for k, v in res.items()}), "ALL", ok_all, "MB", round(size_mb, 1))
