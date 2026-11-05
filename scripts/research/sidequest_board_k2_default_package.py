"""HFSS board package, design (i) of the K2-default board: ONE diode group (`D_*`), default (unpowered, diodes OFF) = beam 1 (25, 0), circuit ON = beam 2 (25, 180), linear polarisation, 10.5 GHz, MODEL_ONLY.
State 1: NO boundary on the D_* sheets (OFF = open). State 2: ONE Impedance boundary on all D_* sheets, Rs = R_on * width / gap = 2.44 * 1.0 / 0.4 = 6.1 ohm/sq (exact, from the real sheet size).
Cells (6.4 mm lattice, integer index, discs fully inside the outline): switched = rim disc (ls 2.8) + two sheets; analog statics = discs with a cross slit of half-length ls (0 = solid, 1.0, 1.8, 2.2, 2.6, 2.8) or strapped footprints (ls 2.2, 2.8).
    python sidequest_board_k2_default_package.py [R_um=150000] [outdir] [syncname]      (needs results/sidequest_hybrid/k2_default_lib.json from sidequest_k2_default_lib.py)"""
import csv
import json
import os
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(Path(__file__).resolve().parent))
os.environ.setdefault("K4_A0", "0.844")  # HFSS-calibrated scalar cell; MUST be set before sidequest_coded_k4 is imported
os.environ.setdefault("K4_BETA_K0", "1.0731")
warnings.simplefilter("ignore")
import cadquery as cq
import numpy as np

import sidequest_coded_k4 as K
import sidequest_k2_default_design as D
from sidequest_board_k2_lp_package import read_names
from src.transmitarray.export import inexact_literals, snap_step_text, verify_snap

R_UM = int(sys.argv[1]) if len(sys.argv) > 1 else 150000
OUTNAME = sys.argv[2] if len(sys.argv) > 2 else "sidequest_board_K2_default"
SYNCNAME = sys.argv[3] if len(sys.argv) > 3 else "board_K2_default"
OUT = ROOT / "hfss_exports" / OUTNAME
SYNC = Path(r"C:\hmsa_sync\data\sidequest") / SYNCNAME
PITCH_UM, HSUB_UM, T_UM, GT_UM = 6400, 1575, 35, 35
R_DISC_UM, G_UM, WD_UM, POS_UM = 3000, 400, 1000, 800
PIN_R_UM, PORT_R_UM, PIN_TOP_UM = 635, 2050, 6000
PAD_XY_UM, AIR_TOP_UM = 10000, 20000
HEX_X_UM = int(round(R_UM * np.sqrt(3) / 2)); HEX_Y_UM = R_UM // 2   # vertices (0, +-R), (+-X, +-R/2); X rounded to 1 um (regular to 1e-6)
F_GHZ, R_ON, C_FF = 10.5, 2.44, 28.2
RS = R_ON * WD_UM / G_UM
mm = lambda u: u / 1000.0  # noqa: E731  (the only division)
LS_UM = {0.0: 0, 1.0: 1000, 1.8: 1800, 2.2: 2200, 2.6: 2600, 2.8: 2800}


# ---------------------------------------------------------------------------------------------- design
def design():
    """Re-run the chosen configuration of design (i) (k2_default_designs.json) for this board size; returns cells (ix, iy, kind, key) and the per-state model numbers."""
    res = json.loads((ROOT / "results/sidequest_hybrid/k2_default_designs.json").read_text(encoding="utf-8"))
    best = res["rows"]["i"]["best"] if R_UM == 150000 else res["rows"]["i_R%d" % (R_UM // 1000)]["best"]
    stat, _pairs, _norm = D.library()
    levels = sorted({round(v[0], 4) for v in stat.values()})
    D.setup(F_GHZ * 1e9, R_UM / 1e6, D.MARGIN)
    bd = K.KBoard([0, 180], [best["pol"], best["pol"]])
    r = D.run_design(bd, D.options(levels, [(1.0, -1.0, 1.0)], [(0, 1)]), [0.0, np.radians(best["eps_B_deg"])], best["lam"], best["g"], mu=best["mu"])
    sm = D.summarise(r, "i", bd)
    sm.update({"pol": best["pol"], "eps_B_deg": best["eps_B_deg"]})
    assert abs(sm["gain_dbi"][0] - best["gain_dbi"][0]) < 0.01 and sm["diodes"] == best["diodes"], "re-run does not reproduce the stored design"
    ix = np.rint(bd.x / 6.4e-3).astype(int)
    iy = np.rint(bd.y / 6.4e-3).astype(int)
    assert np.allclose(ix * 6.4e-3, bd.x, atol=1e-9) and np.allclose(iy * 6.4e-3, bd.y, atol=1e-9)
    keys = list(stat.keys())
    vals = np.array([stat[k][0] for k in keys])
    cells = []
    for i in range(bd.n):
        if r["nd"][i] > 0:
            cells.append((int(ix[i]), int(iy[i]), "SW", ("open", 2.8)))
        else:
            kk = keys[int(np.argmin(np.abs(vals - r["m"][0, 0, i])))]
            cells.append((int(ix[i]), int(iy[i]), "ST", (kk[1], kk[0])))
    return cells, sm, r


# ---------------------------------------------------------------------------------------------- geometry
def box(cx, cy, lx, ly, z0, z1):
    return cq.Workplane("XY").workplane(offset=z0).center(cx, cy).rect(lx, ly).extrude(z1 - z0)


def hex_wp(z0, z1):
    X, Y, R = mm(HEX_X_UM), mm(HEX_Y_UM), mm(R_UM)
    return cq.Workplane("XY").workplane(offset=z0).polyline([(0.0, R), (-X, Y), (-X, -Y), (0.0, -R), (X, -Y), (X, Y)]).close().extrude(z1 - z0)


def sheet(cx, cy, lx, ly, z):
    p = [cq.Vector(cx - lx / 2, cy - ly / 2, z), cq.Vector(cx + lx / 2, cy - ly / 2, z), cq.Vector(cx + lx / 2, cy + ly / 2, z), cq.Vector(cx - lx / 2, cy + ly / 2, z)]
    return cq.Face.makeFromWires(cq.Wire.makePolygon(p, close=True))


def template(kind, ls):
    t, r, g2, pos, wd = mm(T_UM), mm(R_DISC_UM), mm(G_UM), mm(POS_UM), mm(WD_UM)
    disc = cq.Workplane("XY").circle(r).extrude(t)
    if LS_UM[ls] == 0:
        return disc.val()
    h = mm(LS_UM[ls])
    cross = box(0, 0, g2, 2 * h, 0, t).union(box(0, 0, 2 * h, g2, 0, t))
    if kind == "strap":
        cross = cross.cut(box(0, pos, g2, wd, 0, t).union(box(pos, 0, wd, g2, 0, t)))
    return disc.cut(cross).val()


def build(cells):
    tmpl = {k: template(*k) for k in {c[3] for c in cells}}
    H, GT = mm(HSUB_UM), mm(GT_UM)
    metal = [hex_wp(-H - GT, -H).cut(cq.Workplane("XY").workplane(offset=-H - GT).circle(mm(PORT_R_UM)).extrude(GT)).val(),
             cq.Workplane("XY").workplane(offset=-H).circle(mm(PIN_R_UM)).extrude(H + mm(PIN_TOP_UM)).val()]
    sheets = []
    for ix, iy, kind, key in cells:
        cx, cy = mm(ix * PITCH_UM), mm(iy * PITCH_UM)
        metal.append(tmpl[key].moved(cq.Location(cq.Vector(cx, cy, 0.0))))
        if kind == "SW":
            tag = f"{ix:+03d}_{iy:+03d}".replace("+", "p").replace("-", "m")
            sheets.append((f"D_{tag}_DX", sheet(cx, cy + mm(POS_UM), mm(G_UM), mm(WD_UM), mm(T_UM) / 2)))
            sheets.append((f"D_{tag}_DY", sheet(cx + mm(POS_UM), cy, mm(WD_UM), mm(G_UM), mm(T_UM) / 2)))
    sub = hex_wp(-H, 0.0).cut(cq.Workplane("XY").workplane(offset=-H).circle(mm(PIN_R_UM)).extrude(H)).val()
    air = box(0, 0, 2 * mm(HEX_X_UM + PAD_XY_UM), 2 * mm(R_UM + PAD_XY_UM), -H - GT, mm(AIR_TOP_UM)).val()
    port = cq.Face.makeFromWires(cq.Wire.makeCircle(mm(PORT_R_UM), cq.Vector(0, 0, -H), cq.Vector(0, 0, 1)), [cq.Wire.makeCircle(mm(PIN_R_UM), cq.Vector(0, 0, -H), cq.Vector(0, 0, 1))])
    return {"METAL_PEC": cq.Compound.makeCompound(metal), "SUBSTRATE": sub, "AIR": air, "PORT": port}, sheets, tmpl


def intended(cells):
    offs = {0, 200, 300, 1300, 3000} | set(LS_UM.values())
    offs = sorted(offs | {-o for o in offs})
    xs, ys = set(), set()
    for ix, iy, _k, _g in cells:
        for o in offs:
            xs.add(ix * PITCH_UM + o)
            ys.add(iy * PITCH_UM + o)
    glob = [0, PIN_R_UM, PORT_R_UM, HEX_X_UM, R_UM, HEX_Y_UM, HEX_X_UM + PAD_XY_UM, R_UM + PAD_XY_UM]
    glob += [-v for v in glob]
    xy = sorted({mm(v) for v in xs | ys | set(glob)})
    zs = sorted({mm(v) for v in (0, T_UM, T_UM / 2.0, -HSUB_UM, -HSUB_UM - GT_UM, PIN_TOP_UM, AIR_TOP_UM)})
    return dict(xy=np.array(xy), z=np.array(zs), dir=np.array([-1.0, 0.0, 1.0]))


# ---------------------------------------------------------------------------------------------- checks on the re-imported STEP
def check(path, cells, sheets, tmpl):
    parts = read_names(path)
    res, ok_all = {}, True

    def rec(name, ok, **kw):
        nonlocal ok_all
        ok_all &= bool(ok)
        res[name] = {"ok": bool(ok), **kw}

    exp_sh = {n for n, _f in sheets}
    sh = {n for n in parts if n.startswith("D_")}
    rec("bodies_and_names", {"METAL_PEC", "SUBSTRATE", "AIR", "PORT"} <= set(parts) and sh == exp_sh and len(parts) == 4 + len(exp_sh), n_objects=len(parts), n_sheets=len(sh), groups="all sheets D_* (one group)")
    solids = parts["METAL_PEC"].Solids()
    rec("metal_single_part_lumps", len(solids) == 2 + len(cells), lumps=len(solids), expected=2 + len(cells), note="ONE named part; lumps are DC-isolated (patches, ground, probe)")
    vol = {k: v.Volume() for k, v in tmpl.items()}
    found, bad = {}, 0
    for s in solids:
        bb = s.BoundingBox()
        if abs(bb.zmin) < 1e-6 and abs(bb.zmax - mm(T_UM)) < 1e-6:
            cx, cy = (bb.xmin + bb.xmax) / 2, (bb.ymin + bb.ymax) / 2
            ix, iy = round(cx / 6.4), round(cy / 6.4)
            bad += not (abs(cx - mm(ix * PITCH_UM)) < 1e-6 and abs(cy - mm(iy * PITCH_UM)) < 1e-6 and abs(bb.xlen - 6.0) < 1e-6 and abs(bb.ylen - 6.0) < 1e-6)
            found[(ix, iy)] = s.Volume()
    mism = sum(1 for ix, iy, _k, key in cells if (ix, iy) not in found or abs(found[(ix, iy)] - vol[key]) > 1e-6)
    rec("patch_vs_design_map", bad == 0 and mism == 0 and len(found) == len(cells), n_patches=len(found), off_grid=bad, map_mismatches=mism, template_types=len(tmpl))
    ls28 = tmpl[("open", 2.8)] if ("open", 2.8) in tmpl else template("open", 2.8)
    rec("rim_template_dims", abs(ls28.Volume() - (np.pi * 9.0 - (2 * 5.6 * 0.4 - 0.16)) * 0.035) < 1e-6, disc_d=6.0, slit_L=5.6, slit_W=0.4, rim=0.2)
    sub = parts["SUBSTRATE"]
    gp = cq.Compound.makeCompound([s for s in solids if s.BoundingBox().zmin < -1e-6])
    pat = cq.Compound.makeCompound([s for s in solids if s.BoundingBox().zmin > -1e-6 and s.BoundingBox().zmax < 0.04])
    c1, c2 = sub.intersect(gp).Volume(), sub.intersect(pat).Volume()
    rec("no_intersections", c1 < 1e-9 and c2 < 1e-9, substrate_x_ground_probe=c1, substrate_x_patches=c2)
    X, Y, Rr = mm(HEX_X_UM), mm(HEX_Y_UM), mm(R_UM)
    area = 2 * X * 2 * Y + 2 * (0.5 * 2 * X * (Rr - Y))
    rec("substrate", abs(sub.Volume() - (area - np.pi * 0.635**2) * 1.575) < 1e-6 and abs(sub.BoundingBox().zmin + 1.575) < 1e-9 and abs(sub.BoundingBox().zmax) < 1e-9, thickness=1.575, circumradius_mm=Rr, volume_mm3=sub.Volume())
    sb = sub.BoundingBox()
    rec("outline_circumradius", abs(sb.ymax - Rr) < 1e-9 and abs(sb.ymin + Rr) < 1e-9 and abs(sb.xmax - X) < 1e-9, y_extent_mm=[sb.ymin, sb.ymax], x_extent_mm=[sb.xmin, sb.xmax])
    edge_ok = True
    poly = D.hexpoly(Rr / 1000.0) * 1000.0
    dmin = 9.0
    for ix, iy, _k, _key in cells:
        px, py = ix * 6.4, iy * 6.4
        d = min(((b[0] - a[0]) * (py - a[1]) - (b[1] - a[1]) * (px - a[0])) / np.hypot(*(b - a)) for a, b in zip(poly, np.roll(poly, -1, axis=0), strict=True))
        dmin = min(dmin, d - 3.0)
    rec("discs_inside_outline", dmin >= 0.19, min_clearance_disc_to_edge_mm=round(float(dmin), 3))
    del edge_ok
    air = parts["AIR"]
    ab = air.BoundingBox()
    rec("air_box_ground_to_lid_single_face_walls", len(air.Faces()) == 6 and abs(ab.zmin + 1.61) < 1e-9 and abs(ab.zmax - 20.0) < 1e-9 and ab.xmax >= X + 9.99 and ab.ymax >= Rr + 9.99 and ab.zmax - 6.0 >= 7.14,
        faces=len(air.Faces()), z=[ab.zmin, ab.zmax], top_above_probe_mm=ab.zmax - 6.0, lambda_over_4_mm=round(299.792458 / F_GHZ / 4, 2), lateral_pad_mm=round(ab.xmax - X, 3))
    pin = [s for s in solids if abs(s.BoundingBox().zmax - 6.0) < 1e-9]
    pb = pin[0].BoundingBox() if len(pin) == 1 else None
    port = parts["PORT"]
    pbx = port.BoundingBox()
    rec("probe_and_port", pb is not None and abs(pb.xlen - 1.27) < 1e-9 and abs(pb.zmin + 1.575) < 1e-9 and abs(pbx.xlen - 4.1) < 1e-9 and abs(pbx.zmin + 1.575) < 1e-9 and abs(pbx.zmax + 1.575) < 1e-9 and abs(port.Area() - np.pi * (2.05**2 - 0.635**2)) < 1e-6,
        pin_r=0.635, pin_top_above_substrate=6.0, port_r=[0.635, 2.05], port_z=-1.575)
    byname = {n: parts[n] for n in exp_sh}
    bad_sh, gaps, widths, clear, per = 0, set(), set(), 9.0, {}
    for n, _f in sheets:
        s = byname[n]
        bb = s.BoundingBox()
        p_ = n.split("_")
        ix, iy = (int(p_[1].replace("p", "+").replace("m", "-")), int(p_[2].replace("p", "+").replace("m", "-")))
        cx, cy = mm(ix * PITCH_UM), mm(iy * PITCH_UM)
        dx = n.endswith("DX")
        ex = (cx - 0.2, cx + 0.2, cy + 0.3, cy + 1.3) if dx else (cx + 0.3, cx + 1.3, cy - 0.2, cy + 0.2)
        ok = abs(bb.zmin - 0.0175) < 1e-9 and abs(bb.zmax - 0.0175) < 1e-9 and abs(bb.xmin - ex[0]) < 1e-9 and abs(bb.xmax - ex[1]) < 1e-9 and abs(bb.ymin - ex[2]) < 1e-9 and abs(bb.ymax - ex[3]) < 1e-9 and abs(s.Area() - 0.4) < 1e-9
        bad_sh += not ok
        gaps.add(round(bb.xlen if dx else bb.ylen, 6)); widths.add(round(bb.ylen if dx else bb.xlen, 6))
        per.setdefault((ix, iy), []).append(bb)
    for a, b in per.values():
        clear = min(clear, float(np.hypot(max(a.xmin - b.xmax, b.xmin - a.xmax, 0), max(a.ymin - b.ymax, b.ymin - a.ymax, 0))))
    # sheet edges on the metal: distance from the sheet corners to the (re-imported) switched-template patch of a sample cell
    sw0 = next(c for c in cells if c[2] == "SW")
    patch = [s for s in solids if abs((s.BoundingBox().xmin + s.BoundingBox().xmax) / 2 - mm(sw0[0] * PITCH_UM)) < 1e-6 and abs((s.BoundingBox().ymin + s.BoundingBox().ymax) / 2 - mm(sw0[1] * PITCH_UM)) < 1e-6 and s.BoundingBox().zmin > -1e-6][0]
    cx0, cy0 = mm(sw0[0] * PITCH_UM), mm(sw0[1] * PITCH_UM)
    corners = [cq.Vertex.makeVertex(cx0 + sx * 0.2, cy0 + y, 0.0175) for sx in (-1, 1) for y in (0.3, 1.3)] + [cq.Vertex.makeVertex(cx0 + x, cy0 + sy * 0.2, 0.0175) for sy in (-1, 1) for x in (0.3, 1.3)]
    dist = max(patch.distance(v) for v in corners)
    rec("sheet_geometry", bad_sh == 0 and clear >= 0.1 and dist < 1e-6 and gaps == {0.4} and widths == {1.0}, n_checked=len(sheets), size="0.4 (gap, across the slit) x 1.0 (width, along it)", z=0.0175, bad=bad_sh, DX_DY_clearance_min_mm=round(clear, 4),
        sheet_corner_to_patch_max_mm=float(dist), gap_mm=sorted(gaps), width_mm=sorted(widths))
    rec("Rs_from_the_real_sheet_size", abs(R_ON * max(widths) / max(gaps) - RS) < 1e-12, R_on_ohm=R_ON, width_mm=max(widths), gap_mm=max(gaps), Rs_ohm_per_sq=R_ON * max(widths) / max(gaps), formula="Rs = R_on * width / gap")
    return res, ok_all


# ---------------------------------------------------------------------------------------------- package
def main():
    cells, sm, r = design()
    n_sw = sum(1 for c in cells if c[2] == "SW")
    print("R", R_UM, "cells", len(cells), "switched", n_sw, "diodes", 2 * n_sw, sm["gain_dbi"], sm["sll_db"], sm["isolation"], flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    parts, sheets, tmpl = build(cells)
    asm = cq.Assembly()
    cols = {"METAL_PEC": (0.9, 0.7, 0.2), "SUBSTRATE": (0.25, 0.55, 0.25), "AIR": (0.7, 0.8, 1.0), "PORT": (0.9, 0.1, 0.9)}
    for nm, shp in parts.items():
        asm.add(shp, name=nm, color=cq.Color(*cols[nm], 1))
    for nm, f in sheets:
        asm.add(f, name=nm, color=cq.Color(0.9, 0.1, 0.1, 1))
    path = OUT / "board_K2_default.step"
    asm.export(str(path), exportType="STEP")
    raw = path.read_text(encoding="utf-8")
    snapped, stats = snap_step_text(raw, intended(cells), description="sidequest K2 default board, mm, 0.001 mm grid")
    path.write_text(snapped, encoding="utf-8", newline="\n")
    bad = inexact_literals(snapped, intended(cells))
    vs = verify_snap(raw, snapped)
    res, ok_all = check(path, cells, sheets, tmpl)
    exact_ok = len(bad) == 0 and vs.get("collapsed_edges", 1) == 0 and not vs.get("distinct_merge_violations")
    res["EXACT_DIMS"] = {"ok": bool(exact_ok), "inexact_literals": len(bad), "n_snapped": stats["n_snapped"], "max_snap_mm": stats["max_snap_mm"], "grid_mm": 0.001, "positions": "integer index x 6.4 (um integers, one division)",
                         "outline_vertices_mm": f"(0, +-{mm(R_UM)}), (+-{mm(HEX_X_UM)}, +-{mm(HEX_Y_UM)}); circumradius exact, x rounded to 1 um (R*sqrt(3)/2 = {R_UM * np.sqrt(3) / 2000:.5f})"}
    ok_all &= exact_ok
    size_mb = path.stat().st_size / 1e6
    (OUT / "geometry_checks.json").write_text(json.dumps({"label": "MODEL_ONLY geometry checks on the re-imported STEP", "all_pass": bool(ok_all), "step_mb": round(size_mb, 1), "checks": res}, indent=1, default=float), encoding="utf-8", newline="\n")
    kinds = {}
    for _x, _y, kd, key in cells:
        kinds[(kd, *key)] = kinds.get((kd, *key), 0) + 1
    with open(OUT / "cells.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["ix", "iy", "x_mm", "y_mm", "kind", "template", "slit_half_mm"])
        for ix, iy, kd, key in cells:
            w.writerow([ix, iy, f"{mm(ix * PITCH_UM):.3f}", f"{mm(iy * PITCH_UM):.3f}", "switched" if kd == "SW" else "static", key[0], key[1]])
    with open(OUT / "sheets.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["sheet", "cx_mm", "cy_mm", "z_mm", "size_x_mm", "size_y_mm", "current_line_start_x", "current_line_start_y", "current_line_end_x", "current_line_end_y"])
        for n, f in sheets:
            bb = f.BoundingBox()
            cx, cy = (bb.xmin + bb.xmax) / 2, (bb.ymin + bb.ymax) / 2
            a, b = ((bb.xmin, cy), (bb.xmax, cy)) if n.endswith("DX") else ((cx, bb.ymin), (cx, bb.ymax))
            w.writerow([n, f"{cx:.3f}", f"{cy:.3f}", "0.0175", f"{bb.xlen:.3f}", f"{bb.ylen:.3f}", f"{a[0]:.3f}", f"{a[1]:.3f}", f"{b[0]:.3f}", f"{b[1]:.3f}"])
    with open(OUT / "states.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["state", "name", "D_sheets", "diodes", "target_theta_deg", "target_phi_deg", "polarisation"])
        w.writerow(["1", "default (unpowered)", "NO boundary (open)", "OFF", "25", "0", sm["pol"]])
        w.writerow(["2", "circuit ON", f"one Impedance boundary on all D_*: Rs = {RS} ohm/sq (R_on 2.44 ohm x 1.0 / 0.4), Reactance 0", "ON", "25", "180", sm["pol"]])
    return cells, sm, r, ok_all, res, size_mb, kinds


if __name__ == "__main__":
    cells, sm, r, ok_all, res, size_mb, kinds = main()
    (OUT / "design_summary.json").write_text(json.dumps({"summary": sm, "cells_by_type": {"/".join(map(str, k)): v for k, v in kinds.items()}, "R_um": R_UM, "n_cells": len(cells)}, indent=1, default=float), encoding="utf-8", newline="\n")
    print(json.dumps({k: v["ok"] for k, v in res.items()}), "ALL", ok_all, "MB", round(size_mb, 1))
