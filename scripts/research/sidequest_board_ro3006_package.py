"""HFSS board package, RO3006 1.27 mm board (cell r 2.22): ONE diode group (`D_*`), default (unpowered, diodes OFF) = beam 1 (25, 0), circuit ON = beam 2 (25, 180), LP-V, 10.5 GHz, 150 mm hexagon, <= 514 diodes, one control line.  MODEL_ONLY.
State 1: NO boundary on the D_* sheets (OFF = open). State 2: ONE Impedance boundary on all D_* sheets, Rs = R_on * W / L = 2.44 * 1.0 / 0.4 = 6.1 ohm/sq (W/L = 2.5).
Probe: pin r 0.635 mm, PinLength 5.0 mm above the substrate (starting point, tuned in HFSS 4.5-6.5 mm), ground opening r 2.05; launch pad: floating printed annulus r 1.3 .. 3.0 mm on the top layer (no extra part; UNVALIDATED).
    python sidequest_board_ro3006_package.py [design.json] [outdir] [syncname]    (design.json from sidequest_ro3006_design.py final)"""
import csv
import json
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

from sidequest_board_k2_lp_package import read_names
from src.transmitarray.export import inexact_literals, snap_step_text, verify_snap

DESIGN = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results/sidequest_hybrid/ro3006_board_design_r2.22.json"
OUTNAME = sys.argv[2] if len(sys.argv) > 2 else "sidequest_board_RO3006"
SYNCNAME = sys.argv[3] if len(sys.argv) > 3 else "board_RO3006"
OUT = ROOT / "hfss_exports" / OUTNAME
SYNC = Path(r"C:\hmsa_sync\data\sidequest") / SYNCNAME
R_UM = 150000
PITCH_UM, HSUB_UM, T_UM, GT_UM = 6400, 1270, 35, 35
R_DISC_UM, G_UM, WD_UM, POS_UM = 2220, 400, 1000, 800
PIN_R_UM, PORT_R_UM, PIN_TOP_UM = 635, 2050, 5000
PAD_XY_UM, AIR_TOP_UM = 10000, 20000
PADIN_UM, PADOUT_UM = 1300, 3000
HEX_X_UM = int(round(R_UM * np.sqrt(3) / 2))
HEX_Y_UM = R_UM // 2
F_GHZ, R_ON, ER, TAND = 10.5, 2.44, 6.15, 0.002
RS = R_ON * WD_UM / G_UM
LS_UM = {0.0: 0, 0.357: 721, 0.643: 1300, 0.786: 1588, 0.929: 1877, 1.0: 2020}   # slit half-length (rim 0.2 at 1.0); MoM used (r - 0.2) x fraction, differences <= 0.5 um
mm = lambda u: u / 1000.0  # noqa: E731


def box(cx, cy, lx, ly, z0, z1):
    return cq.Workplane("XY").workplane(offset=z0).center(cx, cy).rect(lx, ly).extrude(z1 - z0)


def hex_wp(z0, z1):
    X, Y, R = mm(HEX_X_UM), mm(HEX_Y_UM), mm(R_UM)
    return cq.Workplane("XY").workplane(offset=z0).polyline([(0.0, R), (-X, Y), (-X, -Y), (0.0, -R), (X, -Y), (X, Y)]).close().extrude(z1 - z0)


def sheet(cx, cy, lx, ly, z):
    p = [cq.Vector(cx - lx / 2, cy - ly / 2, z), cq.Vector(cx + lx / 2, cy - ly / 2, z), cq.Vector(cx + lx / 2, cy + ly / 2, z), cq.Vector(cx - lx / 2, cy + ly / 2, z)]
    return cq.Face.makeFromWires(cq.Wire.makePolygon(p, close=True))


def template(kind, fr):
    t, r, g2, pos, wd = mm(T_UM), mm(R_DISC_UM), mm(G_UM), mm(POS_UM), mm(WD_UM)
    disc = cq.Workplane("XY").circle(r).extrude(t)
    if kind == "none":
        return disc.val()
    h = mm(LS_UM[fr])
    cross = box(0, 0, g2, 2 * h, 0, t).union(box(0, 0, 2 * h, g2, 0, t))
    if kind == "strap":
        cross = cross.cut(box(0, pos, g2, wd, 0, t).union(box(pos, 0, wd, g2, 0, t)))
    return disc.cut(cross).val()


def build(cells):
    tmpl = {k: template(*k) for k in {(c[3], c[4]) for c in cells}}
    H, GT = mm(HSUB_UM), mm(GT_UM)
    pad = cq.Workplane("XY").circle(mm(PADOUT_UM)).circle(mm(PADIN_UM)).extrude(mm(T_UM)).val()
    metal = [hex_wp(-H - GT, -H).cut(cq.Workplane("XY").workplane(offset=-H - GT).circle(mm(PORT_R_UM)).extrude(GT)).val(),
             cq.Workplane("XY").workplane(offset=-H).circle(mm(PIN_R_UM)).extrude(H + mm(PIN_TOP_UM)).val(), pad]
    sheets = []
    for ix, iy, kind, tk, fr in cells:
        cx, cy = mm(ix * PITCH_UM), mm(iy * PITCH_UM)
        metal.append(tmpl[(tk, fr)].moved(cq.Location(cq.Vector(cx, cy, 0.0))))
        if kind == "SW":
            tag = f"{ix:+03d}_{iy:+03d}".replace("+", "p").replace("-", "m")
            sheets.append((f"D_{tag}_DX", sheet(cx, cy + mm(POS_UM), mm(G_UM), mm(WD_UM), mm(T_UM) / 2)))
            sheets.append((f"D_{tag}_DY", sheet(cx + mm(POS_UM), cy, mm(WD_UM), mm(G_UM), mm(T_UM) / 2)))
    sub = hex_wp(-H, 0.0).cut(cq.Workplane("XY").workplane(offset=-H).circle(mm(PIN_R_UM)).extrude(H)).val()
    air = box(0, 0, 2 * mm(HEX_X_UM + PAD_XY_UM), 2 * mm(R_UM + PAD_XY_UM), -H - GT, mm(AIR_TOP_UM)).val()
    port = cq.Face.makeFromWires(cq.Wire.makeCircle(mm(PORT_R_UM), cq.Vector(0, 0, -H), cq.Vector(0, 0, 1)), [cq.Wire.makeCircle(mm(PIN_R_UM), cq.Vector(0, 0, -H), cq.Vector(0, 0, 1))])
    return {"METAL_PEC": cq.Compound.makeCompound(metal), "SUBSTRATE": sub, "AIR": air, "PORT": port}, sheets, tmpl


def intended(cells):
    offs = {0, 200, 300, 1300, R_DISC_UM, PADIN_UM, PADOUT_UM} | set(LS_UM.values())
    offs = sorted(offs | {-o for o in offs})
    xs, ys = set(), set()
    for ix, iy, _k, _t, _f in cells:
        for o in offs:
            xs.add(ix * PITCH_UM + o)
            ys.add(iy * PITCH_UM + o)
    glob = [0, PIN_R_UM, PORT_R_UM, PADIN_UM, PADOUT_UM, HEX_X_UM, R_UM, HEX_Y_UM, HEX_X_UM + PAD_XY_UM, R_UM + PAD_XY_UM]
    glob += [-v for v in glob]
    xy = sorted({mm(v) for v in xs | ys | set(glob)})
    zs = sorted({mm(v) for v in (0, T_UM, T_UM / 2.0, -HSUB_UM, -HSUB_UM - GT_UM, PIN_TOP_UM, AIR_TOP_UM)})
    return dict(xy=np.array(xy), z=np.array(zs), dir=np.array([-1.0, 0.0, 1.0]))


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
    rec("metal_single_part_lumps", len(solids) == 3 + len(cells), lumps=len(solids), expected=3 + len(cells), note="ONE named part; lumps are DC-isolated (patches, launch pad, ground, probe)")
    vol = {k: v.Volume() for k, v in tmpl.items()}
    found, bad = {}, 0
    for s in solids:
        bb = s.BoundingBox()
        if abs(bb.zmin) < 1e-6 and abs(bb.zmax - mm(T_UM)) < 1e-6 and abs(bb.xlen - 2 * mm(R_DISC_UM)) < 1e-6:
            cx, cy = (bb.xmin + bb.xmax) / 2, (bb.ymin + bb.ymax) / 2
            ix, iy = round(cx / 6.4), round(cy / 6.4)
            bad += not (abs(cx - mm(ix * PITCH_UM)) < 1e-6 and abs(cy - mm(iy * PITCH_UM)) < 1e-6 and abs(bb.ylen - 2 * mm(R_DISC_UM)) < 1e-6)
            found[(ix, iy)] = s.Volume()
    mism = sum(1 for ix, iy, _k, tk, fr in cells if (ix, iy) not in found or abs(found[(ix, iy)] - vol[(tk, fr)]) > 1e-6)
    rec("patch_vs_design_map", bad == 0 and mism == 0 and len(found) == len(cells), n_patches=len(found), off_grid=bad, map_mismatches=mism, template_types=len(tmpl))
    off_t = tmpl.get(("open", 1.0)) or template("open", 1.0)
    hh = mm(LS_UM[1.0])
    rec("rim_template_dims", abs(off_t.Volume() - (np.pi * mm(R_DISC_UM) ** 2 - (2 * (2 * hh) * mm(G_UM) - mm(G_UM) ** 2)) * mm(T_UM)) < 1e-6, disc_d=2 * mm(R_DISC_UM), slit_L=2 * hh, slit_W=mm(G_UM), rim=mm(R_DISC_UM) - hh)
    pad = [s for s in solids if abs(s.BoundingBox().xlen - 2 * mm(PADOUT_UM)) < 1e-6 and abs(s.BoundingBox().zmin) < 1e-6]
    rec("launch_pad", len(pad) == 1 and abs(pad[0].Volume() - np.pi * (mm(PADOUT_UM) ** 2 - mm(PADIN_UM) ** 2) * mm(T_UM)) < 1e-6, inner_r=mm(PADIN_UM), outer_r=mm(PADOUT_UM), note="floating annulus on the top layer, UNVALIDATED")
    sub = parts["SUBSTRATE"]
    gp = cq.Compound.makeCompound([s for s in solids if s.BoundingBox().zmin < -1e-6])
    pat = cq.Compound.makeCompound([s for s in solids if s.BoundingBox().zmin > -1e-6 and s.BoundingBox().zmax < 0.04])
    c1, c2 = sub.intersect(gp).Volume(), sub.intersect(pat).Volume()
    rec("no_intersections", c1 < 1e-9 and c2 < 1e-9, substrate_x_ground_probe=c1, substrate_x_patches_pad=c2)
    X, Y, Rr = mm(HEX_X_UM), mm(HEX_Y_UM), mm(R_UM)
    area = 2 * X * 2 * Y + 2 * (0.5 * 2 * X * (Rr - Y))
    hs = mm(HSUB_UM)
    rec("substrate", abs(sub.Volume() - (area - np.pi * 0.635**2) * hs) < 1e-6 and abs(sub.BoundingBox().zmin + hs) < 1e-9 and abs(sub.BoundingBox().zmax) < 1e-9, thickness=hs, circumradius_mm=Rr, volume_mm3=sub.Volume())
    sb = sub.BoundingBox()
    rec("outline_circumradius", abs(sb.ymax - Rr) < 1e-9 and abs(sb.ymin + Rr) < 1e-9 and abs(sb.xmax - X) < 1e-9, y_extent_mm=[sb.ymin, sb.ymax], x_extent_mm=[sb.xmin, sb.xmax])
    poly = np.array([(0.0, Rr), (-X, Y), (-X, -Y), (0.0, -Rr), (X, -Y), (X, Y)])
    dmin = 9.0
    for ix, iy, *_r in cells:
        px, py = ix * 6.4, iy * 6.4
        d = min(((b[0] - a[0]) * (py - a[1]) - (b[1] - a[1]) * (px - a[0])) / np.hypot(*(b - a)) for a, b in zip(poly, np.roll(poly, -1, axis=0), strict=True))
        dmin = min(dmin, d - mm(R_DISC_UM))
    rec("discs_inside_outline", dmin >= 0.19, min_clearance_disc_to_edge_mm=round(float(dmin), 3))
    air = parts["AIR"]
    ab = air.BoundingBox()
    top = mm(AIR_TOP_UM)
    rec("air_box_ground_to_lid_single_face_walls", len(air.Faces()) == 6 and abs(ab.zmin + hs + 0.035) < 1e-9 and abs(ab.zmax - top) < 1e-9 and ab.xmax >= X + 9.99 and ab.ymax >= Rr + 9.99 and ab.zmax - mm(PIN_TOP_UM) >= 7.14,
        faces=len(air.Faces()), z=[ab.zmin, ab.zmax], top_above_probe_mm=ab.zmax - mm(PIN_TOP_UM), lambda_over_4_mm=round(299.792458 / F_GHZ / 4, 2), lateral_pad_mm=round(ab.xmax - X, 3))
    pin = [s for s in solids if abs(s.BoundingBox().zmax - mm(PIN_TOP_UM)) < 1e-9]
    pb = pin[0].BoundingBox() if len(pin) == 1 else None
    port = parts["PORT"]
    pbx = port.BoundingBox()
    rec("probe_and_port", pb is not None and abs(pb.xlen - 1.27) < 1e-9 and abs(pb.zmin + hs) < 1e-9 and abs(pbx.xlen - 4.1) < 1e-9 and abs(pbx.zmin + hs) < 1e-9 and abs(pbx.zmax + hs) < 1e-9 and abs(port.Area() - np.pi * (2.05**2 - 0.635**2)) < 1e-6,
        pin_r=0.635, PinLength_above_substrate=mm(PIN_TOP_UM), port_r=[0.635, 2.05], port_z=-hs)
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
        gaps.add(round(bb.xlen if dx else bb.ylen, 6))
        widths.add(round(bb.ylen if dx else bb.xlen, 6))
        per.setdefault((ix, iy), []).append(bb)
    for a, b in per.values():
        clear = min(clear, float(np.hypot(max(a.xmin - b.xmax, b.xmin - a.xmax, 0), max(a.ymin - b.ymax, b.ymin - a.ymax, 0))))
    sw0 = next(c for c in cells if c[2] == "SW" and c[4] == 1.0)
    patch = [s for s in solids if abs((s.BoundingBox().xmin + s.BoundingBox().xmax) / 2 - mm(sw0[0] * PITCH_UM)) < 1e-6 and abs((s.BoundingBox().ymin + s.BoundingBox().ymax) / 2 - mm(sw0[1] * PITCH_UM)) < 1e-6 and s.BoundingBox().zmin > -1e-6][0]
    cx0, cy0 = mm(sw0[0] * PITCH_UM), mm(sw0[1] * PITCH_UM)
    corners = [cq.Vertex.makeVertex(cx0 + sx * 0.2, cy0 + y, 0.0175) for sx in (-1, 1) for y in (0.3, 1.3)] + [cq.Vertex.makeVertex(cx0 + x, cy0 + sy * 0.2, 0.0175) for sy in (-1, 1) for x in (0.3, 1.3)]
    dist = max(patch.distance(v) for v in corners)
    rec("sheet_geometry", bad_sh == 0 and clear >= 0.1 and dist < 1e-6 and gaps == {0.4} and widths == {1.0}, n_checked=len(sheets), size="0.4 (gap, across the slit) x 1.0 (width, along it)", z=0.0175, bad=bad_sh, DX_DY_clearance_min_mm=round(clear, 4),
        sheet_corner_to_patch_max_mm=float(dist), gap_mm=sorted(gaps), width_mm=sorted(widths))
    rec("Rs_from_the_real_sheet_size", abs(R_ON * max(widths) / max(gaps) - RS) < 1e-12, R_on_ohm=R_ON, width_mm=max(widths), gap_mm=max(gaps), Rs_ohm_per_sq=R_ON * max(widths) / max(gaps), W_over_L=max(widths) / max(gaps), formula="Rs = R_on * W / L")
    return res, ok_all


def main():
    ev = json.loads(DESIGN.read_text())
    cells = [tuple(c) for c in ev["cells"]]
    n_sw = sum(1 for c in cells if c[2] == "SW")
    print("cells", len(cells), "switched", n_sw, "diodes", 2 * n_sw, flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    parts, sheets, tmpl = build(cells)
    asm = cq.Assembly()
    cols = {"METAL_PEC": (0.9, 0.7, 0.2), "SUBSTRATE": (0.25, 0.55, 0.25), "AIR": (0.7, 0.8, 1.0), "PORT": (0.9, 0.1, 0.9)}
    for nm, shp in parts.items():
        asm.add(shp, name=nm, color=cq.Color(*cols[nm], 1))
    for nm, f in sheets:
        asm.add(f, name=nm, color=cq.Color(0.9, 0.1, 0.1, 1))
    path = OUT / "board_RO3006.step"
    asm.export(str(path), exportType="STEP")
    raw = path.read_text(encoding="utf-8")
    snapped, stats = snap_step_text(raw, intended(cells), description="sidequest RO3006 K2 board, mm, 0.001 mm grid")
    path.write_text(snapped, encoding="utf-8", newline="\n")
    bad = inexact_literals(snapped, intended(cells))
    vs = verify_snap(raw, snapped)
    res, ok_all = check(path, cells, sheets, tmpl)
    exact_ok = len(bad) == 0 and vs.get("collapsed_edges", 1) == 0 and not vs.get("distinct_merge_violations")
    res["EXACT_DIMS"] = {"ok": bool(exact_ok), "inexact_literals": len(bad), "n_snapped": stats["n_snapped"], "max_snap_mm": stats["max_snap_mm"], "grid_mm": 0.001, "positions": "integer index x 6.4 (um integers)",
                         "outline_vertices_mm": f"(0, +-{mm(R_UM)}), (+-{mm(HEX_X_UM)}, +-{mm(HEX_Y_UM)}); circumradius exact, x rounded to 1 um"}
    ok_all &= exact_ok
    size_mb = path.stat().st_size / 1e6
    (OUT / "geometry_checks.json").write_text(json.dumps({"label": "MODEL_ONLY geometry checks on the re-imported STEP", "all_pass": bool(ok_all), "step_mb": round(size_mb, 1), "checks": res}, indent=1, default=float), encoding="utf-8", newline="\n")
    print({k: v["ok"] for k, v in res.items()}, "ALL", ok_all, "MB", round(size_mb, 1), flush=True)
    kinds = {}
    for _x, _y, kd, tk, fr in cells:
        kinds[(kd, tk, fr)] = kinds.get((kd, tk, fr), 0) + 1
    with open(OUT / "cells.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["ix", "iy", "x_mm", "y_mm", "kind", "template", "slit_half_mm"])
        for ix, iy, kd, tk, fr in cells:
            w.writerow([ix, iy, f"{mm(ix * PITCH_UM):.3f}", f"{mm(iy * PITCH_UM):.3f}", "switched" if kd == "SW" else "static", tk, mm(LS_UM[fr]) if tk != "none" else 0.0])
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
        w.writerow(["1", "default (unpowered)", "NO boundary (open)", "OFF", "25", "0", "LP phi-hat (LP-V)"])
        w.writerow(["2", "circuit ON", f"one Impedance boundary on all D_*: Rs = {RS} ohm/sq (R_on 2.44 ohm x W/L = 1.0 / 0.4), Reactance 0", "ON", "25", "180", "LP phi-hat (LP-V)"])
    names = [n for n, _f in sheets]
    (OUT / "assign_impedance_D.py").write_text(f'''"""AEDT (IronPython) helper - UNTESTED draft: ONE Impedance boundary (Rs = {RS} ohm/sq, reactance 0) on all D_* sheets = state 2. State 1 = delete this boundary (no boundary on the sheets).
Run in HFSS: Tools > Run Script. Check the boundary after the run (Objects list, Rs)."""
oDesign = oProject.GetActiveDesign()
oMod = oDesign.GetModule("BoundarySetup")
names = {json.dumps(names)}
oMod.AssignImpedance(["NAME:Imp_D", "Objects:=", names, "Resistance:=", "{RS}", "Reactance:=", "0", "InfiniteGround:=", False])
''', encoding="utf-8", newline="\n")
    return cells, ev, ok_all, res, size_mb, kinds, n_sw


def docs(cells, ev, ok_all, res, size_mb, kinds, n_sw):
    s = ev["states"]
    n_cells = len(cells)
    n_sheets = 2 * n_sw
    # S11 estimate of the probe (thin-pin FD model, sidequest_probe_trim.py; MODEL_ONLY, X offset a few ohm unknown)
    s11 = {"model": "thin-pin FD, PinLength 5.0 mm, bare launch geometry (cells from 6.5 mm): R 24 ohm, X +5 ohm -> S11 -9 dB; with the launch pad (rho_in 2.6 mm, UNVALIDATED): R 43 ohm, X +9 ohm -> S11 -18 dB",
           "acceptance": "S11 <= -10 dB (goal -15 dB), tuned in HFSS with PinLength 4.5-6.5 mm", "tuning_direction": "+1 mm PinLength: R +11 ohm, X +25 ohm; pad radius / PinLength / pin diameter are the HFSS variables",
           "s11_estimate_db": [-9.1, -9.1], "s11_estimate_with_pad_db": [-18.2, -18.2], "note": "the hologram's loading of the probe is NOT modelled; the model's baseline probe gave -12.7 dB against HFSS -15 dB"}
    states = []
    for k, nm in enumerate(("1", "2")):
        st = s[k]
        states.append({"state": nm, "name": "default (unpowered, diodes OFF)" if nm == "1" else "circuit ON (diodes ON)", "D_sheets": "no boundary" if nm == "1" else f"one Impedance boundary, Rs = {RS} ohm/sq (W/L = 2.5)",
                       "target_theta_deg": 25.0, "target_phi_deg": 0.0 if nm == "1" else 180.0, "peak_dir_theta_phi_deg": st["peak_dir"], "gain_at_design_direction_dbi": round(st["gain_design_dir"], 2), "peak_gain_dbi": round(st["gain_peak"], 2),
                       "sll_aperture_db": round(st["SLL_aperture_db"], 2), "sll_nonbeam_floor_db": round(st["SLL_floor_db"], 2), "sll_total_db": round(st["SLL_total_db"], 2),
                       "other_beam_level_dbi": round(st["other_beam_dbi"], 2), "isolation_db": round(st["isolation_db"], 2), "xpd_db": round(st["XPD_db"], 1), "polarisation": "LP phi-hat (LP-V)",
                       "s11_db": s11["s11_estimate_db"][k], "s11_note": s11["model"]})
    pred = {"labels": ["MODEL_ONLY", "cell: RO3006 r 2.22, MoM mesh 96 with the HFSS calibration (ON fixed, OFF contrast x0.75 from the r 2.20 cell check: ON 10.3617 / OFF 10.0193 GHz at 90 deg); the cell at r 2.22 and the board are NOT HFSS-validated",
                       "launch 0.71 (bare slab, spectral model; the pad would give up to 0.96, UNVALIDATED); gain offsets from board_K2_default HFSS: state 1 -1.78 dB (9.1 dBi), state 2 -2.7 dB (8.3 dBi)"],
            "design": {"id": "RO3006 r 2.22, design (i) one group D_*", "frequency_ghz": F_GHZ, "hexagon_circumradius_mm": R_UM / 1000, "cells": n_cells, "switched_cells": n_sw, "diode_sheets": n_sheets, "diodes": 2 * n_sw,
                       "cells_by_type": {"/".join(map(str, k)): v for k, v in kinds.items()}, "A0_np_per_m": ev["A0"], "beta_over_k0": ev["beta_k0"], "config": ev["cfg"], "lam": ev["lam"], "SW_to_beam_pct": round(100 * ev["conv"], 1), "SW_at_edge_pct": round(100 * ev["edge"], 1), "nonbeam_NB": round(ev["NB"], 3)},
            "states": states, "probe": s11,
            "sidelobe_floor": "SLL_total = max(aperture SLL, 10 log(NB x 7.45) - gain): NB = (1 - launch) + launch x SW_at_edge is the state-independent non-beam power (direct probe + edge); 7.45 = peaking factor from the HFSS sidelobes of board_K2_default (8.3 / 8.2 dBi in both states)",
            "falsifiers": "peak gain more than 4 dB below the prediction, peak more than 3 deg from (25, 0) / (25, 180), isolation below 8 dB, strongest sidelobe above -3 dB, or the two patterns not mirror images",
            "NA": ["board S11 with the hologram", "HFSS gain", "RO3006 r 2.22 cell eigenfrequencies in HFSS (the r 2.20 cell is checked)", "tolerance / pointing over corners (later work, paper notes)"]}
    (OUT / "predictions.json").write_text(json.dumps(pred, indent=1, default=float), encoding="utf-8", newline="\n")
    eqv = lambda c, sh: c + 0.25 * sh  # noqa: E731
    r_t = eqv(n_cells, n_sheets) / eqv(1338, 514)
    cost = {"basis": "HFSS actual for board_K2_default (RT5880, 1338 cells + 514 sheets): 1.28 M tets, 53 GB, 18 min adaptive; scaled by the cell count and x1.2-2.0 for the eps_r 6.15 substrate (ESTIMATE)",
            "tets_M": [round(1.28 * r_t * 1.2, 2), round(1.28 * r_t * 2.0, 2)], "ram_gb": [round(53 * r_t * 1.2), round(53 * r_t * 2.0)], "time_state1_min": [round(18 * r_t * 1.3), round(18 * r_t * 3.0)], "fits_256gb": True, "cells": n_cells, "sheets": n_sheets}
    (OUT / "hfss_cost.json").write_text(json.dumps(cost, indent=1), encoding="utf-8", newline="\n")
    setup = f"""# HFSS_SETUP - RO3006 K2 board: default = beam 1 (25, 0), circuit ON = beam 2 (25, 180), ONE diode group `D_*`, LP-V, exactly {F_GHZ} GHz, {R_UM/1000:g} mm hexagon circumradius. MODEL_ONLY. One page.
**1. Assign material `pec` to `METAL_PEC` (patches, launch pad, ground and probe: ONE named part, {res['metal_single_part_lumps']['lumps']} DC-isolated lumps; the STEP imports as vacuum).** `SUBSTRATE`: RO3006 eps_r {ER}, tan d {TAND}, z -{mm(HSUB_UM):g}..0. `AIR`: vacuum, first in the tree (overlaps the solids on purpose).
**Files** (mm): `board_RO3006.step` ({round(size_mb,1)} MB: `METAL_PEC`, `SUBSTRATE`, `AIR`, `PORT`, {n_sheets} sheets `D_<ix>_<iy>_DX|DY`), `states.csv`, `sheets.csv`, `cells.csv`, `assign_impedance_D.py` (AEDT helper, untested), `predictions.json`, `hfss_cost.json`, `geometry_checks.json` (all PASS on the re-imported STEP, EXACT_DIMS PASS).
**Geometry**: substrate hexagon circumradius {R_UM/1000:g} (vertices (0, +-{R_UM/1000:g}), (+-{mm(HEX_X_UM)}, +-{mm(HEX_Y_UM):g}), X rounded to 1 um), {mm(HSUB_UM):g} thick, probe hole r 0.635; ground 35 um (z {-mm(HSUB_UM)-0.035:.3f}..{-mm(HSUB_UM):g}) with an opening r 2.05; {n_cells} cells on the 6.4 lattice (integer index x 6.4), discs r {mm(R_DISC_UM):g} (35 um) >= 0.2 mm inside the outline, no cell for rho < 6.5; probe r 0.635, z {-mm(HSUB_UM):g}..+{mm(PIN_TOP_UM):g} (**PinLength {mm(PIN_TOP_UM):g} mm above the substrate: starting point, tune 4.5-6.5 mm**); `PORT` annulus r 0.635..2.05 at z {-mm(HSUB_UM):g}.
**Launch pad (no extra part, UNVALIDATED)**: floating printed annulus r {mm(PADIN_UM):g} .. {mm(PADOUT_UM):g} mm on the top layer around the pin (part of `METAL_PEC`); delete or resize it in HFSS if it hurts S11.
**Cells**: {n_sw} switched cells (rim disc + two sheets DX at (0, +0.8), DY at (+0.8, 0), 0.4 x 1.0 at z 0.0175, ALL in the one group `D_*`; slit half-length 2.02 / 1.877 / 1.588 / 1.300 mm: the switched-pair family) and {n_cells - n_sw} analog static cells (disc with a cross slit of half-length 0 / 0.721 / 1.300 / 1.588 / 1.877 / 2.020 mm, or strapped footprints): `cells.csv`. Static cells carry no sheet.
**2. Boundaries**: Radiation on all six faces of `AIR` (z {-mm(HSUB_UM)-0.035:.3f}..+{mm(AIR_TOP_UM):g}, {mm(AIR_TOP_UM)-mm(PIN_TOP_UM):g} mm above the probe top, > lambda/4 = 7.14 mm; x/y pad 10 mm; one face per side). PML only if licence and RAM allow; say which was used.
**3. Port**: `PORT` -> Lumped Port 50 ohm, integration line (0.635, 0, {-mm(HSUB_UM):g}) -> (2.05, 0, {-mm(HSUB_UM):g}), renormalise 50, no de-embedding.
**4. The two states (same mesh, one project)**: **State 1 (default, diodes OFF): NO boundary on any `D_*` sheet.** **State 2 (circuit ON): ONE Impedance boundary on ALL `D_*` sheets, Rs = {RS} ohm/sq, reactance 0** (Rs = R_on x W / L = 2.44 ohm x 1.0 mm / 0.4 mm, **W/L = 2.5**, from the real sheet size).
**5. Solution**: ONE frequency **{F_GHZ} GHz**, Driven Modal, first-order basis, max passes 6, max delta S 0.02, no sweep. Mesh: substrate max length 0.6 mm (eps_r 6.15), patch faces 0.3 mm, `D_*` sheets 0.15 mm, probe and port 0.3 mm, lambda refinement on. Mesh reuse state 1 -> 2 allowed (say so).
**6. Export (per state)**: far-field infinite sphere theta 0..180, phi 0..360 step **2 deg**, **complex rETheta and rEPhi in V** (full sphere), one CSV per state (`rE_state1.csv`, `rE_state2.csv`); complex **S11** at {F_GHZ} GHz; antenna parameters; passes, delta S, tets, peak RAM, solve time per state.
**Probe tuning**: acceptance S11 <= -10 dB (goal -15). Model starting point PinLength 5.0 mm (S11 about -9 dB, R 24 ohm, X +5 ohm); **+1 mm PinLength = R +11 ohm, X +25 ohm**; the pad (r {mm(PADIN_UM):g}..{mm(PADOUT_UM):g}) is expected to raise R (model: R 43 ohm at 5 mm, S11 -18 dB, UNVALIDATED). Sweep PinLength 4.5-6.5 mm and the pin diameter.
**Model predictions (MODEL_ONLY; launch 0.71)**: state 1 gain at (25, 0) {states[0]['gain_at_design_direction_dbi']} dBi (peak {states[0]['peak_gain_dbi']} at {states[0]['peak_dir_theta_phi_deg']}), SLL aperture {states[0]['sll_aperture_db']} dB, non-beam floor {states[0]['sll_nonbeam_floor_db']} dB, **SLL total {states[0]['sll_total_db']} dB**, isolation {states[0]['isolation_db']} dB; state 2 gain at (25, 180) {states[1]['gain_at_design_direction_dbi']} dBi (peak {states[1]['peak_gain_dbi']} at {states[1]['peak_dir_theta_phi_deg']}), SLL aperture {states[1]['sll_aperture_db']} dB, floor {states[1]['sll_nonbeam_floor_db']} dB, **SLL total {states[1]['sll_total_db']} dB**, isolation {states[1]['isolation_db']} dB. The -13 dB SLL target is NOT reached (floor-limited, see `K2_RO3006_BOARD.md`). Falsifiers: `predictions.json`.
**Cost (ESTIMATE)**: {cost['tets_M'][0]}-{cost['tets_M'][1]} M tets, {cost['ram_gb'][0]}-{cost['ram_gb'][1]} GB, {cost['time_state1_min'][0]}-{cost['time_state1_min'][1]} min state 1 (reference: RT5880 board 1.28 M tets, 53 GB, 18 min).
"""
    (OUT / "HFSS_SETUP.md").write_text(setup, encoding="utf-8", newline="\n")
    SYNC.mkdir(parents=True, exist_ok=True)
    for f in OUT.iterdir():
        shutil.copy2(f, SYNC / f.name)
    print("synced", sorted(p.name for p in SYNC.iterdir()), flush=True)


if __name__ == "__main__":
    cells, ev, ok_all, res, size_mb, kinds, n_sw = main()
    if not ok_all:
        sys.exit("checks failed")
    docs(cells, ev, ok_all, res, size_mb, kinds, n_sw)
    (OUT / "design_summary.json").write_text(json.dumps({"design": ev["cfg"], "n_cells": len(cells), "switched_cells": n_sw, "diodes": 2 * n_sw, "cells_by_type": {"/".join(map(str, k)): v for k, v in kinds.items()}}, indent=1, default=float), encoding="utf-8", newline="\n")
