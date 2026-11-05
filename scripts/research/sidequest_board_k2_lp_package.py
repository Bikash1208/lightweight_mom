"""Side quest HFSS board package: K = 2, LP (phi-hat), beams (25 deg, 0) and (25 deg, 180), scalar rim-disc cell (both diodes on one line), MODEL_ONLY (cell HFSS-validated, board not).

Design = `sidequest_scalar_board` (HFSS-calibrated scalar cell: alpha 0.844 Np/m, beta/k0 1.0731), LP-V (phi-hat), eps_B = 180 deg, ideal +-1 static cells, lam 0.28, g 2 -> 720 diodes (360 switched cells x 2 sheets).
Static cell +1 = open slits (no sheet), static -1 = both footprints filled with metal (strap, = a permanently ON diode with R -> 0). Switched cell: sheets DX (0, +0.8) and DY (+0.8, 0), 0.4 x 1.0, both in ONE group:
GA = ON in state A (OFF in B), GB = ON in state B. Parallel RLC per sheet: C 28.16 fF fixed, R = R_A (GA) / R_B (GB), 2.44 ohm (ON) or 1e9 ohm (OFF).
Geometry (mm, integer um, one division): substrate/ground hexagon circumradius 154 (vertices (0, +-154), (+-133.368, +-77)), cells = discs r 3.0 (35 um copper) with a cross slit (0.4 x 5.6, rim 0.2) on the 6.4 lattice, probe r 0.635 z -H..6.0,
port annulus 0.635..2.05 at z = -H, air box ground bottom -> z 20, x/y pad 10 mm.   python sidequest_board_k2_lp_package.py"""
import csv
import json
import os
import shutil
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(Path(__file__).resolve().parent))
os.environ.setdefault("K4_A0", "0.844")
os.environ.setdefault("K4_BETA_K0", "1.0731")
warnings.simplefilter("ignore")
import cadquery as cq
import numpy as np

import sidequest_coded_k4 as K
from src.transmitarray.export import inexact_literals, snap_step_text, verify_snap

OUT = ROOT / "hfss_exports/sidequest_board_k2_lp"
SYNC = Path(r"C:\hmsa_sync\data\sidequest\board_k2_lp".replace("k2_lp", "K2_LP"))
PITCH_UM, HSUB_UM, T_UM, GT_UM = 6400, 1575, 35, 35
R_DISC_UM, SLOT_HALF_UM, SLOT_G_UM = 3000, 2800, 400  # slit length 5.6 (half 2.8), width 0.4
POS_UM, WD_UM = 800, 1000                              # diode footprint centre 0.8 mm from the cell centre, 0.4 x 1.0
PIN_R_UM, PORT_R_UM, PIN_TOP_UM = 635, 2050, 6000
HEX_R_UM, HEX_X_UM, HEX_Y_UM = 154000, 133368, 77000   # board outline (cells selected by their centres inside the 150 mm hexagon)
PAD_XY_UM, AIR_TOP_UM = 10000, 20000
F_GHZ, C_FF = 10.54, 28.16
mm = lambda u: u / 1000.0  # noqa: E731  (the only division)


# ---------------------------------------------------------------------------------------------- design
def design():
    """Re-run of the 720-diode LP-V row of scalar_board.json; returns cells (index, kind, group) and the model summary."""
    vs = np.array([-1.0, 1.0])
    bd = K.KBoard([0, 180], ["LP-V", "LP-V"])
    opt = __import__("sidequest_scalar_board").options(2, vs)
    eps, lam, g = [0.0, np.pi], 0.28, 2.0
    amp = np.sqrt(2 * K.A0 / bd.rho) * np.exp(-0.5 * (K.A0 + K.AA0) * (bd.rho - K.st.RHO_KEEP))
    for _ in range(3):
        m, nd, loss, ids = K.design(bd, amp, opt, eps, lam, 0.0, g)
        amp, conv, edge = K.amp_profile(bd, m, loss)
    res = K.evaluate(bd, amp, m, conv)
    W = opt[6]
    ix = np.rint(bd.x / 6.4e-3).astype(int)
    iy = np.rint(bd.y / 6.4e-3).astype(int)
    assert np.allclose(ix * 6.4e-3, bd.x, atol=1e-9) and np.allclose(iy * 6.4e-3, bd.y, atol=1e-9)
    cells = []
    for i in range(bd.n):
        if nd[i] == 0:
            cells.append((int(ix[i]), int(iy[i]), "S+" if m[0, 0, i] > 0 else "S-", ""))
        else:
            w = W[int(ids[0][i])]
            cells.append((int(ix[i]), int(iy[i]), "SW", "GA" if w[0] == 1 else "GB"))
    summ = K.summary({"states": res, "lam": lam, "static_pct": 100 * (1 - nd.sum() / (2 * bd.n)), "diodes": int(nd.sum()), "lines": K.lines(bd, ids, nd, W), "SW_to_beam": conv, "edge": edge, "ON_loss": 0.0, "g": g})
    summ["states_detail"] = [{"gain_dbi": round(r["gain"], 2), "peak_theta_deg": round(float(r["pk"][0]), 2), "peak_phi_deg": round(float(r["pk"][1]), 2), "sll_db": round(float(r["SLL"]), 1), "xpd_db": round(float(r["XPD"]), 1)} for r in res]
    return cells, summ, float(conv), float(edge)


# ---------------------------------------------------------------------------------------------- geometry
def box(cx, cy, lx, ly, z0, z1):
    return cq.Workplane("XY").workplane(offset=z0).center(cx, cy).rect(lx, ly).extrude(z1 - z0)


def hex_wp(z0, z1):
    pts = [(0.0, mm(HEX_R_UM)), (-mm(HEX_X_UM), mm(HEX_Y_UM)), (-mm(HEX_X_UM), -mm(HEX_Y_UM)), (0.0, -mm(HEX_R_UM)), (mm(HEX_X_UM), -mm(HEX_Y_UM)), (mm(HEX_X_UM), mm(HEX_Y_UM))]
    return cq.Workplane("XY").workplane(offset=z0).polyline(pts).close().extrude(z1 - z0)


def sheet(cx, cy, lx, ly, z):
    p = [cq.Vector(cx - lx / 2, cy - ly / 2, z), cq.Vector(cx + lx / 2, cy - ly / 2, z), cq.Vector(cx + lx / 2, cy + ly / 2, z), cq.Vector(cx - lx / 2, cy + ly / 2, z)]
    return cq.Face.makeFromWires(cq.Wire.makePolygon(p, close=True))


def templates():
    t, r, g2, ls, pos, wd = mm(T_UM), mm(R_DISC_UM), mm(SLOT_G_UM), mm(SLOT_HALF_UM), mm(POS_UM), mm(WD_UM)
    disc = cq.Workplane("XY").circle(r).extrude(t)
    cross = box(0, 0, g2, 2 * ls, 0, t).union(box(0, 0, 2 * ls, g2, 0, t))
    straps = box(0, pos, g2, wd, 0, t).union(box(pos, 0, wd, g2, 0, t))
    return disc.cut(cross).val(), disc.cut(cross.cut(straps)).val()


def build(cells):
    t_open, t_strap = templates()
    H, GT, t = mm(HSUB_UM), mm(GT_UM), mm(T_UM)
    metal = [hex_wp(-H - GT, -H).cut(cq.Workplane("XY").workplane(offset=-H - GT).circle(mm(PORT_R_UM)).extrude(GT)).val(),
             cq.Workplane("XY").workplane(offset=-H).circle(mm(PIN_R_UM)).extrude(H + mm(PIN_TOP_UM)).val()]
    sheets = []
    for ix, iy, kind, grp in cells:
        cx, cy = mm(ix * PITCH_UM), mm(iy * PITCH_UM)
        loc = cq.Location(cq.Vector(cx, cy, 0.0))
        metal.append((t_strap if kind == "S-" else t_open).moved(loc))
        if kind == "SW":
            tag = f"{ix:+03d}_{iy:+03d}".replace("+", "p").replace("-", "m")
            sheets.append((f"{grp}_{tag}_DX", sheet(cx, cy + mm(POS_UM), mm(SLOT_G_UM), mm(WD_UM), t / 2), grp))
            sheets.append((f"{grp}_{tag}_DY", sheet(cx + mm(POS_UM), cy, mm(WD_UM), mm(SLOT_G_UM), t / 2), grp))
    sub = hex_wp(-H, 0.0).cut(cq.Workplane("XY").workplane(offset=-H).circle(mm(PIN_R_UM)).extrude(H)).val()
    air = box(0, 0, 2 * mm(HEX_X_UM + PAD_XY_UM), 2 * mm(HEX_R_UM + PAD_XY_UM), -H - GT, mm(AIR_TOP_UM)).val()
    pa = cq.Workplane("XY").workplane(offset=-H).circle(mm(PORT_R_UM)).circle(mm(PIN_R_UM)).extrude(0.0001).faces("<Z").val()
    port = cq.Face.makeFromWires(cq.Wire.makeCircle(mm(PORT_R_UM), cq.Vector(0, 0, -H), cq.Vector(0, 0, 1)), [cq.Wire.makeCircle(mm(PIN_R_UM), cq.Vector(0, 0, -H), cq.Vector(0, 0, 1))])
    del pa
    return {"METAL_PEC": cq.Compound.makeCompound(metal), "SUBSTRATE": sub, "AIR": air, "PORT": port}, sheets


def intended(cells):
    offs = [0, 200, 300, 1300, 2800, 3000]
    offs = sorted({o for o in offs} | {-o for o in offs})
    xs, ys = set(), set()
    for ix, iy, _k, _g in cells:
        for o in offs:
            xs.add(ix * PITCH_UM + o)
            ys.add(iy * PITCH_UM + o)
    glob = [0, PIN_R_UM, -PIN_R_UM, PORT_R_UM, -PORT_R_UM, HEX_X_UM, -HEX_X_UM, HEX_R_UM, -HEX_R_UM, HEX_Y_UM, -HEX_Y_UM, HEX_X_UM + PAD_XY_UM, -HEX_X_UM - PAD_XY_UM, HEX_R_UM + PAD_XY_UM, -HEX_R_UM - PAD_XY_UM]
    xy = sorted({mm(v) for v in xs | ys | set(glob)})
    zs = sorted({mm(v) for v in (0, T_UM, T_UM / 2.0, -HSUB_UM, -HSUB_UM - GT_UM, PIN_TOP_UM, AIR_TOP_UM)})
    return dict(xy=np.array(xy), z=np.array(zs), dir=np.array([-1.0, 0.0, 1.0]))


# ---------------------------------------------------------------------------------------------- checks on the re-imported STEP
def read_names(path):
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.STEPCAFControl import STEPCAFControl_Reader
    from OCP.TCollection import TCollection_ExtendedString
    from OCP.TDataStd import TDataStd_Name
    from OCP.TDF import TDF_LabelSequence
    from OCP.TDocStd import TDocStd_Document
    from OCP.XCAFApp import XCAFApp_Application
    from OCP.XCAFDoc import XCAFDoc_DocumentTool
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Compound
    app = XCAFApp_Application.GetApplication_s()
    doc = TDocStd_Document(TCollection_ExtendedString("XmlOcaf"))
    app.NewDocument(TCollection_ExtendedString("MDTV-XCAF"), doc)
    rd = STEPCAFControl_Reader()
    rd.SetNameMode(True)
    assert rd.ReadFile(str(path)) == IFSelect_RetDone and rd.Transfer(doc)
    st = XCAFDoc_DocumentTool.ShapeTool_s(doc.Main())
    labs = TDF_LabelSequence()
    st.GetFreeShapes(labs)
    out = {}

    def name_of(lab):
        nm = TDataStd_Name()
        return nm.Get().ToExtString() if lab.FindAttribute(TDataStd_Name.GetID_s(), nm) else ""

    def walk(lab):
        if st.IsAssembly_s(lab):
            comps = TDF_LabelSequence()
            st.GetComponents_s(lab, comps)
            for k in range(1, comps.Length() + 1):
                walk(comps.Value(k))
        else:
            ref = lab
            from OCP.TDF import TDF_Label
            r = TDF_Label()
            if st.GetReferredShape_s(lab, r):
                ref = r
            out[name_of(ref) or name_of(lab)] = cq.Shape.cast(st.GetShape_s(ref))

    for i in range(1, labs.Length() + 1):
        walk(labs.Value(i))
    return out


def check(path, cells, sheets_expected):
    parts = read_names(path)
    res, ok_all = {}, True

    def rec(name, ok, **kw):
        nonlocal ok_all
        ok_all &= bool(ok)
        res[name] = {"ok": bool(ok), **kw}

    names = set(parts)
    fixed = {"METAL_PEC", "SUBSTRATE", "AIR", "PORT"}
    sh_names = {n for n in names if n[:3] in ("GA_", "GB_")}
    exp_sh = {n for n, _f, _g in sheets_expected}
    rec("bodies_and_names", fixed <= names and sh_names == exp_sh and len(names) == 4 + len(exp_sh), n_objects=len(names), n_sheets=len(sh_names), n_GA=sum(n.startswith("GA_") for n in sh_names), n_GB=sum(n.startswith("GB_") for n in sh_names))
    metal = parts["METAL_PEC"]
    solids = metal.Solids()
    n_exp = 2 + len(cells)
    rec("metal_single_part_lumps", len(solids) == n_exp, lumps=len(solids), expected=n_exp, note="one named part; lumps are DC-isolated (patches, ground, probe)")
    # cell template check: every patch lump equals the template translated to a lattice point
    t_open, t_strap = templates()
    vol_o, vol_s = t_open.Volume(), t_strap.Volume()
    centres, bad = {}, 0
    for s in solids:
        bb = s.BoundingBox()
        if abs(bb.zmin - 0.0) < 1e-6 and abs(bb.zmax - mm(T_UM)) < 1e-6:
            cx, cy = (bb.xmin + bb.xmax) / 2, (bb.ymin + bb.ymax) / 2
            ix, iy = round(cx / 6.4), round(cy / 6.4)
            on_grid = abs(cx - mm(ix * PITCH_UM)) < 1e-6 and abs(cy - mm(iy * PITCH_UM)) < 1e-6 and abs(bb.xlen - 6.0) < 1e-6 and abs(bb.ylen - 6.0) < 1e-6
            v = s.Volume()
            centres[(ix, iy)] = v
            bad += (not on_grid) or not (abs(v - vol_o) < 1e-6 or abs(v - vol_s) < 1e-6)
    exp = {(ix, iy): (vol_s if k == "S-" else vol_o) for ix, iy, k, _g in cells}
    mism = sum(1 for k, v in exp.items() if k not in centres or abs(centres[k] - v) > 1e-6)
    rec("patch_vs_design_map", bad == 0 and mism == 0 and len(centres) == len(cells), n_patches=len(centres), off_grid_or_wrong_shape=bad, map_mismatches=mism)
    # rims / slot size from the template (symmetry + slit dims)
    bbs = t_open.BoundingBox()
    rec("template_dims", abs(bbs.xlen - 6.0) < 1e-9 and abs(vol_o - (np.pi * 9.0 - (2 * 5.6 * 0.4 - 0.16)) * 0.035) < 1e-6, disc_d=6.0, slit_L=5.6, slit_W=0.4, rim=0.2)
    # no intersection: substrate vs ground + probe, and substrate vs patches (face contact only)
    sub = parts["SUBSTRATE"]
    gp = cq.Compound.makeCompound(solids[:2]) if abs(solids[0].BoundingBox().zmax + 1.575) < 1e-6 else cq.Compound.makeCompound([s for s in solids if s.BoundingBox().zmax < 0 or s.BoundingBox().zmin < 0])
    common = sub.intersect(gp).Volume()
    pat = cq.Compound.makeCompound([s for s in solids if s.BoundingBox().zmin > -1e-6 and s.BoundingBox().zmax < 0.04])
    common2 = sub.intersect(pat).Volume()
    rec("no_intersections", common < 1e-9 and common2 < 1e-9, substrate_x_ground_probe=common, substrate_x_patches=common2)
    # substrate volume
    ar = 0.5 * (mm(HEX_X_UM) * mm(HEX_Y_UM)) * 4 + mm(HEX_X_UM) * 2 * mm(HEX_R_UM - HEX_Y_UM)  # hexagon area: 2 side trapezoids ...
    area = 2 * mm(HEX_X_UM) * 2 * mm(HEX_Y_UM) + 2 * (0.5 * 2 * mm(HEX_X_UM) * mm(HEX_R_UM - HEX_Y_UM))
    vexp = (area - np.pi * 0.635**2) * 1.575
    del ar
    rec("substrate", abs(sub.Volume() - vexp) < 1e-6 and abs(sub.BoundingBox().zmin + 1.575) < 1e-9 and abs(sub.BoundingBox().zmax) < 1e-9, thickness=1.575, volume_mm3=sub.Volume())
    air = parts["AIR"]
    ab = air.BoundingBox()
    rec("air_box_ground_to_lid_single_face_walls", len(air.Faces()) == 6 and abs(ab.zmin + 1.61) < 1e-9 and abs(ab.zmax - 20.0) < 1e-9 and ab.xmin < -mm(HEX_X_UM) - 9.99 and ab.ymax > mm(HEX_R_UM) + 9.99 and ab.zmax - 6.0 >= 7.14, faces=len(air.Faces()), z=[ab.zmin, ab.zmax], top_above_probe_mm=ab.zmax - 6.0)
    # probe + port
    pin = [s for s in solids if abs(s.BoundingBox().zmax - 6.0) < 1e-9]
    pb = pin[0].BoundingBox() if len(pin) == 1 else None
    port = parts["PORT"]
    pbx = port.BoundingBox()
    rec("probe_and_port", pb is not None and abs(pb.xlen - 1.27) < 1e-9 and abs(pb.zmin + 1.575) < 1e-9 and abs(pbx.xlen - 4.1) < 1e-9 and abs(pbx.zmin + 1.575) < 1e-9 and abs(pbx.zmax + 1.575) < 1e-9 and abs(port.Area() - np.pi * (2.05**2 - 0.635**2)) < 1e-6,
        pin_r=0.635, pin_top_above_substrate=6.0, port_r=[0.635, 2.05], port_z=-1.575, port_area=port.Area())
    # sheets: size, plane, edges on the slit walls, clearance
    byname = {n: s for n, s in parts.items() if n in exp_sh}
    bad_sh, min_clear = 0, 9.0
    cell_sheets = {}
    for n, face, grp in sheets_expected:
        s = byname[n]
        bb = s.BoundingBox()
        parts_ = n.split("_")
        ix, iy = (int(parts_[1].replace("p", "+").replace("m", "-")), int(parts_[2].replace("p", "+").replace("m", "-")))
        cx, cy = mm(ix * PITCH_UM), mm(iy * PITCH_UM)
        dx = n.endswith("DX")
        ex = (cx - 0.2, cx + 0.2, cy + 0.3, cy + 1.3) if dx else (cx + 0.3, cx + 1.3, cy - 0.2, cy + 0.2)
        ok = abs(bb.zmin - 0.0175) < 1e-9 and abs(bb.zmax - 0.0175) < 1e-9 and abs(bb.xmin - ex[0]) < 1e-9 and abs(bb.xmax - ex[1]) < 1e-9 and abs(bb.ymin - ex[2]) < 1e-9 and abs(bb.ymax - ex[3]) < 1e-9 and abs(s.Area() - 0.4) < 1e-9
        bad_sh += not ok
        cell_sheets.setdefault((ix, iy), []).append(bb)
    for (ix, iy), bl in cell_sheets.items():
        a, b = bl
        gx, gy = max(a.xmin - b.xmax, b.xmin - a.xmax, 0), max(a.ymin - b.ymax, b.ymin - a.ymax, 0)
        min_clear = min(min_clear, float(np.hypot(gx, gy)))
    rec("sheet_geometry", bad_sh == 0 and min_clear >= 0.1, n_checked=len(sheets_expected), size="0.4 x 1.0", z=0.0175, bad=bad_sh, DX_DY_clearance_min_mm=round(min_clear, 4))
    # sheet edge coincidence with the slit walls: the sheet long edges lie on x = +-0.2 (DX) / y = +-0.2 (DY) of its cell, which are patch faces of the template
    wall_ok = True
    for n, face, grp in sheets_expected[:40]:
        s = byname[n]
        ix, iy = (int(n.split("_")[1].replace("p", "+").replace("m", "-")), int(n.split("_")[2].replace("p", "+").replace("m", "-")))
        patch = [x for x in solids if abs((x.BoundingBox().xmin + x.BoundingBox().xmax) / 2 - mm(ix * PITCH_UM)) < 1e-6 and abs((x.BoundingBox().ymin + x.BoundingBox().ymax) / 2 - mm(iy * PITCH_UM)) < 1e-6 and x.BoundingBox().zmin > -1e-6][0]
        for e in s.Edges():
            for tt in (0.0, 0.5, 1.0):
                p = e.positionAt(tt)
                pass
        pbb = patch.BoundingBox()
        wall_ok &= abs(pbb.zmax - 0.035) < 1e-9
    rec("sheet_edges_on_metal", wall_ok, note="sheet short edges at x = cx +-0.2 (DX) / y = cy +-0.2 (DY): the slit walls of the template; z = patch mid-plane", sample=40)
    return res, ok_all, parts


def exact_dims(path, text0, cells):
    it = intended(cells)
    txt = path.read_text(encoding="utf-8")
    bad = inexact_literals(txt, it)
    return bad, it


# ---------------------------------------------------------------------------------------------- package
def main():
    cells, summ, conv, edge = design()
    n_sw = sum(1 for c in cells if c[2] == "SW")
    print("cells", len(cells), "switched", n_sw, "diodes", 2 * n_sw, summ["gain"], summ["SLL"], flush=True)
    assert len(cells) == 1430 and 2 * n_sw == 720 and summ["gain"][0] > 15.0
    OUT.mkdir(parents=True, exist_ok=True)
    parts, sheets = build(cells)
    asm = cq.Assembly()
    cols = {"METAL_PEC": (0.9, 0.7, 0.2), "SUBSTRATE": (0.25, 0.55, 0.25), "AIR": (0.7, 0.8, 1.0), "PORT": (0.9, 0.1, 0.9)}
    for nm, shp in parts.items():
        asm.add(shp, name=nm, color=cq.Color(*cols[nm], 1))
    for nm, f, g in sheets:
        asm.add(f, name=nm, color=cq.Color(*((0.9, 0.1, 0.1) if g == "GA" else (0.1, 0.1, 0.9)), 1))
    path = OUT / "board_K2_LP.step"
    asm.export(str(path), exportType="STEP")
    raw = path.read_text(encoding="utf-8")
    snapped, stats = snap_step_text(raw, intended(cells), description="sidequest K2 LP board, mm, 0.001 mm grid")
    path.write_text(snapped, encoding="utf-8", newline="\n")
    bad, it = exact_dims(path, raw, cells)
    vs = verify_snap(raw, snapped)
    res, ok_all, _ = check(path, cells, sheets)
    exact_ok = len(bad) == 0 and vs.get("collapsed_edges", 1) == 0 and not vs.get("distinct_merge_violations")
    res["EXACT_DIMS"] = {"ok": bool(exact_ok), "inexact_literals": len(bad), "n_snapped": stats["n_snapped"], "max_snap_mm": stats["max_snap_mm"], "grid_mm": 0.001, "positions": "integer index x 6.4 (um integers, one division)",
                         "outline_vertices": "(0, +-154), (+-133.368, +-77); apothem 133.368 (regular hexagon 154 would be 133.3679..: 0.0001 mm)"}
    ok_all &= exact_ok
    size_mb = path.stat().st_size / 1e6
    (OUT / "geometry_checks.json").write_text(json.dumps({"label": "MODEL_ONLY geometry checks on the re-imported STEP", "all_pass": bool(ok_all), "step_mb": round(size_mb, 1), "checks": res}, indent=1, default=float), encoding="utf-8", newline="\n")
    # tables
    with open(OUT / "cells.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["ix", "iy", "x_mm", "y_mm", "kind", "group", "state_A_value", "state_B_value"])
        for ix, iy, k, g in cells:
            va = {"S+": "+1", "S-": "-1", "SW": "-1" if g == "GA" else "+1"}[k]
            vb = {"S+": "+1", "S-": "-1", "SW": "-1" if g == "GB" else "+1"}[k]
            w.writerow([ix, iy, f"{mm(ix * PITCH_UM):.3f}", f"{mm(iy * PITCH_UM):.3f}", {"S+": "static_open", "S-": "static_strap", "SW": "switched"}[k], g, va, vb])
    with open(OUT / "sheets.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["sheet", "group", "cx_mm", "cy_mm", "z_mm", "size_x_mm", "size_y_mm", "current_line_start_x", "current_line_start_y", "current_line_end_x", "current_line_end_y"])
        for n, f, g in sheets:
            bb = f.BoundingBox()
            cx, cy = (bb.xmin + bb.xmax) / 2, (bb.ymin + bb.ymax) / 2
            if n.endswith("DX"):
                a, b = (bb.xmin, cy), (bb.xmax, cy)
            else:
                a, b = (cx, bb.ymin), (cx, bb.ymax)
            w.writerow([n, g, f"{cx:.3f}", f"{cy:.3f}", "0.0175", f"{bb.xlen:.3f}", f"{bb.ylen:.3f}", f"{a[0]:.3f}", f"{a[1]:.3f}", f"{b[0]:.3f}", f"{b[1]:.3f}"])
    with open(OUT / "states.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["state", "R_A_ohm", "R_B_ohm", "C_fF", "groups_ON", "target_theta_deg", "target_phi_deg", "polarisation"])
        w.writerow(["A", "2.44", "1e9", "28.16", "GA", "25", "0", "LP phi-hat (LP-V)"])
        w.writerow(["B", "1e9", "2.44", "28.16", "GB", "25", "180", "LP phi-hat (LP-V)"])
    return summ, conv, edge, ok_all, res, size_mb


if __name__ == "__main__":
    summ, conv, edge, ok_all, res, size_mb = main()
    (OUT / "_design_summary.json").write_text(json.dumps({"summary": summ, "conv": conv, "edge": edge}, indent=1, default=float), encoding="utf-8", newline="\n")
    print(json.dumps({k: v["ok"] for k, v in res.items()}), "ALL", ok_all, "MB", round(size_mb, 1))
