"""Geometry-consistency and contact checks of the crossed-diode rim-disc cell (bodies METAL_PEC = patch + ground, SUBSTRATE, AIR, sheets DX_SHEET, DY_SHEET) against the MoM geometry.
Used (a) on the shapes before export and (b) on the RE-IMPORTED STEP file:   python check_step_cell.py <file.step> [--pos 0.8] [--wd 1.0] [--ls 2.8] [--g 0.4]   (mm; ls >= 3.0: no rim)
Checks (reported as GEOM_CHECKS): bodies_and_names, no_intersections (pairwise common volume of all solids), air_is_box_minus_metal (volume), patch_vs_MoM (mask on a 0.05 mm grid, diode windows and the disc boundary band excluded),
rims_at_4_ends, symmetry (mirror x, mirror y), sheet_edge_coincidence (long edges lie on the patch side faces, ends on the faces, centre in the gap, z = patch mid-plane, size), sheet_to_sheet_clearance (>= 0.1 mm),
sheet_to_other_slit_clearance (>= 0.1 mm to the patch faces the other sheet bridges)."""
import argparse
import json
import re
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import cadquery as cq
import numpy as np
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.gp import gp_Pnt
from OCP.TopAbs import TopAbs_IN

R_DISC, T, H, GT, LID, P = 3.0, 0.035, 1.575, 0.035, 30.0, 6.4
NAMES = ("METAL_PEC", "SUBSTRATE", "AIR", "DX_SHEET", "DY_SHEET")
TOL = 1e-7


def _inside(solid, x, y, z):
    return BRepClass3d_SolidClassifier(solid.wrapped, gp_Pnt(float(x), float(y), float(z)), 1e-7).State() == TopAbs_IN


def _pd(edge, t, patch):
    p = edge.positionAt(t)
    return cq.Vertex.makeVertex(p.x, p.y, p.z).distance(patch)


def classify(shape):
    """Split a shape/compound into patch, ground, substrate, air solids and the sheet faces by geometry."""
    out = {"patch": [], "ground": [], "substrate": [], "air": [], "sheets": []}
    for s in shape.Solids():
        b = s.BoundingBox()
        if abs(b.zmax - LID) < 1e-3:
            out["air"].append(s)
        elif abs(b.zmin + H) < 1e-3 and abs(b.zmax) < 1e-3:
            out["substrate"].append(s)
        elif abs(b.zmax + H) < 1e-3:
            out["ground"].append(s)
        elif abs(b.zmax - T) < 1e-4 and abs(b.zmin) < 1e-4:
            out["patch"].append(s)
    for f in shape.Faces():
        b = f.BoundingBox()
        if b.zlen < 1e-9 and b.xlen < 2.0 and b.ylen < 2.0 and abs(b.zmin - T / 2) < 1e-6:
            out["sheets"].append(f)
    return out


def check(shape, pos=0.8, ls=2.8, g=0.4, wd=1.0, step_text=None, min_clearance=0.1, step=0.05):
    c = classify(shape)
    res = {}
    ok_all = True

    def rec(name, ok, **kw):
        nonlocal ok_all
        res[name] = {"ok": bool(ok), **kw}
        ok_all &= bool(ok)

    names_ok = True if step_text is None else all(n in step_text for n in NAMES)
    rec("bodies_and_names", len(c["patch"]) == 1 and len(c["ground"]) == 1 and len(c["substrate"]) == 1 and len(c["air"]) == 1 and len(c["sheets"]) == 2 and names_ok,
        patch=len(c["patch"]), ground=len(c["ground"]), substrate=len(c["substrate"]), air=len(c["air"]), sheets=len(c["sheets"]), names_in_file=names_ok)
    if not (c["patch"] and c["air"] and c["substrate"] and c["ground"]) or len(c["sheets"]) != 2:
        return {"PASS": False, "checks": res}
    patch, ground, sub, air = c["patch"][0], c["ground"][0], c["substrate"][0], c["air"][0]
    solids = {"patch": patch, "ground": ground, "substrate": sub, "air": air}
    inter = {}
    for a_, b_ in [("patch", "ground"), ("patch", "substrate"), ("patch", "air"), ("ground", "substrate"), ("ground", "air"), ("substrate", "air")]:
        v = solids[a_].intersect(solids[b_]).Volume()
        inter[f"{a_}&{b_}"] = float(v)
    rec("no_intersections", max(inter.values()) < 1e-9, common_volume_mm3=inter)
    v_expect = P * P * LID - patch.Volume()
    rec("air_is_box_minus_metal", abs(air.Volume() - v_expect) < 1e-6 * v_expect, air_volume=float(air.Volume()), expected=float(v_expect))
    # MoM mask
    import sidequest_cell_search as C

    mom = C.cell("rt5880", "OFF", "OFF", diam=6.0, g=g, ls=ls, pos=pos, wd=wd)
    xs = np.arange(-3.0, 3.0 + 1e-9, step)
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    keep = np.hypot(X, Y) < R_DISC - 0.12
    win_x = (np.abs(X) <= g / 2 + 1e-9) & (Y >= pos - wd / 2 - 1e-9) & (Y <= pos + wd / 2 + 1e-9)
    win_y = (np.abs(Y) <= g / 2 + 1e-9) & (X >= pos - wd / 2 - 1e-9) & (X <= pos + wd / 2 + 1e-9)
    keep &= ~(win_x | win_y)
    for e in (g / 2, ls, pos - wd / 2, pos + wd / 2):
        keep &= (np.abs(np.abs(X) - e) > 1e-6) & (np.abs(np.abs(Y) - e) > 1e-6)
    mom_mask = mom.is_metal(X * 1e-3, Y * 1e-3)
    idx = np.argwhere(keep)
    zmid = T / 2
    step_mask = np.array([_inside(patch, X[i, j], Y[i, j], zmid) for i, j in idx])
    diff = step_mask != mom_mask[keep]
    rec("patch_vs_MoM", diff.mean() < 0.002, points=int(len(idx)), mismatch_fraction=float(diff.mean()), first_mismatches=[(float(X[tuple(i)]), float(Y[tuple(i)])) for i in idx[diff][:5]])
    want = ls + 0.1 < R_DISC - 1e-6
    rim_pts = [(0, ls + 0.1), (0, -ls - 0.1), (ls + 0.1, 0), (-ls - 0.1, 0)]
    got = [_inside(patch, x, y, zmid) for x, y in rim_pts]
    rec("rims_at_4_ends", all(v == want for v in got), expected_present=want, top_bottom_right_left=[got[0], got[1], got[2], got[3]])
    pts = [(x, y) for x in np.arange(-2.95, 3.0, 0.35) for y in np.arange(-2.95, 3.0, 0.35)]
    base = {p: _inside(patch, p[0], p[1], zmid) for p in pts}
    mx = sum(base[p] != _inside(patch, -p[0], p[1], zmid) for p in pts)
    my = sum(base[p] != _inside(patch, p[0], -p[1], zmid) for p in pts)
    rec("symmetry", mx == 0 and my == 0, mirror_x_asymmetric_points=int(mx), mirror_y_asymmetric_points=int(my), tested=len(pts))
    # sheets
    info, axes = [], []
    for f in c["sheets"]:
        b, cen = f.BoundingBox(), f.Center()
        along_y = b.ylen > b.xlen  # long side along y: bridges the vertical slit (x diode)
        gap_len, width = (b.xlen, b.ylen) if along_y else (b.ylen, b.xlen)
        edges = f.Edges()
        long_edges = [e for e in edges if abs(e.Length() - width) < 1e-6]
        long_ok = len(long_edges) == 2 and all(all(_pd(e, t, patch) < TOL for t in (0, 0.25, 0.5, 0.75, 1.0)) for e in long_edges)
        ends_ok = all(_pd(e, t, patch) < TOL for e in edges for t in (0.0, 1.0))
        centre_gap = not _inside(patch, cen.x, cen.y, cen.z) and cq.Vertex.makeVertex(cen.x, cen.y, cen.z).distance(patch) > 0.01
        size_ok = abs(gap_len - g) < 1e-6 and abs(width - wd) < 1e-6
        axes.append("x" if along_y else "y")
        info.append({"axis": "x" if along_y else "y", "centre": [round(cen.x, 4), round(cen.y, 4), round(cen.z, 4)], "gap_len": round(gap_len, 4), "width": round(width, 4), "long_edges_on_patch": bool(long_ok), "end_points_on_patch": bool(ends_ok), "centre_in_gap": bool(centre_gap), "size_ok": bool(size_ok)})
    rec("sheet_edge_coincidence", all(i["long_edges_on_patch"] and i["end_points_on_patch"] and i["centre_in_gap"] and i["size_ok"] for i in info) and sorted(axes) == ["x", "y"], sheets=info)
    d_ss = c["sheets"][0].distance(c["sheets"][1])
    rec("sheet_to_sheet_clearance", d_ss >= min_clearance - 1e-9, distance_mm=float(d_ss), required=min_clearance)
    # clearance from each sheet to the patch faces bridged by the other sheet
    worst = 1e9
    for k in (0, 1):
        other, this = c["sheets"][1 - k], c["sheets"][k]
        for fc in patch.Faces():
            if fc.distance(other) < TOL and fc.distance(this) > TOL:  # a face the other sheet bridges and this sheet does not touch
                worst = min(worst, fc.distance(this))
    rec("sheet_to_other_slit_clearance", worst >= min_clearance - 1e-9, distance_mm=float(worst), required=min_clearance)
    return {"PASS": bool(ok_all), "checks": res}


def check_file(path, **kw):
    text = Path(path).read_text(encoding="utf-8")
    shp = cq.importers.importStep(str(path)).val()
    comp = cq.Compound.makeCompound(list(cq.importers.importStep(str(path)).vals())) if not hasattr(shp, "Solids") else shp
    return check(comp, step_text=text, **kw)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("step")
    ap.add_argument("--pos", type=float, default=0.8)
    ap.add_argument("--wd", type=float, default=1.0)
    ap.add_argument("--ls", type=float, default=2.8)
    ap.add_argument("--g", type=float, default=0.4)
    a = ap.parse_args()
    r = check_file(a.step, pos=a.pos, ls=a.ls, g=a.g, wd=a.wd)
    print(json.dumps(r, indent=1, default=str))
    sys.exit(0 if r["PASS"] else 1)
