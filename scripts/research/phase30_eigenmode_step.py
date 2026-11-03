"""Phase 30: STEP files for the HFSS eigenmode check (replaces the AEDT IronPython scripts, which do not run on the user's machine).

One STEP per model in hfss_exports/phase30_cell_eigenmode/: E1..E4 unit cells and the S1 leaky 7-cell supercell. mm, origin at the
cell (supercell) centre on the ground plane. Bodies (CadQuery assembly names): SUB (solid, z 0..0.635), GND (zero-thickness sheet
z = 0), PATCH (zero-thickness sheet(s) at z = 0.635: n = 4 super-ellipse, 72-vertex polygon, minus the rectangular slot), AIR (solid
box above the patch, height AIR_H). The slot turns, the outline stays fixed (the tensor model's rotation check holds only so).
The STEP files stay gitignored; copies go to C:\\hmsa_sync\\data\\phase30_cell_eigenmode\\.   python phase30_eigenmode_step.py [--no-sync]
"""

import json
import math
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PKG = ROOT / "hfss_exports" / "phase30_cell_eigenmode"
SYNC = Path(r"C:\hmsa_sync\data\phase30_cell_eigenmode")
P, H, AIR_H, NV, EXPO = 3.2, 0.635, 60.0, 72, 4.0
MAJOR = 3.05  # gap = P - MAJOR = 0.15 mm
MODELS = {  # name: (slot L, slot W, number of cells)
    "E1_strong_M0p64": (1.89, 0.58, 1),
    "E2_longthin_M0p62": (2.29, 0.10, 1),
    "E3_mid_M0p27": (1.89, 0.10, 1),
    "E4_weak_M0p12": (1.53, 0.10, 1),
    "S1_supercell_leaky": (1.89, 0.58, 7),
}


def outline(cx):
    pts = []
    for k in range(NV):
        t = 2 * math.pi * k / NV
        c, s = math.cos(t), math.sin(t)
        pts.append((cx + 0.5 * MAJOR * math.copysign(abs(c) ** (2 / EXPO), c), 0.5 * MAJOR * math.copysign(abs(s) ** (2 / EXPO), s)))
    return pts


def slot(cx, length, width, ang_deg):
    a = math.radians(ang_deg)
    c, s = math.cos(a), math.sin(a)
    return [(cx + u * c - v * s, u * s + v * c) for u, v in ((-length / 2, -width / 2), (length / 2, -width / 2), (length / 2, width / 2), (-length / 2, width / 2))]


def poly_area(pts):
    return 0.5 * abs(sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(pts, pts[1:] + pts[:1], strict=True)))


def cells(n):
    return [(-(n * P) / 2 + (i + 0.5) * P, 180.0 * i / n) for i in range(n)]


def build(name, out):
    import cadquery as cq

    L, W, n = MODELS[name]
    V = cq.Vector
    wx = n * P
    sub = cq.Workplane("XY").box(wx, P, H, centered=(True, True, False)).val()
    air = cq.Workplane("XY").workplane(offset=H).box(wx, P, AIR_H, centered=(True, True, False)).val()
    gnd = cq.Face.makeFromWires(cq.Wire.makePolygon([V(-wx / 2, -P / 2, 0), V(wx / 2, -P / 2, 0), V(wx / 2, P / 2, 0), V(-wx / 2, P / 2, 0)], close=True))
    faces = []
    for cx, ang in cells(n):
        o = cq.Wire.makePolygon([V(x, y, H) for x, y in outline(cx)], close=True)
        h = cq.Wire.makePolygon([V(x, y, H) for x, y in slot(cx, L, W, ang)], close=True)
        faces.append(cq.Face.makeFromWires(o, [h]))
    assy = cq.Assembly(name=name)
    assy.add(sub, name="SUB")
    assy.add(gnd, name="GND")
    assy.add(cq.Compound.makeCompound(faces), name="PATCH")
    assy.add(air, name="AIR")
    assy.save(str(out), exportType="STEP")


def read_back(path):
    """Re-open the STEP with OCP/XCAF: per named body: solids, faces, bbox, area, volume, and for PATCH the hole sizes."""
    from OCP.Bnd import Bnd_Box
    from OCP.BRep import BRep_Tool
    from OCP.BRepBndLib import BRepBndLib
    from OCP.BRepGProp import BRepGProp
    from OCP.BRepTools import BRepTools
    from OCP.GProp import GProp_GProps
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.STEPCAFControl import STEPCAFControl_Reader
    from OCP.TCollection import TCollection_ExtendedString
    from OCP.TDataStd import TDataStd_Name
    from OCP.TDF import TDF_Label, TDF_LabelSequence
    from OCP.TDocStd import TDocStd_Document
    from OCP.TopAbs import TopAbs_FACE, TopAbs_SOLID, TopAbs_VERTEX, TopAbs_WIRE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS
    from OCP.XCAFDoc import XCAFDoc_DocumentTool

    doc = TDocStd_Document(TCollection_ExtendedString("doc"))
    rd = STEPCAFControl_Reader()
    rd.SetNameMode(True)
    assert rd.ReadFile(str(path)) == IFSelect_RetDone and rd.Transfer(doc)
    st = XCAFDoc_DocumentTool.ShapeTool_s(doc.Main())
    labs = TDF_LabelSequence()
    st.GetFreeShapes(labs)
    res = {}

    def explore(shp, kind):
        ex = TopExp_Explorer(shp, kind)
        while ex.More():
            yield ex.Current()
            ex.Next()

    def bbox(shp):
        b = Bnd_Box()
        BRepBndLib.AddOptimal_s(shp, b, False, False)
        x0, y0, z0, x1, y1, z1 = b.Get()
        return [round(v, 4) for v in (x0, y0, z0, x1, y1, z1)]

    def walk(label):
        comps = TDF_LabelSequence()
        if st.IsAssembly_s(label):
            st.GetComponents_s(label, comps)
            for i in range(1, comps.Length() + 1):
                walk(comps.Value(i))
            return
        ref = label
        if st.IsReference_s(label):
            ref = TDF_Label()
            st.GetReferredShape_s(label, ref)
        nm, name = TDataStd_Name(), None
        for lb in (label, ref):
            if lb.FindAttribute(TDataStd_Name.GetID_s(), nm):
                name = nm.Get().ToExtString()
                if name and not name.startswith("="):
                    break
        shp = st.GetShape_s(ref)
        g = GProp_GProps()
        BRepGProp.SurfaceProperties_s(shp, g)
        area = g.Mass()
        g2 = GProp_GProps()
        BRepGProp.VolumeProperties_s(shp, g2)
        d = {"solids": sum(1 for _ in explore(shp, TopAbs_SOLID)), "faces": sum(1 for _ in explore(shp, TopAbs_FACE)),
             "bbox": bbox(shp), "area": round(area, 5), "volume": round(g2.Mass(), 5)}
        if name == "PATCH":
            holes = []
            for fc in explore(shp, TopAbs_FACE):
                f = TopoDS.Face_s(fc)
                outer = BRepTools.OuterWire_s(f)
                for w in explore(f, TopAbs_WIRE):
                    if not w.IsSame(outer):
                        pts = [BRep_Tool.Pnt_s(TopoDS.Vertex_s(v)) for v in explore(w, TopAbs_VERTEX)]
                        xs = sorted({round(p.X(), 4) for p in pts})
                        ys = sorted({round(p.Y(), 4) for p in pts})
                        gw = GProp_GProps()
                        from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace

                        BRepGProp.SurfaceProperties_s(BRepBuilderAPI_MakeFace(TopoDS.Wire_s(w)).Face(), gw)
                        holes.append({"hole_area": round(gw.Mass(), 5), "bbox": bbox(w), "n_vertices": len(set((round(p.X(), 5), round(p.Y(), 5)) for p in pts))})
            d["holes"] = holes
        res[name] = d

    for i in range(1, labs.Length() + 1):
        walk(labs.Value(i))
    return res


def check(name, path):
    L, W, n = MODELS[name]
    r = read_back(path)
    wx = n * P
    tol = 2e-3
    out = {"file": path.name, "bodies": {k: {kk: v[kk] for kk in ("solids", "faces")} for k, v in r.items()}, "checks": {}}
    c = out["checks"]
    c["body_names"] = sorted(r) == ["AIR", "GND", "PATCH", "SUB"]
    c["SUB_box"] = r["SUB"]["solids"] == 1 and all(abs(a - b) < tol for a, b in zip(r["SUB"]["bbox"], (-wx / 2, -P / 2, 0, wx / 2, P / 2, H), strict=True))
    c["AIR_box"] = r["AIR"]["solids"] == 1 and all(abs(a - b) < tol for a, b in zip(r["AIR"]["bbox"], (-wx / 2, -P / 2, H, wx / 2, P / 2, H + AIR_H), strict=True))
    c["GND_sheet"] = r["GND"]["solids"] == 0 and r["GND"]["faces"] == 1 and abs(r["GND"]["area"] - wx * P) < 1e-3 and r["GND"]["bbox"][2] == r["GND"]["bbox"][5] == 0
    pb = r["PATCH"]["bbox"]
    c["PATCH_n_faces"] = r["PATCH"]["faces"] == n and r["PATCH"]["solids"] == 0
    c["PATCH_z"] = abs(pb[2] - H) < tol and abs(pb[5] - H) < tol
    gap = P - (pb[3] - pb[0] - (n - 1) * P)
    out["gap_to_cell_edge_x_mm"] = round(P - (pb[3] - pb[0] - (n - 1) * P), 4)
    out["gap_between_neighbour_patches_mm"] = round(P - MAJOR, 4)
    c["gap_0p15"] = abs(gap - 0.15) < tol and abs(pb[4] - pb[1] - MAJOR) < tol
    exp_area = poly_area(outline(0.0)) - L * W
    c["PATCH_area_per_cell"] = abs(r["PATCH"]["area"] / n - exp_area) < 2e-3
    holes = r["PATCH"]["holes"]
    c["slot_count"] = len(holes) == n
    c["slot_area_LxW"] = all(abs(h["hole_area"] - L * W) < 1e-3 for h in holes)
    h0 = min(holes, key=lambda h: abs(0.5 * (h["bbox"][0] + h["bbox"][3])))  # cell nearest the centre
    out["slot_bbox_first_hole_mm"] = h0["bbox"]
    out["slot_LxW_mm_from_area_and_bbox"] = [round(h0["bbox"][3] - h0["bbox"][0], 4), round(h0["bbox"][4] - h0["bbox"][1], 4)]
    if n == 1:
        c["slot_bbox_LxW"] = abs((h0["bbox"][3] - h0["bbox"][0]) - L) < tol and abs((h0["bbox"][4] - h0["bbox"][1]) - W) < tol
    out["patch_area_mm2"] = r["PATCH"]["area"]
    out["patch_vertices_per_outline"] = NV
    out["all_pass"] = all(c.values())
    return out


def main(sync=True):
    PKG.mkdir(parents=True, exist_ok=True)
    report = {}
    for name in MODELS:
        p = PKG / (name + ".step")
        build(name, p)
        report[name] = check(name, p)
        print(name, "PASS" if report[name]["all_pass"] else "FAIL", report[name]["bodies"], flush=True)
    (PKG / "step_check.json").write_text(json.dumps(report, indent=1), encoding="utf-8", newline="\n")
    if sync and "--no-sync" not in sys.argv:
        SYNC.mkdir(parents=True, exist_ok=True)
        for name in MODELS:
            shutil.copy2(PKG / (name + ".step"), SYNC / (name + ".step"))
        shutil.copy2(PKG / "HFSS_SETUP.md", SYNC / "HFSS_SETUP.md")
    return report


if __name__ == "__main__":
    main()
