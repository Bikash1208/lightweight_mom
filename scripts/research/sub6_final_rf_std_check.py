"""Re-import gate for the RF-only standard STEP: names / counts survive, no overlaps, sheets on the side faces, port, air box, exact dims.  MODEL_ONLY.
    python sub6_final_rf_std_check.py <step> <expected_D_count> <out.json>
"""
import json
import math
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import cadquery as cq
from sidequest_board_k2_lp_package import read_names
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepExtrema import BRepExtrema_DistShapeShape
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps

step, nD, outp = sys.argv[1], int(sys.argv[2]), sys.argv[3]
t0 = time.time()
names = read_names(step)
W = lambda s: s.wrapped if hasattr(s, "wrapped") else s
cqs = lambda s: s if hasattr(s, "Solids") else cq.Shape.cast(s)


def vol(s):
    tot = 0.0
    for so in cqs(s).Solids():
        p = GProp_GProps()
        BRepGProp.VolumeProperties_s(so.wrapped, p)
        tot += p.Mass()
    return tot


def area(s):
    p = GProp_GProps()
    BRepGProp.SurfaceProperties_s(W(s), p)
    return p.Mass()


def bbox(s):
    b = Bnd_Box()
    BRepBndLib.Add_s(W(s), b)
    return b.Get()


res = dict(label="MODEL_ONLY re-import gate", step=step)
D = [n for n in names if n.startswith("D") and "_" in n and n[1:2].isdigit()]
res["names_found"] = dict(METAL_PEC=int("METAL_PEC" in names), SUBSTRATE=int("SUBSTRATE" in names), SUBSTRATE_FR4=int("SUBSTRATE_FR4" in names), AIR=int("AIR" in names), PORT=int("PORT" in names), D_sheets=len(D), total_products=len(names))
res["names_ok"] = bool(res["names_found"]["METAL_PEC"] == 1 and res["names_found"]["PORT"] == 1 and len(D) == nD and res["names_found"]["AIR"] == 1 and res["names_found"]["SUBSTRATE"] == 1 and res["names_found"]["SUBSTRATE_FR4"] == 1 and len(names) == nD + 5)
metal = cqs(names["METAL_PEC"])
lumps = metal.Solids()
res["METAL_PEC_lumps"] = len(lumps)
res["METAL_PEC_all_solids_valid"] = all(BRepCheck_Analyzer(l.wrapped).IsValid() for l in lumps)
res["all_bodies_valid"] = {n: bool(BRepCheck_Analyzer(W(s)).IsValid()) for n, s in names.items() if n not in D}
# NO_OVERLAP: METAL_PEC vs SUBSTRATE / FR4 (common volume), lumps among themselves (bbox-overlapping pairs)
ov = {}
for dn in ("SUBSTRATE", "SUBSTRATE_FR4"):
    c = BRepAlgoAPI_Common(W(names["METAL_PEC"]), W(names[dn]))
    ov[dn] = vol(cq.Shape.cast(c.Shape())) if c.IsDone() else float("nan")
bbs = [bbox(l) for l in lumps]
worst = 0.0
npairs = 0
for i in range(len(lumps)):
    for j in range(i + 1, len(lumps)):
        a, b = bbs[i], bbs[j]
        if min(a[3], b[3]) - max(a[0], b[0]) > 1e-6 and min(a[4], b[4]) - max(a[1], b[1]) > 1e-6 and min(a[5], b[5]) - max(a[2], b[2]) > 1e-6:
            npairs += 1
            c = BRepAlgoAPI_Common(lumps[i].wrapped, lumps[j].wrapped)
            worst = max(worst, vol(cq.Shape.cast(c.Shape())) if c.IsDone() else 1e9)
res["overlap_METAL_vs_dielectrics_mm3"] = ov
res["overlap_lump_pairs_checked"] = npairs
res["overlap_lump_worst_mm3"] = worst
res["NO_OVERLAP"] = "NO_OVERLAP" if all(v <= 1e-9 for v in ov.values()) and worst <= 1e-9 else "OVERLAP"
# AIR: one box, 6 faces
air = cqs(names["AIR"])
res["AIR_faces"] = len(air.Faces())
res["AIR_solids"] = len(air.Solids())
# sheets: touch METAL_PEC, zero overlap area, mid-thickness z, size 1.00 x 0.35, clearance between sheets
fails = []
zs = set()
for dn in D:
    sh = cqs(names[dn])
    bb = bbox(sh)
    zs.add(round(bb[2], 6))
    d = BRepExtrema_DistShapeShape(W(sh), W(names["METAL_PEC"]))
    d.Perform()
    c = BRepAlgoAPI_Common(W(sh), W(names["METAL_PEC"]))
    a = area(cq.Shape.cast(c.Shape())) if c.IsDone() else 1e9
    ar = area(sh)
    if not d.IsDone() or d.Value() > 1e-9 or a > 1e-9 or abs(ar - 0.35) > 1e-3:
        fails.append((dn, d.Value() if d.IsDone() else None, a, ar))
res["sheet_z_values"] = sorted(zs)
res["SHEETS"] = "SHEETS_OK" if not fails and zs == {0.0175} else "SHEETS_FAIL"
res["sheet_failures"] = fails[:20]
# PORT: touches pin and ground only (plus the dielectric faces it lies on), area
port = cqs(names["PORT"])
pbb = bbox(port)
touch = []
for n, s in names.items():
    if n == "PORT" or n in D:
        continue
    b = bbox(s)
    if min(pbb[3], b[3]) - max(pbb[0], b[0]) >= -1e-6 and min(pbb[4], b[4]) - max(pbb[1], b[1]) >= -1e-6 and min(pbb[5], b[5]) - max(pbb[2], b[2]) >= -1e-6:
        dd = BRepExtrema_DistShapeShape(W(port), W(s))
        dd.Perform()
        if dd.IsDone() and dd.Value() <= 1e-9:
            touch.append(n)
# which METAL lumps does the port touch
mt = []
pbbm = pbb
for k, l in enumerate(lumps):
    b = bbs[k]
    if min(pbb[3], b[3]) - max(pbb[0], b[0]) >= -1e-6 and min(pbb[4], b[4]) - max(pbb[1], b[1]) >= -1e-6 and min(pbb[5], b[5]) - max(pbb[2], b[2]) >= -1e-6:
        dd = BRepExtrema_DistShapeShape(port.wrapped, l.wrapped)
        dd.Perform()
        if dd.IsDone() and dd.Value() <= 1e-9:
            mt.append((round(bbox(l)[5] - bbox(l)[2], 3), round(max(abs(bbox(l)[0]), abs(bbox(l)[3])), 2)))
res["PORT"] = dict(touching_bodies=sorted(touch), metal_lumps_touched=mt, area=area(port), expected=math.pi * (2.05 ** 2 - 0.635 ** 2), z=pbb[2])
res["PORT_OK"] = bool("METAL_PEC" in touch and len(mt) == 2 and abs(res["PORT"]["area"] - res["PORT"]["expected"]) < 1e-3)
if len(sys.argv) > 4 and sys.argv[4] == "win":
    # window gate: no FR4 inside the 12 x 12 mm window; list every metal lump and via / pad whose footprint reaches into the window
    box = cq.Workplane("XY").workplane(offset=-5.0).rect(12.0, 12.0).extrude(5.0 - 3.175).val()
    c = BRepAlgoAPI_Common(W(names["SUBSTRATE_FR4"]), box.wrapped)
    res["window_FR4_common_volume_mm3"] = vol(cq.Shape.cast(c.Shape())) if c.IsDone() else float("nan")
    inw = []
    for k, l in enumerate(lumps):
        b = bbs[k]
        if min(6.0, b[3]) - max(-6.0, b[0]) > 1e-6 and min(6.0, b[4]) - max(-6.0, b[1]) > 1e-6:
            inw.append((round(b[2], 3), round(b[5], 3), round(b[3] - b[0], 2)))
    res["window_metal_lumps_reaching_in_xy(zmin,zmax,width)"] = inw
    res["WINDOW"] = "WINDOW_OK" if res["window_FR4_common_volume_mm3"] <= 1e-9 else "WINDOW_FAIL"
res["elapsed_s"] = round(time.time() - t0)
json.dump(res, open(outp, "w"), indent=1, default=str)
print(json.dumps({k: v for k, v in res.items() if k not in ("all_bodies_valid",)}, default=str)[:1800])
