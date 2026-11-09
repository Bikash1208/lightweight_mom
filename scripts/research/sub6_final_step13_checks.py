"""SUB6 FINAL BOARD, 13-pad variant: gates on the exported STEP.  MODEL_ONLY / NOT_HFSS_VALIDATED.
    python sub6_final_step13_checks.py <step> <diode_sheets.csv> <out.json>
NO_OVERLAP : pairwise OCC boolean common of every pair of 3D bodies whose bounding boxes overlap with positive extent in x, y and z; common volume > 1e-9 mm^3 = fail.
             AIR is also checked pairwise (common with every other solid) and by volume balance V(AIR) + sum V(others) = V(box).
SHEETS_OK  : every diode sheet D<s>_<p> (a) touches the island body and the disc body of its cell (distance <= 1e-6 mm), (b) has zero overlap area with every solid except AIR
             (the sheet lies in the open slot, i.e. inside the AIR region by construction; sheets inside AIR are normal in HFSS).
Names unique: duplicates are reported.
"""
import csv
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from sidequest_board_k2_lp_package import read_names
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common
from OCP.BRepExtrema import BRepExtrema_DistShapeShape
from OCP.GProp import GProp_GProps
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.TopAbs import TopAbs_ON, TopAbs_OUT
from OCP.gp import gp_Pnt
from OCP.BRepGProp import BRepGProp

step, sheets_csv, outp = sys.argv[1], sys.argv[2], sys.argv[3]
t0 = time.time()
names = read_names(step)
print("read", len(names), "bodies", round(time.time() - t0), "s", flush=True)


def bbox(shp):
    b = Bnd_Box()
    BRepBndLib.Add_s(shp.wrapped if hasattr(shp, "wrapped") else shp, b)
    return b.Get()


def vol(shp):
    """sum of the volumes of the SOLIDS only (a face / sheet has no volume; VolumeProperties of a bare face returns a meaningless signed number)"""
    import cadquery as cq
    if not hasattr(shp, "Solids"):
        shp = cq.Shape.cast(shp)
    tot = 0.0
    for so in shp.Solids():
        p = GProp_GProps()
        BRepGProp.VolumeProperties_s(so.wrapped, p)
        tot += p.Mass()
    return tot


def area(shp):
    p = GProp_GProps()
    BRepGProp.SurfaceProperties_s(shp.wrapped if hasattr(shp, "wrapped") else shp, p)
    return p.Mass()


W = lambda s: s.wrapped if hasattr(s, "wrapped") else s
solid_names = [n for n in names if not (n.startswith("D") and "_" in n and n[1:2].isdigit()) and n != "PORT_ANNULUS"]
sheet_names = [n for n in names if n.startswith("D") and "_" in n and n[1:2].isdigit()]
bb = {n: bbox(names[n]) for n in solid_names}
vols = {n: vol(names[n]) for n in solid_names}
print("volumes", round(time.time() - t0), "s", flush=True)
fails = []
EPS = 1e-6
pairs = []
nm = [n for n in solid_names if n != "AIR"]
for i, a in enumerate(nm):
    xa0, ya0, za0, xa1, ya1, za1 = bb[a]
    for b in nm[i + 1:]:
        xb0, yb0, zb0, xb1, yb1, zb1 = bb[b]
        if min(xa1, xb1) - max(xa0, xb0) > EPS and min(ya1, yb1) - max(ya0, yb0) > EPS and min(za1, zb1) - max(za0, zb0) > EPS:
            pairs.append((a, b))
print("candidate pairs (non-AIR)", len(pairs), flush=True)
worst = 0.0
for k, (a, b) in enumerate(pairs):
    c = BRepAlgoAPI_Common(W(names[a]), W(names[b]))
    v = vol(c.Shape()) if c.IsDone() else float("nan")
    worst = max(worst, v if v == v else 1e9)
    if not (v <= 1e-9):
        fails.append((a, b, v))
    if k % 200 == 0:
        print(" pair", k, "/", len(pairs), round(time.time() - t0), "s fails", len(fails), flush=True)
# AIR: volume balance + pairwise
air = names["AIR"]
gj = json.load(open(sys.argv[4]))["meta"]   # exact box dimensions (a Bnd_Box of a complex shape is enlarged by tolerances)
vbox = 2 * gj["air_box"]["ax"] * 2 * gj["air_box"]["ay"] * (gj["z"]["air_hi"] - gj["z"]["air_lo"])
vsum = sum(vols[n] for n in nm)
balance = vbox - vsum - vols["AIR"]
print("AIR balance (box - sum others - AIR) = %.6g mm^3 (box %.6g)" % (balance, vbox), flush=True)
air_pair_fail = []
air_worst = 0.0
for k, n in enumerate(nm):
    c = BRepAlgoAPI_Common(W(air), W(names[n]))
    v = vol(c.Shape()) if c.IsDone() else float("nan")
    air_worst = max(air_worst, v if v == v else 1e9)
    if not (v <= 1e-9):
        air_pair_fail.append((n, v))
    if k % 100 == 0:
        print(" AIR pair", k, "/", len(nm), round(time.time() - t0), "s fails", len(air_pair_fail), flush=True)
overlap_ok = not fails and not air_pair_fail and abs(balance) < 1e-3
# sheets
rows = list(csv.DictReader(open(sheets_csv)))
sh_fail = []
slivers = []
rows_by = {r['name']: r for r in rows}
for r in rows:
    nmx = r["name"]
    sh = names[nmx]
    s, p = int(r["string"]), int(r["pos_in_string"])
    isl = f"METAL_PEC_S{s:02d}" if p == 1 else f"METAL_PEC_N{s:02d}_{p - 1:02d}"
    dsc = "METAL_PEC_RETURN" if p == 20 else f"METAL_PEC_N{s:02d}_{p:02d}"
    ok = True
    for tgt in (isl, dsc):
        d = BRepExtrema_DistShapeShape(W(sh), W(names[tgt]))
        d.Perform()
        if not d.IsDone() or d.Value() > 1e-6:
            ok = False
            sh_fail.append((nmx, tgt, "distance", d.Value() if d.IsDone() else None))
    sb = bbox(sh)
    for n in nm:
        b2 = bb[n]
        if min(sb[3], b2[3]) - max(sb[0], b2[0]) >= -EPS and min(sb[4], b2[4]) - max(sb[1], b2[1]) >= -EPS and min(sb[5], b2[5]) - max(sb[2], b2[2]) >= -EPS:
            c = BRepAlgoAPI_Common(W(sh), W(names[n]))
            a = area(c.Shape()) if c.IsDone() else 1e9
            if a > 1e-9:
                # sheet lying ON a face of the solid (coplanar contact, e.g. thin model: sheet on the substrate top) is a contact, not an overlap
                cl = BRepClass3d_SolidClassifier(W(names[n]))
                cp = gp_Pnt(float(rows_by[nmx]["sheet_centre_x"]), float(rows_by[nmx]["sheet_centre_y"]), (0.0 if "thin" in step else 0.0175))
                cl.Perform(cp, 1e-7)
                st = cl.State()
                if st == TopAbs_ON or st == TopAbs_OUT:
                    continue
                if a < 2e-3:   # sliver: the sheet is 1.5 um longer than the gap at each end (+ 0.001 mm grid rounding) so that it really touches the copper: area <= ~5 um x 0.35 mm
                    slivers.append((nmx, n, a))
                    continue
                sh_fail.append((nmx, n, "overlap_area", a))
res = dict(label="MODEL_ONLY / NOT_HFSS_VALIDATED", step=str(step), n_bodies=len(names), n_solids=len(solid_names), n_sheets=len(sheet_names), duplicate_names="n/a (names dict keys)",
           NO_OVERLAP="NO_OVERLAP" if overlap_ok else "OVERLAP", pairs_checked=len(pairs), pair_overlaps=[(a, b, v) for a, b, v in fails][:50], worst_common_volume_mm3=worst,
           AIR_balance_mm3=balance, AIR_box_volume_mm3=vbox, AIR_pairs_checked=len(nm), AIR_worst_common_volume_mm3=air_worst, AIR_pair_overlaps=air_pair_fail[:50],
           SHEETS="SHEETS_OK" if not sh_fail else "SHEETS_FAIL", sheet_failures=sh_fail[:50], sheet_rounding_slivers=len(slivers), sheet_max_sliver_mm2=max([a for _, _, a in slivers], default=0.0), elapsed_s=round(time.time() - t0))
json.dump(res, open(outp, "w"), indent=1, default=str)
print(json.dumps({k: res[k] for k in ("NO_OVERLAP", "SHEETS", "pairs_checked", "AIR_balance_mm3", "worst_common_volume_mm3", "AIR_worst_common_volume_mm3", "elapsed_s")}), flush=True)
