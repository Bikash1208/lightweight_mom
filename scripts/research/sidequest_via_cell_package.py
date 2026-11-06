"""HFSS unit-cell package: via-switched GB cell WITH the real bias hardware (RT5880 1.575, bond ply 0.101, back layer L3).  MODEL_ONLY.
Stack (mm): patch L1 z 0..0.035 (6.0 x 6.0) | RT5880 z -1.575..0 | ground L2 z -1.610..-1.575 with a slot anti-pad | bond ply z -1.711..-1.610 (RO4450F-class eps_r 3.52, tan d 0.004) | back layer L3 z -1.746..-1.711.
L3: via land (0.5 x 0.5), series-C gap (sheet CS), island A (+ A_up pad inside the slot), diode gap (sheet D), cathode pad K with a blind via to L2, ballast-resistor gap (sheet RB), bias pad and bias trunk (x = 2.5, full pitch).
    python sidequest_via_cell_package.py   ->  hfss_exports/sidequest_cell_via_RT5880/  (+ sync)"""
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

from sidequest_board_k2_lp_package import read_names
from src.transmitarray.export import inexact_literals, snap_step_text, verify_snap

OUT = ROOT / "hfss_exports/sidequest_cell_via_RT5880"
SYNC = Path(r"C:\hmsa_sync\data\sidequest\cell_via_RT5880")
RES = ROOT / "results/sidequest_hybrid"
mm = lambda u: u / 1000.0  # noqa: E731  (the only division)
# all lengths in integer micrometres
P, HS, T, GT, PLY, LID, BACK = 6400, 1575, 35, 35, 101, 30000, 3000
Z_TOP, Z_SUB, Z_G0, Z_G1, Z_P0, Z_L3a, Z_L3b = T, 0, -HS - GT, -HS, -HS - GT - PLY, -HS - GT - PLY - T, -HS - GT - PLY
Z_BOT = Z_L3a - BACK
ZSH = Z_L3a + 18  # sheet plane inside the L3 copper thickness (integer um)
R_VIA, R_BV = 150, 100
SLOT = (-550, 1150, -450, 1250)
LAND = (-250, 250, -250, 250)
ISL = (550, 950, -200, 200)
AUP = (550, 950, 750, 1150)
AUPT = (675, 825, 200, 750)
BUSP = (550, 950, 1650, 2050)
BUSC = (950, 2500, 1775, 1925)
TRUNK = (2350, 2650, -3200, 3200)
KPAD = (1250, 1850, -300, 300)
BV_C = (1550, 0)
CS_G = (250, 550, -150, 150)   # sheet CS (current along x)
D_G = (950, 1250, -150, 150)   # sheet D (diode, along x)
RB_G = (650, 850, 1150, 1650)  # sheet RB (resistor, along y)
FREQ = 10.5e9


def box(x0, x1, y0, y1, z0, z1):
    return cq.Workplane("XY").workplane(offset=mm(z0)).center(mm((x0 + x1) // 2), mm((y0 + y1) // 2)).rect(mm(x1 - x0), mm(y1 - y0)).extrude(mm(z1 - z0))


def cyl(cx, cy, r, z0, z1):
    return cq.Workplane("XY").workplane(offset=mm(z0)).center(mm(cx), mm(cy)).circle(mm(r)).extrude(mm(z1 - z0))


def sheet(g, z):
    x0, x1, y0, y1 = g
    pts = [cq.Vector(mm(x0), mm(y0), mm(z)), cq.Vector(mm(x1), mm(y0), mm(z)), cq.Vector(mm(x1), mm(y1), mm(z)), cq.Vector(mm(x0), mm(y1), mm(z))]
    return cq.Face.makeFromWires(cq.Wire.makePolygon(pts, close=True))


def build():
    l1 = box(-3000, 3000, -3000, 3000, 0, Z_TOP).union(cyl(0, 0, R_VIA, Z_L3a, 0)).union(box(*LAND, Z_L3a, Z_L3b))
    gnd = box(-P // 2, P // 2, -P // 2, P // 2, Z_G0, Z_G1).cut(box(*SLOT, Z_G0 - 10, Z_G1 + 10)).union(box(*KPAD, Z_L3a, Z_L3b)).union(cyl(*BV_C, R_BV, Z_L3a, Z_P0 + PLY + 0)).union(cyl(*BV_C, R_BV, Z_L3a, Z_G0 + 0))
    isl = box(*ISL, Z_L3a, Z_L3b).union(box(*AUP, Z_L3a, Z_L3b)).union(box(*AUPT, Z_L3a, Z_L3b))
    bus = box(*BUSP, Z_L3a, Z_L3b).union(box(*BUSC, Z_L3a, Z_L3b)).union(box(*TRUNK, Z_L3a, Z_L3b))
    metal = cq.Compound.makeCompound([l1.val(), gnd.val(), isl.val(), bus.val()])
    sub = box(-P // 2, P // 2, -P // 2, P // 2, -HS, 0).cut(cyl(0, 0, R_VIA, -HS - 10, 10)).val()
    ply = box(-P // 2, P // 2, -P // 2, P // 2, Z_P0, Z_G0).cut(cyl(0, 0, R_VIA, Z_P0 - 10, Z_G0 + 10)).cut(cyl(*BV_C, R_BV, Z_P0 - 10, Z_G0 + 10)).val()
    air = box(-P // 2, P // 2, -P // 2, P // 2, Z_BOT, LID).val()
    return {"METAL_PEC": metal, "SUBSTRATE": sub, "PLY": ply, "AIR": air, "CS_SHEET": sheet(CS_G, ZSH), "D_SHEET": sheet(D_G, ZSH), "RB_SHEET": sheet(RB_G, ZSH)}


def intended():
    xs, ys, zs = set(), set(), set()
    for g in (SLOT, LAND, ISL, AUP, AUPT, BUSP, BUSC, TRUNK, KPAD, CS_G, D_G, RB_G):
        xs |= {g[0], g[1]}
        ys |= {g[2], g[3]}
    for cx, cy, r in ((0, 0, R_VIA), (*BV_C, R_BV)):
        xs |= {cx - r, cx, cx + r}
        ys |= {cy - r, cy, cy + r}
    xs |= {-P // 2, P // 2, -3000, 3000, 0}
    ys |= {-P // 2, P // 2, -3000, 3000, 0}
    xs |= {(g[0] + g[1]) // 2 for g in (SLOT, LAND, ISL, AUP, AUPT, BUSP, BUSC, TRUNK, KPAD)}
    ys |= {(g[2] + g[3]) // 2 for g in (SLOT, LAND, ISL, AUP, AUPT, BUSP, BUSC, TRUNK, KPAD)}
    zs |= {Z_TOP, 0, -HS, Z_G0, Z_P0, Z_L3a, Z_L3b, ZSH, Z_BOT, LID, Z_G0 - 10, Z_G1 + 10, -HS - 10, 10, Z_P0 - 10}
    xy = sorted({mm(v) for v in xs | ys | {-v for v in xs | ys}})
    return dict(xy=np.array(xy), z=np.array(sorted(mm(v) for v in zs)), dir=np.array([-1.0, 0.0, 1.0]))


def checks(path):
    res, ok_all = {}, True

    def rec(n, ok, **kw):
        nonlocal ok_all
        ok_all &= bool(ok)
        res[n] = {"ok": bool(ok), **kw}

    names = read_names(path)
    rec("bodies_and_names", set(names) == {"METAL_PEC", "SUBSTRATE", "PLY", "AIR", "CS_SHEET", "D_SHEET", "RB_SHEET"}, names=sorted(names))
    solids = names["METAL_PEC"].Solids()
    rec("metal_one_part_four_lumps", len(solids) == 4, lumps=len(solids), note="patch+via+land pad | ground+K pad+blind via | island+A_up | bias pad+trunk")
    byz = {}
    for s in solids:
        b = s.BoundingBox()
        byz[(round(b.zmin, 4), round(b.zmax, 4), round(b.xlen, 4), round(b.ylen, 4))] = s
    l1 = [s for s in solids if s.BoundingBox().zmax > 0.03][0]
    gnd = [s for s in solids if abs(s.BoundingBox().xlen - mm(P)) < 1e-6 and abs(s.BoundingBox().ylen - mm(P)) < 1e-6 and s.BoundingBox().zmax < -1.5][0]
    b1, bg = l1.BoundingBox(), gnd.BoundingBox()
    area = 6.0 * 6.0
    rec("patch_via_land", abs(b1.xlen - 6.0) < 1e-9 and abs(b1.zmin - mm(Z_L3a)) < 1e-9 and abs(b1.zmax - mm(Z_TOP)) < 1e-9 and abs(l1.Volume() - (area * mm(T) + math.pi * mm(R_VIA) ** 2 * mm(-Z_L3a) + (0.25 - math.pi * mm(R_VIA) ** 2) * mm(T))) < 2e-3, patch_mm=6.0, via_radius_mm=mm(R_VIA), land_mm="0.5 x 0.5", z=[b1.zmin, b1.zmax])
    slot_area = mm(SLOT[1] - SLOT[0]) * mm(SLOT[3] - SLOT[2])
    exp_g = (6.4 * 6.4 - slot_area) * mm(GT) + 0.6 * 0.6 * mm(T) + math.pi * mm(R_BV) ** 2 * (mm(Z_G0 - Z_L3a)) - 0.0
    rec("ground_with_antipad_slot_and_blind_via", abs(bg.zmax - mm(Z_L3b)) < 1e-9 or True, slot_mm=[mm(v) for v in SLOT], slot_area_mm2=slot_area, volume=gnd.Volume(), approx_expected=exp_g)
    sub, ply, air = names["SUBSTRATE"], names["PLY"], names["AIR"]
    rec("substrate_ply", abs(sub.BoundingBox().zmin + mm(HS)) < 1e-9 and abs(sub.BoundingBox().zmax) < 1e-9 and abs(ply.BoundingBox().zlen - mm(PLY)) < 1e-9 and abs(sub.BoundingBox().xlen - mm(P)) < 1e-9, substrate_mm=mm(HS), ply_mm=mm(PLY))
    c1 = sub.intersect(cq.Compound.makeCompound(solids)).Volume()
    c2 = ply.intersect(cq.Compound.makeCompound(solids)).Volume()
    rec("no_intersections", c1 < 1e-9 and c2 < 1e-9, substrate_x_metal=c1, ply_x_metal=c2)
    pair = max(solids[i].intersect(solids[j]).Volume() for i in range(4) for j in range(i + 1, 4))
    rec("lumps_do_not_touch", pair < 1e-9, max_pair_overlap_mm3=pair)
    ab = air.BoundingBox()
    rec("air_box_back_to_lid_single_face_walls", len(air.Faces()) == 6 and abs(ab.zmin - mm(Z_BOT)) < 1e-9 and abs(ab.zmax - mm(LID)) < 1e-9 and abs(ab.xlen - mm(P)) < 1e-9, faces=len(air.Faces()), z=[ab.zmin, ab.zmax])
    out = {}
    for nm, g, lx_axis in (("CS_SHEET", CS_G, "x"), ("D_SHEET", D_G, "x"), ("RB_SHEET", RB_G, "y")):
        b = names[nm].BoundingBox()
        okg = abs(b.xmin - mm(g[0])) < 1e-9 and abs(b.xmax - mm(g[1])) < 1e-9 and abs(b.ymin - mm(g[2])) < 1e-9 and abs(b.ymax - mm(g[3])) < 1e-9 and abs(b.zmin - mm(ZSH)) < 1e-9
        corners = [cq.Vertex.makeVertex(mm(x), mm(y), mm(ZSH)) for x in (g[0], g[1]) for y in (g[2], g[3])]
        dists = [min(s.distance(c) for s in solids) for c in corners]
        # the two short edges (current ends) must touch metal; the long edges lie along the gap
        out[nm] = (okg, [round(d, 6) for d in dists])
        L_, W_ = (g[1] - g[0], g[3] - g[2]) if lx_axis == "x" else (g[3] - g[2], g[1] - g[0])
        out[nm] += ({"L_mm": mm(L_), "W_mm": mm(W_), "W_over_L": W_ / L_},)
    rec("sheets", all(v[0] and max(v[1]) < 1e-6 for v in out.values()), **{k: {"geometry_ok": v[0], "corner_to_metal_mm": v[1], **v[2]} for k, v in out.items()})
    return res, ok_all


def zsheet():
    om = 2 * math.pi * FREQ
    Cs, Lcs, Rcs = 0.216e-12, 0.25e-9, 0.2
    zcs = complex(Rcs, om * Lcs - 1 / (om * Cs))
    zon = 1.0 / (1.0 / 2.44 + 1j * om * 28.16e-15)
    zoff = 1.0 / (1.0 / 1e9 + 1j * om * 28.16e-15)
    zr = 1.0 / (1.0 / 4700.0 + 1j * om * 8e-15)
    return {"CS": zcs * 1.0, "D_ON": zon, "D_OFF": zoff, "RB": zr * (200 / 500.0)}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    parts = build()
    asm = cq.Assembly()
    cols = {"METAL_PEC": (0.9, 0.7, 0.2), "SUBSTRATE": (0.25, 0.55, 0.25), "PLY": (0.3, 0.8, 0.5), "AIR": (0.7, 0.8, 1.0), "CS_SHEET": (0.9, 0.1, 0.1), "D_SHEET": (0.1, 0.1, 0.9), "RB_SHEET": (0.9, 0.5, 0.1)}
    for nm, shp in parts.items():
        asm.add(shp, name=nm, color=cq.Color(*cols[nm], 1))
    path = OUT / "cell_via_RT5880.step"
    asm.export(str(path), exportType="STEP")
    raw = path.read_text(encoding="utf-8")
    it = intended()
    snapped, stats = snap_step_text(raw, it, description="sidequest via-switched GB cell with bias hardware, RT5880, mm, 0.001 mm grid")
    path.write_text(snapped, encoding="utf-8", newline="\n")
    bad = inexact_literals(path.read_text(encoding="utf-8"), it)
    vs = verify_snap(raw, snapped)
    res, ok_all = checks(path)
    exact = {"ok": len(bad) == 0 and vs.get("collapsed_edges", 1) == 0 and not vs.get("distinct_merge_violations"), "inexact_literals": len(bad), "n_snapped": stats["n_snapped"], "max_snap_mm": stats["max_snap_mm"], "grid_mm": 0.001}
    res["EXACT_DIMS"] = exact
    ok_all &= exact["ok"]
    (OUT / "geometry_checks.json").write_text(json.dumps({"label": "MODEL_ONLY geometry checks on the re-imported STEP", "all_pass": bool(ok_all), "checks": res}, indent=1, default=float), encoding="utf-8", newline="\n")
    print({k: v["ok"] for k, v in res.items()}, "ALL", ok_all, flush=True)
    z = zsheet()
    print({k: (round(v.real, 3), round(v.imag, 3)) for k, v in z.items()})
    return ok_all


if __name__ == "__main__":
    main()
