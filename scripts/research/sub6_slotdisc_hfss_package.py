"""SUB6 slot-disc HFSS EIGENMODE package (MODEL_ONLY predictions, written BEFORE any HFSS run).  EXPORT STANDARD (hfss_exports/EXPORT_STANDARD.md).

Cells (RT5880 3.175 mm, p 11.6, slot along x, notch 1.00 mm at the +x slot end through the rim to the disc edge):
  CELL1 gap 0.73, L/2r 0.8, W 2.9  (r 5.435, L 8.696)      CELL2 gap 1.00, L/2r 0.9, W 2.9  (r 5.300, L 9.540)
Geometry files (4): CELL{1,2}_notch.step (OFF and ON: same geometry, different RLC on DIODE_SHEET) and CELL{1,2}_nocut.step (closed slot, no notch, no sheet).
Objects: METAL_PEC (patch + ground), SUBSTRATE, AIR (ground bottom to lid, one box, single face per side wall, overlaps on purpose), DIODE_SHEET (notch files only).
    python sub6_slotdisc_hfss_package.py geometry      -> STEP files + geometry_checks.json
    python sub6_slotdisc_hfss_package.py docs          -> states.csv, predictions.json, HFSS_SETUP.md, results_template.csv (needs hfss_pred/pred_table.json)
"""
import csv
import json
import math
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np

OUT = ROOT / "hfss_exports" / "sub6_slotdisc_cell"
RES = ROOT / "results" / "sub6_slotdisc"
P, H, T, GT, LID = 11.6, 3.175, 0.035, 0.035, 150.0
EPS_R, TAND = 2.2, 0.0009
GC, ALONG = 1.00, 0.35
F0 = 5.8e9
CELLS = {
    "CELL1": dict(gap=0.73, ratio=0.8, W=2.9),
    "CELL2": dict(gap=1.0, ratio=0.9, W=2.9),
}
for _n, _c in CELLS.items():
    _c["r"] = round((P - _c["gap"]) / 2, 3)
    _c["L"] = round(_c["ratio"] * 2 * _c["r"], 3)
    _c["diode_x"] = round(0.5 * (_c["L"] / 2 + _c["r"]), 3)  # MoM: 4.8915 for CELL1 -> 4.891 (0.5 um shift, on the 0.001 mm grid)
C0 = 299792458.0
EPS0 = 8.8541878128e-12


# ------------------------------------------------------------------------------------------------ RLC equivalents
def z_off(f, C=0.15e-12, R=5e3, L=0.6e-9):
    w = 2 * math.pi * f
    return 1j * w * L + 1.0 / (1.0 / R + 1j * w * C)


def z_on(f, R=3.0, L=0.6e-9):
    return R + 1j * 2 * math.pi * f * L


def parallel_equiv(z, f):
    """Y = 1/Z = G + jB -> R_p = 1/G; B > 0: C = B/w ; B < 0: L = -1/(w B)."""
    y = 1.0 / z
    w = 2 * math.pi * f
    out = dict(R_ohm=1.0 / y.real, Z_series=[z.real, z.imag])
    if y.imag > 0:
        out.update(C_pF=y.imag / w * 1e12, L_nH=None)
    else:
        out.update(C_pF=None, L_nH=-1.0 / (w * y.imag) * 1e9)
    return out


# ------------------------------------------------------------------------------------------------ geometry
def make(name, notch):
    import cadquery as cq

    c = CELLS[name]
    r, L, W = c["r"], c["L"], c["W"]

    def box(x0, x1, y0, y1, z0, z1):
        return cq.Workplane("XY").workplane(offset=z0).center((x0 + x1) / 2, (y0 + y1) / 2).rect(x1 - x0, y1 - y0).extrude(z1 - z0)

    patch = cq.Workplane("XY").circle(r).extrude(T).cut(box(-L / 2, L / 2, -W / 2, W / 2, -0.01, T + 0.01))
    if notch:
        patch = patch.cut(box(L / 2 - 0.0, r + 1.0, -GC / 2, GC / 2, -0.01, T + 0.01))
        # the cut starts at the slot end x = L/2 exactly, so the notch opens the slot to the rim
    ground = box(-P / 2, P / 2, -P / 2, P / 2, -H - GT, -H)
    metal = cq.Compound.makeCompound([patch.val(), ground.val()])
    sub = box(-P / 2, P / 2, -P / 2, P / 2, -H, 0.0)
    air = box(-P / 2, P / 2, -P / 2, P / 2, -H - GT, LID)
    parts = {"METAL_PEC": metal, "SUBSTRATE": sub.val(), "AIR": air.val()}
    if notch:
        x0, x1, z = c["diode_x"] - ALONG / 2, c["diode_x"] + ALONG / 2, T / 2
        pts = [cq.Vector(x0, -GC / 2, z), cq.Vector(x1, -GC / 2, z), cq.Vector(x1, GC / 2, z), cq.Vector(x0, GC / 2, z)]
        parts["DIODE_SHEET"] = cq.Face.makeFromWires(cq.Wire.makePolygon(pts, close=True))
    return parts


def intended(name):
    c = CELLS[name]
    xy = {0.0, P / 2, c["r"], c["L"] / 2, c["W"] / 2, GC / 2, c["diode_x"] - ALONG / 2, c["diode_x"] + ALONG / 2, c["diode_x"]}
    xy |= {-v for v in xy}
    z = {0.0, T, T / 2, -H, -H - GT, LID, -0.01, T + 0.01}
    return dict(xy=np.array(sorted(round(v, 6) for v in xy)), z=np.array(sorted(round(v, 6) for v in z)), dir=np.array([-1.0, 0.0, 1.0]))


def inside(solid, x, y, z):
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.gp import gp_Pnt
    from OCP.TopAbs import TopAbs_IN

    return BRepClass3d_SolidClassifier(solid.wrapped, gp_Pnt(float(x), float(y), float(z)), 1e-7).State() == TopAbs_IN


def checks(path, name, notch):
    import cadquery as cq
    from sidequest_board_k2_lp_package import read_names
    from sub6_slotdisc_lib import SlotDiscNotch

    c = CELLS[name]
    r, L, W, dx = c["r"], c["L"], c["W"], c["diode_x"]
    res, ok_all = {}, True

    def rec(n, ok, **kw):
        nonlocal ok_all
        ok_all &= bool(ok)
        res[n] = {"ok": bool(ok), **kw}

    names = read_names(path)
    want = {"METAL_PEC", "SUBSTRATE", "AIR"} | ({"DIODE_SHEET"} if notch else set())
    rec("bodies_and_names", set(names) == want, names=sorted(names))
    solids = names["METAL_PEC"].Solids()
    rec("metal_two_parts_patch_ground", len(solids) == 2, parts=len(solids))
    patch = [s for s in solids if s.BoundingBox().zmax > -1.0][0]
    gnd = [s for s in solids if s.BoundingBox().zmax < -1.0][0]
    sub, air = names["SUBSTRATE"], names["AIR"]
    # patch volume against the analytic area
    yy = np.linspace(0, GC / 2, 20001)
    notch_area = 2 * np.trapezoid(np.sqrt(r * r - yy * yy) - L / 2, yy) if notch else 0.0
    area = math.pi * r * r - L * W - notch_area
    rec("patch_volume_vs_analytic", abs(patch.Volume() - area * T) / (area * T) < 2e-4, volume_mm3=patch.Volume(), analytic_mm3=area * T, rel_err=patch.Volume() / (area * T) - 1)
    # MoM mask (SlotDiscNotch, OPEN notch, i.e. the geometry HFSS sees) vs the solid on a 0.1 mm grid; exclude +-0.06 mm bands at the disc edge, slot walls, notch walls
    geo = SlotDiscNotch(r * 1e-3, L * 1e-3, W * 1e-3, GC * 1e-3 if notch else 0.0, bridge=False)
    bad = tot = 0
    for x in np.arange(-r - 0.5, r + 0.5, 0.1):
        for y in np.arange(-r - 0.5, r + 0.5, 0.1):
            rad = math.hypot(x, y)
            if abs(rad - r) < 0.06:
                continue
            if abs(abs(x) - L / 2) < 0.06 and abs(y) <= W / 2 + 0.06:
                continue
            if abs(abs(y) - W / 2) < 0.06 and abs(x) <= L / 2 + 0.06:
                continue
            if notch and (abs(abs(y) - GC / 2) < 0.06 and x > L / 2 - 0.06 or abs(x - L / 2) < 0.06 and abs(y) <= GC / 2 + 0.06):
                continue
            tot += 1
            if bool(geo.in_feature(x * 1e-3, y * 1e-3)) != inside(patch, x, y, T / 2):
                bad += 1
    rec("patch_vs_MoM_mask", bad == 0, points=tot, mismatches=bad)
    # symmetry y -> -y
    rng = np.random.default_rng(1)
    pts = rng.uniform(-r, r, size=(400, 2))
    asym = sum(inside(patch, x, y, T / 2) != inside(patch, x, -y, T / 2) for x, y in pts)
    rec("mirror_symmetry_y", asym == 0, asymmetric_points=int(asym), note="mirror x is broken on purpose by the notch" if notch else "")
    ov = {"patch_x_ground": patch.intersect(gnd).Volume(), "metal_x_substrate": sub.intersect(cq.Compound.makeCompound(solids)).Volume()}
    rec("no_intersections", max(ov.values()) < 1e-9, **ov)
    ab = air.BoundingBox()
    allin = all(ab.xmin - 1e-9 <= s.BoundingBox().xmin and s.BoundingBox().xmax <= ab.xmax + 1e-9 and ab.zmin - 1e-9 <= s.BoundingBox().zmin and s.BoundingBox().zmax <= ab.zmax + 1e-9 for s in list(solids) + [sub])
    rec("air_box_ground_to_lid_single_face_walls", len(air.Faces()) == 6 and abs(ab.zmin - (-H - GT)) < 1e-9 and abs(ab.zmax - LID) < 1e-9 and allin, faces=len(air.Faces()), z=[ab.zmin, ab.zmax], contains_all_solids=allin)
    rec("substrate_stack", abs(sub.BoundingBox().zmin + H) < 1e-9 and abs(sub.BoundingBox().zmax) < 1e-9 and abs(patch.BoundingBox().zmax - T) < 1e-9 and abs(gnd.BoundingBox().zmax + H) < 1e-9, substrate_mm=H, patch_t_mm=T, ground_t_mm=GT)
    rim_end = r - L / 2
    rim_corner = r - math.hypot(L / 2, W / 2)
    rec("rims", rim_end >= 0.35 and rim_corner >= 0.2, slot_end_rim_mm=rim_end, slot_corner_rim_mm=rim_corner, note="grid rules: rim >= 0.35 (r - L/2), corner rim >= 0.2")
    if notch:
        b = names["DIODE_SHEET"].BoundingBox()
        x0, x1 = dx - ALONG / 2, dx + ALONG / 2
        geom_ok = abs(b.xmin - x0) < 1e-9 and abs(b.xmax - x1) < 1e-9 and abs(b.ymin + GC / 2) < 1e-9 and abs(b.ymax - GC / 2) < 1e-9 and abs(b.zmin - T / 2) < 1e-9 and abs(b.zmax - T / 2) < 1e-9
        corners = [cq.Vertex.makeVertex(x, y, T / 2) for x in (x0, x1) for y in (-GC / 2, GC / 2)]
        mids = [cq.Vertex.makeVertex(dx, y, T / 2) for y in (-GC / 2, GC / 2)]
        d_corner = [round(patch.distance(v), 9) for v in corners]
        d_mid = [round(patch.distance(v), 9) for v in mids]
        centre_in_metal = inside(patch, dx, 0.0, T / 2)
        rec("sheet_size_position", geom_ok, across_y_mm=GC, along_x_mm=ALONG, centre_x_mm=dx, z_mm=T / 2, note="centred in the rim between slot end and disc edge; CELL1 dx rounded 4.8915 -> 4.891 (0.5 um)")
        rec("sheet_edges_on_notch_walls", max(d_corner + d_mid) < 1e-6 and not centre_in_metal, corner_to_metal_mm=d_corner, long_edge_mid_to_metal_mm=d_mid, centre_in_metal=bool(centre_in_metal))
        xr = math.sqrt(r * r - (GC / 2) ** 2)
        rec("sheet_inside_notch_clearances", (x0 - L / 2) > 0 and (xr - x1) > 0, to_slot_end_mm=x0 - L / 2, to_disc_edge_at_y_half_gc_mm=xr - x1, to_disc_edge_at_y0_mm=r - x1,
            note="informational: the sheet ends face the open slot / open air, not metal; metal touches only the long edges (notch walls). Under 0.1 mm when the rim is narrow (CELL2)")
        rec("notch_width", abs(GC - 1.0) < 1e-12, notch_width_mm=GC, notch_x_from_to_mm=[L / 2, xr], bar_each_side_of_notch_mm=(W - GC) / 2)
    return res, ok_all


def geometry():
    import cadquery as cq
    from src.transmitarray.export import inexact_literals, snap_step_text, verify_snap

    OUT.mkdir(parents=True, exist_ok=True)
    allres, allok = {}, True
    for name in CELLS:
        for tag, notch in (("notch", True), ("nocut", False)):
            parts = make(name, notch)
            asm = cq.Assembly()
            cols = {"METAL_PEC": (0.9, 0.7, 0.2), "SUBSTRATE": (0.25, 0.55, 0.25), "AIR": (0.7, 0.8, 1.0), "DIODE_SHEET": (0.9, 0.1, 0.1)}
            for nm, shp in parts.items():
                asm.add(shp, name=nm, color=cq.Color(*cols[nm], 1))
            path = OUT / f"{name}_{tag}.step"
            asm.export(str(path), exportType="STEP")
            raw = path.read_text(encoding="utf-8")
            it = intended(name)
            snapped, stats = snap_step_text(raw, it, description=f"sub6 slot-disc {name} {tag}, RT5880 3.175, mm, 0.001 mm grid")
            path.write_text(snapped, encoding="utf-8", newline="\n")
            badl = inexact_literals(path.read_text(encoding="utf-8"), it)
            vs = verify_snap(raw, snapped)
            res, ok = checks(path, name, notch)
            exact = {"ok": len(badl) == 0 and vs.get("collapsed_edges", 1) == 0 and not vs.get("distinct_merge_violations"), "inexact_literals": len(badl), "n_snapped": stats["n_snapped"], "max_snap_mm": stats["max_snap_mm"], "grid_mm": 0.001}
            res["EXACT_DIMS"] = exact
            ok &= exact["ok"]
            allres[f"{name}_{tag}"] = {"all_pass": bool(ok), "checks": res}
            allok &= ok
            print(name, tag, ok, {k: v["ok"] for k, v in res.items()}, flush=True)
    meta = dict(P_mm=P, substrate_mm=H, patch_t_mm=T, ground_t_mm=GT, lid_mm=LID, notch_mm=GC, sheet_across_mm=GC, sheet_along_mm=ALONG, cells=CELLS)
    (OUT / "geometry_checks.json").write_text(json.dumps({"label": "MODEL_ONLY geometry checks on the re-imported STEP files", "all_pass": bool(allok), "parameters": meta, "files": allres}, indent=1, default=float), encoding="utf-8", newline="\n")
    print("ALL", allok)
    return allok


if __name__ == "__main__":
    if sys.argv[1] == "geometry":
        sys.exit(0 if geometry() else 1)
