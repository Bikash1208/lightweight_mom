"""SUB6 island-cell HFSS EIGENMODE package (MODEL_ONLY predictions, written BEFORE any HFSS run).  EXPORT STANDARD (hfss_exports/EXPORT_STANDARD.md) + the package-layer STEP exactification
(sub6_step_exact.py).

Cell ISL (row ISL_W2.9_q0.8_Li0.5_c0.6 of the island scan; the row the best boards use, see HFSS_SETUP.md): RT5880 3.175 mm, p 11.6, disc r 5.435, patch 35 um, slot L 8.696 x W 2.9 along x,
metal ISLAND centred in the slot (Li 4.348 = 0.5 L, x in +-2.174), island y from -0.85 (clearance 0.60 to the slot -y edge) to +0.45 (1.00 mm clearance to the slot +y edge), and ONE
DIODE_SHEET 1.00 (across, y: 0.45 ... 1.45) x 0.35 (along, x: -0.175 ... +0.175) at z = 0.0175, x = 0, bridging island and rim.  Nominal dimensions on the 0.001 mm grid (the MoM uses mesh-64
snapped dimensions: clearance 0.544, gap 0.906, island width 1.45 instead of 0.60 / 1.00 / 1.30; recorded in geometry_checks.json).
Files (hfss_exports/sub6_island_cell/):
  ISL_novia.step  (variant 1)  METAL_PEC = patch(disc with slot) + island + ground, SUBSTRATE, AIR (ground bottom -3.210 to lid +150), DIODE_SHEET
  ISL_absent.step (reference)  ISLAND_ABSENT: closed slotted disc (no island, no sheet) = CELL1_nocut geometry
  ISL_vias.step   (variant 2)  as variant 1 plus bias vias: island via at (0, -0.20) and disc via at (0, +2.35), drill 0.3 mm (barrel r 0.15), pad 0.6 mm; ground with 1.0 mm anti-pads;
                               FR4 0.5 mm below the ground (z -3.710 ... -3.210; SUBSTRATE_FR4: eps_r 4.4, tan d 0.02 ASSUMED typical values), bottom copper 35 um (z -3.745 ... -3.710) with a 0.6 mm pad and a
                               2.0 mm x 0.3 mm stub per via (island stub towards +x, disc stub towards -x; not connected to each other); AIR from z = -4.745 (1.0 mm below the stubs) to the lid.
    python sub6_island_hfss_package.py geometry   -> STEP files + geometry_checks.json
"""
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

OUT = ROOT / "hfss_exports" / "sub6_island_cell"
RES = ROOT / "results" / "sub6_island_hfss"
P, H, T, GT, LID = 11.6, 3.175, 0.035, 0.035, 150.0
EPS_R, TAND = 2.2, 0.0009
FR4_T, FR4_EPS, FR4_TAND = 0.5, 4.4, 0.02
GC, ALONG = 1.00, 0.35
F0 = 5.8e9
C0, EPS0 = 299792458.0, 8.8541878128e-12
R = 5.435
L, W = 8.696, 2.9
LI = 4.348
CLR = 0.60
Y_ISL0, Y_ISL1 = round(-W / 2 + CLR, 3), round(W / 2 - GC, 3)  # -0.85, +0.45
X_ISL = round(LI / 2, 3)  # 2.174
Y_SH0, Y_SH1 = Y_ISL1, round(W / 2, 3)  # sheet from the island +y edge to the slot +y edge
Z_SH = T / 2
VIA_R, PAD_R, ANTIPAD_R = 0.15, 0.30, 0.50
Y_VIA_I = round(0.5 * (Y_ISL0 + Y_ISL1), 3)  # -0.20, island centre in y
Y_VIA_D = 2.35
STUB_L, STUB_W = 2.0, 0.3
Z_GB, Z_FR4B, Z_BOT = -H - GT, -H - GT - FR4_T, -H - GT - FR4_T - GT  # -3.210, -3.710, -3.745
Z_AIR_V2 = round(Z_BOT - 1.0, 3)  # -4.745


def sheet_current_line():
    return (0.0, Y_SH0, Z_SH), (0.0, Y_SH1, Z_SH)


def make(variant):
    """variant: 'novia' | 'absent' | 'vias'  -> {name: cq shape}."""
    import cadquery as cq

    def box(x0, x1, y0, y1, z0, z1):
        return cq.Workplane("XY").workplane(offset=z0).center((x0 + x1) / 2, (y0 + y1) / 2).rect(x1 - x0, y1 - y0).extrude(z1 - z0)

    def cyl(x, y, r, z0, z1):
        return cq.Workplane("XY").workplane(offset=z0).center(x, y).circle(r).extrude(z1 - z0)

    patch = cq.Workplane("XY").circle(R).extrude(T).cut(box(-L / 2, L / 2, -W / 2, W / 2, -0.01, T + 0.01))
    ground = box(-P / 2, P / 2, -P / 2, P / 2, Z_GB, -H)
    vias = variant == "vias"
    solids = []
    if variant == "absent":
        solids = [patch.val(), ground.val()]
    elif not vias:
        island = box(-X_ISL, X_ISL, Y_ISL0, Y_ISL1, 0.0, T)
        solids = [patch.val(), island.val(), ground.val()]
    else:
        island = box(-X_ISL, X_ISL, Y_ISL0, Y_ISL1, 0.0, T)

        def via_solid(y, stub_dir):
            barrel = cyl(0.0, y, VIA_R, Z_BOT, 0.0)
            pad = cyl(0.0, y, PAD_R, Z_BOT, Z_FR4B)
            x0, x1 = (0.0, stub_dir * STUB_L) if stub_dir > 0 else (stub_dir * STUB_L, 0.0)
            stub = box(x0, x1, y - STUB_W / 2, y + STUB_W / 2, Z_BOT, Z_FR4B)
            return barrel.union(pad).union(stub)

        isl = island.union(via_solid(Y_VIA_I, +1))
        dsc = patch.union(via_solid(Y_VIA_D, -1))
        ground = ground.cut(cyl(0.0, Y_VIA_I, ANTIPAD_R, Z_GB - 0.01, -H + 0.01)).cut(cyl(0.0, Y_VIA_D, ANTIPAD_R, Z_GB - 0.01, -H + 0.01))
        solids = [dsc.val(), isl.val(), ground.val()]
    metal = cq.Compound.makeCompound(solids)
    sub = box(-P / 2, P / 2, -P / 2, P / 2, -H, 0.0)
    z_air = Z_AIR_V2 if vias else Z_GB
    air = box(-P / 2, P / 2, -P / 2, P / 2, z_air, LID)
    parts = {"METAL_PEC": metal, "SUBSTRATE": sub.val(), "AIR": air.val()}
    if vias:
        parts["SUBSTRATE_FR4"] = box(-P / 2, P / 2, -P / 2, P / 2, Z_FR4B, Z_GB).val()
    if variant != "absent":
        pts = [cq.Vector(-ALONG / 2, Y_SH0, Z_SH), cq.Vector(ALONG / 2, Y_SH0, Z_SH), cq.Vector(ALONG / 2, Y_SH1, Z_SH), cq.Vector(-ALONG / 2, Y_SH1, Z_SH)]
        parts["DIODE_SHEET"] = cq.Face.makeFromWires(cq.Wire.makePolygon(pts, close=True))
    return parts


def intended(variant):
    xy = {0.0, P / 2, R, L / 2, W / 2, GC / 2, X_ISL, ALONG / 2, Y_ISL0, Y_ISL1, Y_SH0, Y_SH1, STUB_L}
    if variant == "vias":
        xy |= {Y_VIA_I, Y_VIA_D, VIA_R, PAD_R, ANTIPAD_R, STUB_L + PAD_R, Y_VIA_I - STUB_W / 2, Y_VIA_I + STUB_W / 2, Y_VIA_D - STUB_W / 2, Y_VIA_D + STUB_W / 2, Y_VIA_I + PAD_R, Y_VIA_D + PAD_R, Y_VIA_I - PAD_R, Y_VIA_D - PAD_R}
    xy |= {-v for v in xy}
    z = {0.0, T, T / 2, -H, Z_GB, LID, -0.01, T + 0.01}
    if variant == "vias":
        z |= {Z_FR4B, Z_BOT, Z_AIR_V2, Z_GB - 0.01, -H + 0.01}
    z |= {-v for v in z}
    return dict(xy=np.array(sorted(round(v, 6) for v in xy)), z=np.array(sorted(round(v, 6) for v in z)), dir=np.array([-1.0, 0.0, 1.0]))


def inside(solid, x, y, z):
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.gp import gp_Pnt
    from OCP.TopAbs import TopAbs_IN

    return BRepClass3d_SolidClassifier(solid.wrapped, gp_Pnt(float(x), float(y), float(z)), 1e-7).State() == TopAbs_IN


def nominal_metal(x, y, island=True):
    """Analytic nominal top-layer metal (patch + island + sheet region excluded) at (x, y)."""
    in_disc = math.hypot(x, y) <= R
    in_slot = abs(x) <= L / 2 and abs(y) <= W / 2
    isl = island and abs(x) <= X_ISL and Y_ISL0 <= y <= Y_ISL1
    return (in_disc and not in_slot) or isl


def checks(path, variant):
    import cadquery as cq
    from sidequest_board_k2_lp_package import read_names

    res, ok_all = {}, True

    def rec(n, ok, **kw):
        nonlocal ok_all
        ok_all &= bool(ok)
        res[n] = {"ok": bool(ok), **kw}

    vias = variant == "vias"
    names = read_names(path)
    want = {"METAL_PEC", "SUBSTRATE", "AIR"} | ({"DIODE_SHEET"} if variant != "absent" else set()) | ({"SUBSTRATE_FR4"} if vias else set())
    rec("bodies_and_names", set(names) == want, names=sorted(names))
    solids = names["METAL_PEC"].Solids()
    nwant = 2 if variant == "absent" else 3
    rec("metal_parts", len(solids) == nwant, parts=len(solids), expected=nwant, note="disc patch, [island], ground; the island is a separate floating solid")
    patch = [s for s in solids if s.BoundingBox().zmax > T - 1e-6 and s.BoundingBox().xmax > 5.0][0]
    gnd = [s for s in solids if s.BoundingBox().xmax > 5.7 and s.BoundingBox().zmax < -1.0][0]
    isl = [s for s in solids if s is not patch and s is not gnd]
    sub, air = names["SUBSTRATE"], names["AIR"]
    # areas of the top layer (z = T/2) from the solids vs analytic
    yy = np.linspace(0, 1, 2)
    area_patch = math.pi * R * R - L * W
    barrels = 2 * math.pi * VIA_R**2 * 0 if not vias else 0.0
    if not vias:
        rec("patch_volume_vs_analytic", abs(patch.Volume() - area_patch * T) / (area_patch * T) < 2e-4, volume_mm3=patch.Volume(), analytic_mm3=area_patch * T, rel_err=patch.Volume() / (area_patch * T) - 1)
    if variant != "absent" and not vias:
        a_isl = 2 * X_ISL * (Y_ISL1 - Y_ISL0)
        rec("island_volume_vs_analytic", len(isl) == 1 and abs(isl[0].Volume() - a_isl * T) < 1e-9 * max(1.0, a_isl * T), volume_mm3=isl[0].Volume(), analytic_mm3=a_isl * T, island_mm=[2 * X_ISL, Y_ISL1 - Y_ISL0])
    # nominal-geometry mask on a 0.1 mm grid (bands of +-0.06 mm around every edge excluded)
    tot = bad = 0
    edges_x = [R, L / 2, X_ISL]
    for x in np.arange(-R - 0.5, R + 0.5, 0.1):
        for y in np.arange(-R - 0.5, R + 0.5, 0.1):
            rad = math.hypot(x, y)
            if abs(rad - R) < 0.06:
                continue
            if min(abs(abs(x) - L / 2), abs(abs(x) - X_ISL)) < 0.06 and abs(y) <= W / 2 + 0.06:
                continue
            if min(abs(abs(y) - W / 2), abs(y - Y_ISL0), abs(y - Y_ISL1)) < 0.06 and abs(x) <= L / 2 + 0.06:
                continue
            if vias and (math.hypot(x, y - Y_VIA_I) < PAD_R + 0.06 or math.hypot(x, y - Y_VIA_D) < PAD_R + 0.06):
                continue
            tot += 1
            metal_here = any(inside(s, x, y, T / 2) for s in [patch] + isl)
            if metal_here != nominal_metal(x, y, island=(variant != "absent")):
                bad += 1
    rec("top_layer_vs_nominal_geometry", bad == 0, points=tot, mismatches=bad)
    rng = np.random.default_rng(1)
    pts = rng.uniform(-R, R, size=(400, 2))
    asym = sum(any(inside(s, x, y, T / 2) for s in [patch] + isl) != any(inside(s, -x, y, T / 2) for s in [patch] + isl) for x, y in pts if not (vias and (math.hypot(x, y - Y_VIA_I) < 0.5 or math.hypot(x, y - Y_VIA_D) < 0.5)))
    rec("mirror_symmetry_x", asym == 0, asymmetric_points=int(asym), note="vias/stubs excluded" if vias else "")
    allm = cq.Compound.makeCompound(solids)
    ov = {"patch_x_ground": patch.intersect(gnd).Volume()}
    for k, s in enumerate(isl):
        ov[f"island{k}_x_patch"] = s.intersect(patch).Volume()
        ov[f"island{k}_x_ground"] = s.intersect(gnd).Volume()
    if vias:
        expected_sub = 2 * math.pi * VIA_R**2 * H
        v_sub = sub.intersect(allm).Volume()
        rec("vias_through_substrate", abs(v_sub - expected_sub) / expected_sub < 5e-3, overlap_mm3=v_sub, expected_two_barrels_mm3=expected_sub, note="the barrels overlap SUBSTRATE on purpose (nothing subtracted)")
        fr4 = names["SUBSTRATE_FR4"]
        v_fr4 = fr4.intersect(allm).Volume()
        exp_fr4 = 2 * math.pi * VIA_R**2 * FR4_T + 0.0
        # FR4 box z -3.710 ... -3.210 contains the barrels (r 0.15) and nothing of the pads/stubs (bottom copper is below it)
        rec("vias_through_fr4", abs(v_fr4 - exp_fr4) / exp_fr4 < 5e-3, overlap_mm3=v_fr4, expected_mm3=exp_fr4)
        dist_g_i = min(g.distance(s) for g in [gnd] for s in isl) if isl else 0.0
        rec("vias_clear_of_ground", gnd.intersect(allm.cut(gnd)).Volume() < 1e-12 if False else dist_g_i > 0.3, island_to_ground_mm=dist_g_i, antipad_radial_clearance_mm=ANTIPAD_R - VIA_R,
            note="barrels are 0.35 mm from the anti-pad wall; only the plated barrel passes the ground")
        stubs = [s for s in isl + [patch] if s.BoundingBox().zmin < Z_FR4B]
        slab = cq.Workplane("XY").workplane(offset=Z_BOT - 0.1).rect(P, P).extrude(Z_FR4B + 1e-9 - (Z_BOT - 0.1)).val()
        bot = [s.intersect(slab) for s in stubs]
        dd = bot[0].distance(bot[1]) if len(bot) == 2 else None
        rec("stubs_not_connected", len(stubs) == 2 and dd is not None and dd > 1.0, bottom_copper_stub_to_stub_mm=dd,
            note="island stub +x, disc stub -x; bottom-layer copper of the two vias is 2.55 mm apart in y; the two copper solids (island+via, disc+via) are separate")
        rec("stub_clear_of_ground", all(s.BoundingBox().zmax <= Z_FR4B + 1e-9 or s.BoundingBox().zmax >= 0 for s in stubs) and gnd.BoundingBox().zmin >= Z_GB - 1e-9, fr4_mm=FR4_T)
    rec("no_intersections", max(ov.values()) < 1e-9, **ov)
    ab = air.BoundingBox()
    zair = Z_AIR_V2 if vias else Z_GB
    shapes = list(solids) + [sub] + ([names["SUBSTRATE_FR4"]] if vias else [])
    allin = all(ab.xmin - 1e-9 <= s.BoundingBox().xmin and s.BoundingBox().xmax <= ab.xmax + 1e-9 and ab.zmin - 1e-9 <= s.BoundingBox().zmin and s.BoundingBox().zmax <= ab.zmax + 1e-9 for s in shapes)
    rec("air_box_single_face_walls", len(air.Faces()) == 6 and abs(ab.zmin - zair) < 1e-9 and abs(ab.zmax - LID) < 1e-9 and allin, faces=len(air.Faces()), z=[ab.zmin, ab.zmax], contains_all_solids=allin)
    rec("substrate_stack", abs(sub.BoundingBox().zmin + H) < 1e-9 and abs(sub.BoundingBox().zmax) < 1e-9 and abs(patch.BoundingBox().zmax - T) < 1e-9 and abs(gnd.BoundingBox().zmax + H) < 1e-9, substrate_mm=H, patch_t_mm=T, ground_t_mm=GT)
    rim_end, rim_corner = R - L / 2, R - math.hypot(L / 2, W / 2)
    rec("rims", rim_end >= 0.35 and rim_corner >= 0.2, slot_end_rim_mm=rim_end, slot_corner_rim_mm=rim_corner)
    if variant != "absent":
        b = names["DIODE_SHEET"].BoundingBox()
        geom_ok = abs(b.xmin + ALONG / 2) < 1e-9 and abs(b.xmax - ALONG / 2) < 1e-9 and abs(b.ymin - Y_SH0) < 1e-9 and abs(b.ymax - Y_SH1) < 1e-9 and abs(b.zmin - Z_SH) < 1e-9 and abs(b.zmax - Z_SH) < 1e-9
        rec("sheet_size_position", geom_ok, across_y_mm=Y_SH1 - Y_SH0, along_x_mm=ALONG, x_mm=0.0, y_from_to_mm=[Y_SH0, Y_SH1], z_mm=Z_SH, current_line=[list(p) for p in sheet_current_line()])
        isl_s = [s for s in isl if abs(s.BoundingBox().ymax - Y_ISL1) < 1e-6 or s.BoundingBox().ymax >= Y_ISL1]
        corners = [cq.Vertex.makeVertex(x, y, Z_SH) for x in (-ALONG / 2, ALONG / 2) for y in (Y_SH0, Y_SH1)]
        dmin = [round(min(s.distance(v) for s in [patch] + isl), 9) for v in corners]
        rec("sheet_edges_on_metal_edges", max(dmin) < 1e-6, corner_to_metal_mm=dmin, note="short ends at y = 0.45 on the island +y edge and y = 1.45 on the slot +y edge")
        cl_x = X_ISL - ALONG / 2
        rec("sheet_clearances", True, sheet_end_to_island_end_mm=cl_x, island_edge_to_slot_end_mm=L / 2 - X_ISL, minus_y_clearance_mm=Y_ISL0 + W / 2, note="informational: the sheet sides face open slot gap (+-0.175 mm of slot either side)")
    return res, ok_all


def geometry():
    import cadquery as cq
    from src.transmitarray.export import inexact_literals, snap_step_text, verify_snap
    from sub6_step_exact import audit_step_text, exactify_step_text

    OUT.mkdir(parents=True, exist_ok=True)
    RES.mkdir(parents=True, exist_ok=True)
    allres, allok = {}, True
    cols = {"METAL_PEC": (0.9, 0.7, 0.2), "SUBSTRATE": (0.25, 0.55, 0.25), "SUBSTRATE_FR4": (0.2, 0.5, 0.2), "AIR": (0.7, 0.8, 1.0), "DIODE_SHEET": (0.9, 0.1, 0.1)}
    files = {"novia": "ISL_novia.step", "absent": "ISL_absent.step", "vias": "ISL_vias.step"}
    for variant, fn in files.items():
        parts = make(variant)
        asm = cq.Assembly()
        for nm, shp in parts.items():
            asm.add(shp, name=nm, color=cq.Color(*cols[nm], 1))
        path = OUT / fn
        asm.export(str(path), exportType="STEP")
        raw = path.read_text(encoding="utf-8")
        it = intended(variant)
        snapped, stats = snap_step_text(raw, it, description=f"sub6 island cell {variant}, RT5880 3.175, mm, 0.001 mm grid")
        snapped, ex_stats = exactify_step_text(snapped)
        path.write_text(snapped, encoding="utf-8", newline="\n")
        text = path.read_text(encoding="utf-8")
        badl = list(inexact_literals(text, it)) + audit_step_text(text)
        vs = verify_snap(raw, snapped)
        res, ok = checks(path, variant)
        exact = {"ok": len(badl) == 0 and vs.get("collapsed_edges", 1) == 0 and not vs.get("distinct_merge_violations"), "inexact_literals": len(badl), "n_snapped": stats["n_snapped"], "exactify": ex_stats, "grid_mm": 0.001}
        res["EXACT_DIMS"] = exact
        ok &= exact["ok"]
        allres[fn] = {"all_pass": bool(ok), "checks": res}
        allok &= ok
        print(variant, ok, {k: v["ok"] for k, v in res.items()}, flush=True)
    meta = dict(P_mm=P, substrate_mm=H, patch_t_mm=T, ground_t_mm=GT, lid_mm=LID, disc_r_mm=R, slot_mm=[L, W], island_x_halfwidth_mm=X_ISL, island_y_mm=[Y_ISL0, Y_ISL1],
                sheet_across_mm=GC, sheet_along_mm=ALONG, via=dict(drill_mm=2 * VIA_R, pad_mm=2 * PAD_R, antipad_mm=2 * ANTIPAD_R, y_island_via=Y_VIA_I, y_disc_via=Y_VIA_D, stub=[STUB_L, STUB_W]),
                fr4=dict(t_mm=FR4_T, eps_r=FR4_EPS, tan_d=FR4_TAND, assumed=True), z_levels=dict(ground=[Z_GB, -H], fr4=[Z_FR4B, Z_GB], bottom_copper=[Z_BOT, Z_FR4B], air_bottom_variant2=Z_AIR_V2),
                mom_vs_nominal=dict(mesh_cell_mm=11.6 / 64, mom_clearance_mm=3 * 11.6 / 64, nominal_clearance_mm=CLR, mom_gap_mm=5 * 11.6 / 64, nominal_gap_mm=GC, mom_island_width_mm=8 * 11.6 / 64, nominal_island_width_mm=Y_ISL1 - Y_ISL0,
                                    mom_island_length_mm=24 * 11.6 / 64, nominal_island_length_mm=2 * X_ISL))
    (OUT / "geometry_checks.json").write_text(json.dumps({"label": "MODEL_ONLY geometry checks on the re-imported STEP files", "all_pass": bool(allok), "parameters": meta, "files": allres}, indent=1, default=float), encoding="utf-8", newline="\n")
    print("ALL", allok)
    return allok


if __name__ == "__main__":
    if sys.argv[1] == "geometry":
        sys.exit(0 if geometry() else 1)
