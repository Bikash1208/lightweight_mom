"""SUB6 island cell: HFSS eigenmode package with bias vias, variant 2 (vias, no choke) and variant 3 (vias + single 120-degree radial-stub chokes).  MODEL_ONLY.
Takes over the generator of commit 32f7c8a (compute/sub6-island-chokes) with the user's decisions of 2026-11-09:
  * FR4 0.8 mm, choke line 0.15 mm (Z0 121 ohm), lambda_g/4 = 7.571 mm meandered inside the cell, ONE 120-deg radial stub r 3.148 mm at each line end, bias trace 0.3 mm (via pad to choke line);
  * vias at the both-direction |J| minima (island frame, slot along x): island (-1.78, -0.53), disc (-3.63, -3.63); pads r 0.30, drill 0.30 (barrel r 0.15), anti-pads r 0.50 in the ground;
  * variant 2 and variant 3 have IDENTICAL top metal, vias, pads and anti-pads and differ only below the ground: v2 = 2.0 x 0.3 mm bias stubs (toward +x of the island frame, ends inside the cell), v3 = the chokes;
  * frame ROTATED by +90 deg about z: slot along y, ACROSS the slot = x = $ph_x, ALONG = $ph_y;
  * stack (z, mm): patch 0..0.035, RT5880 -3.175..0, ground -3.210..-3.175, FR4 0.8 -4.010..-3.210, bottom copper 35 um -4.045..-4.010, AIR from -34.045 (30 mm below the bottom copper) to +150 (PEC lid);
    RADIATION boundary on the bottom face of AIR, lattice pairs on the four AIR side faces.
    python sub6_island_choke_package.py geometry
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

import sub6_island_choke_layout2 as lay
import sub6_island_hfss_package as base

OUT = ROOT / "hfss_exports" / "sub6_island_choke"
RES = ROOT / "results" / "sub6_island_choke"
P, H, T, GT, LID = base.P, base.H, base.T, base.GT, base.LID
FR4_T, FR4_EPS, FR4_TAND = 0.8, 4.4, 0.02
AIR_BELOW = 30.0
VIA_I, VIA_D = lay.VIA_I, lay.VIA_D
VIA_R, PAD_R, ANTIPAD_R = base.VIA_R, base.PAD_R, base.ANTIPAD_R
STUB_L, STUB_W = base.STUB_L, base.STUB_W
GRID = 0.001
Z_GB = -H - GT
Z_FR4B = round(Z_GB - FR4_T, 6)
Z_BOT = round(Z_FR4B - GT, 6)
Z_AIR = round(Z_BOT - AIR_BELOW, 6)


def g(v):
    return round(round(v / GRID) * GRID, 6)


def load_chokes():
    fit = json.loads((RES / "layout_fit_single.json").read_text())
    out = []
    for c in fit["chokes"]:
        out.append(lay.build_single(tuple(c["via"]), c["dir0"], c["dir"], c["flip"], c["spacing"], c["n_runs"], c["link"], c["lead"], fit["design"]["l_mm"], fit["design"]["r2_mm"]))
    return fit, out


def make(variant, chokes=None):
    """variant 'vias' (2) | 'choke' (3) -> {name: cq shape} in the rotated frame (ACROSS = x)."""
    import cadquery as cq

    def box(x0, x1, y0, y1, za, zb):
        return cq.Workplane("XY").workplane(offset=za).center((x0 + x1) / 2, (y0 + y1) / 2).rect(x1 - x0, y1 - y0).extrude(zb - za)

    def cyl(x, y, r, za, zb):
        return cq.Workplane("XY").workplane(offset=za).center(x, y).circle(r).extrude(zb - za)

    patch = cq.Workplane("XY").circle(base.R).extrude(T).cut(box(-base.L / 2, base.L / 2, -base.W / 2, base.W / 2, -0.01, T + 0.01))
    island = box(-base.X_ISL, base.X_ISL, base.Y_ISL0, base.Y_ISL1, 0.0, T)
    ground = box(-P / 2, P / 2, -P / 2, P / 2, Z_GB, -H)

    def barrel_pad(v):
        return cyl(v[0], v[1], VIA_R, Z_BOT, 0.0).union(cyl(v[0], v[1], PAD_R, Z_BOT, Z_FR4B))

    isl = island.union(barrel_pad(VIA_I))
    dsc = patch.union(barrel_pad(VIA_D))
    if variant == "vias":
        for v, which in ((VIA_I, "i"), (VIA_D, "d")):
            st = box(v[0], v[0] + STUB_L, v[1] - STUB_W / 2, v[1] + STUB_W / 2, Z_BOT, Z_FR4B)
            if which == "i":
                isl = isl.union(st)
            else:
                dsc = dsc.union(st)
    else:
        for which, c in (("i", chokes[0]), ("d", chokes[1])):
            solid = None
            for k, poly in enumerate(c["parts"]):
                if k == 0:
                    continue  # the pad
                if k == len(c["parts"]) - 1:  # the 120-degree radial stub: apex + exact-centre arc
                    apex = (g(c["apex"][0]), g(c["apex"][1]))
                    pts = [(g(x), g(y)) for x, y in poly[1:]]
                    mid_ang = c["axis_angle"]
                    mid = (g(apex[0] + c["r2"] * math.cos(mid_ang)), g(apex[1] + c["r2"] * math.sin(mid_ang)))
                    part = cq.Workplane("XY").workplane(offset=Z_BOT).moveTo(*apex).lineTo(*pts[0]).threePointArc(mid, pts[-1]).close().extrude(Z_FR4B - Z_BOT)
                else:
                    xs = [g(x) for x, _ in poly]
                    ys = [g(y) for _, y in poly]
                    part = box(min(xs), max(xs), min(ys), max(ys), Z_BOT, Z_FR4B)
                solid = part if solid is None else solid.union(part)
            if which == "i":
                isl = isl.union(solid)
            else:
                dsc = dsc.union(solid)
    ground = ground.cut(cyl(VIA_I[0], VIA_I[1], ANTIPAD_R, Z_GB - 0.01, -H + 0.01)).cut(cyl(VIA_D[0], VIA_D[1], ANTIPAD_R, Z_GB - 0.01, -H + 0.01))
    sub = box(-P / 2, P / 2, -P / 2, P / 2, -H, 0.0)
    fr4 = box(-P / 2, P / 2, -P / 2, P / 2, Z_FR4B, Z_GB)
    air = box(-P / 2, P / 2, -P / 2, P / 2, Z_AIR, LID)
    pts = [cq.Vector(-base.ALONG / 2, base.Y_SH0, base.Z_SH), cq.Vector(base.ALONG / 2, base.Y_SH0, base.Z_SH), cq.Vector(base.ALONG / 2, base.Y_SH1, base.Z_SH), cq.Vector(-base.ALONG / 2, base.Y_SH1, base.Z_SH)]
    parts = {"METAL_PEC": cq.Compound.makeCompound([dsc.val(), isl.val(), ground.val()]), "SUBSTRATE": sub.val(), "SUBSTRATE_FR4": fr4.val(), "AIR": air.val(),
             "DIODE_SHEET": cq.Face.makeFromWires(cq.Wire.makePolygon(pts, close=True))}
    return {k: v.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), 90.0) for k, v in parts.items()}


def intended(chokes):
    xy = {0.0, P / 2, base.R, base.L / 2, base.W / 2, base.GC / 2, base.X_ISL, base.ALONG / 2, base.Y_ISL0, base.Y_ISL1, base.Y_SH0, base.Y_SH1, STUB_L, VIA_R, PAD_R, ANTIPAD_R, STUB_W / 2}
    for v in (VIA_I, VIA_D):
        xy |= {v[0], v[1], v[0] + STUB_L, v[0] + PAD_R, v[0] - PAD_R, v[1] + PAD_R, v[1] - PAD_R, v[1] + STUB_W / 2, v[1] - STUB_W / 2}
    for c in chokes or []:
        for poly in c["parts"][1:]:
            xy |= {g(v) for pt in poly for v in pt}
        xy |= {g(c["apex"][0]), g(c["apex"][1])}
    xy |= {-v for v in xy}
    z = {0.0, T, T / 2, -H, LID, -0.01, T + 0.01, Z_GB, Z_FR4B, Z_BOT, Z_AIR, Z_GB - 0.01, -H + 0.01}
    z |= {-v for v in z}
    return dict(xy=np.array(sorted(round(v, 6) for v in xy)), z=np.array(sorted(round(v, 6) for v in z)), dir=np.array([-1.0, 0.0, 1.0]))


def inside(solid, x, y, z):
    return base.inside(solid, x, y, z)


def checks(path, variant, chokes, fit):
    import cadquery as cq
    from sidequest_board_k2_lp_package import read_names

    res, ok_all = {}, True

    def rec(n, ok, **kw):
        nonlocal ok_all
        ok_all &= bool(ok)
        res[n] = {"ok": bool(ok), **kw}

    names = read_names(path)
    rec("bodies_and_names", set(names) == {"METAL_PEC", "SUBSTRATE", "SUBSTRATE_FR4", "AIR", "DIODE_SHEET"}, names=sorted(names))
    solids = names["METAL_PEC"].Solids()
    rec("metal_parts", len(solids) == 3, parts=len(solids), note="disc patch + via(+choke), floating island + via(+choke), ground with two anti-pads")
    bb = [s.BoundingBox() for s in solids]
    gnd = [s for s, b in zip(solids, bb) if b.xmax > 5.7 and b.zmax < -1.0][0]
    rest = [s for s in solids if s is not gnd]
    patch = [s for s in rest if s.BoundingBox().xmax > 4.5 and s.BoundingBox().ymax > 4.5 or s.BoundingBox().xmin < -4.5][0] if False else max(rest, key=lambda s: s.BoundingBox().xlen * s.BoundingBox().ylen if s.BoundingBox().zmax > 0 else 0)
    isl = [s for s in rest if s is not patch][0]
    sub, fr4, air = names["SUBSTRATE"], names["SUBSTRATE_FR4"], names["AIR"]
    # top layer (z = T/2) against the nominal geometry in the rotated frame: rotate the query point back (x, y) -> (y, -x)
    tot = bad = 0
    for x in np.arange(-base.R - 0.5, base.R + 0.5, 0.1):
        for y in np.arange(-base.R - 0.5, base.R + 0.5, 0.1):
            xo, yo = y, -x
            rad = math.hypot(xo, yo)
            if abs(rad - base.R) < 0.06 or (min(abs(abs(xo) - base.L / 2), abs(abs(xo) - base.X_ISL)) < 0.06 and abs(yo) <= base.W / 2 + 0.06) or (
                    min(abs(abs(yo) - base.W / 2), abs(yo - base.Y_ISL0), abs(yo - base.Y_ISL1)) < 0.06 and abs(xo) <= base.L / 2 + 0.06):
                continue
            if any(math.hypot(xo - v[0], yo - v[1]) < PAD_R + 0.06 for v in (VIA_I, VIA_D)):
                continue
            tot += 1
            if (inside(patch, x, y, T / 2) or inside(isl, x, y, T / 2)) != base.nominal_metal(xo, yo, True):
                bad += 1
    rec("top_layer_vs_nominal_geometry", bad == 0, points=tot, mismatches=bad)
    ov_sub = sub.intersect(cq.Compound.makeCompound(solids)).Volume()
    exp_sub = 2 * math.pi * VIA_R**2 * H
    rec("vias_through_substrate", abs(ov_sub - exp_sub) / exp_sub < 5e-3, overlap_mm3=ov_sub, expected_two_barrels_mm3=exp_sub)
    ov_fr4 = fr4.intersect(cq.Compound.makeCompound(solids)).Volume()
    exp_fr4 = 2 * math.pi * VIA_R**2 * FR4_T
    rec("vias_through_fr4", abs(ov_fr4 - exp_fr4) / exp_fr4 < 5e-3, overlap_mm3=ov_fr4, expected_mm3=exp_fr4)
    d_gi = min(gnd.distance(patch), gnd.distance(isl))
    rec("vias_clear_of_ground", d_gi > ANTIPAD_R - VIA_R - 0.01, copper_to_ground_mm=d_gi, expected_radial_clearance_mm=ANTIPAD_R - VIA_R)
    ov = {"patch_x_ground": patch.intersect(gnd).Volume(), "island_x_ground": isl.intersect(gnd).Volume(), "island_x_patch": isl.intersect(patch).Volume()}
    rec("no_intersections", max(ov.values()) < 1e-9, **ov)
    slab = cq.Workplane("XY").workplane(offset=Z_BOT - 0.1).rect(2 * P, 2 * P).extrude(Z_FR4B + 1e-9 - (Z_BOT - 0.1)).val()
    bot = [s.intersect(slab) for s in (patch, isl)]
    d_between = bot[0].distance(bot[1])
    d_img = min(bot[a].distance(bot[b].translate(cq.Vector(ix * P, iy * P, 0))) for a in (0, 1) for b in (0, 1) for ix in (-1, 0, 1) for iy in (-1, 0, 1) if (ix, iy) != (0, 0))
    bbs = [b.BoundingBox() for b in bot]
    inside_cell = all(-P / 2 <= bx.xmin and bx.xmax <= P / 2 and -P / 2 <= bx.ymin and bx.ymax <= P / 2 for bx in bbs)
    if variant == "choke":
        rec("bottom_layer_gaps", d_between >= 0.15 and d_img >= 0.15 and inside_cell, between_chokes_mm=d_between, to_periodic_images_mm=d_img, inside_cell=inside_cell, layout_gaps_polygon_model=fit["gaps"],
            note="nothing crosses the cell boundary; trace / gap >= 0.15 mm")
    else:
        rec("bottom_layer_gaps", d_between >= 0.15 and d_img >= 0.15 and inside_cell, between_stubs_mm=d_between, to_periodic_images_mm=d_img, inside_cell=inside_cell)
    ab = air.BoundingBox()
    shapes = list(solids) + [sub, fr4]
    allin = all(ab.xmin - 1e-9 <= s.BoundingBox().xmin and s.BoundingBox().xmax <= ab.xmax + 1e-9 and ab.zmin - 1e-9 <= s.BoundingBox().zmin and s.BoundingBox().zmax <= ab.zmax + 1e-9 for s in shapes)
    rec("air_box_single_face_walls", len(air.Faces()) == 6 and abs(ab.zmin - Z_AIR) < 1e-9 and abs(ab.zmax - LID) < 1e-9 and allin, faces=len(air.Faces()), z=[ab.zmin, ab.zmax], below_bottom_copper_mm=Z_BOT - ab.zmin, contains_all_solids=allin)
    rec("stack", abs(sub.BoundingBox().zmin + H) < 1e-9 and abs(patch.BoundingBox().zmax - T) < 1e-9 and abs(gnd.BoundingBox().zmax + H) < 1e-9 and abs(fr4.BoundingBox().zmin - Z_FR4B) < 1e-9
        and abs(min(patch.BoundingBox().zmin, isl.BoundingBox().zmin) - Z_BOT) < 1e-9, substrate_mm=H, fr4_mm=FR4_T, bottom_copper_z=[Z_BOT, Z_FR4B])
    b = names["DIODE_SHEET"].BoundingBox()
    # rotated sheet: x -> -y_old, y -> x_old : x in [-1.45, -0.45], y in [-0.175, 0.175]
    rec("sheet_size_position", abs(b.xmin + base.Y_SH1) < 1e-9 and abs(b.xmax + base.Y_SH0) < 1e-9 and abs(b.ymin + base.ALONG / 2) < 1e-9 and abs(b.ymax - base.ALONG / 2) < 1e-9 and abs(b.zmin - base.Z_SH) < 1e-9,
        across_x_mm=[b.xmin, b.xmax], along_y_mm=[b.ymin, b.ymax], z_mm=base.Z_SH, current_line=[[-base.Y_SH0, 0.0, base.Z_SH], [-base.Y_SH1, 0.0, base.Z_SH]])
    corners = [cq.Vertex.makeVertex(-y, x, base.Z_SH) for x in (-base.ALONG / 2, base.ALONG / 2) for y in (base.Y_SH0, base.Y_SH1)]
    dmin = [round(min(s.distance(v) for s in (patch, isl)), 9) for v in corners]
    rec("sheet_edges_on_metal_edges", max(dmin) < 1e-6, corner_to_metal_mm=dmin)
    return res, ok_all


def geometry():
    import cadquery as cq
    from src.transmitarray.export import inexact_literals, snap_step_text, verify_snap
    from sub6_step_exact import audit_step_text, exactify_step_text

    OUT.mkdir(parents=True, exist_ok=True)
    fit, chokes = load_chokes()
    cols = {"METAL_PEC": (0.9, 0.7, 0.2), "SUBSTRATE": (0.25, 0.55, 0.25), "SUBSTRATE_FR4": (0.2, 0.5, 0.2), "AIR": (0.7, 0.8, 1.0), "DIODE_SHEET": (0.9, 0.1, 0.1)}
    allres, allok = {}, True
    for variant, fn in (("vias", "ISL_vias2.step"), ("choke", "ISL_choke3.step")):
        parts = make(variant, chokes if variant == "choke" else None)
        asm = cq.Assembly()
        for nm, shp in parts.items():
            asm.add(shp, name=nm, color=cq.Color(*cols[nm], 1))
        path = OUT / fn
        asm.export(str(path), exportType="STEP")
        raw = path.read_text(encoding="utf-8")
        it = intended(chokes if variant == "choke" else None)
        snapped, stats = snap_step_text(raw, it, description=f"sub6 island cell {variant}, FR4 {FR4_T} mm, ACROSS = x, mm, 0.001 mm grid")
        snapped, ex = exactify_step_text(snapped)
        path.write_text(snapped, encoding="utf-8", newline="\n")
        text = path.read_text(encoding="utf-8")
        bad = list(inexact_literals(text, it)) + audit_step_text(text)
        vs = verify_snap(raw, snapped)
        res, ok = checks(path, variant, chokes, fit)
        res["EXACT_DIMS"] = {"ok": len(audit_step_text(text)) == 0, "inexact_literals_audit": len(audit_step_text(text)), "inexact_vs_design_table": len(list(inexact_literals(text, it))), "n_snapped": stats["n_snapped"], "exactify": ex, "grid_mm": GRID,
                             "verify_snap": {k: (v if isinstance(v, (int, float, bool, str)) else str(v)) for k, v in vs.items()}}
        ok &= res["EXACT_DIMS"]["ok"] and vs.get("collapsed_edges", 1) == 0
        allres[fn] = {"all_pass": bool(ok), "checks": res}
        allok &= ok
        print(variant, ok, {k: v["ok"] for k, v in res.items()}, flush=True)
    meta = dict(P_mm=P, substrate_mm=H, patch_t_mm=T, ground_t_mm=GT, fr4=dict(t_mm=FR4_T, eps_r=FR4_EPS, tan_d=FR4_TAND, assumed=True), z_levels=dict(ground=[Z_GB, -H], fr4=[Z_FR4B, Z_GB], bottom_copper=[Z_BOT, Z_FR4B], air=[Z_AIR, LID]),
                vias_island_frame=dict(island=VIA_I, disc=VIA_D, drill_mm=2 * VIA_R, pad_mm=2 * PAD_R, antipad_mm=2 * ANTIPAD_R), frame="rotated +90 deg about z: slot along y, ACROSS = x = $ph_x",
                chokes=fit, v2_stubs=dict(length_mm=STUB_L, width_mm=STUB_W, direction="+x of the island frame"))
    (OUT / "geometry_checks.json").write_text(json.dumps({"label": "MODEL_ONLY geometry checks on the re-imported STEP files", "all_pass": bool(allok), "parameters": meta, "files": allres}, indent=1, default=float), encoding="utf-8", newline="\n")
    print("ALL", allok)
    return allok


if __name__ == "__main__":
    if sys.argv[1] == "geometry":
        sys.exit(0 if geometry() else 1)
