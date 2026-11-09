"""SUB6 FINAL BOARD, RF-only board STEP exactly to the 10.5 GHz export standard (hfss_exports/EXPORT_STANDARD.md, sidequest_board_CP_LP_hybrid/HFSS_SETUP.md).  MODEL_ONLY / NOT_HFSS_VALIDATED.
    python sub6_final_rf_std.py full|test9
Writes hfss_exports/sub6_final_board_m96_13pad/BOARD_N240_RF_full.step / BOARD_N240_RF_9cell.step, rf_std_<variant>.json (re-import gate), diode_sheets_rf_<variant>.csv
Bodies (names that must survive the re-import): METAL_PEC (ONE part: one compound of disjoint solids, all copper 35 um SOLIDS), SUBSTRATE (RT5880), SUBSTRATE_FR4, AIR (one box overlapping everything, nothing subtracted, 6 faces),
  PORT (annulus r 0.635 ... 2.05 at z -3.175, in the ground opening), D<string>_<pos> (diode sheets).  No other sheets.
Stack (thick): top Cu z 0 ... 0.035, RT5880 -3.175 ... 0, ground Cu -3.210 ... -3.175 (opening r 2.05), FR4 -4.010 ... -3.175 minus ground plate / barrels (the 0.035 band of the ground is FR4 outside the plate),
  via pads (r 0.30) -4.045 ... -4.010, barrels r 0.15 from the pad top (z -4.010) to the disc / island bottom (z 0), pin r 0.635 z -3.175 ... +8.925.  METAL_PEC lumps: each switched disc / island fused with its own via (barrel + pad),
  518 static discs, the ground plate, the pin.  SUBSTRATE has holes for the barrels and the pin, SUBSTRATE_FR4 holes for barrels and the ground plate (no pin hole: the FR4 is solid under the port).
Diode sheets: 1.00 x 0.35 mm at z 0.0175 (copper mid-thickness), in the gap between the island edge and the slot edge.  To make the sheet edges lie EXACTLY on the copper side faces on the 0.001 mm grid for rotated cells,
  the four sheet corner points (snapped to 0.001) are also vertices of the island polygon (anode edge) and of the slot polygon of the disc (cathode edge): the island top face and the slot side face carry the sheet edge as one of
  their edges, no overlap, no gap (the copper outline differs from the ideal one by <= 0.7 um).
test9: the 9 switched cells nearest to the feed + all static cells within 17 mm of them, same stack, same feed (pin, port, opening), hexagon sized to the included cells.
"""
import csv
import json
import math
import os
import sys
import time
import warnings
from pathlib import Path

os.environ.setdefault("SUB6_MESH", "96")
sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
SRC = ROOT / "results" / "sub6_final_board" / "mesh96"
OUT = ROOT / "hfss_exports" / "sub6_final_board_m96_13pad"
R3 = lambda v: round(float(v) + 0.0, 3)
H, T = 3.175, 0.035
Z_GND_LO, Z_FR4_LO, Z_PAD_LO = -H - T, -H - T - 0.8, -H - T - 0.8 - T
Z_FR4_BOT = Z_FR4_LO            # -4.010 (top of the via pads)
R_PIN, R_OUT, R_BAR, R_PAD, R_ANTI = 0.635, 2.05, 0.15, 0.30, 0.50
PIN_TOP = round(-H + 12.1, 3)
ISL = (2.174, 0.65, 0.0, -0.20)
SH_W, GAP0, GAP1, Z_SHEET = 0.35, 0.45, 1.45, T / 2
AIR_XY, AIR_BOT, AIR_TOP = 15.0, 15.0, 40.0


def rot(v, a):
    c, s = math.cos(a), math.sin(a)
    return (c * v[0] - s * v[1], s * v[0] + c * v[1])


def lp(x, y, th, u, v):
    """cell-local (u, v) -> board point rounded to the 0.001 mm grid"""
    r = rot((u, v), th)
    return (R3(x + r[0]), R3(y + r[1]))


def main():
    variant = sys.argv[1] if len(sys.argv) > 1 else "full"
    import cadquery as cq
    from sub6_step_exact import audit_step_text, exactify_step_text
    from sidequest_board_k2_lp_package import read_names
    cells = list(csv.DictReader(open(SRC / "fab" / "cell_table.csv")))
    rt = json.load(open(SRC / "routing.json"))
    rout = rt["r_out"]
    t0 = time.time()
    if variant == "test9":
        sw = sorted([c for c in cells if c["class"] == "SWITCHED"], key=lambda c: math.hypot(float(c["x_mm"]), float(c["y_mm"])))[:9]
        keep = {c["site"] for c in sw}
        for c in cells:
            if c["class"] != "SWITCHED" and any(math.hypot(float(c["x_mm"]) - float(s["x_mm"]), float(c["y_mm"]) - float(s["y_mm"])) <= 17.0 for s in sw):
                keep.add(c["site"])
        cells = [c for c in cells if c["site"] in keep]
        dmax = max(math.hypot(float(c["x_mm"]), float(c["y_mm"])) + float(c["disc_r_mm"]) for c in cells)
        rout = (dmax + 6.0) / math.cos(math.radians(30))
    hexpts = [(R3(rout * math.cos(math.radians(90 + 60 * i))), R3(rout * math.sin(math.radians(90 + 60 * i)))) for i in range(6)]
    lumps, barrels, sheets, rows = [], [], {}, []
    via_xy = []
    for c in cells:
        x, y, th = float(c["x_mm"]), float(c["y_mm"]), math.radians(float(c["slot_axis_deg"]))
        R, Lh, Wh = float(c["disc_r_mm"]), float(c["slot_L_mm"]) / 2, float(c["slot_W_mm"]) / 2
        sw = c["class"] == "SWITCHED"
        # slot polygon (counter-clockwise from the lower left), with the two cathode-edge sheet corners on its top edge
        slot = [lp(x, y, th, -Lh, -Wh), lp(x, y, th, Lh, -Wh), lp(x, y, th, Lh, Wh)]
        if sw:
            k2, k1 = lp(x, y, th, SH_W / 2, GAP1), lp(x, y, th, -SH_W / 2, GAP1)
            slot += [k2, k1]
        slot += [lp(x, y, th, -Lh, Wh)]
        wp = cq.Workplane("XY").center(R3(x), R3(y)).circle(R).polyline([(p[0] - R3(x), p[1] - R3(y)) for p in slot]).close()
        disc = wp.extrude(T).val()
        if not sw:
            lumps.append(disc)
            continue
        s, p = int(c["string"]), int(c["pos_in_string"])
        a1, a2 = lp(x, y, th, -SH_W / 2, GAP0), lp(x, y, th, SH_W / 2, GAP0)
        ipts = [lp(x, y, th, ISL[2] - ISL[0], ISL[3] - ISL[1]), lp(x, y, th, ISL[2] + ISL[0], ISL[3] - ISL[1]), lp(x, y, th, ISL[2] + ISL[0], ISL[3] + ISL[1]), a2, a1, lp(x, y, th, ISL[2] - ISL[0], ISL[3] + ISL[1])]
        isl = cq.Workplane("XY").polyline(ipts).close().extrude(T).val()
        for key, body in (("disc", disc), ("island", isl)):
            vx, vy = R3(float(c[key + "_via_x"])), R3(float(c[key + "_via_y"]))
            via_xy.append((vx, vy))
            barrel = cq.Workplane("XY").workplane(offset=Z_FR4_BOT).center(vx, vy).circle(R_BAR).extrude(0.0 - Z_FR4_BOT).val()
            pad = cq.Workplane("XY").workplane(offset=Z_PAD_LO).center(vx, vy).circle(R_PAD).extrude(Z_FR4_BOT - Z_PAD_LO).val()
            barrels.append(barrel)
            lumps.append(body.fuse(barrel, pad).clean())
        sheet = cq.Face.makeFromWires(cq.Wire.makePolygon([cq.Vector(a1[0], a1[1], Z_SHEET), cq.Vector(a2[0], a2[1], Z_SHEET), cq.Vector(k2[0], k2[1], Z_SHEET), cq.Vector(k1[0], k1[1], Z_SHEET)], close=True))
        sheets[f"D{s}_{p}"] = sheet
        d = rot((0.0, 1.0), th)
        rows.append([f"D{s}_{p}", c["site"], a1[0], a1[1], a2[0], a2[1], k1[0], k1[1], k2[0], k2[1], round(d[0], 6), round(d[1], 6), s, p])
    print(f"cells done {time.time() - t0:.0f}s lumps {len(lumps)} sheets {len(sheets)}", flush=True)
    gw = cq.Workplane("XY").workplane(offset=Z_GND_LO).polyline(hexpts).close()
    for vx, vy in via_xy:
        gw = gw.moveTo(vx, vy).circle(R_ANTI)
    gw = gw.moveTo(0, 0).circle(R_OUT)
    ground = gw.extrude(T).val()
    lumps.append(ground)
    pin = cq.Workplane("XY").workplane(offset=-H).circle(R_PIN).extrude(PIN_TOP + H).val()
    lumps.append(pin)
    metal = cq.Compound.makeCompound(lumps)
    sub = cq.Workplane("XY").workplane(offset=-H).polyline(hexpts).close().extrude(H).val().cut(*([b for b in barrels] + [pin])).clean()
    fr4 = cq.Workplane("XY").workplane(offset=Z_FR4_LO).polyline(hexpts).close().extrude((-H) - Z_FR4_LO).val()
    fr4 = fr4.cut(*(barrels + [ground])).clean()
    ax, ay = rout * math.cos(math.radians(30)) + AIR_XY, rout + AIR_XY
    z0, z1 = round(Z_PAD_LO - AIR_BOT, 3), round(T + AIR_TOP, 3)
    air = cq.Workplane("XY").workplane(offset=z0).rect(2 * R3(ax), 2 * R3(ay)).extrude(z1 - z0).val()
    port = cq.Face.makeFromWires(cq.Wire.makeCircle(R_OUT, cq.Vector(0, 0, -H), cq.Vector(0, 0, 1)), [cq.Wire.makeCircle(R_PIN, cq.Vector(0, 0, -H), cq.Vector(0, 0, 1))])
    print(f"bodies built {time.time() - t0:.0f}s", flush=True)
    asm = cq.Assembly(None, name="BOARD_N240_RF")
    asm.add(metal, name="METAL_PEC", color=cq.Color(0.9, 0.7, 0.2, 1))
    asm.add(sub, name="SUBSTRATE", color=cq.Color(0.25, 0.55, 0.25, 1))
    asm.add(fr4, name="SUBSTRATE_FR4", color=cq.Color(0.2, 0.5, 0.2, 1))
    asm.add(air, name="AIR", color=cq.Color(0.7, 0.8, 1.0, 1))
    asm.add(port, name="PORT", color=cq.Color(0.1, 0.1, 0.9, 1))
    for nm, f in sheets.items():
        asm.add(f, name=nm, color=cq.Color(0.9, 0.1, 0.1, 1))
    path = OUT / ("BOARD_N240_RF_full.step" if variant == "full" else "BOARD_N240_RF_9cell.step")
    asm.export(str(path), exportType="STEP")
    raw = path.read_text(encoding="utf-8")
    snapped, ex = exactify_step_text(raw)
    path.write_text(snapped, encoding="utf-8", newline="\n")
    bad = audit_step_text(path.read_text(encoding="utf-8"))
    with open(OUT / f"diode_sheets_rf_{variant}.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["name", "site", "anode_edge_x1", "anode_edge_y1", "anode_edge_x2", "anode_edge_y2", "cathode_edge_x1", "cathode_edge_y1", "cathode_edge_x2", "cathode_edge_y2", "current_dir_x", "current_dir_y", "string", "pos_in_string"])
        w.writerows(rows)
    print(f"exported {path.name} {path.stat().st_size / 1e6:.1f} MB, inexact literals {len(bad)}, {time.time() - t0:.0f}s", flush=True)
    json.dump(dict(variant=variant, path=str(path), n_cells=len(cells), n_switched=len(sheets), hex_circumradius=rout, audit_inexact_literals=len(bad), n_lumps=len(lumps)), open(OUT / f"rf_std_build_{variant}.json", "w"), indent=1)


if __name__ == "__main__":
    main()
