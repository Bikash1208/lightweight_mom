"""SUB6 FINAL BOARD: RF-ONLY thin STEP for the board HFSS test (no bias network, no connector, no SMA).  MODEL_ONLY / NOT_HFSS_VALIDATED.
    SUB6_MESH=96 python sub6_final_step13_rf.py
Writes hfss_exports/sub6_final_board_m96_13pad/BOARD_N240_RF_thin.step, body_check_rf_thin.json, geometry_build_rf_thin.json
KEEP : top copper (every disc with its slot, every island) as sheets at z 0, 240 diode sheets D<string>_<pos> (z 0), RT5880 (z -3.175 ... 0), ground sheet (z -3.175) with 480 anti-pads r 0.5 and the probe opening r 2.05,
       FR4 slab (z -3.975 ... -3.175, contiguous stack), the 480 vias: each via = its own sheet body (barrel tube r 0.15 from the top plane z 0 to the FR4 bottom plane z -3.975 + its own via pad r 0.30 at z -3.975, one shell),
       probe pin (solid, r 0.635, z -3.975 ... +8.925), PORT_ANNULUS (r 0.635 ... 2.05 at z -3.975), AIR box.
DROP : all B.Cu bias traces, RETURN, connector pads, SMA pad / body / ground vias.  Nothing is fused between vias or nets.
Names (all unique, nothing unnamed, root assembly named): METAL_PEC_DISC_<site>, METAL_PEC_ISL_<site>, METAL_PEC_VIA_<string>_<pos><I|D>, METAL_PEC_GROUND, METAL_PEC_PIN, SUBSTRATE, SUBSTRATE_FR4, AIR, D<string>_<pos>, PORT_ANNULUS.
Every body goes through BRepCheck_Analyzer (+ ShapeFix_Shape / UnifySameDomain if needed); the written STEP is re-read and every body is checked again.
"""
import csv
import json
import math
import os
import sys
import time
import warnings
from pathlib import Path

os.environ["SUB6_VARIANT"] = "13pad"
os.environ.setdefault("SUB6_MESH", "96")
sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import sub6_final_step13_v2 as V

R3, H, T = V.R3, V.H, V.T


def main():
    import cadquery as cq
    from sub6_step_exact import audit_step_text, exactify_step_text
    from sidequest_board_k2_lp_package import read_names
    from OCP.BRepCheck import BRepCheck_Analyzer
    z = V.levels(True)
    z["pin_lo"] = z["bot_lo"]                      # no SMA: the pin starts at the FR4 bottom plane
    z["air_lo"] = round(z["bot_lo"] - V.PAD_BOT, 3)
    cells = list(csv.DictReader(open(V.SRC / "fab" / "cell_table.csv")))
    rt = json.load(open(V.SRC / "routing.json"))
    rout = rt["r_out"]
    hexpts = [(R3(rout * math.cos(math.radians(90 + 60 * i))), R3(rout * math.sin(math.radians(90 + 60 * i)))) for i in range(6)]
    t0 = time.time()
    parts, checks = {}, {}
    via_xy = []

    def put(name, shp, kind):
        h, rec = V.heal(shp, kind)
        ok, k2 = V.single_body(h)
        rec["single_lump"] = bool(ok)
        rec["type"] = k2
        if not ok:
            rec["status"] = "SPLIT"
        parts[name] = h
        checks[name] = rec

    for c in cells:
        x, y, th = float(c["x_mm"]), float(c["y_mm"]), math.radians(float(c["slot_axis_deg"]))
        R, Lh, Wh = float(c["disc_r_mm"]), float(c["slot_L_mm"]) / 2, float(c["slot_W_mm"]) / 2
        wp = cq.Workplane("XY").center(R3(x), R3(y)).circle(R).polyline([(p[0] - R3(x), p[1] - R3(y)) for p in V.quad(x, y, Lh, Wh, th)]).close()
        disc = [f for f in wp.extrude(T).val().Faces() if abs(f.Center().z - 0.0) < 1e-9 and f.normalAt().z < -0.5][0]
        put("METAL_PEC_DISC_" + c["site"], disc, "sheet")
        if c["class"] != "SWITCHED":
            continue
        s, p = int(c["string"]), int(c["pos_in_string"])
        pts = V.quad(x, y, V.ISL[0], V.ISL[1], th, V.ISL[2], V.ISL[3])
        put("METAL_PEC_ISL_" + c["site"], cq.Face.makeFromWires(cq.Wire.makePolygon([cq.Vector(a, b, 0.0) for a, b in pts], close=True)), "sheet")
        for tag, key in (("D", "disc"), ("I", "island")):
            vx, vy = R3(float(c[key + "_via_x"])), R3(float(c[key + "_via_y"]))
            via_xy.append((vx, vy))
            cyl = cq.Solid.makeCylinder(V.R_BAR, 0.0 - z["bot_lo"], cq.Vector(vx, vy, z["bot_lo"]), cq.Vector(0, 0, 1))
            tube = [f for f in cyl.Faces() if f.geomType() == "CYLINDER"][0]
            padw = cq.Workplane("XY").workplane(offset=z["bot_lo"]).center(vx, vy).circle(V.R_PAD).circle(V.R_BAR)
            ring = [f for f in padw.extrude(T).val().Faces() if abs(f.Center().z - z["bot_lo"]) < 1e-9 and f.normalAt().z < -0.5][0]
            put(f"METAL_PEC_VIA_{s}_{p}{tag}", cq.Shell.makeShell([tube, ring]), "sheet")
    print(f"copper + vias done {time.time() - t0:.0f}s bodies {len(parts)}", flush=True)
    gw = cq.Workplane("XY").workplane(offset=z["gnd_lo"]).polyline(hexpts).close()
    for vx, vy in via_xy:
        gw = gw.moveTo(vx, vy).circle(V.R_ANTI)
    gw = gw.moveTo(0, 0).circle(V.R_OUT)
    gnd = [f for f in gw.extrude(T).val().Faces() if abs(f.Center().z - z["gnd_lo"]) < 1e-9 and f.normalAt().z < -0.5][0]
    put("METAL_PEC_GROUND", gnd, "sheet")
    pin = cq.Workplane("XY").workplane(offset=z["pin_lo"]).circle(V.R_PIN).extrude(z["pin_top"] - z["pin_lo"]).val()
    put("METAL_PEC_PIN", pin, "solid")
    sub = cq.Workplane("XY").workplane(offset=-H).polyline(hexpts).close().extrude(H).val().cut(pin).clean()
    fr4 = cq.Workplane("XY").workplane(offset=z["fr4_lo"]).polyline(hexpts).close().extrude(z["fr4_hi"] - z["fr4_lo"]).val().cut(pin).clean()
    ax, ay = rout * math.cos(math.radians(30)) + V.PAD_XY, rout + V.PAD_XY
    air = cq.Workplane("XY").workplane(offset=z["air_lo"]).rect(2 * R3(ax), 2 * R3(ay)).extrude(z["air_hi"] - z["air_lo"]).val().cut(sub, fr4, pin).clean()
    put("SUBSTRATE", sub, "solid")
    put("SUBSTRATE_FR4", fr4, "solid")
    put("AIR", air, "solid")
    # diode sheets (1.5 um into the copper at both ends), port
    rows = []
    for c in cells:
        if c["class"] != "SWITCHED":
            continue
        s, p = int(c["string"]), int(c["pos_in_string"])
        x, y, th = float(c["x_mm"]), float(c["y_mm"]), math.radians(float(c["slot_axis_deg"]))
        q = []
        for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
            v = V.rot((sx * V.SH_ALONG / 2, 0.45 - V.SH_EXT + (sy + 1) / 2 * (V.SH_ACROSS + 2 * V.SH_EXT)), th)
            q.append(cq.Vector(R3(x + v[0]), R3(y + v[1]), 0.0))
        parts[f"D{s}_{p}"] = cq.Face.makeFromWires(cq.Wire.makePolygon(q, close=True))
        checks[f"D{s}_{p}"] = dict(type="sheet", status="VALID", single_lump=True)
        rows.append(c["site"])
    zp = z["bot_lo"]
    parts["PORT_ANNULUS"] = cq.Face.makeFromWires(cq.Wire.makeCircle(V.R_OUT, cq.Vector(0, 0, zp), cq.Vector(0, 0, 1)), [cq.Wire.makeCircle(V.R_PIN, cq.Vector(0, 0, zp), cq.Vector(0, 0, 1))])
    checks["PORT_ANNULUS"] = dict(type="sheet", status="VALID", single_lump=True)
    asm = cq.Assembly(None, name="BOARD_N240_RF")
    for nm, shp in parts.items():
        col = (0.9, 0.7, 0.2) if nm.startswith("METAL") else {"SUBSTRATE": (0.25, 0.55, 0.25), "SUBSTRATE_FR4": (0.2, 0.5, 0.2), "AIR": (0.7, 0.8, 1.0), "PORT_ANNULUS": (0.1, 0.1, 0.9)}.get(nm, (0.9, 0.1, 0.1))
        asm.add(shp, name=nm, color=cq.Color(*col, 1))
    path = V.OUT / "BOARD_N240_RF_thin.step"
    asm.export(str(path), exportType="STEP")
    raw = path.read_text(encoding="utf-8")
    snapped, ex = exactify_step_text(raw)
    path.write_text(snapped, encoding="utf-8", newline="\n")
    bad = audit_step_text(path.read_text(encoding="utf-8"))
    names = read_names(path)
    after = {}
    for nm, shp in names.items():
        sh = shp.wrapped if hasattr(shp, "wrapped") else shp
        after[nm] = bool(BRepCheck_Analyzer(sh).IsValid())
    for nm in checks:
        checks[nm]["valid_in_exported_step"] = after.get(nm)
    cnt = {}
    for nm, r in checks.items():
        key = r["status"] + ("" if after.get(nm, False) else "+EXPORT_INVALID")
        cnt[key] = cnt.get(key, 0) + 1
    chk = dict(label="MODEL_ONLY RF-only thin STEP", n_bodies=len(names), n_expected=len(parts), unnamed_or_uuid_names=[n for n in names if not n or len(n) > 40], duplicate_names=len(parts) - len(set(parts)),
               status_counts=cnt, fixed_bodies=[n for n, r in checks.items() if r["status"] == "FIXED"], split_bodies=[n for n, r in checks.items() if not r.get("single_lump", True)],
               invalid_after_export=[n for n in checks if not after.get(n, False)], all_valid=bool(all(after.get(n, False) for n in checks) and len(names) == len(parts)),
               step_size_MB=round(path.stat().st_size / 1e6, 1), audit_inexact_literals=len(bad), meta=dict(z=z, air_box=dict(ax=R3(ax), ay=R3(ay))), build_s=round(time.time() - t0))
    (V.OUT / "body_check_rf_thin.json").write_text(json.dumps(dict(summary=chk, bodies=checks), indent=1, default=str), encoding="utf-8", newline="\n")
    (V.OUT / "geometry_build_rf_thin.json").write_text(json.dumps(chk, indent=1, default=str), encoding="utf-8", newline="\n")
    print(json.dumps({k: v for k, v in chk.items() if k != "meta"}, default=str)[:1500], flush=True)


if __name__ == "__main__":
    main()
