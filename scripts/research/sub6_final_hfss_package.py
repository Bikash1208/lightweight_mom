"""SUB6 FINAL BOARD step 5: full-board HFSS package (N 240, calibrated design).  PRELIMINARY_EXTRAPOLATED_CELL / MODEL_ONLY / NOT_HFSS_VALIDATED.
Reads results/sub6_final_board/fab/{traces.json,cell_table.csv,layout_summary.json}; writes hfss_exports/sub6_final_board/:
  BOARD_N240.step (METAL_PEC, SUBSTRATE, SUBSTRATE_FR4, AIR, DIODE_SHEET, PORT_ANNULUS), diode_sheets.csv, geometry_checks.json, predictions.json, results_template.csv, HFSS_SETUP.md
    python sub6_final_hfss_package.py
Frame: board frame of the layout (feed at the origin, +z up).  z levels: patch 0 ... 0.035; substrate RT5880 -3.175 ... 0; ground -3.210 ... -3.175; FR4 0.8 -4.010 ... -3.210 (eps_r 4.4, tan d 0.02 ASSUMED);
bottom copper -4.045 ... -4.010.  Probe: PEC pin r 0.635 from the bottom copper plane (z -4.045) to 12.1 mm above the ground top (z +8.925), through a ground opening r 2.05 (port annulus r 0.635 ... 2.05 at z -4.010).
"""
import csv
import json
import os
import math
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np

import sub6_final_layout as L

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
FAB = L.OUT / "fab"
OUT = ROOT / "hfss_exports" / (("sub6_final_board_m96" if os.environ.get("SUB6_MESH") == "96" else "sub6_final_board") + ("_13pad" if os.environ.get("SUB6_VARIANT") == "13pad" else ""))
H, T, GT = 3.175, 0.035, 0.035
Z_GB, Z_FR4B, Z_BOT = -H - GT, -H - GT - 0.8, -H - GT - 0.8 - 0.035
R_PIN, R_OUT, H_PIN = 0.635, 2.05, 12.1
Z_PIN_TOP = round(-H + H_PIN, 3)
PAD_XY, PAD_TOP, PAD_BOT = 15.0, 40.0, 15.0
Z_AIR0, Z_LID = round(Z_BOT - PAD_BOT, 3), round(T + PAD_TOP, 3)
SH_ALONG, SH_ACROSS, Z_SH = 0.35, 1.00, T / 2
ISL = 2.174
OFF = dict(R_ohm=3876.0, C_pF=0.1703)
ON = dict(R_ohm=162.4, L_nH=0.6113)
R3 = lambda v: round(float(v) + 0.0, 3)


def load():
    cells = list(csv.DictReader(open(FAB / "cell_table.csv")))
    traces = json.load(open(FAB / "traces.json"))
    meta = json.load(open(FAB / "layout_summary.json"))
    return cells, traces, meta


def seg_box(cq, a, b, w, z0, z1):
    dx, dy = b[0] - a[0], b[1] - a[1]
    ln = math.hypot(dx, dy)
    if ln < 1e-9:
        return None
    ux, uy = dx / ln, dy / ln
    ext = w / 2
    pts = [(a[0] - ux * ext - uy * w / 2, a[1] - uy * ext + ux * w / 2), (b[0] + ux * ext - uy * w / 2, b[1] + uy * ext + ux * w / 2),
           (b[0] + ux * ext + uy * w / 2, b[1] + uy * ext - ux * w / 2), (a[0] - ux * ext + uy * w / 2, a[1] - uy * ext - ux * w / 2)]
    pts = [(R3(x), R3(y)) for x, y in pts]
    return cq.Workplane("XY").workplane(offset=z0).polyline(pts).close().extrude(z1 - z0)


def build():
    import cadquery as cq
    cells, traces, meta = load()
    rout = meta["outline_hex_circumradius_mm"]
    hexpts = [(R3(rout * math.cos(math.radians(90 + 60 * i))), R3(rout * math.sin(math.radians(90 + 60 * i)))) for i in range(6)]
    solids = []
    # top copper: discs with slots (wire with a hole), islands, vias
    for c in cells:
        x, y, th = float(c["x_mm"]), float(c["y_mm"]), math.radians(float(c["slot_axis_deg"]))
        Lc, Wc, r = float(c["slot_L_mm"]), float(c["slot_W_mm"]), float(c["disc_r_mm"])

        def pts(hx, hy, ox=0.0, oy=0.0):
            out = []
            for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
                v = L.rot(np.array([ox + sx * hx, oy + sy * hy]), th)
                out.append((R3(x + v[0]), R3(y + v[1])))
            return out

        wp = cq.Workplane("XY").center(R3(x), R3(y)).circle(r).polyline([(p[0] - R3(x), p[1] - R3(y)) for p in pts(Lc / 2, Wc / 2)]).close().extrude(T)
        solids.append(wp.val())
        if c["class"] == "SWITCHED":
            solids.append(cq.Workplane("XY").polyline(pts(ISL, 0.65, 0.0, -0.20)).close().extrude(T).val())
            for key in (("island_via_x", "island_via_y"), ("disc_via_x", "disc_via_y")):
                vx, vy = float(c[key[0]]), float(c[key[1]])
                solids.append(cq.Workplane("XY").workplane(offset=Z_BOT).center(vx, vy).circle(L.PAD_R).extrude(Z_FR4B - Z_BOT).val())  # bottom pad
                solids.append(cq.Workplane("XY").workplane(offset=Z_BOT).center(vx, vy).circle(0.15).extrude(0.0 - Z_BOT).val())  # barrel (touches the top copper at z=0)
    top = cq.Compound.makeCompound(solids)
    # ground with anti-pads and the probe opening
    gw = cq.Workplane("XY").workplane(offset=Z_GB).polyline(hexpts).close()
    holes = [(0.0, 0.0, R_OUT)]
    for c in cells:
        if c["class"] == "SWITCHED":
            holes += [(float(c["island_via_x"]), float(c["island_via_y"]), 0.5), (float(c["disc_via_x"]), float(c["disc_via_y"]), 0.5)]
    for hx, hy, hr in holes:
        gw = gw.moveTo(hx, hy).circle(hr)
    ground = gw.extrude(GT)
    # bottom copper traces (one union per net) + connector pads
    bt = []
    for tr in traces:
        poly = tr["poly"]
        u = None
        for a, b in zip(poly[:-1], poly[1:]):
            s = seg_box(cq, a, b, L.W_TR, Z_BOT, Z_FR4B)
            if s is None:
                continue
            u = s if u is None else u.union(s)
        if u is not None:
            bt.append(u.val())
    pads = []
    for k in range(24):
        py = (11.5 - k) * 2.54
        pads.append(cq.Workplane("XY").workplane(offset=Z_BOT).center(meta["connector_x_mm"], R3(py)).rect(1.0, 1.0).extrude(Z_FR4B - Z_BOT).val())
    pin = cq.Workplane("XY").workplane(offset=Z_BOT).circle(R_PIN).extrude(Z_PIN_TOP - Z_BOT).val()
    metal = cq.Compound.makeCompound([top, ground.val(), pin] + bt + pads)
    sub = cq.Workplane("XY").workplane(offset=-H).polyline(hexpts).close().extrude(H).val()
    fr4 = cq.Workplane("XY").workplane(offset=Z_FR4B).polyline(hexpts).close().extrude(Z_GB - Z_FR4B).val()
    ax = rout * math.cos(math.radians(30)) + PAD_XY
    ay = rout + PAD_XY
    air = cq.Workplane("XY").workplane(offset=Z_AIR0).center(0, 0).rect(2 * R3(ax), 2 * R3(ay)).extrude(Z_LID - Z_AIR0).val()
    parts = {"METAL_PEC": metal, "SUBSTRATE": sub, "SUBSTRATE_FR4": fr4, "AIR": air}
    shs, rows = [], []
    for c in cells:
        if c["class"] != "SWITCHED":
            continue
        x, y, th = float(c["x_mm"]), float(c["y_mm"]), math.radians(float(c["slot_axis_deg"]))
        q = []
        for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
            v = L.rot(np.array([sx * SH_ALONG / 2, 0.45 + (sy + 1) / 2 * SH_ACROSS]), th)
            q.append(cq.Vector(R3(x + v[0]), R3(y + v[1]), Z_SH))
        shs.append(cq.Face.makeFromWires(cq.Wire.makePolygon(q, close=True)))
        d = L.rot(np.array([0.0, 1.0]), th)
        cen = np.array([x, y]) + L.rot(np.array([0.0, 0.95]), th)
        rows.append([c["site"], R3(cen[0]), R3(cen[1]), Z_SH, round(d[0], 6), round(d[1], 6), c["string"], c["pos_in_string"]])
    parts["DIODE_SHEET"] = cq.Compound.makeCompound(shs)
    ang = cq.Workplane("XY").workplane(offset=Z_FR4B).circle(R_OUT).circle(R_PIN).extrude(0.0).val() if False else None
    # port annulus sheet at z = Z_FR4B (bottom of the FR4, the coax opening r 0.635 ... 2.05)
    ring = cq.Face.makeFromWires(cq.Wire.makeCircle(R_OUT, cq.Vector(0, 0, Z_FR4B), cq.Vector(0, 0, 1)), [cq.Wire.makeCircle(R_PIN, cq.Vector(0, 0, Z_FR4B), cq.Vector(0, 0, 1))])
    parts["PORT_ANNULUS"] = ring
    return parts, rows, hexpts


def main():
    import cadquery as cq
    from src.transmitarray.export import inexact_literals, snap_step_text, verify_snap
    from sub6_step_exact import audit_step_text, exactify_step_text
    from sidequest_board_k2_lp_package import read_names
    OUT.mkdir(parents=True, exist_ok=True)
    parts, rows, hexpts = build()
    with open(OUT / "diode_sheets.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["site", "sheet_centre_x", "sheet_centre_y", "sheet_z", "current_dir_x", "current_dir_y", "string", "pos_in_string"])
        w.writerows(rows)
    cols = {"METAL_PEC": (0.9, 0.7, 0.2), "SUBSTRATE": (0.25, 0.55, 0.25), "SUBSTRATE_FR4": (0.2, 0.5, 0.2), "AIR": (0.7, 0.8, 1.0), "DIODE_SHEET": (0.9, 0.1, 0.1), "PORT_ANNULUS": (0.1, 0.1, 0.9)}
    asm = cq.Assembly()
    for nm, shp in parts.items():
        asm.add(shp, name=nm, color=cq.Color(*cols[nm], 1))
    path = OUT / "BOARD_N240.step"
    asm.export(str(path), exportType="STEP")
    raw = path.read_text(encoding="utf-8")
    snapped, ex = exactify_step_text(raw)
    path.write_text(snapped, encoding="utf-8", newline="\n")
    text = path.read_text(encoding="utf-8")
    bad = audit_step_text(text)
    names = read_names(path)
    chk = dict(bodies=sorted(names), n_METAL_solids=len(names["METAL_PEC"].Solids()), n_diode_sheets=len(rows), exactify=ex, audit_inexact_literals=len(bad), audit_examples=[str(b) for b in bad[:5]],
               step_size_MB=round(path.stat().st_size / 1e6, 2), all_pass=bool(len(bad) == 0 and {"METAL_PEC", "SUBSTRATE", "SUBSTRATE_FR4", "AIR", "DIODE_SHEET", "PORT_ANNULUS"} <= set(names)))
    (OUT / "geometry_checks.json").write_text(json.dumps({"label": "MODEL_ONLY geometry checks", **chk}, indent=1, default=str))
    print(chk, flush=True)


if __name__ == "__main__":
    main()
