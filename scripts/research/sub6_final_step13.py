"""SUB6 FINAL BOARD, 13-pad variant: board STEP builder (thick and thin) with NO overlapping bodies.  MODEL_ONLY / NOT_HFSS_VALIDATED.
    SUB6_MESH=96 python sub6_final_step13.py thick|thin  [tracks.json]
Reads  : fab/cell_table.csv, fab/assembly/diode_placement.csv (24-pad dir, RF unchanged), kicad/tracks.json (routed B.Cu segments: net, width, a, b), kicad/nets.json (connector pads)
Writes : hfss_exports/sub6_final_board_m96_13pad/BOARD_N240.step (thick) or BOARD_N240_thin.step, diode_sheets.csv, geometry_checks_<variant>.json

Bodies (every name unique):
  METAL_PEC_<id>   one fused solid per electrical net (cell copper + via barrels + via pads + B.Cu traces + connector pad), one per static cell (METAL_PEC_C<site>), METAL_PEC_GROUND, METAL_PEC_PIN
  SUBSTRATE        RT5880 prism with holes for all barrels and the probe pin
  SUBSTRATE_FR4    FR4 prism (including the 0.035 bond band of the ground plane) minus the ground plate, barrels and pin  -> anti-pad voids are filled with FR4, no enclosed air pockets
  AIR              box minus every other solid
  D<string>_<pos>  240 diode sheets (1.00 x 0.35 mm, between island edge and disc edge), PORT_ANNULUS  coax port sheet
Thick: copper 35 um (z 0..0.035 top, ground -3.210..-3.175, bottom -4.045..-4.010), FR4 -4.010..-3.175.
Thin : copper/ground as zero-thickness sheets (top z 0, ground z -3.175, bottom z -3.975); RT5880 -3.175..0, FR4 -3.975..-3.175 (contiguous stack, FR4 0.8); barrels remain 0.3 mm rods.
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
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
SRC = ROOT / "results" / "sub6_final_board" / "mesh96"
KDIR = ROOT / "results" / "sub6_final_board" / "mesh96_13pad" / "kicad"
OUT = ROOT / "hfss_exports" / "sub6_final_board_m96_13pad"
R3 = lambda v: round(float(v) + 0.0, 3)
H, T = 3.175, 0.035
R_PIN, R_OUT, H_PIN = 0.635, 2.05, 12.1
R_BAR, R_PAD, R_ANTI = 0.15, 0.30, 0.50
ISL = (2.174, 0.65, 0.0, -0.20)
SH_ALONG, SH_ACROSS, SH_EXT = 0.35, 1.00, 0.0015
PAD_XY, PAD_TOP, PAD_BOT = 15.0, 40.0, 15.0
SMA_BODY, SMA_PAD = 6.35, 7.35          # SMA body square and B.Cu ground pad (body + 0.5 mm per side); pad opening r 2.05 (= R_OUT)
H_SMA = 6.0                             # SMA body height below the board: ASSUMED (no part number)
N_GV, R_GV = 12, 2.7                    # ring of ground vias (pitch 1.41 mm), via r 0.15
SMA_THROUGH = os.environ.get("SMA_VIA", "blind") == "through"


def levels(thin):
    if thin:
        z = dict(top=0.0, ttop=0.0, gnd_lo=-H, gnd_hi=-H, fr4_hi=-H, fr4_lo=-H - 0.8, bot_lo=-H - 0.8, bot_hi=-H - 0.8, sheet=0.0)
    else:
        z = dict(top=T, ttop=T, gnd_lo=-H - T, gnd_hi=-H, fr4_hi=-H, fr4_lo=-H - T - 0.8, bot_lo=-H - T - 0.8 - T, bot_hi=-H - T - 0.8, sheet=T / 2)
    z["pin_top"] = round(-H + H_PIN, 3)
    z["air_lo"] = round(z["bot_lo"] - PAD_BOT, 3)
    z["air_hi"] = round(T + PAD_TOP, 3)
    z["pin_lo"] = z["bot_lo"] - H_SMA   # the centre pin runs through the SMA bore down to the bottom of the body block
    return z


def rot(v, a):
    c, s = math.cos(a), math.sin(a)
    return (c * v[0] - s * v[1], s * v[0] + c * v[1])


def quad(cx, cy, hx, hy, th, ox=0.0, oy=0.0):
    out = []
    for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
        v = rot((ox + sx * hx, oy + sy * hy), th)
        out.append((R3(cx + v[0]), R3(cy + v[1])))
    return out


def seg_poly(a, b, w):
    dx, dy = b[0] - a[0], b[1] - a[1]
    ln = math.hypot(dx, dy)
    if ln < 1e-9:
        return None
    ux, uy = dx / ln, dy / ln
    e = w / 2
    return [(R3(a[0] - ux * e - uy * w / 2), R3(a[1] - uy * e + ux * w / 2)), (R3(b[0] + ux * e - uy * w / 2), R3(b[1] + uy * e + ux * w / 2)),
            (R3(b[0] + ux * e + uy * w / 2), R3(b[1] + uy * e - ux * w / 2)), (R3(a[0] - ux * e + uy * w / 2), R3(a[1] - uy * e - ux * w / 2))]


def build(thin, tracks_file=None):
    import cadquery as cq
    z = levels(thin)
    cells = list(csv.DictReader(open(SRC / "fab" / "cell_table.csv")))
    dio = {r["ref"]: r for r in csv.DictReader(open(SRC / "fab" / "assembly" / "diode_placement.csv"))}
    rt = json.load(open(SRC / "routing.json"))
    rout = rt["r_out"]
    hexpts = [(R3(rout * math.cos(math.radians(90 + 60 * i))), R3(rout * math.sin(math.radians(90 + 60 * i)))) for i in range(6)]
    tracks = json.load(open(tracks_file or (KDIR / "tracks.json")))
    nets_js = json.load(open(KDIR / "nets.json"))
    isl_net = lambda s, p: f"S{s:02d}" if p == 1 else f"N{s:02d}_{p - 1:02d}"
    dsc_net = lambda s, p: "RETURN" if p == 20 else f"N{s:02d}_{p:02d}"
    t0 = time.time()
    # ---------------- metal pieces per net
    pieces = {}      # net -> list of cq.Shape (solids; or faces in thin mode)
    static = {}
    barrels = []     # barrel solids (for the dielectric holes)
    via_xy = []

    def add(net, shp):
        pieces.setdefault(net, []).append(shp)

    for c in cells:
        x, y, th = float(c["x_mm"]), float(c["y_mm"]), math.radians(float(c["slot_axis_deg"]))
        R, Lh, Wh = float(c["disc_r_mm"]), float(c["slot_L_mm"]) / 2, float(c["slot_W_mm"]) / 2
        z0 = z["top"] - (0.0 if thin else T)
        base = cq.Workplane("XY").workplane(offset=z0).center(R3(x), R3(y)).circle(R).polyline([(p[0] - R3(x), p[1] - R3(y)) for p in quad(x, y, Lh, Wh, th)]).close()
        disc_solid = base.extrude(T).val()
        if thin:
            # zero-thickness copper: the bottom face of the solid at z = top (z0 = 0)
            fcs = [f for f in disc_solid.Faces() if abs(f.Center().z - z0) < 1e-9 and f.normalAt().z < -0.5]
            disc_shape = fcs[0]
        else:
            disc_shape = disc_solid
        if c["class"] != "SWITCHED":
            static[c["site"]] = disc_shape
            continue
        s, p = int(c["string"]), int(c["pos_in_string"])
        add(dsc_net(s, p), disc_shape)
        pts = quad(x, y, ISL[0], ISL[1], th, ISL[2], ISL[3])
        isl = cq.Workplane("XY").workplane(offset=z0).polyline(pts).close()
        isl_shape = isl.extrude(T).val() if not thin else cq.Face.makeFromWires(cq.Wire.makePolygon([cq.Vector(a, b, z0) for a, b in pts], close=True))
        add(isl_net(s, p), isl_shape)
        for (vx, vy), nm in (((float(c["disc_via_x"]), float(c["disc_via_y"])), dsc_net(s, p)), ((float(c["island_via_x"]), float(c["island_via_y"])), isl_net(s, p))):
            vx, vy = R3(vx), R3(vy)
            via_xy.append((vx, vy))
            barrel = cq.Workplane("XY").workplane(offset=z["bot_hi"]).center(vx, vy).circle(R_BAR).extrude(0.0 - z["bot_hi"]).val()   # barrel up to z = 0 (bottom face of the top copper)
            barrels.append(barrel)
            add(nm, barrel)
            if not thin:
                add(nm, cq.Workplane("XY").workplane(offset=z["bot_lo"]).center(vx, vy).circle(R_PAD).extrude(z["bot_hi"] - z["bot_lo"]).val())
            else:
                add(nm, cq.Face.makeFromWires(cq.Wire.makeCircle(R_PAD, cq.Vector(vx, vy, z["bot_lo"]), cq.Vector(0, 0, 1))))
    print(f"cells done {time.time() - t0:.0f}s  nets {len(pieces)} static {len(static)}", flush=True)
    # ---------------- bottom copper traces
    for t in tracks:
        poly = seg_poly(t["a"], t["b"], t["w"])
        if poly is None:
            continue
        if thin:
            add(t["net"], cq.Face.makeFromWires(cq.Wire.makePolygon([cq.Vector(a, b, z["bot_lo"]) for a, b in poly], close=True)))
        else:
            add(t["net"], cq.Workplane("XY").workplane(offset=z["bot_lo"]).polyline(poly).close().extrude(z["bot_hi"] - z["bot_lo"]).val())
        # round the joint between segments with a pad-like disc of the trace width (keeps width, closes the polyline gaps)
        for q in (t["a"], t["b"]):
            if thin:
                add(t["net"], cq.Face.makeFromWires(cq.Wire.makeCircle(t["w"] / 2, cq.Vector(R3(q[0]), R3(q[1]), z["bot_lo"]), cq.Vector(0, 0, 1))))
            else:
                add(t["net"], cq.Workplane("XY").workplane(offset=z["bot_lo"]).center(R3(q[0]), R3(q[1])).circle(t["w"] / 2).extrude(z["bot_hi"] - z["bot_lo"]).val())
    for pd in nets_js["pads"]:
        if thin:
            add(pd["net"], cq.Face.makeFromWires(cq.Wire.makePolygon([cq.Vector(R3(pd["x"]) + dx, R3(pd["y"]) + dy, z["bot_lo"]) for dx, dy in ((-.5, -.5), (.5, -.5), (.5, .5), (-.5, .5))], close=True)))
        else:
            add(pd["net"], cq.Workplane("XY").workplane(offset=z["bot_lo"]).center(R3(pd["x"]), R3(pd["y"])).rect(1.0, 1.0).extrude(z["bot_hi"] - z["bot_lo"]).val())
    print(f"traces done {time.time() - t0:.0f}s", flush=True)
    # ---------------- fuse per net
    metal = {}
    for i, (net, lst) in enumerate(sorted(pieces.items())):
        if len(lst) == 1:
            metal["METAL_PEC_" + net] = lst[0]
        else:
            shp = lst[0].fuse(*lst[1:]).clean() if not thin else cq.Compound.makeCompound(lst)
            metal["METAL_PEC_" + net] = shp
    for site, shp in static.items():
        metal["METAL_PEC_C" + site] = shp
    print(f"fuse done {time.time() - t0:.0f}s  metal bodies {len(metal)}", flush=True)
    # ---------------- ground plate (anti-pads + probe opening)
    gw = cq.Workplane("XY").workplane(offset=z["gnd_lo"]).polyline(hexpts).close()
    for vx, vy in via_xy:
        gw = gw.moveTo(vx, vy).circle(R_ANTI)
    gw = gw.moveTo(0, 0).circle(R_OUT)
    if thin:
        gnd = [f for f in gw.extrude(T).val().Faces() if abs(f.Center().z - z["gnd_lo"]) < 1e-9 and f.normalAt().z < -0.5][0]
        metal["METAL_PEC_GROUND"] = gnd
        ground_solid = None
    else:
        ground_solid = gw.extrude(z["gnd_hi"] - z["gnd_lo"]).val()
        metal["METAL_PEC_GROUND"] = ground_solid
    # ---------------- SMA (flat, no legs): ground pad on B.Cu with opening r 2.05, PEC body block with the coax bore, ring of ground vias (blind B.Cu-In1.Cu = FR4 only; env SMA_VIA=through: to the top)
    sma_rods = []
    z_rod_top = 0.0 if SMA_THROUGH else z["gnd_lo"]
    for k in range(N_GV):
        ang_ = 2 * math.pi * k / N_GV + math.pi / N_GV
        sma_rods.append(cq.Workplane("XY").workplane(offset=z["bot_lo"]).center(R3(R_GV * math.cos(ang_)), R3(R_GV * math.sin(ang_))).circle(R_BAR).extrude(z_rod_top - z["bot_lo"]).val())
    body = cq.Workplane("XY").workplane(offset=z["bot_lo"] - H_SMA).rect(SMA_BODY, SMA_BODY).circle(R_OUT).extrude(H_SMA).val()
    if thin:
        pad_f = [f for f in cq.Workplane("XY").workplane(offset=z["bot_lo"]).rect(SMA_PAD, SMA_PAD).circle(R_OUT).extrude(T).val().Faces() if abs(f.Center().z - z["bot_lo"]) < 1e-9 and f.normalAt().z < -0.5][0]
        metal["METAL_PEC_GROUND"] = cq.Compound.makeCompound([metal["METAL_PEC_GROUND"], pad_f, body] + sma_rods)
    else:
        pad_s = cq.Workplane("XY").workplane(offset=z["bot_lo"]).rect(SMA_PAD, SMA_PAD).circle(R_OUT).extrude(z["bot_hi"] - z["bot_lo"]).val()
        ground_solid = ground_solid.fuse(pad_s, body, *sma_rods).clean()
        metal["METAL_PEC_GROUND"] = ground_solid
    pin = cq.Workplane("XY").workplane(offset=z["pin_lo"]).circle(R_PIN).extrude(z["pin_top"] - z["pin_lo"]).val()
    metal["METAL_PEC_PIN"] = pin
    # ---------------- dielectrics
    holes = barrels + [pin]
    sub = cq.Workplane("XY").workplane(offset=-H).polyline(hexpts).close().extrude(H).val()
    sub = sub.cut(*(holes + (sma_rods if SMA_THROUGH else []))).clean()
    fr4 = cq.Workplane("XY").workplane(offset=z["fr4_lo"]).polyline(hexpts).close().extrude(z["fr4_hi"] - z["fr4_lo"]).val()
    cut_tools = holes + ([ground_solid] if ground_solid is not None else sma_rods)
    fr4 = fr4.cut(*cut_tools).clean()
    print(f"dielectrics done {time.time() - t0:.0f}s", flush=True)
    # ---------------- air = box minus everything solid
    ax, ay = rout * math.cos(math.radians(30)) + PAD_XY, rout + PAD_XY
    air = cq.Workplane("XY").workplane(offset=z["air_lo"]).rect(2 * R3(ax), 2 * R3(ay)).extrude(z["air_hi"] - z["air_lo"]).val()
    solids = [sub, fr4]
    for nm, shp in metal.items():
        if not thin or nm in ("METAL_PEC_PIN",):
            solids.append(shp)
    if thin:
        solids += barrels + sma_rods + [body]
    air = air.cut(*solids).clean()
    print(f"air done {time.time() - t0:.0f}s", flush=True)
    parts = {}
    parts.update(metal)
    parts["SUBSTRATE"] = sub
    parts["SUBSTRATE_FR4"] = fr4
    parts["AIR"] = air
    # ---------------- diode sheets D<string>_<pos>, port annulus
    rows = []
    for c in cells:
        if c["class"] != "SWITCHED":
            continue
        s, p = int(c["string"]), int(c["pos_in_string"])
        x, y, th = float(c["x_mm"]), float(c["y_mm"]), math.radians(float(c["slot_axis_deg"]))
        q = []
        for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
            v = rot((sx * SH_ALONG / 2, 0.45 - SH_EXT + (sy + 1) / 2 * (SH_ACROSS + 2 * SH_EXT)), th)   # 1.5 um into the copper at both ends: the 0.001 mm grid rounding of slanted edges would otherwise leave sub-micron gaps
            q.append(cq.Vector(R3(x + v[0]), R3(y + v[1]), z["sheet"]))
        parts[f"D{s}_{p}"] = cq.Face.makeFromWires(cq.Wire.makePolygon(q, close=True))
        d = rot((0.0, 1.0), th)
        a_pt = (R3(x + rot((0.0, 0.45), th)[0]), R3(y + rot((0.0, 0.45), th)[1]))
        k_pt = (R3(x + rot((0.0, 1.45), th)[0]), R3(y + rot((0.0, 1.45), th)[1]))
        cen = (R3(x + rot((0.0, 0.95), th)[0]), R3(y + rot((0.0, 0.95), th)[1]))
        rows.append([f"D{s}_{p}", cen[0], cen[1], z["sheet"], round(d[0], 6), round(d[1], 6), a_pt[0], a_pt[1], k_pt[0], k_pt[1], s, p, c["site"], round(math.degrees(th) % 360, 3)])
    zp = z["fr4_lo"]
    parts["PORT_ANNULUS"] = cq.Face.makeFromWires(cq.Wire.makeCircle(R_OUT, cq.Vector(0, 0, zp), cq.Vector(0, 0, 1)), [cq.Wire.makeCircle(R_PIN, cq.Vector(0, 0, zp), cq.Vector(0, 0, 1))])
    meta = dict(z=z, thin=thin, n_metal=len(metal), n_sheets=len(rows), air_box=dict(ax=R3(ax), ay=R3(ay)))
    return parts, rows, meta


def main():
    variant = sys.argv[1] if len(sys.argv) > 1 else "thick"
    thin = variant == "thin"
    import cadquery as cq
    from sub6_step_exact import audit_step_text, exactify_step_text
    from sidequest_board_k2_lp_package import read_names
    OUT.mkdir(parents=True, exist_ok=True)
    parts, rows, meta = build(thin, sys.argv[2] if len(sys.argv) > 2 else None)
    with open(OUT / "diode_sheets.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["name", "sheet_centre_x", "sheet_centre_y", "sheet_z", "current_dir_x", "current_dir_y", "anode_edge_x", "anode_edge_y", "cathode_edge_x", "cathode_edge_y", "string", "pos_in_string", "site", "slot_axis_deg"])
        w.writerows(rows)
    asm = cq.Assembly()
    for nm, shp in parts.items():
        col = (0.9, 0.7, 0.2) if nm.startswith("METAL") else {"SUBSTRATE": (0.25, 0.55, 0.25), "SUBSTRATE_FR4": (0.2, 0.5, 0.2), "AIR": (0.7, 0.8, 1.0), "PORT_ANNULUS": (0.1, 0.1, 0.9)}.get(nm, (0.9, 0.1, 0.1))
        asm.add(shp, name=nm, color=cq.Color(*col, 1))
    path = OUT / ("BOARD_N240_thin_clean.step" if thin else "BOARD_N240_clean.step")
    t0 = time.time()
    asm.export(str(path), exportType="STEP")
    raw = path.read_text(encoding="utf-8")
    snapped, ex = exactify_step_text(raw)
    path.write_text(snapped, encoding="utf-8", newline="\n")
    bad = audit_step_text(path.read_text(encoding="utf-8"))
    names = read_names(path)
    chk = dict(label="MODEL_ONLY geometry build", variant=variant, n_bodies=len(names), duplicate_names=len(parts) - len(set(parts)), step_size_MB=round(path.stat().st_size / 1e6, 1), audit_inexact_literals=len(bad), meta=meta, export_s=round(time.time() - t0))
    (OUT / f"geometry_build_{variant}.json").write_text(json.dumps(chk, indent=1, default=str), encoding="utf-8", newline="\n")
    print(chk, flush=True)


if __name__ == "__main__":
    main()
