"""SUB6 FINAL BOARD, 13-pad variant: build the KiCad board (.kicad_pcb + .kicad_pro + .kicad_dru) from the script-generated RF geometry.  MODEL_ONLY / NOT_HFSS_VALIDATED.
    SUB6_MESH=96 python sub6_final_kicad_build.py            (any python with numpy)
Inputs : results/sub6_final_board/mesh96/fab/{cell_table.csv, assembly/diode_placement.csv}, routing.json (outline, connector x)
Output : results/sub6_final_board/mesh96_13pad/kicad/sub6_board.{kicad_pcb,kicad_pro,kicad_dru}, nets.json (net table, pad assignment)
Layers : F.Cu cells (custom pads, locked), In1.Cu ground (zone, anti-pad = via land 0.3 + 0.2 clearance = 0.5), In2.Cu unused (KiCad needs an even copper count), B.Cu bias (routed by Freerouting).
Frame  : KiCad y axis points down, so y_kicad = -y_model (the Gerber export flips it back; the physical board is identical).
Nets   : S<string> (drive, anode end = island of position 1), N<string>_<i> (disc of position i + island of position i+1, i = 1..19), RETURN (disc of position 20 of all 12 strings + RETURN pad), GND.
Connector: 13 pads 1.0 x 1.0, pitch 2.54, RETURN in the middle (index 7 from the top); S pads ordered by the angle of the string start seen from the connector (planar fan).
"""
import csv
import json
import math
import os
import sys
import uuid
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "results" / "sub6_final_board" / "mesh96"
OUT = ROOT / "results" / "sub6_final_board" / "mesh96_13pad" / "kicad"
PITCH, N_PAD, RET_K = 2.54, 13, 6
VIA_D, VIA_DRILL = 0.6, 0.3
ISL = (2.174, 0.65, 0.0, -0.20)  # half length, half width, local offset x, y
UID = lambda: str(uuid.uuid4())
F = lambda v: f"{v:.6f}".rstrip("0").rstrip(".") if abs(v) > 1e-9 else "0"


def rot(v, a):
    c, s = math.cos(a), math.sin(a)
    return (c * v[0] - s * v[1], s * v[0] + c * v[1])


def kx(p):
    return (p[0], -p[1])


def arc_pts(R, a0, a1, step_deg=2.0):
    n = max(2, int(math.ceil(abs(a1 - a0) / math.radians(step_deg))))
    return [(R * math.cos(a0 + (a1 - a0) * i / n), R * math.sin(a0 + (a1 - a0) * i / n)) for i in range(n + 1)]


def disc_pieces(R, a, b):
    """disc minus the slot |x|<=a, |y|<=b as four convex polygons (local, model frame)."""
    t0 = math.asin(b / R)
    top = arc_pts(R, t0, math.pi - t0)  # right -> left along the upper arc
    bot = arc_pts(R, -t0, -(math.pi - t0))
    b2 = b + 0.05  # side pieces overlap the top/bottom caps by 0.05 so that the union is a single polygon
    t2 = math.asin(b2 / R)
    right = [(a, -b2), (a, b2)] + arc_pts(R, t2, -t2)
    left = [(-a, b2), (-a, -b2)] + arc_pts(R, math.pi + t2, math.pi - t2)
    return [top, bot, right, left]


def poly(pts, anchor, th, flip=True):
    """local model points -> pad-relative kicad points (pad anchor given in model absolute coords, th = rotation)."""
    out = []
    for p in pts:
        q = rot(p, th)
        out.append((q[0], -q[1]))
    return out


def gr_poly(pts):
    return "(gr_poly (pts " + " ".join(f"(xy {F(x)} {F(y)})" for x, y in pts) + ") (width 0) (fill yes))"


class Nets:
    def __init__(self):
        self.names = {"": 0}

    def get(self, n):
        if n not in self.names:
            self.names[n] = len(self.names)
        return self.names[n]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cells = list(csv.DictReader(open(SRC / "fab" / "cell_table.csv")))
    dio = {r["ref"]: r for r in csv.DictReader(open(SRC / "fab" / "assembly" / "diode_placement.csv"))}
    rt = json.load(open(SRC / "routing.json"))
    outline = rt["outline"]
    conn_x = rt["conn"][0]
    nets = Nets()
    gnd = nets.get("GND")
    # ---- strings and net names per cell
    by = {}
    for c in cells:
        if c["class"] == "SWITCHED":
            by[(int(c["string"]), int(c["pos_in_string"]))] = c
    sids = sorted({s for s, _ in by})
    start = {s: (float(by[(s, 1)]["island_via_x"]), float(by[(s, 1)]["island_via_y"])) for s in sids}
    ang = sorted(sids, key=lambda s: -math.atan2(start[s][1], start[s][0] - conn_x))
    spads = [k for k in range(N_PAD) if k != RET_K]
    pad_string = {k: ang[i] for i, k in enumerate(spads)}
    isl_net = lambda s, p: f"S{s:02d}" if p == 1 else f"N{s:02d}_{p - 1:02d}"
    dsc_net = lambda s, p: "RETURN" if p == 20 else f"N{s:02d}_{p:02d}"
    lines = []
    A = lines.append
    # ---- footprints: cells
    fps = []
    for c in cells:
        x, y, th = float(c["x_mm"]), float(c["y_mm"]), math.radians(float(c["slot_axis_deg"]))
        R, Lh, Wh = float(c["disc_r_mm"]), float(c["slot_L_mm"]) / 2, float(c["slot_W_mm"]) / 2
        sw = c["class"] == "SWITCHED"
        pads = []
        if sw:
            s, p = int(c["string"]), int(c["pos_in_string"])
            nd, ni = dsc_net(s, p), isl_net(s, p)
            vd = (float(c["disc_via_x"]), float(c["disc_via_y"]))
            vi = (float(c["island_via_x"]), float(c["island_via_y"]))
        else:
            nd = ni = ""
            vd = None
        # disc pad anchored at the disc via (switched) or at a point inside the upper cap (static)
        anc = vd if vd else (x + rot((0.0, (Wh + R) / 2), th)[0], y + rot((0.0, (Wh + R) / 2), th)[1])
        prims = []
        for pc in disc_pieces(R, Lh, Wh):
            pts = []
            for q in pc:
                m = rot(q, th)
                pts.append((x + m[0] - anc[0], -(y + m[1] - anc[1])))
            prims.append(gr_poly(pts))
        net = f' (net {nets.get(nd)} "{nd}")' if nd else ""
        pads.append(f'(pad "1" smd custom (at {F(anc[0])} {F(-anc[1])}) (size 0.3 0.3) (layers "F.Cu"){net} (options (clearance outline) (anchor circle)) (primitives {" ".join(prims)}) (uuid "{UID()}"))')
        if sw:
            pts = []
            for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
                m = rot((ISL[2] + sx * ISL[0], ISL[3] + sy * ISL[1]), th)
                pts.append((x + m[0] - vi[0], -(y + m[1] - vi[1])))
            pads.append(f'(pad "2" smd custom (at {F(vi[0])} {F(-vi[1])}) (size 0.3 0.3) (layers "F.Cu") (net {nets.get(ni)} "{ni}") (options (clearance outline) (anchor circle)) (primitives {gr_poly(pts)}) (uuid "{UID()}"))')
        fps.append(f'(footprint "SUB6:CELL" (locked yes) (layer "F.Cu") (uuid "{UID()}") (at 0 0) (property "Reference" "C{c["site"]}" (at 0 0) (layer "F.Fab") (hide yes) (effects (font (size 0.5 0.5)))) (attr smd exclude_from_pos_files exclude_from_bom) {" ".join(pads)})')
    # ---- diodes (SC-79: 0.35 x 0.35 lands, 1.35 pitch; A = anode on island, K = cathode on disc)
    for (s, p), c in sorted(by.items()):
        ref = f"D{s}-{p:02d}"
        r = dio[ref]
        cx, cy = float(r["x_mm"]), float(r["y_mm"])
        th = float(r["slot_axis_deg"])
        pa = f'(pad "A" smd rect (at {F(float(r["anode_x"]) - cx)} {F(-(float(r["anode_y"]) - cy))} {F((th) % 360)}) (size 0.35 0.35) (layers "F.Cu" "F.Paste") (net {nets.get(isl_net(s, p))} "{isl_net(s, p)}") (uuid "{UID()}"))'
        pk = f'(pad "K" smd rect (at {F(float(r["cathode_x"]) - cx)} {F(-(float(r["cathode_y"]) - cy))} {F((th) % 360)}) (size 0.35 0.35) (layers "F.Cu" "F.Paste") (net {nets.get(dsc_net(s, p))} "{dsc_net(s, p)}") (uuid "{UID()}"))'
        fps.append(f'(footprint "SUB6:SC79" (locked yes) (layer "F.Cu") (uuid "{UID()}") (at {F(cx)} {F(-cy)}) (property "Reference" "{ref}" (at 0 0) (layer "F.Fab") (hide yes) (effects (font (size 0.5 0.5)))) (attr smd) {pa} {pk})')
    # ---- connector J1 (B.Cu pads)
    cp = []
    pad_xy = []
    for k in range(N_PAD):
        py = (RET_K - k) * PITCH
        nm = "RETURN" if k == RET_K else f"S{pad_string[k]:02d}"
        pad_xy.append((conn_x, py, nm))
        cp.append(f'(pad "{k + 1}" smd rect (at {F(conn_x)} {F(-py)}) (size 1 1) (layers "B.Cu") (net {nets.get(nm)} "{nm}") (uuid "{UID()}"))')
    conn_fp = f'(footprint "SUB6:EDGE_CONN13" (locked yes) (layer "B.Cu") (uuid "{UID()}") (at 0 0) (property "Reference" "J1" (at 0 0) (layer "B.Fab") (hide yes) (effects (font (size 0.5 0.5)) (justify mirror))) (attr smd) {" ".join(cp)})'
    fps.append(conn_fp)
    probe_fp = f'(footprint "SUB6:PROBE" (locked yes) (layer "F.Cu") (uuid "{UID()}") (at 0 0) (property "Reference" "P1" (at 0 0) (layer "F.Fab") (hide yes) (effects (font (size 0.5 0.5)))) (attr exclude_from_pos_files exclude_from_bom) (pad "" np_thru_hole circle (at 0 0) (size 1.27 1.27) (drill 1.27) (layers "*.Cu") (uuid "{UID()}")))'
    fps.append(probe_fp)
    # SMA jack, NO legs, flat on the bottom (FR4) side: body 6.35 x 6.35 mm soldered to a B.Cu ground pad = body + 0.5 mm per side (7.35 x 7.35) with a circular clearance r 2.05 around the centre pin
    # (the pin hole 1.27 mm is the existing probe NPTH); ring of 12 ground vias (R 2.7 mm, pitch 1.41 mm) tie the pad to the In1.Cu ground (default: BLIND B.Cu-In1.Cu, FR4 only).
    PAD_H, RO = 3.675, 2.05
    def quadrant(sgx, sgy):
        arc = arc_pts(RO, 0.0, math.pi / 2)
        pts = [(RO, 0.0), (PAD_H, 0.0), (PAD_H, PAD_H), (0.0, PAD_H), (0.0, RO)] + [(x, y) for x, y in arc[::-1]][1:-1]
        return [(sgx * x, sgy * y) for x, y in pts]
    prims = []
    for sgx, sgy in ((1, 1), (-1, 1), (-1, -1), (1, -1)):
        q = quadrant(sgx, sgy)
        cx_, cy_ = sgx * 0.05, sgy * 0.05    # overlap the neighbours by 0.05 over copper so that the union is one polygon
        q2 = [(x + (0 if abs(x) > 1e-9 else -sgx * 0.05) * 1, y + (0 if abs(y) > 1e-9 else -sgy * 0.05)) for x, y in q]
        prims.append(gr_poly([(x - 3.0, -y) for x, y in q2]))   # primitives are relative to the pad anchor at (3, 0)
    sma_pad = f'(pad "G" smd custom (at 3 0) (size 0.4 0.4) (layers "B.Cu") (net {gnd} "GND") (options (clearance outline) (anchor circle)) (primitives {" ".join(prims)}) (uuid "{UID()}"))'
    sma_pad = sma_pad.replace("(at 3 0)", "(at 3 0)")
    body_outline = f'(fp_rect (start -3.175 -3.175) (end 3.175 3.175) (stroke (width 0.1) (type default)) (fill no) (layer "B.Fab") (uuid "{UID()}"))'
    sma_fp = f'(footprint "SUB6:SMA_FLAT_PAD_6.35" (locked yes) (layer "B.Cu") (uuid "{UID()}") (at 0 0) (property "Reference" "J2" (at 0 0) (layer "B.Fab") (hide yes) (effects (font (size 0.5 0.5)) (justify mirror))) (attr smd) {body_outline} {sma_pad})'
    fps.append(sma_fp)
    SMA_VIA_TYPE = os.environ.get("SMA_VIA", "blind")
    sma_vias = []
    for k in range(12):
        ang_ = 2 * math.pi * k / 12 + math.pi / 12
        vx_, vy_ = 2.7 * math.cos(ang_), 2.7 * math.sin(ang_)
        lay = '"In1.Cu" "B.Cu"' if SMA_VIA_TYPE == "blind" else '"F.Cu" "B.Cu"'
        sma_vias.append(f'(via {"blind " if SMA_VIA_TYPE == "blind" else ""}(at {F(vx_)} {F(-vy_)}) (size {VIA_D}) (drill {VIA_DRILL}) (layers {lay}) (locked yes) (net {gnd}) (uuid "{UID()}"))')
    # ---- vias
    vias = []
    for c in cells:
        if c["class"] != "SWITCHED":
            continue
        s, p = int(c["string"]), int(c["pos_in_string"])
        for (vx, vy), nm in (((float(c["island_via_x"]), float(c["island_via_y"])), isl_net(s, p)), ((float(c["disc_via_x"]), float(c["disc_via_y"])), dsc_net(s, p))):
            vias.append(f'(via (at {F(vx)} {F(-vy)}) (size {VIA_D}) (drill {VIA_DRILL}) (layers "F.Cu" "B.Cu") (locked yes) (net {nets.get(nm)}) (uuid "{UID()}"))')
    vias += sma_vias
    # ---- zones: ground (In1.Cu) hexagon; keepouts: probe (B.Cu r 4.5, In1.Cu r 2.05)
    hexpts = " ".join(f"(xy {F(px)} {F(-py)})" for px, py in outline)
    zones = [f'(zone (net {gnd}) (net_name "GND") (layer "In1.Cu") (uuid "{UID()}") (name "GND_PLANE") (hatch edge 0.5) (connect_pads yes (clearance 0.2)) (min_thickness 0.2) (filled_areas_thickness no) (fill yes (thermal_gap 0.5) (thermal_bridge_width 0.5)) (polygon (pts {hexpts})))']

    def circle(r, n=64):
        return " ".join(f"(xy {F(r / math.cos(math.pi / n) * math.cos(2 * math.pi * i / n))} {F(r / math.cos(math.pi / n) * math.sin(2 * math.pi * i / n))})" for i in range(n))

    for layer, r, nm, keep in (("B.Cu", None, "SMA_KEEPOUT_B", "(tracks not_allowed) (vias allowed) (pads allowed) (copperpour allowed) (footprints allowed)"), ("In1.Cu", 2.05, "PROBE_OPENING_GND", "(tracks not_allowed) (vias not_allowed) (pads allowed) (copperpour not_allowed) (footprints allowed)")):
        zones.append(f'(zone (net 0) (net_name "") (layer "{layer}") (uuid "{UID()}") (name "{nm}") (hatch edge 0.5) (connect_pads (clearance 0)) (min_thickness 0.2) (keepout {keep}) (fill (thermal_gap 0.5) (thermal_bridge_width 0.5)) (polygon (pts {circle(r) if r else "(xy -5.175 -5.175) (xy 5.175 -5.175) (xy 5.175 5.175) (xy -5.175 5.175)"})))')
    edge = "\n".join(f'(gr_line (start {F(a[0])} {F(-a[1])}) (end {F(b[0])} {F(-b[1])}) (stroke (width 0.1) (type default)) (layer "Edge.Cuts") (uuid "{UID()}"))' for a, b in zip(outline, outline[1:] + outline[:1]))
    nl = "\n".join(f'(net {i} "{n}")' for n, i in sorted(nets.names.items(), key=lambda kv: kv[1]))
    hdr = '''(kicad_pcb (version 20240108) (generator "sub6_final_kicad_build") (generator_version "8.0")
(general (thickness 4.08) (legacy_teardrops no))
(paper "A0")
(layers (0 "F.Cu" signal) (1 "In1.Cu" signal) (2 "In2.Cu" signal) (31 "B.Cu" signal) (32 "B.Adhes" user "B.Adhesive") (33 "F.Adhes" user "F.Adhesive") (34 "B.Paste" user) (35 "F.Paste" user) (36 "B.SilkS" user "B.Silkscreen") (37 "F.SilkS" user "F.Silkscreen") (38 "B.Mask" user) (39 "F.Mask" user) (40 "Dwgs.User" user "User.Drawings") (41 "Cmts.User" user "User.Comments") (44 "Edge.Cuts" user) (45 "Margin" user) (46 "B.CrtYd" user "B.Courtyard") (47 "F.CrtYd" user "F.Courtyard") (48 "B.Fab" user) (49 "F.Fab" user))
(setup (pad_to_mask_clearance 0) (allow_soldermask_bridges_in_footprints no) (pcbplotparams (layerselection 0x00010fc_ffffffff) (plot_on_all_layers_selection 0x0000000_00000000) (disableapertmacros no) (usegerberextensions no) (usegerberattributes yes) (usegerberadvancedattributes yes) (creategerberjobfile no) (dashed_line_dash_ratio 12.0) (dashed_line_gap_ratio 3.0) (svgprecision 4) (plotframeref no) (viasonmask no) (mode 1) (useauxorigin no) (hpglpennumber 1) (hpglpenspeed 20) (hpglpendiameter 15.0) (pdf_front_fp_property_popups yes) (pdf_back_fp_property_popups yes) (dxfpolygonmode yes) (dxfimperialunits no) (dxfusepcbnewfont yes) (psnegative no) (psa4output no) (plotreference no) (plotvalue no) (plotfptext no) (plotinvisibletext no) (sketchpadsonfab no) (subtractmaskfromsilk no) (outputformat 1) (mirror no) (drillshape 0) (scaleselection 1) (outputdirectory "")))
'''
    txt = hdr + nl + "\n" + "\n".join(fps) + "\n" + "\n".join(vias) + "\n" + "\n".join(zones) + "\n" + edge + "\n)\n"
    (OUT / "sub6_board.kicad_pcb").write_text(txt, encoding="utf-8", newline="\n")
    vp = []
    for c in [c for c in cells if c["class"] == "SWITCHED"]:
        s_, p_ = int(c["string"]), int(c["pos_in_string"])
        for tag, (vx, vy), nm in (("I", (float(c["island_via_x"]), float(c["island_via_y"])), isl_net(s_, p_)), ("D", (float(c["disc_via_x"]), float(c["disc_via_y"])), dsc_net(s_, p_))):
            vp.append(f'(pad "{s_}-{p_}{tag}" thru_hole circle (at {F(vx)} {F(-vy)}) (size {VIA_D + 2 * 0.222:.3f} {VIA_D + 2 * 0.222:.3f}) (drill {VIA_DRILL}) (layers "*.Cu") (net {nets.get(nm)} "{nm}") (uuid "{UID()}"))')
    vias_fp = f'(footprint "SUB6:VIAPINS" (locked yes) (layer "F.Cu") (uuid "{UID()}") (at 0 0) (property "Reference" "V1" (at 0 0) (layer "F.Fab") (hide yes) (effects (font (size 0.5 0.5)))) (attr through_hole) {" ".join(vp)})'
    # staged routing: ACTIVE letters (N links, S drive leads, R return) decide which nets keep their pins' net; previous tracks (json) are added as locked B.Cu segments
    active = sys.argv[1] if len(sys.argv) > 1 else "NSR"
    prev = json.load(open(sys.argv[2])) if len(sys.argv) > 2 else []
    pins = conn_fp + "\n" + vias_fp
    for nm, i in nets.names.items():
        if not nm or nm == "GND":
            continue
        kind = "R" if nm == "RETURN" else ("S" if nm.startswith("S") else "N")
        if kind not in active:
            pins = pins.replace(f' (net {i} "{nm}")', "")
    seg = "\n".join(f'(segment (start {F(t["a"][0])} {F(-t["a"][1])}) (end {F(t["b"][0])} {F(-t["b"][1])}) (width {F(t["w"])}) (layer "B.Cu") (locked yes) (net {nets.names[t["net"]]}) (uuid "{UID()}"))' for t in prev)
    txt_r = hdr + nl + "\n" + pins + "\n" + probe_fp + "\n" + zones[1] + "\n" + seg + "\n" + edge + "\n)\n"
    (OUT / "sub6_board_route.kicad_pcb").write_text(txt_r, encoding="utf-8", newline="\n")
    # ---- project (net classes) + custom rules
    def cls(name, clr, w):
        return dict(name=name, clearance=clr, track_width=w, via_diameter=VIA_D, via_drill=VIA_DRILL, microvia_diameter=0.3, microvia_drill=0.1, diff_pair_gap=0.25, diff_pair_via_gap=0.25, diff_pair_width=0.2, bus_width=12, wire_width=6, line_style=0, pcb_color="rgba(0, 0, 0, 0.000)", schematic_color="rgba(0, 0, 0, 0.000)", priority=0)

    pro = dict(meta=dict(filename="sub6_board.kicad_pro", version=3),
               board=dict(design_settings=dict(rules=dict(min_clearance=0.15, min_copper_edge_clearance=0.9, min_track_width=0.15, min_via_diameter=0.5, min_through_hole_diameter=0.3, min_hole_clearance=0.25, min_hole_to_hole=0.25, min_connection=0.15, min_via_annular_width=0.1),
                                           defaults=dict(copper_line_width=0.25, board_outline_line_width=0.1), rule_severities=dict(silk_over_copper="ignore", silk_overlap="ignore", silk_edge_clearance="ignore", solder_mask_bridge="ignore", via_dangling="ignore", track_dangling="warning", lib_footprint_issues="ignore", lib_footprint_mismatch="ignore", footprint_symbol_mismatch="ignore", footprint_type_mismatch="ignore", courtyards_overlap="ignore", missing_courtyard="ignore", footprint="ignore", extra_footprint="ignore", missing_footprint="ignore", duplicate_footprints="ignore", isolated_copper="warning", starved_thermal="ignore", copper_sliver="ignore", solder_mask_bridge_=""))),
               net_settings=dict(classes=[cls("Default", 0.15, 0.25), cls("BIAS", 0.15, 0.25), cls("RETURN", 0.15, 0.5), cls("GND", 0.2, 0.25)],
                                 netclass_patterns=[dict(netclass="BIAS", pattern="S*"), dict(netclass="BIAS", pattern="N*"), dict(netclass="RETURN", pattern="RETURN"), dict(netclass="GND", pattern="GND")], meta=dict(version=3)))
    (OUT / "sub6_board.kicad_pro").write_text(json.dumps(pro, indent=2), encoding="utf-8", newline="\n")
    dru = '''(version 1)
(rule "gnd_antipad_0p2" (condition "A.Type == 'Via' && B.Type == 'Zone'") (constraint clearance (min 0.2)))
(rule "gnd_antipad_0p2b" (condition "A.Type == 'Zone' && B.Type == 'Via'") (constraint clearance (min 0.2)))
(rule "trace_to_foreign_via_antipad_0p172" (condition "A.Type == 'Track' && B.Type == 'Via'") (constraint clearance (min 0.372)))
'''
    (OUT / "sub6_board.kicad_dru").write_text(dru, encoding="utf-8", newline="\n")
    js = dict(nets=nets.names, pad_string={str(k): v for k, v in pad_string.items()}, pads=[dict(index_from_top=k + 1, net=n, x=px, y=py) for k, (px, py, n) in enumerate(pad_xy)], return_index_from_top=RET_K + 1)
    (OUT / "nets.json").write_text(json.dumps(js, indent=1), encoding="utf-8", newline="\n")
    print("wrote", OUT / "sub6_board.kicad_pcb", len(txt) / 1e6, "MB; nets", len(nets.names), "vias", len(vias), "footprints", len(fps))


if __name__ == "__main__":
    main()
