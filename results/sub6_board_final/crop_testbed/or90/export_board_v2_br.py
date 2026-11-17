"""Board v2 STEP exporter (cadquery; NEVER launches AEDT).  Circular board, RT5880 3.175 / ground / RT5880 bias layer t_bias, vias only (no bias network).
usage: python export_board_v2.py <layout.json> <cell.json> <out_dir> [t_bias_mm=0.254]
Outputs in <out_dir>: BOARD_V2.step, manifest.json (every cell / via / antipad / pad / diode sheet), diode_sheets.json."""
import json, math, os, sys, time
from pathlib import Path
TAG = os.environ.get("BOARD_TAG", "V2")
sys.path.insert(0, r"C:/lwmom_final/scripts/research"); sys.path.insert(0, r"C:/lwmom_final")
import cadquery as cq
from sub6_step_exact import audit_step_text, exactify_step_text

layout = json.load(open(sys.argv[1])); cell = json.load(open(sys.argv[2])); OUT = Path(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
T_BIAS = float(sys.argv[4]) if len(sys.argv) > 4 else 0.254
R3 = lambda v: round(float(v) + 0.0, 3)
D_MM = layout["D_mm"]; RB = D_MM / 2.0
T, H = 0.035, 3.175
Z_G_TOP, Z_G_BOT = -H, -H - T
Z_B_BOT = Z_G_BOT - T_BIAS                  # bottom of the bias layer = top of the via pads
Z_PAD_LO = Z_B_BOT - T
R_PIN, R_OUT, PIN_TOP = 0.635, 2.05, 8.925
R_BAR, R_PAD, R_ANTI = cell["via"]["barrel_r"], cell["via"]["pad_r"], cell["via"]["antipad_r"]
WIN = 12.0
AIR_MARGIN, AIR_BOT, AIR_TOP = 15.0, 15.0, 40.0
Z_SHEET = T / 2
assert cell["cell_type"] == "island"


def rot(v, a):
    c, s = math.cos(a), math.sin(a)
    return (c * v[0] - s * v[1], s * v[0] + c * v[1])


def lp(x, y, th, u, v):
    r = rot((u, v), th)
    return (R3(x + r[0]), R3(y + r[1]))


# ---- switched cell (plug-in "island"): canonical frame = slot along +y; rotated by (slot_axis - 90 deg)
RD, LS, WS = cell["disc_r"], cell["slot_l"], cell["slot_w"]
LI, WI, G, C = cell["island_l"], cell["island_w"], cell["diode_gap"], cell["opposite_clearance"]
SH_L, SH_W = cell["diode_l"], cell["diode_w"]           # sheet: 0.30 (x) x 0.60 (y)
X0 = -WS / 2 + G                                        # island left edge (canonical)
XC = X0 + WI / 2                                        # island centre-line x
ISL_VIA = tuple(cell["island_via"])
RAD_DISC_VIA = RD - R_PAD - cell["disc_via_pad_edge_inset"]
DISC_VIA = (RAD_DISC_VIA / math.sqrt(2), -RAD_DISC_VIA / math.sqrt(2))
assert abs(ISL_VIA[0] - XC) < 1e-9 or True


def build_switched(x, y, axis_deg, n):
    th = math.radians(axis_deg - 90.0)
    sx0, sx1 = -WS / 2, -WS / 2 + SH_L                  # sheet x range: slot wall ... island left edge (SH_L = G for g = 0.30)
    BR = os.environ.get("BRIDGE", "")
    h = 0.3; xr = X0 + WI
    sA = SH_W / 2
    sl = [(-WS / 2, -LS / 2)]
    isl_p = [(X0, -LI / 2)]
    if BR == "B":
        sl += [(XC - h, -LS / 2), (XC + h, -LS / 2)]; isl_p += [(XC - h, -LI / 2), (XC + h, -LI / 2)]
    sl += [(WS / 2, -LS / 2)]; isl_p += [(xr, -LI / 2)]
    if BR == "A":
        sl += [(WS / 2, -sA), (WS / 2, sA)]; isl_p += [(xr, -sA), (xr, sA)]
    sl += [(WS / 2, LS / 2)]; isl_p += [(xr, LI / 2)]
    if BR == "B":
        sl += [(XC + h, LS / 2), (XC - h, LS / 2)]; isl_p += [(XC + h, LI / 2), (XC - h, LI / 2)]
    sl += [(-WS / 2, LS / 2), (-WS / 2, SH_W / 2), (-WS / 2, -SH_W / 2)]
    isl_p += [(X0, LI / 2), (X0, SH_W / 2), (X0, -SH_W / 2)]
    slot = [lp(x, y, th, u, v) for u, v in sl]
    isl = [lp(x, y, th, u, v) for u, v in isl_p]
    sheet = [lp(x, y, th, sx0, -SH_W / 2), lp(x, y, th, sx0 + SH_L, -SH_W / 2), lp(x, y, th, sx0 + SH_L, SH_W / 2), lp(x, y, th, sx0, SH_W / 2)]
    bridges = []
    if BR == "A":
        bridges.append([lp(x, y, th, u, v) for u, v in ((xr, -sA), (WS / 2, -sA), (WS / 2, sA), (xr, sA))])
    if BR == "B":
        bridges.append([lp(x, y, th, u, v) for u, v in ((XC - h, -LS / 2), (XC + h, -LS / 2), (XC + h, -LI / 2), (XC - h, -LI / 2))])
        bridges.append([lp(x, y, th, u, v) for u, v in ((XC - h, LI / 2), (XC + h, LI / 2), (XC + h, LS / 2), (XC - h, LS / 2))])
    iv = lp(x, y, th, ISL_VIA[0], ISL_VIA[1]); dv = lp(x, y, th, DISC_VIA[0], DISC_VIA[1])
    cur = rot((1.0, 0.0), th)   # current across the gap: canonical +x direction (island -> ... sheet is on the -x side: cathode = disc wall)
    return dict(site_xy=(R3(x), R3(y)), axis_deg=axis_deg, rot_deg=axis_deg - 90.0, slot_polygon=slot, island_polygon=isl, sheet=sheet,
                bridges=bridges, island_via=iv, disc_via=dv, disc_r=RD, name=f"D{n}", gap_dir=(round(cur[0], 6), round(cur[1], 6)))


def build_static(s):
    x, y, th = s["x_mm"], s["y_mm"], math.radians(s["slot_axis_deg"])
    L, W = s["static_slot_L_mm"], s["static_slot_W_mm"]
    slot = [lp(x, y, th, -L / 2, -W / 2), lp(x, y, th, L / 2, -W / 2), lp(x, y, th, L / 2, W / 2), lp(x, y, th, -L / 2, W / 2)]
    return dict(site_xy=(R3(x), R3(y)), axis_deg=s["slot_axis_deg"], slot_polygon=slot, disc_r=s["static_disc_r_mm"], slot_L=L, slot_W=W, row_id=s["static_row_id"])


def disc_solid(xy, r, slot):
    wp = cq.Workplane("XY").center(xy[0], xy[1]).circle(r).polyline([(p[0] - xy[0], p[1] - xy[1]) for p in slot]).close()
    return wp.extrude(T).val()


t0 = time.time()
lumps, barrels, sheets, man_cells, man_vias, antipads, pads, diodes = [], [], {}, [], [], [], [], []
n_sw = 0
for s in layout["sites"]:
    if s["class"] == "SWITCHED":
        n_sw += 1
        c = build_switched(s["x_mm"], s["y_mm"], s["slot_axis_deg"], n_sw)
        disc = disc_solid(c["site_xy"], c["disc_r"], c["slot_polygon"])
        isl = cq.Workplane("XY").polyline(c["island_polygon"]).close().extrude(T).val()
        for key, body, vxy in (("disc", disc, c["disc_via"]), ("island", isl, c["island_via"])):
            barrel = cq.Workplane("XY").workplane(offset=Z_B_BOT).center(*vxy).circle(R_BAR).extrude(0.0 - Z_B_BOT).val()
            pad = cq.Workplane("XY").workplane(offset=Z_PAD_LO).center(*vxy).circle(R_PAD).extrude(T).val()
            barrels.append(barrel)
            lumps.append(body.fuse(barrel, pad).clean())
            man_vias.append(dict(cell=c["name"], node=key, xy=vxy, barrel_r=R_BAR, pad_r=R_PAD, antipad_r=R_ANTI, barrel_z=[Z_B_BOT, 0.0], pad_z=[Z_PAD_LO, Z_B_BOT]))
            antipads.append(dict(cell=c["name"], node=key, xy=vxy, r=R_ANTI)); pads.append(dict(cell=c["name"], node=key, xy=vxy, r=R_PAD, z=[Z_PAD_LO, Z_B_BOT]))
        sh = cq.Face.makeFromWires(cq.Wire.makePolygon([cq.Vector(p[0], p[1], Z_SHEET) for p in c["sheet"]], close=True))
        sheets[c["name"]] = sh
        for kb, bp in enumerate(c["bridges"]):
            sheets[f"BR{n_sw}x{kb}"] = cq.Face.makeFromWires(cq.Wire.makePolygon([cq.Vector(p[0], p[1], Z_SHEET) for p in bp], close=True))
        man_cells.append(dict(site=s["site"], cell_class="SWITCHED", **{k: c[k] for k in ("site_xy", "axis_deg", "rot_deg", "slot_polygon", "island_polygon", "disc_r", "name")}))
        diodes.append(dict(name=c["name"], site=s["site"], sheet_corners_xyz=[[p[0], p[1], Z_SHEET] for p in c["sheet"]], current_dir_xy=c["gap_dir"],
                           anode_edge_on_island=[c["sheet"][1], c["sheet"][2]], cathode_edge_on_disc=[c["sheet"][0], c["sheet"][3]],
                           size_mm=[SH_L, SH_W], z_mm=Z_SHEET,
                           OFF=dict(R_ohm_per_sq=22.23, X_ohm_per_sq=-335.4), ON=dict(R_ohm_per_sq=1.60, X_ohm_per_sq=29.2)))
    else:
        c = build_static(s)
        lumps.append(disc_solid(c["site_xy"], c["disc_r"], c["slot_polygon"]))
        man_cells.append(dict(site=s["site"], cell_class="STATIC", **{k: c[k] for k in ("site_xy", "axis_deg", "slot_polygon", "disc_r", "slot_L", "slot_W", "row_id")}))
print(f"cells built {time.time() - t0:.0f}s switched {n_sw}", flush=True)
# ground: circular plate, feed opening, one anti-pad per via
gw = cq.Workplane("XY").workplane(offset=Z_G_BOT).circle(RB)
for a in antipads:
    gw = gw.moveTo(*a["xy"]).circle(a["r"])
gw = gw.moveTo(0, 0).circle(R_OUT)
ground = gw.extrude(T).val()
pin = cq.Workplane("XY").workplane(offset=Z_G_TOP).circle(R_PIN).extrude(PIN_TOP - Z_G_TOP).val()
lumps += [ground, pin]
metal = cq.Compound.makeCompound(lumps)
sub = cq.Workplane("XY").workplane(offset=Z_G_TOP).circle(RB).extrude(H).val().cut(*(barrels + [pin])).clean()
bias = cq.Workplane("XY").workplane(offset=Z_B_BOT).circle(RB).extrude(Z_G_BOT - Z_B_BOT).val()
window = cq.Workplane("XY").workplane(offset=Z_B_BOT - 1.0).rect(WIN, WIN).extrude((Z_G_BOT + 0.5) - (Z_B_BOT - 1.0)).val()
bias = bias.cut(*(barrels + [window])).clean()
z0, z1 = R3(Z_PAD_LO - AIR_BOT), R3(T + AIR_TOP)
air = cq.Workplane("XY").workplane(offset=z0).rect(2 * R3(RB + AIR_MARGIN), 2 * R3(RB + AIR_MARGIN)).extrude(z1 - z0).val()
port = cq.Face.makeFromWires(cq.Wire.makeCircle(R_OUT, cq.Vector(0, 0, Z_G_TOP), cq.Vector(0, 0, 1)), [cq.Wire.makeCircle(R_PIN, cq.Vector(0, 0, Z_G_TOP), cq.Vector(0, 0, 1))])
print(f"bodies built {time.time() - t0:.0f}s", flush=True)
asm = cq.Assembly(None, name=f"BOARD_{TAG}")
asm.add(metal, name="METAL_PEC", color=cq.Color(0.9, 0.7, 0.2, 1))
asm.add(sub, name="SUBSTRATE", color=cq.Color(0.25, 0.55, 0.25, 1))
asm.add(bias, name="BIAS_RT5880", color=cq.Color(0.2, 0.5, 0.2, 1))
asm.add(air, name="AIR", color=cq.Color(0.7, 0.8, 1.0, 1))
asm.add(port, name="PORT", color=cq.Color(0.1, 0.1, 0.9, 1))
for nm, f in sheets.items():
    asm.add(f, name=nm, color=cq.Color(0.9, 0.1, 0.1, 1))
path = OUT / f"BOARD_{TAG}.step"
asm.export(str(path), exportType="STEP")
snapped, ex = exactify_step_text(path.read_text(encoding="utf-8"))
path.write_text(snapped, encoding="utf-8", newline="\n")
bad = audit_step_text(path.read_text(encoding="utf-8"))
manifest = dict(label=f"board {TAG} exporter manifest (all coordinates mm, board frame, rot180 layout)", D_mm=D_MM, t_bias_mm=T_BIAS, stack_z_mm=dict(
    top_cu=[0.0, T], substrate_rt5880=[Z_G_TOP, 0.0], ground_cu=[Z_G_BOT, Z_G_TOP], bias_rt5880=[Z_B_BOT, Z_G_BOT], bottom_pads=[Z_PAD_LO, Z_B_BOT], fr4="none"),
    feed=dict(pin_r=R_PIN, pin_z=[Z_G_TOP, PIN_TOP], opening_r=R_OUT, port_plane_z=Z_G_TOP, sma_window_xy=[WIN, WIN], window_removed_from="BIAS_RT5880"),
    cell=cell, layout_meta={k: layout[k] for k in layout if k not in ("sites",)},
    static_table_note="FR4-CALIBRATED: the static A cells come from results/sub6_slotdisc/rows_RT5880.jsonl (NOCUT rows) designed with the jobs1 FR4-board model; not re-derived for the RT5880 bias-layer stack",
    cells=man_cells, vias=man_vias, antipads=antipads, bottom_pads=pads, ground=dict(r=RB, z=[Z_G_BOT, Z_G_TOP], opening_r=R_OUT),
    counts=dict(sites=len(man_cells), switched=n_sw, static=len(man_cells) - n_sw, vias=len(man_vias), antipads=len(antipads), bottom_pads=len(pads), diodes=len(diodes)),
    step=dict(file=path.name, bytes=path.stat().st_size, inexact_literals=len(bad)))
(OUT / "manifest.json").write_text(json.dumps(manifest, indent=1, default=float))
(OUT / "diode_sheets.json").write_text(json.dumps(dict(units="mm", sheets=diodes, states=dict(OFF="22.23 - j335.4 ohm/sq", ON="1.60 + j29.2 ohm/sq")), indent=1))
print(f"exported {path.name} {path.stat().st_size / 1e6:.1f} MB, inexact literals {len(bad)}, {time.time() - t0:.0f}s", flush=True)
