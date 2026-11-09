"""SUB6 FINAL BOARD, 13-pad + SMA: board STEP v2 (thick and thin), healed, one body per net, every body named.  MODEL_ONLY / NOT_HFSS_VALIDATED.
    SUB6_MESH=96 python sub6_final_step13_v2.py thick|thin
Writes hfss_exports/sub6_final_board_m96_13pad/BOARD_N240_clean_v2.step / BOARD_N240_thin_clean_v2.step, body_check_<variant>_v2.json, geometry_build_<variant>_v2.json, diode_sheets.csv (unchanged data)

Changes with respect to v1 (geometry is not changed beyond healing):
 * bias traces are built as ONE stroke face per polyline (OCC offset2D of the centre line, round caps and joins) instead of one rectangle per segment + a disc per end; collinear segments of a chain are merged first
   (no sub-micron slivers at the joints);
 * thick: every net is fused into ONE solid (cell copper + barrels + via pads + stroke solids + connector pad);
 * thin  : every net is ONE connected SHEET body (a single sewn shell): cell faces with a hole at every via, the via barrels as cylindrical faces (tubes) whose end circles are the hole edges, and one merged B.Cu face
           (stroke faces + via pads + connector pad, fused, with the via holes cut).  Barrels are sheets, so no solid is mixed with sheets in one body and no body is a compound of loose pieces.
           The dielectrics keep no holes for the barrel tubes (sheets inside a dielectric are legal); only the probe pin has its hole;
 * the ground (thin) is one sewn shell: ground face + SMA pad face + 12 tube faces; the SMA body block is its own solid body METAL_PEC_SMA_BODY;
 * every body goes through BRepCheck_Analyzer + ShapeFix_Shape + UnifySameDomain (merges collinear edges / coplanar faces) and small-edge / small-face (< 1 um) detection; per-body result in body_check_*.json;
 * the root assembly is named BOARD_N240 (v1 exported an unnamed root product with a generated UUID name, which AEDT showed as the unnamed body 'ffd89c01-..._246').
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

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
SRC = ROOT / "results" / "sub6_final_board" / "mesh96"
KDIR = ROOT / "results" / "sub6_final_board" / "mesh96_13pad" / "kicad"
OUT = ROOT / "hfss_exports" / "sub6_final_board_m96_13pad"
R3 = lambda v: round(float(v) + 0.0, 3)
H, T = 3.175, 0.035
R_PIN, R_OUT = 0.635, 2.05
R_BAR, R_PAD, R_ANTI = 0.15, 0.30, 0.50
ISL = (2.174, 0.65, 0.0, -0.20)
SH_ALONG, SH_ACROSS, SH_EXT = 0.35, 1.00, 0.0015
PAD_XY, PAD_TOP, PAD_BOT = 15.0, 40.0, 15.0
SMA_BODY, SMA_PAD, H_SMA, N_GV, R_GV = 6.35, 7.35, 6.0, 12, 2.7
SMALL = 1e-3   # 1 um in mm


def levels(thin):
    if thin:
        z = dict(top=0.0, gnd_lo=-H, gnd_hi=-H, fr4_hi=-H, fr4_lo=-H - 0.8, bot_lo=-H - 0.8, bot_hi=-H - 0.8, sheet=0.0)
    else:
        z = dict(top=T, gnd_lo=-H - T, gnd_hi=-H, fr4_hi=-H, fr4_lo=-H - T - 0.8, bot_lo=-H - T - 0.8 - T, bot_hi=-H - T - 0.8, sheet=T / 2)
    z["pin_top"] = round(-H + 12.1, 3)
    z["air_lo"] = round(z["bot_lo"] - PAD_BOT, 3)
    z["air_hi"] = round(T + PAD_TOP, 3)
    z["pin_lo"] = z["bot_lo"] - H_SMA
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


# ------------------------------------------------------------------ track chains (merge collinear segments, one polyline per chain)
def chains_from_tracks(tracks):
    by = {}
    for t in tracks:
        a = (R3(t["a"][0]), R3(t["a"][1]))
        b = (R3(t["b"][0]), R3(t["b"][1]))
        if a != b:
            by.setdefault(t["net"], []).append((a, b, round(t["w"], 4)))
    out = {}
    for net, segs in by.items():
        adj = {}
        for i, (a, b, w) in enumerate(segs):
            adj.setdefault(a, []).append(i)
            adj.setdefault(b, []).append(i)
        used = set()

        def passthru(n, i):
            lst = adj[n]
            return len(lst) == 2 and segs[lst[0]][2] == segs[lst[1]][2]

        def far(i, n):
            a, b, _ = segs[i]
            return b if a == n else a

        chains = []
        order = sorted(range(len(segs)), key=lambda i: (0 if (not passthru(segs[i][0], i) or not passthru(segs[i][1], i)) else 1, i))
        for i in order:
            if i in used:
                continue
            used.add(i)
            a, b, w = segs[i]
            pts = [a, b]
            for end in (1, 0):   # extend at the end, then at the start
                while True:
                    n = pts[-1] if end == 1 else pts[0]
                    if not passthru(n, i):
                        break
                    nxt = [j for j in adj[n] if j != i and j not in used]
                    if not nxt:
                        break
                    j = nxt[0]
                    used.add(j)
                    m = far(j, n)
                    if end == 1:
                        pts.append(m)
                    else:
                        pts.insert(0, m)
                    i = j
            # merge collinear interior points
            res = [pts[0]]
            for k in range(1, len(pts) - 1):
                p0, p1, p2 = res[-1], pts[k], pts[k + 1]
                cr = (p1[0] - p0[0]) * (p2[1] - p1[1]) - (p1[1] - p0[1]) * (p2[0] - p1[0])
                dt = (p1[0] - p0[0]) * (p2[0] - p1[0]) + (p1[1] - p0[1]) * (p2[1] - p1[1])
                if abs(cr) < 1e-9 and dt > 0:
                    continue
                res.append(p1)
            res.append(pts[-1])
            chains.append((res, w))
        out[net] = chains
    return out


FALLBACKS = []


def stroke_face(pts, w, zz, name=""):
    """one planar face for a polyline of half width w/2 with round caps/joins (OCC offset2D); fallback: fused rectangles + discs."""
    import cadquery as cq
    es = [cq.Edge.makeLine(cq.Vector(a[0], a[1], zz), cq.Vector(b[0], b[1], zz)) for a, b in zip(pts[:-1], pts[1:])]
    try:
        o = cq.Wire.assembleEdges(es).offset2D(w / 2, "arc")
        return cq.Face.makeFromWires(o[0])
    except Exception:
        FALLBACKS.append((name, len(pts)))
        pieces = []
        for a, b in zip(pts[:-1], pts[1:]):
            dx, dy = b[0] - a[0], b[1] - a[1]
            ln = math.hypot(dx, dy)
            ux, uy = dx / ln, dy / ln
            nx, ny = -uy * w / 2, ux * w / 2
            pieces.append(cq.Face.makeFromWires(cq.Wire.makePolygon([cq.Vector(a[0] + nx, a[1] + ny, zz), cq.Vector(b[0] + nx, b[1] + ny, zz), cq.Vector(b[0] - nx, b[1] - ny, zz), cq.Vector(a[0] - nx, a[1] - ny, zz)], close=True)))
        for q_ in pts:
            pieces.append(cq.Face.makeFromWires(cq.Wire.makeCircle(w / 2, cq.Vector(q_[0], q_[1], zz), cq.Vector(0, 0, 1))))
        return pieces[0].fuse(*pieces[1:]).clean()


# ------------------------------------------------------------------ healing and per-body check
def heal(shape, kind):
    """returns (healed shape, record).  kind: 'solid' or 'sheet'."""
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.ShapeFix import ShapeFix_Shape
    from OCP.ShapeUpgrade import ShapeUpgrade_UnifySameDomain
    from OCP.GProp import GProp_GProps
    from OCP.BRepGProp import BRepGProp
    import cadquery as cq

    def edges_small(s):
        n = 0
        for e in s.Edges():
            if e.Length() < SMALL:
                n += 1
        return n

    def faces_small(s):
        return sum(1 for f in s.Faces() if f.Area() < SMALL * SMALL)

    rec = dict(valid_before=bool(BRepCheck_Analyzer(shape.wrapped).IsValid()), small_edges_before=edges_small(shape), small_faces_before=faces_small(shape), actions=[])
    s = shape
    if not rec["valid_before"] or rec["small_edges_before"] or rec["small_faces_before"]:
        sf = ShapeFix_Shape(s.wrapped)
        sf.SetPrecision(1e-6)
        sf.SetMinTolerance(1e-7)
        sf.Perform()
        s = cq.Shape.cast(sf.Shape())
        rec["actions"].append("ShapeFix_Shape")
    u = ShapeUpgrade_UnifySameDomain(s.wrapped, True, True, True)
    u.Build()
    s2 = cq.Shape.cast(u.Shape())
    if len(s2.Edges()) < len(s.Edges()) or len(s2.Faces()) < len(s.Faces()):
        rec["actions"].append("UnifySameDomain(edges %d->%d, faces %d->%d)" % (len(s.Edges()), len(s2.Edges()), len(s.Faces()), len(s2.Faces())))
    if BRepCheck_Analyzer(s2.wrapped).IsValid() or not BRepCheck_Analyzer(s.wrapped).IsValid():
        s = s2
    rec["valid_after"] = bool(BRepCheck_Analyzer(s.wrapped).IsValid())
    rec["small_edges_after"] = edges_small(s)
    rec["small_faces_after"] = faces_small(s)
    rec["n_solids"] = len(s.Solids())
    rec["n_shells"] = len(s.Shells())
    rec["n_faces"] = len(s.Faces())
    rec["status"] = "VALID" if rec["valid_before"] and not rec["actions"] else ("FIXED" if rec["valid_after"] else "INVALID")
    return s, rec


def single_body(shape):
    """a body must be ONE lump: one solid (thick) or one shell (sheet)."""
    import cadquery as cq
    n_sol = len(shape.Solids())
    if n_sol >= 1:
        return n_sol == 1, "solid"
    if len(shape.Shells()) == 0:
        return len(shape.Faces()) == 1, "sheet"   # a single bare face
    return len(shape.Shells()) == 1, "sheet"


def build(thin, tracks_file=None):
    import cadquery as cq
    z = levels(thin)
    cells = list(csv.DictReader(open(SRC / "fab" / "cell_table.csv")))
    rt = json.load(open(SRC / "routing.json"))
    rout = rt["r_out"]
    hexpts = [(R3(rout * math.cos(math.radians(90 + 60 * i))), R3(rout * math.sin(math.radians(90 + 60 * i)))) for i in range(6)]
    tracks = json.load(open(tracks_file or (KDIR / "tracks.json")))
    nets_js = json.load(open(KDIR / "nets.json"))
    isl_net = lambda s, p: f"S{s:02d}" if p == 1 else f"N{s:02d}_{p - 1:02d}"
    dsc_net = lambda s, p: "RETURN" if p == 20 else f"N{s:02d}_{p:02d}"
    t0 = time.time()
    chains = chains_from_tracks(tracks)
    n_seg_in = len(tracks)
    n_chain = sum(len(v) for v in chains.values())
    top = {}      # net -> list of top faces (thin) / solids (thick)
    vias = {}     # net -> list of (x, y)
    static = {}
    barrel_solids = []   # solid rods for the thick dielectric holes
    via_xy = []

    def add(d, net, shp):
        d.setdefault(net, []).append(shp)

    def disc_slot_wire_face(x, y, th, R, Lh, Wh, zz):
        wp = cq.Workplane("XY").workplane(offset=zz).center(R3(x), R3(y)).circle(R).polyline([(p[0] - R3(x), p[1] - R3(y)) for p in quad(x, y, Lh, Wh, th)]).close()
        return wp

    def circle_face(cx, cy, r, zz):
        return cq.Face.makeFromWires(cq.Wire.makeCircle(r, cq.Vector(cx, cy, zz), cq.Vector(0, 0, 1)))

    def poly_face(pts, zz):
        return cq.Face.makeFromWires(cq.Wire.makePolygon([cq.Vector(a, b, zz) for a, b in pts], close=True))

    z0_top = z["top"] - (0.0 if thin else T)
    for c in cells:
        x, y, th = float(c["x_mm"]), float(c["y_mm"]), math.radians(float(c["slot_axis_deg"]))
        R, Lh, Wh = float(c["disc_r_mm"]), float(c["slot_L_mm"]) / 2, float(c["slot_W_mm"]) / 2
        base = disc_slot_wire_face(x, y, th, R, Lh, Wh, z0_top)
        sol = base.extrude(T).val()
        if thin:
            disc = [f for f in sol.Faces() if abs(f.Center().z - z0_top) < 1e-9 and f.normalAt().z < -0.5][0]
        else:
            disc = sol
        if c["class"] != "SWITCHED":
            static[c["site"]] = disc
            continue
        s, p = int(c["string"]), int(c["pos_in_string"])
        vd = (R3(float(c["disc_via_x"])), R3(float(c["disc_via_y"])))
        vi = (R3(float(c["island_via_x"])), R3(float(c["island_via_y"])))
        ipts = quad(x, y, ISL[0], ISL[1], th, ISL[2], ISL[3])
        if thin:
            disc = disc.cut(circle_face(vd[0], vd[1], R_BAR, z0_top))
            isl = poly_face(ipts, z0_top).cut(circle_face(vi[0], vi[1], R_BAR, z0_top))
        else:
            isl = cq.Workplane("XY").workplane(offset=z0_top).polyline(ipts).close().extrude(T).val()
        add(top, dsc_net(s, p), disc)
        add(top, isl_net(s, p), isl)
        for (vx, vy), nm in ((vd, dsc_net(s, p)), (vi, isl_net(s, p))):
            via_xy.append((vx, vy))
            add(vias, nm, (vx, vy))
            rod = cq.Workplane("XY").workplane(offset=z["bot_hi"]).center(vx, vy).circle(R_BAR).extrude(0.0 - z["bot_hi"]).val()
            barrel_solids.append(rod)
    print(f"cells done {time.time() - t0:.0f}s nets {len(top)} static {len(static)} chains {n_chain} (from {n_seg_in} segments)", flush=True)
    # ---------------- per-net bodies
    metal = {}
    conn_pads = {pd["net"]: (R3(pd["x"]), R3(pd["y"])) for pd in nets_js["pads"]}
    for net in sorted(top):
        vlist = vias.get(net, [])
        zb0 = z["bot_lo"]
        if thin:
            faces = []
            for pts, w in chains.get(net, []):
                faces.append(stroke_face(pts, w, zb0, net))
            for vx, vy in vlist:
                faces.append(circle_face(vx, vy, R_PAD, zb0))
            if net in conn_pads:
                cx_, cy_ = conn_pads[net]
                faces.append(poly_face([(cx_ - .5, cy_ - .5), (cx_ + .5, cy_ - .5), (cx_ + .5, cy_ + .5), (cx_ - .5, cy_ + .5)], zb0))
            bf = faces[0].fuse(*faces[1:]).clean() if len(faces) > 1 else faces[0]
            bf = bf.cut(*[circle_face(vx, vy, R_BAR, zb0) for vx, vy in vlist])
            tubes = []
            for vx, vy in vlist:
                cyl = cq.Solid.makeCylinder(R_BAR, 0.0 - zb0, cq.Vector(vx, vy, zb0), cq.Vector(0, 0, 1))
                tubes.append([f for f in cyl.Faces() if f.geomType() == "CYLINDER"][0])
            allf = []
            for sh in top[net] + [bf]:
                allf += list(sh.Faces())
            allf += tubes
            metal["METAL_PEC_" + net] = cq.Shell.makeShell(allf)
        else:
            sols = list(top[net])
            for pts, w in chains.get(net, []):
                f = stroke_face(pts, w, zb0, net)
                sols.append(cq.Solid.extrudeLinear(f, cq.Vector(0, 0, z["bot_hi"] - zb0)))
            for vx, vy in vlist:
                sols.append(cq.Workplane("XY").workplane(offset=zb0).center(vx, vy).circle(R_PAD).extrude(z["bot_hi"] - zb0).val())
                sols.append(cq.Workplane("XY").workplane(offset=z["bot_hi"]).center(vx, vy).circle(R_BAR).extrude(0.0 - z["bot_hi"]).val())
            if net in conn_pads:
                cx_, cy_ = conn_pads[net]
                sols.append(cq.Workplane("XY").workplane(offset=zb0).center(cx_, cy_).rect(1.0, 1.0).extrude(z["bot_hi"] - zb0).val())
            metal["METAL_PEC_" + net] = sols[0].fuse(*sols[1:]).clean()
    for site, shp in static.items():
        metal["METAL_PEC_C" + site] = shp
    print(f"nets done {time.time() - t0:.0f}s  metal bodies {len(metal)}", flush=True)
    # ---------------- ground (+ SMA pad, ring vias, body)
    gv = [(R3(R_GV * math.cos(2 * math.pi * k / N_GV + math.pi / N_GV)), R3(R_GV * math.sin(2 * math.pi * k / N_GV + math.pi / N_GV))) for k in range(N_GV)]
    gw = cq.Workplane("XY").workplane(offset=z["gnd_lo"]).polyline(hexpts).close()
    for vx, vy in via_xy:
        gw = gw.moveTo(vx, vy).circle(R_ANTI)
    gw = gw.moveTo(0, 0).circle(R_OUT)
    body = cq.Workplane("XY").workplane(offset=z["bot_lo"] - H_SMA).rect(SMA_BODY, SMA_BODY).circle(R_OUT).extrude(H_SMA).val()
    sma_rods = [cq.Workplane("XY").workplane(offset=z["bot_lo"]).center(vx, vy).circle(R_BAR).extrude(z["gnd_lo"] - z["bot_lo"]).val() for vx, vy in gv]
    if thin:
        gwx = gw
        for vx, vy in gv:
            gwx = gwx.moveTo(vx, vy).circle(R_BAR)
        gnd_face = [f for f in gwx.extrude(T).val().Faces() if abs(f.Center().z - z["gnd_lo"]) < 1e-9 and f.normalAt().z < -0.5][0]
        padw = cq.Workplane("XY").workplane(offset=z["bot_lo"]).rect(SMA_PAD, SMA_PAD).circle(R_OUT)
        for vx, vy in gv:
            padw = padw.moveTo(vx, vy).circle(R_BAR)
        pad_face = [f for f in padw.extrude(T).val().Faces() if abs(f.Center().z - z["bot_lo"]) < 1e-9 and f.normalAt().z < -0.5][0]
        tubes = []
        for vx, vy in gv:
            cyl = cq.Solid.makeCylinder(R_BAR, z["gnd_lo"] - z["bot_lo"], cq.Vector(vx, vy, z["bot_lo"]), cq.Vector(0, 0, 1))
            tubes.append([f for f in cyl.Faces() if f.geomType() == "CYLINDER"][0])
        metal["METAL_PEC_GROUND"] = cq.Shell.makeShell([gnd_face, pad_face] + tubes)
        metal["METAL_PEC_SMA_BODY"] = body
        ground_solid = None
    else:
        ground_solid = gw.extrude(z["gnd_hi"] - z["gnd_lo"]).val()
        pad_s = cq.Workplane("XY").workplane(offset=z["bot_lo"]).rect(SMA_PAD, SMA_PAD).circle(R_OUT).extrude(z["bot_hi"] - z["bot_lo"]).val()
        ground_solid = ground_solid.fuse(pad_s, body, *sma_rods).clean()
        metal["METAL_PEC_GROUND"] = ground_solid
    pin = cq.Workplane("XY").workplane(offset=z["pin_lo"]).circle(R_PIN).extrude(z["pin_top"] - z["pin_lo"]).val()
    metal["METAL_PEC_PIN"] = pin
    # ---------------- healing + per-body checks
    checks = {}
    healed = {}
    for nm, shp in metal.items():
        h, rec = heal(shp, "sheet" if thin else "solid")
        ok1, kind = single_body(h)
        rec["single_lump"] = bool(ok1)
        rec["type"] = kind
        if not ok1:
            rec["status"] = "SPLIT"
        healed[nm] = h
        checks[nm] = rec
    print(f"healing done {time.time() - t0:.0f}s", flush=True)
    # ---------------- dielectrics
    sub = cq.Workplane("XY").workplane(offset=-H).polyline(hexpts).close().extrude(H).val()
    fr4 = cq.Workplane("XY").workplane(offset=z["fr4_lo"]).polyline(hexpts).close().extrude(z["fr4_hi"] - z["fr4_lo"]).val()
    if thin:
        sub = sub.cut(pin).clean()
        fr4 = fr4.cut(pin).clean()
    else:
        sub = sub.cut(*(barrel_solids + [pin])).clean()
        fr4 = fr4.cut(*(barrel_solids + [ground_solid, pin])).clean()
    ax, ay = rout * math.cos(math.radians(30)) + PAD_XY, rout + PAD_XY
    air = cq.Workplane("XY").workplane(offset=z["air_lo"]).rect(2 * R3(ax), 2 * R3(ay)).extrude(z["air_hi"] - z["air_lo"]).val()
    tools = [sub, fr4]
    for nm, shp in healed.items():
        if not thin or nm in ("METAL_PEC_PIN", "METAL_PEC_SMA_BODY"):
            tools.append(shp)
    air = air.cut(*tools).clean()
    print(f"air done {time.time() - t0:.0f}s", flush=True)
    parts = dict(healed)
    for nm, shp in (("SUBSTRATE", sub), ("SUBSTRATE_FR4", fr4), ("AIR", air)):
        h, rec = heal(shp, "solid")
        rec["single_lump"] = bool(len(h.Solids()) == 1)
        rec["type"] = "solid"
        if not rec["single_lump"]:
            rec["status"] = "SPLIT"
        parts[nm] = h
        checks[nm] = rec
    # ---------------- diode sheets, port
    rows = []
    for c in cells:
        if c["class"] != "SWITCHED":
            continue
        s, p = int(c["string"]), int(c["pos_in_string"])
        x, y, th = float(c["x_mm"]), float(c["y_mm"]), math.radians(float(c["slot_axis_deg"]))
        q = []
        for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
            v = rot((sx * SH_ALONG / 2, 0.45 - SH_EXT + (sy + 1) / 2 * (SH_ACROSS + 2 * SH_EXT)), th)
            q.append(cq.Vector(R3(x + v[0]), R3(y + v[1]), z["sheet"]))
        parts[f"D{s}_{p}"] = cq.Face.makeFromWires(cq.Wire.makePolygon(q, close=True))
        d = rot((0.0, 1.0), th)
        a_pt = (R3(x + rot((0.0, 0.45), th)[0]), R3(y + rot((0.0, 0.45), th)[1]))
        k_pt = (R3(x + rot((0.0, 1.45), th)[0]), R3(y + rot((0.0, 1.45), th)[1]))
        cen = (R3(x + rot((0.0, 0.95), th)[0]), R3(y + rot((0.0, 0.95), th)[1]))
        rows.append([f"D{s}_{p}", cen[0], cen[1], z["sheet"], round(d[0], 6), round(d[1], 6), a_pt[0], a_pt[1], k_pt[0], k_pt[1], s, p, c["site"], round(math.degrees(th) % 360, 3)])
        checks[f"D{s}_{p}"] = dict(type="sheet", n_faces=1, status="VALID", single_lump=True)
    zp = z["fr4_lo"]
    parts["PORT_ANNULUS"] = cq.Face.makeFromWires(cq.Wire.makeCircle(R_OUT, cq.Vector(0, 0, zp), cq.Vector(0, 0, 1)), [cq.Wire.makeCircle(R_PIN, cq.Vector(0, 0, zp), cq.Vector(0, 0, 1))])
    checks["PORT_ANNULUS"] = dict(type="sheet", n_faces=1, status="VALID", single_lump=True)
    meta = dict(stroke_fallbacks=FALLBACKS, z=z, thin=thin, n_metal=len(metal), n_sheets=len(rows), air_box=dict(ax=R3(ax), ay=R3(ay)), n_segments_in=n_seg_in, n_chains=n_chain)
    return parts, rows, meta, checks


def main():
    variant = sys.argv[1] if len(sys.argv) > 1 else "thick"
    thin = variant == "thin"
    import cadquery as cq
    from sub6_step_exact import audit_step_text, exactify_step_text
    from sidequest_board_k2_lp_package import read_names
    from OCP.BRepCheck import BRepCheck_Analyzer
    OUT.mkdir(parents=True, exist_ok=True)
    parts, rows, meta, checks = build(thin, sys.argv[2] if len(sys.argv) > 2 else None)
    with open(OUT / "diode_sheets.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["name", "sheet_centre_x", "sheet_centre_y", "sheet_z_thick", "current_dir_x", "current_dir_y", "anode_edge_x", "anode_edge_y", "cathode_edge_x", "cathode_edge_y", "string", "pos_in_string", "site", "slot_axis_deg", "sheet_z_thin"])
        for r in rows:
            w.writerow(r[:3] + [0.0175] + r[4:] + [0.0])
    asm = cq.Assembly(None, name="BOARD_N240")
    for nm, shp in parts.items():
        col = (0.9, 0.7, 0.2) if nm.startswith("METAL") else {"SUBSTRATE": (0.25, 0.55, 0.25), "SUBSTRATE_FR4": (0.2, 0.5, 0.2), "AIR": (0.7, 0.8, 1.0), "PORT_ANNULUS": (0.1, 0.1, 0.9)}.get(nm, (0.9, 0.1, 0.1))
        asm.add(shp, name=nm, color=cq.Color(*col, 1))
    path = OUT / ("BOARD_N240_thin_clean_v2.step" if thin else "BOARD_N240_clean_v2.step")
    t0 = time.time()
    asm.export(str(path), exportType="STEP")
    raw = path.read_text(encoding="utf-8")
    snapped, ex = exactify_step_text(raw)
    path.write_text(snapped, encoding="utf-8", newline="\n")
    bad = audit_step_text(path.read_text(encoding="utf-8"))
    names = read_names(path)
    # re-read the written STEP (after the 0.001 mm literal snapping) and check every body again
    after = {}
    for nm, shp in names.items():
        sh = shp.wrapped if hasattr(shp, "wrapped") else shp
        after[nm] = bool(BRepCheck_Analyzer(sh).IsValid())
    for nm in checks:
        checks[nm]["valid_in_exported_step"] = after.get(nm)
    cnt = {}
    for nm, r in checks.items():
        key = r["status"] if after.get(nm, False) else (r["status"] + "+EXPORT_INVALID")
        cnt[key] = cnt.get(key, 0) + 1
    unnamed = [n for n in names if not n or n.startswith("ffd8") or len(n) > 40]
    chk = dict(label="MODEL_ONLY geometry build v2", variant=variant, n_bodies=len(names), n_expected=len(parts), unnamed_or_uuid_names=unnamed, duplicate_names=len(parts) - len(set(parts)), status_counts=cnt,
               split_bodies=[n for n, r in checks.items() if not r.get("single_lump", True)], invalid_after_export=[n for n in checks if not after.get(n, False)],
               fixed_bodies=[n for n, r in checks.items() if r.get("status") == "FIXED"], step_size_MB=round(path.stat().st_size / 1e6, 1), audit_inexact_literals=len(bad), meta=meta, export_s=round(time.time() - t0))
    (OUT / f"body_check_{variant}_v2.json").write_text(json.dumps(dict(summary=chk, bodies=checks), indent=1, default=str), encoding="utf-8", newline="\n")
    (OUT / f"geometry_build_{variant}_v2.json").write_text(json.dumps(chk, indent=1, default=str), encoding="utf-8", newline="\n")
    print(json.dumps({k: v for k, v in chk.items() if k != "meta"}, default=str)[:1500], flush=True)


if __name__ == "__main__":
    main()
