"""SUB6 FINAL BOARD: Gerber RS-274X (+ Excellon drill) from the layout of sub6_final_export.py (fab/cell_table.csv, traces.json, layout_summary.json).  MODEL_ONLY / NOT_HFSS_VALIDATED.
    SUB6_MESH=96 python sub6_final_gerber.py
Layers (mm, FSLAX46Y46, 1e-6 mm resolution, all coordinates on the 0.001 mm grid):
  L1 top copper     : discs (dark), slots (clear), islands (dark) of all 758 cells.  The 1.00 mm diode gap is left open (no diode footprint is defined; see README).
  L2 ground (mid)   : hexagon (dark) with 1.0 mm anti-pads at all 480 vias and the 4.1 mm probe opening (clear)
  L3 bottom copper  : bias traces 0.15 mm, via pads 0.6 mm, 24 connector pads 1.0 x 1.0, probe pad
  Outline           : hexagon, 0.1 mm line
  Drill (Excellon)  : T1 0.30 mm plated vias (480), T2 1.27 mm probe hole
"""
import csv
import json
import math
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sub6_final_layout as L

FAB = L.OUT / "fab"
GB = FAB / "gerber"
Q = lambda v: int(round(float(v) * 1e6))


class G:
    def __init__(self, fn, function, polarity="Positive"):
        self.lines = ["%TF.GenerationSoftware,sub6_final_gerber.py,MODEL_ONLY*%", f"%TF.FileFunction,{function}*%", f"%TF.FilePolarity,{polarity}*%", "%FSLAX46Y46*%", "%MOMM*%", "G75*"]
        self.ap = {}
        self.fn = fn
        self.cur = None
        self.body = []

    def aperture(self, kind, *dims):
        key = (kind,) + tuple(round(d, 4) for d in dims)
        if key not in self.ap:
            code = 10 + len(self.ap)
            self.ap[key] = code
            self.lines.append(f"%ADD{code}{kind},{'X'.join(f'{d:.4f}' for d in dims)}*%")
        return self.ap[key]

    def select(self, code):
        if self.cur != code:
            self.body.append(f"D{code}*")
            self.cur = code

    def flash(self, kind, dims, x, y):
        self.select(self.aperture(kind, *dims))
        self.body.append(f"X{Q(x)}Y{Q(y)}D03*")

    def draw(self, d, pts, close=False):
        self.select(self.aperture("C", d))
        p = list(pts) + ([pts[0]] if close else [])
        self.body.append(f"X{Q(p[0][0])}Y{Q(p[0][1])}D02*")
        for q in p[1:]:
            self.body.append(f"X{Q(q[0])}Y{Q(q[1])}D01*")

    def region(self, pts):
        p = list(pts) + [pts[0]]
        self.body.append("G36*")
        self.body.append(f"X{Q(p[0][0])}Y{Q(p[0][1])}D02*")
        for q in p[1:]:
            self.body.append(f"G01X{Q(q[0])}Y{Q(q[1])}D01*")
        self.body.append("G37*")

    def polarity(self, dark):
        self.body.append("%LPD*%" if dark else "%LPC*%")

    def save(self):
        GB.mkdir(parents=True, exist_ok=True)
        txt = "\n".join(self.lines + ["%LPD*%"] + self.body + ["M02*"]) + "\n"
        (GB / self.fn).write_text(txt, encoding="ascii", newline="\n")
        return dict(file=self.fn, n_lines=len(self.body), n_apertures=len(self.ap))


def rect(x, y, hx, hy, th, ox=0.0, oy=0.0):
    out = []
    for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
        v = L.rot(np.array([ox + sx * hx, oy + sy * hy]), th)
        out.append((round(x + v[0], 3), round(y + v[1], 3)))
    return out


def main():
    cells = list(csv.DictReader(open(FAB / "cell_table.csv")))
    traces = json.load(open(FAB / "traces.json"))
    meta = json.load(open(FAB / "layout_summary.json"))
    rout = meta["outline_hex_circumradius_mm"]
    hexp = [(round(rout * math.cos(math.radians(90 + 60 * i)), 3), round(rout * math.sin(math.radians(90 + 60 * i)), 3)) for i in range(6)]
    vias = []
    # ---- L1 top copper
    top = G("sub6_final-F_Cu.gtl", "Copper,L1,Top")
    for c in cells:
        top.flash("C", (2 * float(c["disc_r_mm"]),), float(c["x_mm"]), float(c["y_mm"]))
    top.polarity(False)
    for c in cells:
        th = math.radians(float(c["slot_axis_deg"]))
        top.region(rect(float(c["x_mm"]), float(c["y_mm"]), float(c["slot_L_mm"]) / 2, float(c["slot_W_mm"]) / 2, th))
    top.polarity(True)
    for c in cells:
        if c["class"] == "SWITCHED":
            th = math.radians(float(c["slot_axis_deg"]))
            top.region(rect(float(c["x_mm"]), float(c["y_mm"]), 2.174, 0.65, th, 0.0, -0.20))
            vias += [(float(c["island_via_x"]), float(c["island_via_y"])), (float(c["disc_via_x"]), float(c["disc_via_y"]))]
    r1 = top.save()
    # ---- L2 ground
    gnd = G("sub6_final-In1_Cu.g2", "Copper,L2,Inr")
    gnd.region(hexp)
    gnd.polarity(False)
    for vx, vy in vias:
        gnd.flash("C", (1.0,), vx, vy)
    gnd.flash("C", (4.1,), 0.0, 0.0)
    r2 = gnd.save()
    if os.environ.get("GERBER_RF_ONLY") != "1":
        # ---- L3 bottom copper
        bot = G("sub6_final-B_Cu.gbl", "Copper,L3,Bot")
        for t in traces:
            bot.draw(L.W_TR, [tuple(p) for p in t["poly"]])
        for vx, vy in vias:
            bot.flash("C", (0.6,), vx, vy)
        for k in range(24):
            bot.flash("R", (1.0, 1.0), meta["connector_x_mm"], round((11.5 - k) * L.PAD_PITCH, 3))
        bot.flash("C", (2.0,), 0.0, 0.0)
        r3 = bot.save()
    # ---- outline
    out = G("sub6_final-Edge_Cuts.gm1", "Profile,NP")
    out.draw(0.1, hexp, close=True)
    r4 = out.save()
    # ---- Excellon
    ex = ["M48", "; sub6_final_gerber.py, MODEL_ONLY", ";TF.FileFunction,Plated,1,3,PTH", "METRIC,TZ", "FMAT,2", "T1C0.300", "T2C1.270", "%", "T1"]
    ex += [f"X{vx:.3f}Y{vy:.3f}" for vx, vy in vias]
    ex += ["T2", "X0.000Y0.000", "M30"]
    (GB / "sub6_final-PTH.drl").write_text("\n".join(ex) + "\n", encoding="ascii", newline="\n")
    if os.environ.get("GERBER_RF_ONLY") == "1":
        print("RF-only Gerbers written (F_Cu, In1_Cu, Edge_Cuts, PTH)")
        return
    # ---- self-check: re-parse flash/region counts
    chk = {}
    for fn in ("sub6_final-F_Cu.gtl", "sub6_final-In1_Cu.g2", "sub6_final-B_Cu.gbl", "sub6_final-Edge_Cuts.gm1"):
        t = (GB / fn).read_text()
        chk[fn] = dict(flashes=t.count("D03*"), regions=t.count("G36*"), regions_end=t.count("G37*"), draws=t.count("D01*"), ends_with_M02=t.rstrip().endswith("M02*"), lpc=t.count("%LPC*%"))
    chk["n_cells"], chk["n_vias"] = len(cells), len(vias)
    chk["ok"] = bool(chk["sub6_final-F_Cu.gtl"]["flashes"] == len(cells) and chk["sub6_final-F_Cu.gtl"]["regions"] == len(cells) + len(vias) // 2 and chk["sub6_final-In1_Cu.g2"]["flashes"] == len(vias) + 1)
    pts = [(float(c["x_mm"]), float(c["y_mm"]), float(c["disc_r_mm"])) for c in cells]
    chk["min_cell_centre_radius_mm"] = min(math.hypot(x, y) for x, y, _ in pts)
    (GB / "gerber_selfcheck.json").write_text(json.dumps(chk, indent=1), encoding="utf-8", newline="\n")
    print(json.dumps(chk, indent=1))
    readme = f"""# Fabrication README: SUB6 final board (hexagon, N 240, mesh-96 switched-cell design)
Status: MODEL_ONLY / NOT_HFSS_VALIDATED; design status in ../layout_summary.json (FINAL_CANDIDATE).  Files are script-generated (scripts/research/sub6_final_gerber.py); Gerber RS-274X, format 4.6 mm, Excellon METRIC.

## Stack-up (top to bottom)
| layer | material | thickness |
|---|---|---|
| L1 top copper | 35 um Cu (cells: discs, slots, islands) | 0.035 |
| core | Rogers RT5880 (eps_r 2.2, tan d 0.0009) | 3.175 |
| L2 ground | 35 um Cu, anti-pads 1.0 mm at every via, 4.1 mm opening at the probe | 0.035 |
| dielectric | FR4 (eps_r 4.4, tan d 0.02 ASSUMED) | 0.8 |
| L3 bottom copper | 35 um Cu (bias traces, via pads, connector pads) | 0.035 |
Total about 4.08 mm.  Board outline: regular hexagon, circumradius {rout:.3f} mm (flat-to-flat {2 * rout * math.cos(math.radians(30)):.3f} mm), vertices up (+y), origin = probe centre.

## Files
sub6_final-F_Cu.gtl (L1), sub6_final-In1_Cu.g2 (L2 ground), sub6_final-B_Cu.gbl (L3), sub6_final-Edge_Cuts.gm1 (outline), sub6_final-PTH.drl (Excellon); ../cell_table.csv, ../string_table.csv, ../drill.csv, ../sub6_final_board.dxf.
Polarity: dark = copper; the top layer uses clear polarity for the slots, so apply the file in order.

## Drill and rules
Vias: {len(vias)} plated holes, drill 0.30 mm, pad 0.6 mm (top land = the disc or island itself), anti-pad 1.0 mm in L2.  Probe hole 1.27 mm (pin r 0.635) at the origin, 4.1 mm opening in L2, 2.0 mm pad on L3.
Min trace 0.15 mm, min copper-to-copper gap on L3 0.274 mm (all trace/pad/anti-pad/outline clearances >= 0.15 mm, see ../layout_summary.json).  Connector: 24 pads 1.0 x 1.0 mm, pitch 2.54 mm (one column near the left flat edge, 3.5 mm from the pad centres to the outline).

## Points to check with the manufacturer (not resolved here)
1. Via aspect ratio: 0.3 mm drill through about 4.08 mm = 13.6 : 1, above the usual 8-10 : 1 plating limit; mixed RT5880 / FR4 stack with a bonded ground layer.  Ask for feasibility or a larger drill.
2. The BAR50 diode footprint is NOT defined: the Gerbers leave the 1.00 mm gap between island and rim open (the HFSS model uses a 1.00 x 0.35 mm RLC sheet).  Pads for the real diode package must be added; diode orientation: anode on the island via, cathode on the disc via of the same cell (series string island(i) - diode - disc(i) - trace - island(i+1)).
3. No solder mask, silkscreen or paste layers are generated; no copper thickness or finish is specified beyond 35 um.
4. FR4 properties are assumed typical values.  The bias strings (12 x 20 diodes, 19.0 V at 10 mA each) end on the 24 connector pads (see ../string_table.csv for the pad of each string start/end).
"""
    (GB / "README_FAB.md").write_text(readme, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
