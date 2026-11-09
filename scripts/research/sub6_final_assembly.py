"""SUB6 FINAL BOARD: diode assembly data (SC-79 pads, paste layer, assembly drawing DXF/PDF, pick-and-place CSV).  MODEL_ONLY / NOT_HFSS_VALIDATED.
    SUB6_MESH=96 python sub6_final_assembly.py
Reads fab/{cell_table.csv,string_table.csv,traces.json,layout_summary.json}; writes fab/assembly/:
  sub6_final-F_Paste.gtp   paste layer, 480 SC-79 pads 0.35 x 0.35 mm (Gerber RS-274X)
  sub6_final-F_Assy.gto    assembly layer: package body outline + cathode bar (Gerber, not for fabrication)
  diode_placement.csv      240 rows: ref, string, position, centre, rotation, anode/cathode pad centres
  sub6_final_assembly.dxf  assembly drawing, layers ASM_*
  sub6_final_assembly.pdf  assembly drawing (overview, string order, per-string pages, detail); NOT committed (PDFs are git-ignored by rule)
  assembly_check.json      pad-on-copper, pitch, gap, polarity-vs-string checks
  README_ASSEMBLY.md
Convention: the diode of cell k sits in the 1.00 mm gap between the island (anode side) and the disc rim (cathode side).
  Anode pad on the ISLAND (the island via is the same copper), cathode pad on the DISC (the disc via is the same copper).
  Current (conventional) flows along the string: connector pad_start -> island(1) -> diode(1) -> disc(1) -> bottom trace -> island(2) -> ... -> disc(20) -> connector pad_end.
  Rotation = angle of the anode->cathode vector (the cell-local +y axis) from board +x, counter-clockwise, degrees.
No solder-mask layer is generated for the top side (pads sit on bare, unmasked copper); no top mask opening exists in any file.
"""
import csv
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sub6_final_layout as L

FAB = L.OUT / "fab"
AS = FAB / "assembly"
Q = lambda v: int(round(float(v) * 1e6))
PAD = 0.35          # SC-79 land, square
PITCH = 1.35        # land centre pitch across the 1.00 mm gap
GAP_C = 0.95        # sheet / gap centre, cell-local y
BODY_L, BODY_W = 1.25, 0.85   # nominal SC-79 body (assumed typical, not from a datasheet here)
ISL = (2.174, 0.65, 0.0, -0.20)  # island half-length, half-width, local offset x, y


def lrot(v, th):
    return L.rot(np.asarray(v, float), th)


def rect_pts(cx, cy, hx, hy, th, ox=0.0, oy=0.0):
    return [(round(cx + lrot([ox + sx * hx, oy + sy * hy], th)[0], 3), round(cy + lrot([ox + sx * hx, oy + sy * hy], th)[1], 3)) for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1))]


def gerber(fn, function, regions, extra=None):
    out = ["%TF.GenerationSoftware,sub6_final_assembly.py,MODEL_ONLY*%", f"%TF.FileFunction,{function}*%", "%TF.FilePolarity,Positive*%", "%FSLAX46Y46*%", "%MOMM*%", "G75*", "%ADD10C,0.0500*%", "%LPD*%"]
    for pts in regions:
        p = list(pts) + [pts[0]]
        out += ["G36*", f"X{Q(p[0][0])}Y{Q(p[0][1])}D02*"] + [f"G01X{Q(q[0])}Y{Q(q[1])}D01*" for q in p[1:]] + ["G37*"]
    if extra:
        out += ["D10*"]
        for pts in extra:
            p = list(pts)
            out.append(f"X{Q(p[0][0])}Y{Q(p[0][1])}D02*")
            out += [f"X{Q(q[0])}Y{Q(q[1])}D01*" for q in p[1:]]
    out.append("M02*")
    (AS / fn).write_text("\n".join(out) + "\n", encoding="ascii", newline="\n")
    return len(regions)


def in_copper(pt, c):
    """pt in board frame -> inside (disc minus slot) or inside island of cell c (tolerance 2e-3 mm)."""
    th = math.radians(float(c["slot_axis_deg"]))
    d = np.array([pt[0] - float(c["x_mm"]), pt[1] - float(c["y_mm"])])
    lx, ly = lrot(d, -th)
    e = 2e-3   # pad corners are rounded to the 0.001 mm grid and the pad edges are flush with the copper edges
    isl = abs(lx - ISL[2]) <= ISL[0] + e and abs(ly - ISL[3]) <= ISL[1] + e
    in_slot = abs(lx) < float(c["slot_L_mm"]) / 2 - e and abs(ly) < float(c["slot_W_mm"]) / 2 - e
    disc = math.hypot(lx, ly) <= float(c["disc_r_mm"]) + e and not in_slot
    return isl, disc


def main():
    AS.mkdir(parents=True, exist_ok=True)
    cells = [c for c in csv.DictReader(open(FAB / "cell_table.csv")) if c["class"] == "SWITCHED"]
    strings = {int(r["string"]): r for r in csv.DictReader(open(FAB / "string_table.csv"))}
    traces = json.load(open(FAB / "traces.json"))
    meta = json.load(open(FAB / "layout_summary.json"))
    rout = meta["outline_hex_circumradius_mm"]
    hexp = [(round(rout * math.cos(math.radians(90 + 60 * i)), 3), round(rout * math.sin(math.radians(90 + 60 * i)), 3)) for i in range(6)]
    rows, pads, bodies, bars, chk_bad = [], [], [], [], []
    for c in cells:
        x, y, ang = float(c["x_mm"]), float(c["y_mm"]), float(c["slot_axis_deg"])
        th = math.radians(ang)
        a = np.array([x, y]) + lrot([0.0, GAP_C - PITCH / 2], th)
        k = np.array([x, y]) + lrot([0.0, GAP_C + PITCH / 2], th)
        cen = np.array([x, y]) + lrot([0.0, GAP_C], th)
        rot = (ang + 90.0) % 360.0
        s, p = int(c["string"]), int(c["pos_in_string"])
        ref = f"D{s}-{p:02d}"
        pa = rect_pts(a[0], a[1], PAD / 2, PAD / 2, th)
        pk = rect_pts(k[0], k[1], PAD / 2, PAD / 2, th)
        pads += [pa, pk]
        bodies.append(rect_pts(cen[0], cen[1], BODY_W / 2, BODY_L / 2, th))   # body long axis along the anode->cathode (local y)
        bars.append(rect_pts(cen[0], cen[1], BODY_W / 2, 0.06, th, 0.0, BODY_L / 2 - 0.12))   # cathode bar at the cathode end
        ia = [in_copper(q, c) for q in pa]
        ik = [in_copper(q, c) for q in pk]
        if not (all(v[0] for v in ia) and all(v[1] for v in ik)):
            chk_bad.append(ref)
        rows.append(dict(ref=ref, site=c["site"], string=s, pos=p, x_mm=round(cen[0], 3), y_mm=round(cen[1], 3), rotation_deg=round(rot, 3), anode_x=round(a[0], 3), anode_y=round(a[1], 3), cathode_x=round(k[0], 3), cathode_y=round(k[1], 3),
                         anode_on="island", cathode_on="disc", cell_x=x, cell_y=y, slot_axis_deg=ang))
    rows.sort(key=lambda r: (r["string"], r["pos"]))
    with open(AS / "diode_placement.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    n_paste = gerber("sub6_final-F_Paste.gtp", "Paste,Top", pads)
    gerber("sub6_final-F_Assy.gto", "Other,Assembly,Top", [], [b + [b[0]] for b in bodies] + [b for b in bars])
    # polarity vs string order: every link trace must go from the DISC via of pos i to the ISLAND via of pos i+1
    byk = {(int(c["string"]), int(c["pos_in_string"])): c for c in cells}
    bad_link = 0
    d = lambda p, q: math.hypot(p[0] - q[0], p[1] - q[1])
    for t in traces:
        if t["kind"] != "link":
            continue
        s, i = int(t["string"]), int(t["pos"]) + 1   # link pos 0 = cell pos 1 -> pos 2
        c0, c1 = byk[(s, i)], byk[(s, i + 1)]
        if d(t["poly"][0], (float(c0["disc_via_x"]), float(c0["disc_via_y"]))) > 0.01 or d(t["poly"][-1], (float(c1["island_via_x"]), float(c1["island_via_y"]))) > 0.01:
            bad_link += 1
    lead_ok = True
    for s in strings:
        st = [t for t in traces if t["kind"] == "lead" and int(t["string"]) == s]
        for t in st:
            c = byk[(s, 1)] if t["pos"] == "start" else byk[(s, 20)]
            via = (float(c["island_via_x"]), float(c["island_via_y"])) if t["pos"] == "start" else (float(c["disc_via_x"]), float(c["disc_via_y"]))
            if d(t["poly"][0], via) > 0.01:
                lead_ok = False
    chk = dict(label="MODEL_ONLY / NOT_HFSS_VALIDATED", n_diodes=len(rows), n_paste_pads=n_paste, pad_mm=[PAD, PAD], pitch_mm=PITCH, pad_to_pad_gap_mm=round(PITCH - PAD, 3),
               pads_on_copper_violations=chk_bad, link_traces_checked=sum(t["kind"] == "link" for t in traces), link_polarity_violations=bad_link, lead_polarity_ok=lead_ok,
               strings=len(strings), diodes_per_string=sorted({sum(1 for r in rows if r["string"] == s) for s in strings}),
               note="anode pad edge is flush with the island edge (local y 0.45) and cathode pad edge flush with the slot edge (local y 1.45); land 0.35 x 0.35 on copper, 1.00 mm gap between lands")
    chk["ok"] = bool(not chk_bad and bad_link == 0 and lead_ok and len(rows) == 240 and n_paste == 480)
    (AS / "assembly_check.json").write_text(json.dumps(chk, indent=1), encoding="utf-8", newline="\n")
    print(json.dumps(chk, indent=1), flush=True)
    dxf(cells, rows, pads, bodies, bars, hexp, traces)
    pdf(cells, rows, pads, bodies, bars, hexp, strings, meta)
    readme(rout, strings)


def dxf(cells, rows, pads, bodies, bars, hexp, traces):
    import ezdxf
    doc = ezdxf.new("R2010", setup=True)
    doc.units = 4
    msp = doc.modelspace()
    for name, col in (("OUTLINE", 7), ("CELL_DISC", 8), ("CELL_ISLAND", 3), ("ASM_PAD_ANODE", 1), ("ASM_PAD_CATHODE", 5), ("ASM_BODY", 2), ("ASM_CATHODE_BAR", 6), ("ASM_ARROW", 4), ("ASM_STRING_PATH", 30), ("ASM_REF", 7), ("CONNECTOR", 50)):
        doc.layers.add(name, color=col)
    msp.add_lwpolyline(hexp, close=True, dxfattribs={"layer": "OUTLINE"})
    for c in cells:
        th = math.radians(float(c["slot_axis_deg"]))
        x, y = float(c["x_mm"]), float(c["y_mm"])
        msp.add_circle((x, y), float(c["disc_r_mm"]), dxfattribs={"layer": "CELL_DISC"})
        msp.add_lwpolyline(rect_pts(x, y, ISL[0], ISL[1], th, ISL[2], ISL[3]), close=True, dxfattribs={"layer": "CELL_ISLAND"})
    byref = {}
    for c in cells:
        byref[(int(c["string"]), int(c["pos_in_string"]))] = c
    cen_of = {r["ref"]: r for r in rows}
    # pads/bodies were appended in cell_table order; rebuild by index
    for idx, c in enumerate(cells):
        s, p = int(c["string"]), int(c["pos_in_string"])
        r = cen_of[f"D{s}-{p:02d}"]
        msp.add_lwpolyline(pads[2 * idx], close=True, dxfattribs={"layer": "ASM_PAD_ANODE"})
        msp.add_lwpolyline(pads[2 * idx + 1], close=True, dxfattribs={"layer": "ASM_PAD_CATHODE"})
        msp.add_lwpolyline(bodies[idx], close=True, dxfattribs={"layer": "ASM_BODY"})
        msp.add_lwpolyline(bars[idx], close=True, dxfattribs={"layer": "ASM_CATHODE_BAR"})
        a, k = (r["anode_x"], r["anode_y"]), (r["cathode_x"], r["cathode_y"])
        msp.add_line(a, k, dxfattribs={"layer": "ASM_ARROW"})
        ux, uy = (k[0] - a[0]) / PITCH, (k[1] - a[1]) / PITCH
        tip = (k[0] - ux * 0.0, k[1] - uy * 0.0)
        for sgn in (-1, 1):
            msp.add_line(tip, (tip[0] - ux * 0.3 - sgn * uy * 0.15, tip[1] - uy * 0.3 + sgn * ux * 0.15), dxfattribs={"layer": "ASM_ARROW"})
        t = msp.add_text(r["ref"], dxfattribs={"layer": "ASM_REF", "height": 0.9, "rotation": 0.0})
        t.set_placement((r["x_mm"] + 1.4, r["y_mm"] + 1.0))
    # string paths: diode centres in order
    for s in sorted({r["string"] for r in rows}):
        pl = [(r["x_mm"], r["y_mm"]) for r in sorted(rows, key=lambda r: r["pos"]) if r["string"] == s]
        msp.add_lwpolyline(pl, dxfattribs={"layer": "ASM_STRING_PATH"})
        msp.add_text(f"STRING {s} start (+)", dxfattribs={"layer": "ASM_REF", "height": 2.0}).set_placement((pl[0][0] - 8, pl[0][1] + 4))
    doc.saveas(AS / "sub6_final_assembly.dxf")


def pdf(cells, rows, pads, bodies, bars, hexp, strings, meta):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.patches import Circle, Polygon
    cm = plt.get_cmap("tab20")
    idx = {(int(c["string"]), int(c["pos_in_string"])): i for i, c in enumerate(cells)}
    conv = ("MODEL_ONLY / NOT_HFSS_VALIDATED. SC-79 diodes, 240 positions.\nAnode pad on the ISLAND, cathode pad on the DISC of the same cell. Rotation = angle of the anode->cathode vector from board +x (ccw).\n"
            "Current along the string: pad_start (+) -> island(1) -> diode(1) -> disc(1) -> bottom trace -> island(2) -> ... -> disc(20) -> pad_end (-). Land 0.35 x 0.35 mm, pitch 1.35 mm, gap 1.00 mm.")

    def base(ax, cells_sel=None):
        ax.add_patch(Polygon(hexp, closed=True, fill=False, ec="k", lw=0.8))
        for c in (cells_sel or cells):
            th = math.radians(float(c["slot_axis_deg"]))
            ax.add_patch(Circle((float(c["x_mm"]), float(c["y_mm"])), float(c["disc_r_mm"]), fc="#f1e7c4", ec="#b9ad7e", lw=0.3))
            ax.add_patch(Polygon(rect_pts(float(c["x_mm"]), float(c["y_mm"]), ISL[0], ISL[1], th, ISL[2], ISL[3]), fc="#e0c15a", ec="#8a7330", lw=0.3))
        ax.set_aspect("equal")

    def diode(ax, i, col, lw=0.4, label=None, fs=4):
        r = rows_by[i]
        ax.add_patch(Polygon(bodies[i], fc="#222", ec="k", lw=0.2))
        ax.add_patch(Polygon(pads[2 * i], fc="#d33", ec="none"))
        ax.add_patch(Polygon(pads[2 * i + 1], fc="#25b", ec="none"))
        ax.annotate("", xy=(r["cathode_x"], r["cathode_y"]), xytext=(r["anode_x"], r["anode_y"]), arrowprops=dict(arrowstyle="->", color=col, lw=lw))
        if label:
            ax.text(r["x_mm"], r["y_mm"] + 1.8, label, fontsize=fs, ha="center", va="bottom", color=col)

    rows_by = {}
    for c_i, c in enumerate(cells):
        rows_by[c_i] = next(r for r in rows if r["ref"] == f"D{int(c['string'])}-{int(c['pos_in_string']):02d}")
    with PdfPages(AS / "sub6_final_assembly.pdf") as pp:
        fig, ax = plt.subplots(figsize=(11.7, 8.3))
        base(ax)
        for s in strings:
            for p in range(1, 21):
                diode(ax, idx[(s, p)], cm((s - 1) % 20), lw=0.6)
            xs = [rows_by[idx[(s, p)]]["x_mm"] for p in range(1, 21)]
            ys = [rows_by[idx[(s, p)]]["y_mm"] for p in range(1, 21)]
            ax.plot(xs, ys, color=cm((s - 1) % 20), lw=0.5, alpha=0.6)
            ax.text(xs[0], ys[0] + 4, f"S{s} +", fontsize=6, color=cm((s - 1) % 20), ha="center")
        ax.autoscale_view()
        ax.set_title("SUB6 final board: diode assembly overview (240 x SC-79)", fontsize=10)
        fig.text(0.04, 0.02, conv, fontsize=6.5, va="bottom")
        pp.savefig(fig)
        plt.close(fig)
        for s in strings:
            fig, ax = plt.subplots(figsize=(11.7, 8.3))
            sel = [cells[idx[(s, p)]] for p in range(1, 21)]
            base(ax, sel)
            for p in range(1, 21):
                diode(ax, idx[(s, p)], "k", 0.7, f"D{s}-{p:02d}\n{rows_by[idx[(s, p)]]['rotation_deg']:.0f}°", 5)
            xs = [rows_by[idx[(s, p)]]["x_mm"] for p in range(1, 21)]
            ys = [rows_by[idx[(s, p)]]["y_mm"] for p in range(1, 21)]
            ax.plot(xs, ys, color="g", lw=0.6, ls="--")
            ax.autoscale_view()
            ax.set_title(f"String {s}: 20 diodes, pos 1 (+, pad {strings[s]['pad_start']}) -> pos 20 (-, pad {strings[s]['pad_end']}); red = anode land (island), blue = cathode land (disc)", fontsize=8)
            fig.text(0.04, 0.02, conv, fontsize=6.5, va="bottom")
            pp.savefig(fig)
            plt.close(fig)
        c0 = cells[idx[(1, 1)]]
        fig, ax = plt.subplots(figsize=(8.3, 8.3))
        th = math.radians(float(c0["slot_axis_deg"]))
        x, y = float(c0["x_mm"]), float(c0["y_mm"])
        base(ax, [c0])
        i = idx[(1, 1)]
        diode(ax, i, "k", 1.0, rows_by[i]["ref"], 9)
        r = rows_by[i]
        ax.set_xlim(r["x_mm"] - 4, r["x_mm"] + 4)
        ax.set_ylim(r["y_mm"] - 4, r["y_mm"] + 4)
        ax.set_title(f"Detail {r['ref']}: land 0.35 x 0.35, pitch 1.35, gap 1.00; anode (red) on island, cathode (blue) on disc; body outline nominal", fontsize=8)
        fig.text(0.04, 0.02, conv, fontsize=6.5, va="bottom")
        pp.savefig(fig)
        plt.close(fig)


def readme(rout, strings):
    t = f"""# Diode assembly README (MODEL_ONLY / NOT_HFSS_VALIDATED)
Generated by scripts/research/sub6_final_assembly.py from the same layout as the Gerbers (../cell_table.csv, ../string_table.csv, ../traces.json).

## Footprint
SC-79 (SOD-523 class). Two square lands 0.35 x 0.35 mm, centre pitch 1.35 mm, i.e. a 1.00 mm gap between the lands, centred on the 1.00 mm diode gap of each switched cell (cell-local y 0.95; lands at local y 0.275 and 1.625, x 0). Each land edge is flush with the copper edge it sits on: the anode land lies inside the island (island edge at local y 0.45), the cathode land lies on the disc rim (slot edge at local y 1.45). The lands are on the existing top copper; no extra copper is added.
The package body outline in the assembly layer ({BODY_L} x {BODY_W} mm, with a cathode bar) is a NOMINAL drawing aid, not taken from a datasheet: use the chosen part's datasheet for body size and land recommendation.

## Polarity convention
**Anode on the island, cathode on the disc** of the same cell. The island via is the anode node, the disc via the cathode node. Conventional current along a string: connector pad_start (+) -> island(1) -> diode(1) -> disc(1) -> bottom trace -> island(2) -> diode(2) -> ... -> disc(20) -> pad_end (-). Position numbers (pos 1..20) increase along the current direction; the (+) pad of string s is `pad_start` and the (-) pad is `pad_end` in ../string_table.csv. Each diode is forward biased (ON state) with 10 mA at about 0.95 V; a string needs about 19.0 V. The script checks that every link trace joins the disc via of position i to the island via of position i+1 and that the start/end leads attach to the island via of position 1 and the disc via of position 20.
Rotation: angle of the anode->cathode vector (the cell-local +y axis) from board +x, counter-clockwise, degrees, in diode_placement.csv (also the explicit anode and cathode land centres). Because each cell is rotated by its own slot axis (plus 180 deg where the string optimisation flipped it), rotations differ from cell to cell; the CSV gives the true value for every diode.

## Files
| file | content |
|---|---|
| sub6_final-F_Paste.gtp | paste layer: 480 lands 0.35 x 0.35 (Gerber RS-274X, FSLAX46Y46, MOMM) |
| sub6_final-F_Assy.gto | assembly layer: body outline and cathode bar (not for fabrication) |
| diode_placement.csv | 240 diodes: ref D<string>-<pos>, centre, rotation, anode/cathode land centres |
| sub6_final_assembly.dxf | drawing, layers ASM_PAD_ANODE, ASM_PAD_CATHODE, ASM_BODY, ASM_CATHODE_BAR, ASM_ARROW (anode->cathode), ASM_STRING_PATH, ASM_REF, plus outline and cell copper |
| sub6_final_assembly.pdf | overview, one page per string (12), one detail page; generated locally, not committed |
| assembly_check.json | pads on copper, pitch/gap, polarity-vs-string checks |

## Solder mask
No top solder mask layer is generated and none is intended: the lands are on bare top copper (the cells are radiating copper). If the manufacturer applies a default mask, order it without top mask (or with the mask clear of the whole top side). Bottom: no mask file either.

## Open points
1. Part choice: a BAR50-type diode is assumed only for the model; check the SC-79 land pattern of the chosen part (this README uses the requested 0.35 x 0.35 mm at 1.35 mm pitch).
2. Solder paste on the thin radiating copper is a reflow thermal-mass question for the assembler (islands are small, discs large); an assembler may prefer conductive epoxy or a stencil adjustment.
3. 240 diodes x 0.95 V per string = about 19 V and 10 mA: the supply and current limit are outside this package.
"""
    (AS / "README_ASSEMBLY.md").write_text(t, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
