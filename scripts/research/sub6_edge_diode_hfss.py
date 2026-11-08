"""SUB6 edge-diode cell: HFSS EIGENMODE package (MODEL_ONLY / IDEAL_BIAS / NOT_HFSS_VALIDATED; predictions written BEFORE any HFSS run).

Cell window = the edge-diode model shifted by p/2 in x: x, y in [-p/2, p/2], the diode gap (1.00 mm) and the diode sheets at the CENTRE, the cell
edges x = +-p/2 run through the disc centres (window = half of two discs + the gap + the diode).  No vias.  RT5880 3.175 mm on a PEC ground at z = 0;
AIR starts at the substrate top (no air below the ground); PEC lid at 60 mm.  Slot (rim intact) cut along the slot axis through each disc centre
(theta 90 deg: along y, so the slot is a half-slot on each window edge).

Diode = HFSS Lumped RLC on two series sheets (each 0.50 x 0.35 mm, current line along x), which realises (C || R) + L EXACTLY
(a single parallel-RLC boundary cannot make a series L):  OFF = sheet A {R 5 kohm || C 0.15 pF} + sheet B {L 0.6 nH};  ON = sheet A {R 3 ohm} + sheet B {L 0.6 nH}.

Eigenmode runs (crucial first):  {OFF,ON} x {along (Bloch phase along the slot axis y), across (phase along x)}.  Tier-1 phase = beta(state, direction) * p, with beta
from the MoM Y_uu / Y_vv at 5.8 GHz (self-consistent, mesh 64).  Prediction: the bound eigenmode sits at 5.8 GHz in each run.
STEP: exact-coordinate (all 3-D vertices on the 0.001 mm grid, verified), arcs are analytic circles.

    python sub6_edge_diode_hfss.py [row_id]       -> hfss_exports/sub6_edge_diode_eigenmode/
"""
import json
import sys
from decimal import Decimal
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import phase26_eigenmode_package as p26
import sub6_edge_diode_cell as E
import sub6_slotdisc_grid as G
import sub6_slotdisc_lib as S
from src.dispersion.sheet_model import find_beta_from_sheet_admittance
from src.frequency import free_space_wavenumber

OUT = G.ROOT / "hfss_exports/sub6_edge_diode_eigenmode"
P, H = E.P, E.SUB[3]
RR = E.R
TAB = 0.02  # mm: PEC tab straddling each disc edge inside the sheet width, so the sheets touch metal (the disc edge deviates 0.003 mm over +-0.175 mm)

TEMPLATE_PATCH = '''
def rect_sheet(name, x0, y0, x1, y1, z):
    oEditor.CreateRectangle(["NAME:RectangleParameters", "IsCovered:=", True, "XStart:=", mm(x0), "YStart:=", mm(y0), "ZStart:=", mm(z),
        "Width:=", mm(x1 - x0), "Height:=", mm(y1 - y0), "WhichAxis:=", "Z"], attrs(name, "vacuum", True))

# two half discs (centres on the cell edges x = -+P/2, radius RD), each minus the half slot (slot along y: length SLOT_L, width SLOT_W)
for nm, cx, sgn in (("PatchL", -P / 2, 1.0), ("PatchR", P / 2, -1.0)):
    oEditor.CreateCircle(["NAME:CircleParameters", "IsCovered:=", True, "XCenter:=", mm(cx), "YCenter:=", "0mm",
        "ZCenter:=", mm(H), "Radius:=", mm(RD), "WhichAxis:=", "Z", "NumSegments:=", "0"], attrs(nm, "vacuum", True))
    rect_sheet(nm + "_slot", min(cx - SLOT_W / 2, cx + SLOT_W / 2), -SLOT_L / 2, max(cx - SLOT_W / 2, cx + SLOT_W / 2), SLOT_L / 2, H)
    oEditor.Subtract(["NAME:Selections", "Blank Parts:=", nm, "Tool Parts:=", nm + "_slot"], ["NAME:SubtractParameters", "KeepOriginals:=", False])
    # cut away the half of the disc outside the window
    rect_sheet(nm + "_out", (cx - RD - 1.0) if sgn > 0 else cx, -RD - 1.0, cx if sgn > 0 else (cx + RD + 1.0), RD + 1.0, H)
    oEditor.Subtract(["NAME:Selections", "Blank Parts:=", nm, "Tool Parts:=", nm + "_out"], ["NAME:SubtractParameters", "KeepOriginals:=", False])
# PEC tabs under the diode-sheet ends
rect_sheet("TabL", -GAP / 2 - TAB, -SW / 2, -GAP / 2, SW / 2, H)
rect_sheet("TabR", GAP / 2, -SW / 2, GAP / 2 + TAB, SW / 2, H)
# diode sheets (series): A (x -0.5..0.0) R||C, B (x 0.0..0.5) L
rect_sheet("DiodeA", -GAP / 2, -SW / 2, 0.0, SW / 2, H)
rect_sheet("DiodeB", 0.0, -SW / 2, GAP / 2, SW / 2, H)
'''


def mkscript(name, phase_x, phase_y, state):
    t = p26.TEMPLATE
    a = t.index("# patch: circle minus a rotated rectangular slot")
    b = t.index('oEditor.CreateRectangle(["NAME:RectangleParameters", "IsCovered:=", True, "XStart:=", mm(-P / 2), "YStart:=", mm(-P / 2),\n    "ZStart:=", "0mm"')
    t = t[:a] + TEMPLATE_PATCH + "\n" + t[b:]
    t = t.replace("P = 6.4; H = 1.575; LID = @@LID@@; D = @@D@@", f"P = {P:g}; H = {H:g}; LID = @@LID@@; RD = {RR:g}; GAP = {E.GAP:g}; SW = 0.35; TAB = {TAB:g}")
    t = t.replace("SLOT_L = @@SLOT_L@@; SLOT_W = @@SLOT_W@@; SLOT_ANGLE = @@ANGLE@@", "SLOT_L = @@SL@@; SLOT_W = @@SW@@")
    t = t.replace('PHASE_Y = "0deg"', f'PHASE_Y = "{phase_y:.3f}deg"')
    t = t.replace('for nm in ("Patch", "Ground", "Lid"):', 'for nm in ("PatchL", "PatchR", "TabL", "TabR", "Ground", "Lid"):')
    rc = {
        "OFF": '"UseResist:=", True, "Resistance:=", "5000ohm", "UseInduct:=", False, "UseCap:=", True, "Capacitance:=", "0.15pF"',
        "ON": '"UseResist:=", True, "Resistance:=", "3ohm", "UseInduct:=", False, "UseCap:=", False',
    }[state]
    rlc = f"""
for nm, obj, xs, xe, opts in (("RLC_A", "DiodeA", -GAP / 2, 0.0, [{rc}]),
                              ("RLC_B", "DiodeB", 0.0, GAP / 2, ["UseResist:=", False, "UseInduct:=", True, "Inductance:=", "0.6nH", "UseCap:=", False])):
    oBnd.AssignLumpedRLC(["NAME:" + nm, "Objects:=", [obj],
        ["NAME:CurrentLine", "Start:=", [mm(xs), "0mm", mm(H)], "End:=", [mm(xe), "0mm", mm(H)]]] + opts)
"""
    t = t.replace("def faces_at(x, y):", rlc + "\ndef faces_at(x, y):", 1)
    hdr = [ln for ln in t.split("\n") if ln.startswith("# HFSS (AEDT) IronPython script")][0]
    t = t.replace(hdr, f'# HFSS (AEDT) IronPython script -- sub6 edge-diode cell EIGENMODE "{name}" (generated by lwmom). NOT run in AEDT by the author. MODEL_ONLY predictions: README / predictions.json.')
    vals = {"NAME": name, "EPS": repr(E.SUB[1]), "TAND": repr(E.SUB[2]), "LID": "60.0", "SLOTDESC": "(see body)", "D": "0", "PX": f"{phase_x:.3f}", "PD": "0"}
    for k, v in vals.items():
        t = t.replace("@@%s@@" % k, v)
    return t


def tier1(row_id):
    """beta(state, direction) from the MoM Y_uu (kt along the slot axis) / Y_vv (kt across) at 5.8 GHz, mesh-64 nodes."""
    nodes = G.load_nodes(E.ckpt(row_id))
    sub = S.make_sub(E.SUB)
    k0 = free_space_wavenumber(S.F)
    out = {}
    for st in ("OFF_BAR50", "ON_BAR50", "NO_DIODE"):
        out[st] = {}
        for d, idx, nm in (("x", 0, "along"), ("y", 1, "across")):
            ks = sorted(nodes)
            kt = 1.06
            for _ in range(60):
                i = int(np.clip(np.searchsorted(ks, kt) - 1, 0, len(ks) - 2))
                w = (kt - ks[i]) / (ks[i + 1] - ks[i])
                y = (1 - w) * nodes[ks[i]][d][st][idx, idx] + w * nodes[ks[i + 1]][d][st][idx, idx]
                b = find_beta_from_sheet_admittance(1j * y.imag, sub, S.F).beta_rad_per_m / k0
                if abs(b - kt) < 1e-6:
                    break
                kt = 0.5 * (kt + b)
            out[st][nm] = dict(beta_over_k0=float(b), phase_deg_per_cell=float(np.degrees(b * k0 * P * 1e-3)), Y_mS_imag=float(y.imag * 1e3))
    return out


def step(row, path):
    import cadquery as cq
    from src.export.exact_board import verify_step_exact

    def d(v):
        return float(Decimal(str(round(v, 6))).quantize(Decimal("0.001")))

    L2, W2, r, p, g = d(row["L"] / 2), d(row["W"] / 2), d(RR), d(P / 2), d(E.GAP / 2)
    sw = 0.175
    t = TAB

    def vec(x, y, z):
        return cq.Vector(x, y, z)

    def half(sgn):  # sgn=+1: PatchL (centre -p, window side x > -p); sgn=-1: PatchR
        cx = -sgn * p
        ex = cx + sgn * r
        wp = (cq.Workplane("XY").workplane(offset=H).moveTo(cx, -r).threePointArc((ex, 0.0), (cx, r)).lineTo(cx, L2).lineTo(cx + sgn * W2, L2)
              .lineTo(cx + sgn * W2, -L2).lineTo(cx, -L2).close())
        return cq.Face.makeFromWires(wp.wires().val())

    def rect(x0, y0, x1, y1, z):
        return cq.Face.makeFromWires(cq.Wire.makePolygon([vec(x0, y0, z), vec(x1, y0, z), vec(x1, y1, z), vec(x0, y1, z)], close=True))

    patches = [half(1), half(-1), rect(-g - t, -sw, -g, sw, H), rect(g, -sw, g + t, sw, H)]
    asm = cq.Assembly(name="sub6_edge_diode_cell")
    asm.add(cq.Solid.makeBox(2 * p, 2 * p, H, vec(-p, -p, 0)), name="Substrate")
    asm.add(cq.Face.makePlane(2 * p, 2 * p, cq.Vector(0, 0, 0), cq.Vector(0, 0, 1)), name="Ground")
    asm.add(cq.Compound.makeCompound(patches), name="Patches")
    asm.add(rect(-g, -sw, 0.0, sw, H), name="DiodeA")
    asm.add(rect(0.0, -sw, g, sw, H), name="DiodeB")
    asm.save(str(path))
    return verify_step_exact(path, strict_points=False)


def main(row_id):
    OUT.mkdir(parents=True, exist_ok=True)
    row = [r for r in E.rows() if r["id"] == row_id][0]
    t1 = tier1(row_id)
    runs = {}
    for st, key in (("OFF", "OFF_BAR50"), ("ON", "ON_BAR50")):
        for dirn in ("along", "across"):
            ph = t1[key][dirn]["phase_deg_per_cell"]
            phx, phy = (0.0, ph) if dirn == "along" else (ph, 0.0)  # theta 90: the slot axis is y
            nm = f"EDGE_{st}_{dirn}"
            s = mkscript(nm, phx, phy, st).replace("@@SL@@", repr(round(row["L"], 3))).replace("@@SW@@", repr(round(row["W"], 3)))
            (OUT / f"build_{nm}.py").write_text(s, encoding="utf-8", newline="\n")
            runs[nm] = dict(state=st, direction=dirn, phase_x_deg=phx, phase_y_deg=phy, beta_over_k0=t1[key][dirn]["beta_over_k0"], predicted_eigenfrequency_ghz=5.8)
    n_pts = step(row, OUT / "cell_edge_diode.step")
    preds = dict(labels=["MODEL_ONLY", "IDEAL_BIAS", "NOT_HFSS_VALIDATED", "PREDICTED_BEFORE_HFSS"], row=row, tier1_beta_phase=t1, runs=runs, step_vertices_checked=n_pts)
    (OUT / "predictions.json").write_text(json.dumps(preds, indent=1))
    print({k: (round(v["phase_x_deg"], 2), round(v["phase_y_deg"], 2)) for k, v in runs.items()}, "STEP coords checked", n_pts)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "EDGE_th90_q0.7_W2.4")
