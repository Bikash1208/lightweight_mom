"""Side quest step 5 (gate test): STEP package of an 83.2 x 70.4 mm RT5880 test board, centre probe, uniform reference surface, 3-row gate ring (+1 transition row each
side) on the -x half only, diode sites as named sheets D_ON_k / D_OFF_k.   python scripts/research/sidequest_gate_package.py   -> hfss_exports/sidequest_gate_test/
MODEL_ONLY: reference/transition cells are square patches (equivalent sheets of the hologram mean, NOT the slotted circles of the hologram)."""
import csv
import json
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import cadquery as cq
import numpy as np

import sidequest_gate_cell as gc
import sidequest_gate_final as gf

OUT = ROOT / "hfss_exports/sidequest_gate_test"
P = 6.4
NX, NY = 13, 11  # cells -> 83.2 x 70.4 mm
H, T = 1.575, 0.035  # substrate / copper slab thickness [mm]
A_REF, A_TRANS = 5.8, 5.4  # square patch sides [mm] (mesh-64 b_xx = 7.49 / 5.44 mS)
R_PIN, R_OUT, PIN_TOP = 0.635, 2.05, 6.0
ZONES = {"inner": (9.6, 16.0), "gate": (16.0, 35.2), "outer": (35.2, 42.0)}


def cells():
    out = []
    for i in range(-(NX // 2), NX // 2 + 1):
        for j in range(-(NY // 2), NY // 2 + 1):
            x, y = i * P, j * P
            rho = float(np.hypot(x, y))
            if rho < 6.5:  # probe keepout (centre cell + 4 neighbours), as Phase 26
                continue
            kind, orient = "ref", ""
            if x < -1e-9:  # gate half
                for k, (lo, hi) in ZONES.items():
                    if lo <= rho < hi:
                        kind = "trans" if k in ("inner", "outer") else "gate"
                if kind == "gate":
                    orient = "x" if abs(y) <= abs(x) + 1e-9 else "y"
            out.append((i, j, x, y, rho, kind, orient))
    return out


def rect(cx, cy, lx, ly, z0, t):
    return cq.Workplane("XY").workplane(offset=z0).center(cx, cy).rect(lx, ly).extrude(t)


def dogbone_solids(x, y, orient):
    c = gf.CAND
    uc = gc.make_dogbone(c["e"] * 1e-3, c["s"] * 1e-3, c["h"] * 1e-3, c["g"] * 1e-3, "OFF", lt=c["lt"] * 1e-3, wd=c["wd"] * 1e-3, direction=orient)
    rects = uc.geometry.metal[:-1]  # last rectangle is the diode bridge (a gap in CAD)
    solids = [rect(x + r.cx * 1e3, y + r.cy * 1e3, r.length * 1e3, r.width * 1e3, 0.0, T) for r in rects]
    b = uc.geometry.metal[-1]
    site = (x + b.cx * 1e3, y + b.cy * 1e3, b.length * 1e3, b.width * 1e3)
    return solids, site


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cl = cells()
    ref, trans, gate = [], [], []
    sites = []
    for (i, j, x, y, rho, kind, orient) in cl:
        if kind == "ref":
            ref.append(rect(x, y, A_REF, A_REF, 0.0, T))
        elif kind == "trans":
            trans.append(rect(x, y, A_TRANS, A_TRANS, 0.0, T))
        else:
            sol, site = dogbone_solids(x, y, orient)
            m = sol[0]
            for s in sol[1:]:
                m = m.union(s)
            gate.append(m)
            sites.append((x, y, orient, rho, site))
    sx, sy = NX * P, NY * P
    sub = cq.Workplane("XY").workplane(offset=-H).rect(sx, sy).extrude(H).cut(cq.Workplane("XY").workplane(offset=-H).circle(R_PIN).extrude(H))
    gnd = cq.Workplane("XY").workplane(offset=-H - T).rect(sx, sy).extrude(T).cut(cq.Workplane("XY").workplane(offset=-H - T).circle(R_OUT).extrude(T))
    pin = cq.Workplane("XY").workplane(offset=-H - T).circle(R_PIN).extrude(H + T + PIN_TOP)
    port = cq.Workplane("XY").workplane(offset=-H - T).circle(R_OUT).circle(R_PIN).extrude(0.0001) if False else None
    ring = cq.Face.makeFromWires(cq.Wire.makeCircle(R_OUT, cq.Vector(0, 0, -H - T), cq.Vector(0, 0, 1)), [cq.Wire.makeCircle(R_PIN, cq.Vector(0, 0, -H - T), cq.Vector(0, 0, 1))])
    port = cq.Workplane("XY").add(ring)

    def comp(lst):
        return cq.Workplane("XY").add(cq.Compound.makeCompound([s.val() for s in lst]))

    asm = cq.Assembly()
    asm.add(sub, name="Substrate", color=cq.Color(0.25, 0.55, 0.25, 1))
    asm.add(gnd, name="Ground", color=cq.Color(0.6, 0.6, 0.65, 1))
    asm.add(pin, name="Probe_Pin", color=cq.Color(1, 1, 0, 1))
    asm.add(port, name="Port1", color=cq.Color(0, 0, 1, 1))
    asm.add(comp(ref), name="Metal_REF", color=cq.Color(0.8, 0.55, 0.1, 1))
    asm.add(comp(trans), name="Metal_TRANS", color=cq.Color(0.9, 0.7, 0.2, 1))
    asm.add(comp(gate), name="Metal_GATE", color=cq.Color(0.9, 0.4, 0.1, 1))
    zs = T / 2
    for grp in ("D_ON", "D_OFF"):
        for k, (x, y, orient, rho, (cx, cy, lx, ly)) in enumerate(sites, 1):
            face = cq.Face.makeFromWires(cq.Wire.makePolygon([cq.Vector(cx - lx / 2, cy - ly / 2, zs), cq.Vector(cx + lx / 2, cy - ly / 2, zs), cq.Vector(cx + lx / 2, cy + ly / 2, zs), cq.Vector(cx - lx / 2, cy + ly / 2, zs)], close=True))
            asm.add(cq.Workplane("XY").add(face), name=f"{grp}_{k:02d}", color=cq.Color(0.9, 0.1, 0.1, 1) if grp == "D_ON" else cq.Color(0.1, 0.1, 0.9, 1))
    asm.export(str(OUT / "gate_test_board.step"), exportType="STEP")
    with open(OUT / "sites.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["k", "x_mm", "y_mm", "current_axis", "rho_mm", "sheet_x_mm", "sheet_y_mm", "sheet_len_along_axis_mm", "sheet_width_mm", "D_ON_name", "D_OFF_name"])
        for k, (x, y, orient, rho, (cx, cy, lx, ly)) in enumerate(sites, 1):
            w.writerow([k, x, y, orient, round(rho, 2), round(cx, 3), round(cy, 3), lx if orient == "x" else ly, ly if orient == "x" else lx, f"D_ON_{k:02d}", f"D_OFF_{k:02d}"])
    with open(OUT / "layout.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["i", "j", "x_mm", "y_mm", "rho_mm", "kind", "axis"])
        for r in cl:
            w.writerow([r[0], r[1], r[2], r[3], round(r[4], 2), r[5], r[6]])
    return {"n_ref": len(ref), "n_trans": len(trans), "n_gate": len(gate), "n_sites": len(sites), "n_x": sum(1 for s in sites if s[2] == "x"), "n_y": sum(1 for s in sites if s[2] == "y")}


if __name__ == "__main__":
    print(main())
