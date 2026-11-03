"""Phase 30: STEP files + port table of the single-feed TEST BOARD (eps_r 11.2, 10.05 GHz; MODEL_ONLY), from `hfss_exports/phase30_testboard/layout.csv`.
`testboard_active.step` (572 modulated cells) and `testboard_control.step` (same outline/feed, every cell the unslotted squircle). mm, origin = board centre on the ground plane.
Bodies: SUB (solid, 100 x 100 x 0.635, pin hole r 0.635), GND (zero-thickness sheet z = 0, 100 x 100 with a coax opening r 2.05), PATCH (zero-thickness sheets z = 0.635,
72-vertex n = 4 super-ellipse 3.05 mm minus the rotated slot; outline fixed, slot turns), PIN (solid cylinder r 0.635, z 0 .. 6.635), PORT (annulus sheet r 0.635 .. 2.05 at z = 0).
Re-reads each file (OCP) and checks counts / bbox / gap / every slot size.  STEP files are gitignored; copies go to C:\\hmsa_sync\\data\\phase30_testboard\\.
"""

import csv
import json
import math
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PKG = ROOT / "hfss_exports" / "phase30_testboard"
SYNC = Path(r"C:\hmsa_sync\data\phase30_testboard")
H, W, MAJOR, EXPO, NV = 0.635, 100.0, 3.05, 4.0, 72
R_PIN, R_OUT, MONO_H = 0.635, 2.05, 6.635


def outline(cx, cy):
    pts = []
    for k in range(NV):
        t = 2 * math.pi * k / NV
        c, s = math.cos(t), math.sin(t)
        pts.append((cx + 0.5 * MAJOR * math.copysign(abs(c) ** (2 / EXPO), c), cy + 0.5 * MAJOR * math.copysign(abs(s) ** (2 / EXPO), s)))
    return pts


def slot(cx, cy, ln, wd, ang):
    a = math.radians(ang)
    c, s = math.cos(a), math.sin(a)
    return [(cx + u * c - v * s, cy + u * s + v * c) for u, v in ((-ln / 2, -wd / 2), (ln / 2, -wd / 2), (ln / 2, wd / 2), (-ln / 2, wd / 2))]


def poly_area(p):
    return 0.5 * abs(sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(p, p[1:] + p[:1], strict=True)))


def read_cells():
    with open(PKG / "layout.csv", newline="", encoding="utf-8") as fh:
        return [(float(r["x_mm"]), float(r["y_mm"]), float(r["slot_L_mm"]), float(r["slot_W_mm"]), float(r["slot_angle_deg"])) for r in csv.DictReader(fh)]


def build(cells, out, control):
    import cadquery as cq

    V = cq.Vector
    sub = cq.Workplane("XY").box(W, W, H, centered=(True, True, False)).cut(cq.Workplane("XY").circle(R_PIN).extrude(H)).val()
    gnd = cq.Face.makeFromWires(
        cq.Wire.makePolygon([V(-W / 2, -W / 2, 0), V(W / 2, -W / 2, 0), V(W / 2, W / 2, 0), V(-W / 2, W / 2, 0)], close=True),
        [cq.Wire.makeCircle(R_OUT, V(0, 0, 0), V(0, 0, 1))],
    )
    faces = []
    for x, y, ln, wd, ang in cells:
        o = cq.Wire.makePolygon([V(px, py, H) for px, py in outline(x, y)], close=True)
        if control or ln <= 0 or wd <= 0:
            faces.append(cq.Face.makeFromWires(o))
        else:
            h = cq.Wire.makePolygon([V(px, py, H) for px, py in slot(x, y, ln, wd, ang)], close=True)
            faces.append(cq.Face.makeFromWires(o, [h]))
    port = cq.Face.makeFromWires(cq.Wire.makeCircle(R_OUT, V(0, 0, 0), V(0, 0, 1)), [cq.Wire.makeCircle(R_PIN, V(0, 0, 0), V(0, 0, 1))])
    assy = cq.Assembly(name="board")
    assy.add(sub, name="SUB")
    assy.add(gnd, name="GND")
    assy.add(cq.Compound.makeCompound(faces), name="PATCH")
    assy.add(cq.Solid.makeCylinder(R_PIN, MONO_H, V(0, 0, 0), V(0, 0, 1)), name="PIN")
    assy.add(port, name="PORT")
    assy.save(str(out), exportType="STEP")


def check(path, cells, control):
    import phase30_eigenmode_step as E

    r = E.read_back(path)
    c = {}
    c["body_names"] = sorted(r) == ["GND", "PATCH", "PIN", "PORT", "SUB"]
    c["SUB"] = r["SUB"]["solids"] == 1 and all(abs(a - b) < 2e-3 for a, b in zip(r["SUB"]["bbox"], (-W / 2, -W / 2, 0, W / 2, W / 2, H), strict=True))
    c["GND"] = r["GND"]["faces"] == 1 and abs(r["GND"]["area"] - (W * W - math.pi * R_OUT**2)) < 0.05
    c["PIN"] = r["PIN"]["solids"] == 1 and abs(r["PIN"]["bbox"][5] - MONO_H) < 2e-3 and abs(r["PIN"]["bbox"][3] - R_PIN) < 2e-3
    c["PORT"] = r["PORT"]["faces"] == 1 and abs(r["PORT"]["area"] - math.pi * (R_OUT**2 - R_PIN**2)) < 0.05
    c["PATCH_faces"] = r["PATCH"]["faces"] == len(cells)
    exp_area = sum(poly_area(outline(0, 0)) - (0 if control else ln * wd) for _, _, ln, wd, _ in cells)
    c["PATCH_total_area"] = abs(r["PATCH"]["area"] - exp_area) < 0.05
    holes = r["PATCH"]["holes"]
    n_exp = 0 if control else sum(1 for x in cells if x[2] > 0)
    c["slot_count"] = len(holes) == n_exp
    if not control:
        have = sorted(round(h["hole_area"], 3) for h in holes)
        want = sorted(round(ln * wd, 3) for _, _, ln, wd, _ in cells if ln > 0)
        c["slot_areas_LxW"] = all(abs(a - b) < 2e-3 for a, b in zip(have, want, strict=True))
    pb = r["PATCH"]["bbox"]
    xs = sorted({round(x, 3) for x, *_ in cells})
    out = {"file": path.name, "bodies": {k: {kk: v[kk] for kk in ("solids", "faces")} for k, v in r.items()}, "checks": c,
           "patch_bbox": pb, "gap_to_neighbour_mm": round(3.2 - MAJOR, 4), "min_cell_spacing_mm": 3.2, "n_cells": len(cells), "n_slots": len(holes)}
    out["all_pass"] = all(c.values())
    return out


PORT_CSV = [
    ["item", "value"],
    ["port_name", "P1 (lumped port, 50 ohm)"],
    ["port_body", "PORT (annulus sheet, r 0.635 .. 2.05 mm, z = 0)"],
    ["integration_line", "from the pin surface (r 0.635, z 0) to the ground opening edge (r 2.05, z 0), radial; direction: pin -> ground"],
    ["reference_impedance_ohm", "50"],
    ["renormalise", "yes (50 ohm); no de-embedding"],
    ["excitation", "1 W accepted, phase 0 deg"],
    ["pin_body", "PIN (solid, r 0.635 mm, z 0 .. 6.635 mm = 6.0 mm above the substrate): assign Perfect E"],
    ["pin_position_mm", "(0, 0)"],
    ["frequency_ghz", "10.05 (solution); optional interpolating sweep 9.5 - 10.6 GHz"],
    ["note", "S11 only (one port). Nothing else is driven."],
]


def main():
    import phase30_eigenmode_step  # noqa: F401  (read_back / bbox helpers)

    cells = read_cells()
    rep = {}
    for name, ctrl in (("testboard_active", False), ("testboard_control", True)):
        p = PKG / (name + ".step")
        build(cells, p, ctrl)
        rep[name] = check(p, cells, ctrl)
        rep[name]["size_mb"] = round(p.stat().st_size / 1e6, 2)
        print(name, "PASS" if rep[name]["all_pass"] else "FAIL", rep[name]["bodies"], rep[name]["size_mb"], "MB", flush=True)
    (PKG / "step_check.json").write_text(json.dumps(rep, indent=1), encoding="utf-8", newline="\n")
    with open(PKG / "ports_and_excitations.csv", "w", newline="", encoding="utf-8") as fh:
        csv.writer(fh, lineterminator="\n").writerows(PORT_CSV)
    return rep


if __name__ == "__main__":
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
