"""HFSS package (NOT SENT): the chosen 1-bit cell pair of the grouped coded hologram, ONE unit cell, two STEP variants (ON: sheet D_ON, OFF: sheet D_OFF), eigenmode at 2 phases, PEC lid.
Cell: split patch 6.0 x 4.8 mm (two plates 2.8 x 4.8, diode gap 0.4 mm, current along x), MA4AGP907 footprint 0.4 x 1.2 mm, RT5880 1.575 mm, 6.4 mm lattice. MODEL_ONLY predictions.
    python scripts/research/sidequest_coded_package.py  -> hfss_exports/sidequest_coded_cell/"""
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
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

import sidequest_gate_cell as gc
from src.dispersion import sheet_model as sm
from src.frequency import free_space_wavenumber
from src.hologram.coded_cell_eval import sheet_tensor
from src.transmitarray.export import inexact_literals, snap_step_text, verify_snap

OUT = ROOT / "hfss_exports/sidequest_coded_cell"
# all lengths in mm, multiples of 0.001
P, H, T, LID = 6.4, 1.575, 0.035, 30.0
A, HP, G, WD = 6.0, 4.8, 0.4, 1.2  # patch length (x), height (y), diode gap, diode footprint width
PLATE = (A - G) / 2  # 2.8
XC = G / 2 + PLATE / 2  # plate centre 1.6
PHASES = (90, 120)  # periodic phase advance along x [deg]; y faces 0 deg


def box(cx, cy, lx, ly, z0, z1):
    return cq.Workplane("XY").workplane(offset=z0).center(cx, cy).rect(lx, ly).extrude(z1 - z0)


def parts():
    pl, pr = box(-XC, 0, PLATE, HP, 0, T), box(XC, 0, PLATE, HP, 0, T)
    sub, gnd = box(0, 0, P, P, -H, 0), box(0, 0, P, P, -H - T, -H)
    air = box(0, 0, P, P, 0, LID).cut(pl).cut(pr)
    return {"Substrate": sub, "Ground": gnd, "Patch_L": pl, "Patch_R": pr, "Air": air}


def sheet(zs=T / 2):
    pts = [cq.Vector(-G / 2, -WD / 2, zs), cq.Vector(G / 2, -WD / 2, zs), cq.Vector(G / 2, WD / 2, zs), cq.Vector(-G / 2, WD / 2, zs)]
    return cq.Workplane("XY").add(cq.Face.makeFromWires(cq.Wire.makePolygon(pts, close=True)))


def intended():
    xy = {0.0, P / 2, A / 2, G / 2, XC, XC - PLATE / 2, XC + PLATE / 2, HP / 2, WD / 2}
    xy |= {-v for v in xy}
    z = {0.0, T, T / 2, -H, -H - T, LID}
    return dict(xy=np.array(sorted(round(v, 6) for v in xy)), z=np.array(sorted(round(v, 6) for v in z)), dir=np.array([-1.0, 0.0, 1.0]))


def build(state):
    asm = cq.Assembly()
    cols = {"Substrate": (0.25, 0.55, 0.25), "Ground": (0.6, 0.6, 0.65), "Patch_L": (0.9, 0.7, 0.2), "Patch_R": (0.9, 0.7, 0.2), "Air": (0.7, 0.8, 1.0)}
    for nm, o in parts().items():
        asm.add(o, name=nm, color=cq.Color(*cols[nm], 1))
    asm.add(sheet(), name="D_ON" if state == "ON" else "D_OFF", color=cq.Color(0.9, 0.1, 0.1, 1))
    path = OUT / f"cell_{state}.step"
    asm.export(str(path), exportType="STEP")
    raw = path.read_text(encoding="utf-8")
    snapped, stats = snap_step_text(raw, intended(), description=f"sidequest coded 1-bit cell {state}, mm, 0.001 mm grid")
    path.write_text(snapped, encoding="utf-8", newline="\n")
    return raw, snapped, stats


def check(state, raw, snapped, stats):
    """EXACT_DIMS: re-open the written STEP; body count, per-body bounding box and volume against the intended values; no inexact literal left; snap safe."""
    shp = cq.importers.importStep(str(OUT / f"cell_{state}.step"))
    solids = shp.solids().vals()
    faces = shp.faces().vals()
    want = {  # (dx, dy, dz, volume)
        "Substrate": (P, P, H), "Ground": (P, P, T), "Patch_L": (PLATE, HP, T), "Patch_R": (PLATE, HP, T), "Air": (P, P, LID)}
    vol = {"Substrate": P * P * H, "Ground": P * P * T, "Patch_L": PLATE * HP * T, "Patch_R": PLATE * HP * T, "Air": P * P * LID - 2 * PLATE * HP * T}
    got = sorted(((round(s.BoundingBox().xlen, 6), round(s.BoundingBox().ylen, 6), round(s.BoundingBox().zlen, 6), s.Volume()) for s in solids))
    exp = sorted(((round(w[0], 6), round(w[1], 6), round(w[2], 6), vol[k]) for k, w in want.items()))
    ok_dims = len(got) == len(exp) and all(max(abs(g[i] - e[i]) for i in range(3)) < 1e-6 and abs(g[3] - e[3]) < 1e-6 * max(1.0, e[3]) for g, e in zip(got, exp, strict=True))
    sheet_faces = [f for f in faces if f.BoundingBox().zlen < 1e-9 and abs(f.BoundingBox().xlen - G) < 1e-6 and abs(f.BoundingBox().ylen - WD) < 1e-6 and abs(f.BoundingBox().zmin - T / 2) < 1e-9]
    text = (OUT / f"cell_{state}.step").read_text(encoding="utf-8")
    bad = inexact_literals(text, intended())
    snap_ok = verify_snap(raw, snapped) if True else None
    return {"state": state, "bodies": len(solids), "sheets": len(sheet_faces), "EXACT_DIMS": "PASS" if (ok_dims and not bad and len(solids) == 5 and len(sheet_faces) == 1) else "FAIL",
            "inexact_literals": len(bad), "n_snapped": stats["n_snapped"], "max_snap_mm": stats["max_snap_mm"], "snap_check": str(snap_ok)[:200], "bbox_mm": [round(v, 6) for v in (shp.val().BoundingBox().xmin, shp.val().BoundingBox().xmax, shp.val().BoundingBox().zmin, shp.val().BoundingBox().zmax)]}


def dispersion():
    """Uniform-sheet TM0 band of the cell: Y_xx(f) of the MoM cell (mesh 32), phase advance beta*a = phase -> f. Surface band below the light line."""
    fg = np.array([7.0, 8.5, 9.5, 10.0, 10.54, 11.0, 11.5, 12.5, 13.5, 14.5]) * 1e9
    ys = {}
    for st in ("OFF", "ON"):
        ys[st] = np.array([sheet_tensor(gc.make_cell(6.0e-3, 4.8e-3, 0.4e-3, 1.2e-3, st), f, 32)[0, 0] for f in fg])
    out = {"freq_GHz": (fg / 1e9).tolist(), "Yxx_mS": {s: [[y.real * 1e3, y.imag * 1e3] for y in v] for s, v in ys.items()}}
    pred = []
    for ph in PHASES:
        row = {"phase_deg": ph, "light_line_GHz": round(299792458.0 * np.radians(ph) / (P * 1e-3) / (2 * np.pi) / 1e9, 3)}
        for st in ("OFF", "ON"):
            cs = CubicSpline(fg, ys[st].imag)

            def f_res(f):
                b = sm.find_beta_from_sheet_admittance(1j * float(cs(f)), gc.SUB, f).beta_rad_per_m
                return b * P * 1e-3 - np.radians(ph)

            try:
                fr = brentq(f_res, fg[0], fg[-1], xtol=1e3)
                if fr > fg[-1] - 5e7:  # at the end of the Y(f) table: not a prediction
                    row[f"{st}_f_GHz"] = f">{fg[-1] / 1e9:.1f} (outside the Y(f) table)"
                    continue
                row[f"{st}_f_GHz"] = round(fr / 1e9, 3)
                row[f"{st}_beta_over_k0"] = round(float(np.radians(ph) / (P * 1e-3) / free_space_wavenumber(fr)), 4)
            except ValueError:
                row[f"{st}_f_GHz"] = None
        pred.append(row)
    out["predictions"] = pred
    return out


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    checks = []
    for st in ("ON", "OFF"):
        raw, snapped, stats = build(st)
        checks.append(check(st, raw, snapped, stats))
        print(checks[-1], flush=True)
    d = dispersion()
    d["label"] = "MODEL_ONLY: MoM Y_xx(f) of the cell as a uniform sheet (mesh 32), TM0 transverse resonance; the cell's gaps and the band edge near 180 deg are not modelled; uncertainty about +-5 % in frequency"
    d["cell_mm"] = {"lattice": P, "patch_x": A, "patch_y": HP, "diode_gap": G, "footprint_y": WD, "substrate": H, "metal": T, "lid_above_patch_plane": LID}
    d["diode_values_10p54"] = {"D_ON_R_ohm": 2.44, "D_OFF_C_fF": 28.16}
    d["exact_dims"] = checks
    (OUT / "predictions.json").write_text(json.dumps(d, indent=1, default=float), encoding="utf-8", newline="\n")
    for r in d["predictions"]:
        print(r)
