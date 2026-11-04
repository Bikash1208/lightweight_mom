"""HFSS package (NOT SENT): ONE unit cell of the crossed two-diode rim cell, four STEP files (states OFF/OFF, ON/OFF, OFF/ON, ON/ON; x diode / y diode) with named lumped sheets, eigenmode 90 and 120 deg
along x and along y, PEC lid. Cell: 6.0 mm disc, two 0.4 mm slits (length 2*(3.0-rim)), diode footprint 0.4 x 1.2 mm at (0, +1.6) [x diode, current along x] and (+1.6, 0) [y diode]. MODEL_ONLY.
    python sidequest_coded_cross_package.py <rim_mm>   -> hfss_exports/sidequest_coded_cross_cell/"""
import itertools
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

import sidequest_coded_cross_cell as X
import sidequest_gate_cell as gc
from src.dispersion import sheet_model as sm
from src.frequency import free_space_wavenumber
from src.hologram.coded_cell_eval import sheet_tensor
from src.transmitarray.export import inexact_literals, snap_step_text, verify_snap

OUT = ROOT / "hfss_exports/sidequest_coded_cross_cell"
P, H, T, LID = 6.4, 1.575, 0.035, 30.0
R, G, WD, POS = 3.0, 0.4, 1.2, 1.6
RIM = float(sys.argv[1]) if len(sys.argv) > 1 else 0.2
LS = round(R - RIM, 4)  # slit half-length
PHASES = (90, 120)
STATES = [f"{a}{b}" for a, b in itertools.product(("OFF", "ON"), repeat=2)]  # x-diode state first


def box(cx, cy, lx, ly, z0, z1):
    return cq.Workplane("XY").workplane(offset=z0).center(cx, cy).rect(lx, ly).extrude(z1 - z0)


def patch():
    disc = cq.Workplane("XY").circle(R).extrude(T)
    return disc.cut(box(0, 0, G, 2 * LS, 0, T)).cut(box(0, 0, 2 * LS, G, 0, T))


def sheet(cx, cy, lx, ly):
    zs = T / 2
    pts = [cq.Vector(cx - lx / 2, cy - ly / 2, zs), cq.Vector(cx + lx / 2, cy - ly / 2, zs), cq.Vector(cx + lx / 2, cy + ly / 2, zs), cq.Vector(cx - lx / 2, cy + ly / 2, zs)]
    return cq.Workplane("XY").add(cq.Face.makeFromWires(cq.Wire.makePolygon(pts, close=True)))


def intended():
    xy = {0.0, P / 2, R, G / 2, LS, POS - WD / 2, POS + WD / 2, WD / 2}
    xy |= {-v for v in xy}
    z = {0.0, T, T / 2, -H, -H - T, LID}
    return dict(xy=np.array(sorted(round(v, 6) for v in xy)), z=np.array(sorted(round(v, 6) for v in z)), dir=np.array([-1.0, 0.0, 1.0]))


def build(state):
    sx, sy = state.split("_")
    pa = patch()
    sub, gnd = box(0, 0, P, P, -H, 0), box(0, 0, P, P, -H - T, -H)
    air = box(0, 0, P, P, 0, LID).cut(pa)
    asm = cq.Assembly()
    for nm, o, c in (("Substrate", sub, (0.25, 0.55, 0.25)), ("Ground", gnd, (0.6, 0.6, 0.65)), ("Patch", pa, (0.9, 0.7, 0.2)), ("Air", air, (0.7, 0.8, 1.0))):
        asm.add(o, name=nm, color=cq.Color(*c, 1))
    asm.add(sheet(0, POS, G, WD), name=f"Dx_{sx}", color=cq.Color(0.9, 0.1, 0.1, 1))
    asm.add(sheet(POS, 0, WD, G), name=f"Dy_{sy}", color=cq.Color(0.1, 0.1, 0.9, 1))
    path = OUT / f"cell_x{sx}_y{sy}.step"
    asm.export(str(path), exportType="STEP")
    raw = path.read_text(encoding="utf-8")
    snapped, stats = snap_step_text(raw, intended(), description=f"sidequest crossed rim cell x{sx} y{sy}, mm, 0.001 mm grid")
    path.write_text(snapped, encoding="utf-8", newline="\n")
    return path, raw, snapped, stats


def check(path, raw, snapped, stats):
    shp = cq.importers.importStep(str(path))
    solids = shp.solids().vals()
    faces = shp.faces().vals()
    area_cut = 2 * (2 * LS * G) - G * G
    vol_patch = (np.pi * R * R - area_cut) * T
    exp = sorted([(P, P, H, P * P * H), (P, P, T, P * P * T), (2 * R, 2 * R, T, vol_patch), (P, P, LID, P * P * LID - vol_patch)])
    got = sorted((round(s.BoundingBox().xlen, 6), round(s.BoundingBox().ylen, 6), round(s.BoundingBox().zlen, 6), s.Volume()) for s in solids)
    ok = len(got) == 4 and all(max(abs(g[i] - e[i]) for i in range(3)) < 1e-6 and abs(g[3] - e[3]) < 1e-5 * max(1.0, e[3]) for g, e in zip(got, exp, strict=True))
    sh = [f for f in faces if f.BoundingBox().zlen < 1e-9 and abs(f.BoundingBox().zmin - T / 2) < 1e-9]
    dims = sorted((round(f.BoundingBox().xlen, 6), round(f.BoundingBox().ylen, 6)) for f in sh)
    sheets_ok = dims == sorted([(G, WD), (WD, G)])
    bad = inexact_literals(path.read_text(encoding="utf-8"), intended())
    vs = verify_snap(raw, snapped)
    return {"file": path.name, "bodies": len(solids), "sheets": len(sh), "EXACT_DIMS": "PASS" if (ok and sheets_ok and not bad and len(sh) == 2) else "FAIL", "inexact_literals": len(bad),
            "n_snapped": stats["n_snapped"], "max_snap_mm": stats["max_snap_mm"], "collapsed_edges": vs.get("collapsed_edges"), "merge_violations": len(vs.get("distinct_merge_violations", []))}


def dispersion():
    """Uniform-sheet TM0: Y_xx(f) (x propagation) and Y_yy(f) (y propagation) of each state from the MoM (mesh 64); OFF/ON = transpose of ON/OFF."""
    fg = np.array([8.5, 9.5, 10.54, 11.5, 12.5, 13.5, 14.5]) * 1e9
    tens = {}
    for key, (sx, sy) in {"OFF_OFF": ("OFF", "OFF"), "ON_OFF": ("ON", "OFF"), "ON_ON": ("ON", "ON")}.items():
        tens[key] = [sheet_tensor(X.cross_cell(sx, sy, ls=LS), f, 64) for f in fg]
        print("solved", key, flush=True)
    yxx = {"OFF_OFF": [t[0, 0] for t in tens["OFF_OFF"]], "ON_OFF": [t[0, 0] for t in tens["ON_OFF"]], "OFF_ON": [t[1, 1] for t in tens["ON_OFF"]], "ON_ON": [t[0, 0] for t in tens["ON_ON"]]}
    yyy = {"OFF_OFF": [t[1, 1] for t in tens["OFF_OFF"]], "ON_OFF": [t[1, 1] for t in tens["ON_OFF"]], "OFF_ON": [t[0, 0] for t in tens["ON_OFF"]], "ON_ON": [t[1, 1] for t in tens["ON_ON"]]}
    out = {"freq_GHz": (fg / 1e9).tolist(), "Yxx_mS": {k: [[y.real * 1e3, y.imag * 1e3] for y in v] for k, v in yxx.items()}, "Yyy_mS": {k: [[y.real * 1e3, y.imag * 1e3] for y in v] for k, v in yyy.items()}, "predictions": []}
    for axis, table in (("x", yxx), ("y", yyy)):
        for ph in PHASES:
            row = {"axis": axis, "phase_deg": ph, "light_line_GHz": round(299792458.0 * np.radians(ph) / (P * 1e-3) / (2 * np.pi) / 1e9, 3)}
            for st_, ys in table.items():
                cs = CubicSpline(fg, np.array(ys).imag)

                def f_res(f):
                    return sm.find_beta_from_sheet_admittance(1j * float(cs(f)), gc.SUB, f).beta_rad_per_m * P * 1e-3 - np.radians(ph)

                try:
                    fr = brentq(f_res, fg[0], fg[-1], xtol=1e3)
                    row[st_] = f">{fg[-1] / 1e9:.1f}" if fr > fg[-1] - 5e7 else round(fr / 1e9, 3)
                except ValueError:
                    row[st_] = None
            out["predictions"].append(row)
    return out


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    checks = []
    for a, b in itertools.product(("OFF", "ON"), repeat=2):
        res = build(f"{a}_{b}")
        checks.append(check(*res))
        print(checks[-1], flush=True)
    d = dispersion()
    d["rim_mm"] = RIM
    d["label"] = "MODEL_ONLY: MoM Y(f) of the cell as a uniform sheet (mesh 64), TM0 transverse resonance; +-5 % in frequency. Keys: x-diode state_y-diode state. Disc is a true circle in STEP (the MoM uses a 48-gon, area -0.7 %)."
    d["exact_dims"] = checks
    (OUT / "predictions.json").write_text(json.dumps(d, indent=1, default=float), encoding="utf-8", newline="\n")
    for r in d["predictions"]:
        print(r)
