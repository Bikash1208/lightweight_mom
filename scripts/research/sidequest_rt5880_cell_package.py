"""HFSS package, RT5880 chosen cell (rim disc, diodes 0.8 mm from the centre): ONE geometry STEP (METAL_PEC, SUBSTRATE, AIR, DX_SHEET, DY_SHEET) + state table (states differ only in RLC values) + predictions.
Geometry = the MoM model: disc r 3.000, slits 0.4 wide, 5.6 long (rim 0.200 at all four ends), x diode across the vertical slit at y = +0.8, y diode across the horizontal slit at x = +0.8, footprint 0.4 x 1.0 mm
(inner edge 0.3 mm from the slit crossing: the two sheets are 0.14 mm apart; footprint 1.2 mm would make their corners touch). Checks run before export and again on the re-imported STEP.
    python sidequest_rt5880_cell_package.py   -> hfss_exports/sidequest_rt5880_cell/"""
import itertools
import json
import sys
import warnings
from multiprocessing import Pool
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import cadquery as cq
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

import check_step_cell as CK
import sidequest_cell_search as C
import sidequest_gate_cell as gc
from src.dispersion import sheet_model as sm
from src.transmitarray.export import inexact_literals, snap_step_text, verify_snap

OUT = ROOT / "hfss_exports/sidequest_rt5880_cell"
P, H, T, GT, LID = 6.4, 1.575, 0.035, 0.035, 30.0
R, G, WD, POS, LS = 3.0, 0.4, 1.0, 0.8, 2.8
C0 = 299792458.0
EPS0 = 8.8541878128e-12
BETA = np.radians(90.0) / (P * 1e-3)
FG = np.array([9.6, 10.0, 10.4, 10.8, 11.2, 11.6]) * 1e9


def box(cx, cy, lx, ly, z0, z1):
    return cq.Workplane("XY").workplane(offset=z0).center(cx, cy).rect(lx, ly).extrude(z1 - z0)


def sheet(cx, cy, lx, ly):
    z = T / 2
    pts = [cq.Vector(cx - lx / 2, cy - ly / 2, z), cq.Vector(cx + lx / 2, cy - ly / 2, z), cq.Vector(cx + lx / 2, cy + ly / 2, z), cq.Vector(cx - lx / 2, cy + ly / 2, z)]
    return cq.Face.makeFromWires(cq.Wire.makePolygon(pts, close=True))


def build(ls=LS, pos=POS, wd=WD):
    patch = cq.Workplane("XY").circle(R).extrude(T).cut(box(0, 0, G, 2 * ls, 0, T)).cut(box(0, 0, 2 * ls, G, 0, T))
    ground = box(0, 0, P, P, -H - GT, -H)
    metal = cq.Compound.makeCompound([patch.val(), ground.val()])
    sub = box(0, 0, P, P, -H, 0)
    air = box(0, 0, P, P, 0, LID).cut(patch)
    return {"METAL_PEC": metal, "SUBSTRATE": sub.val(), "AIR": air.val(), "DX_SHEET": sheet(0, pos, G, wd), "DY_SHEET": sheet(pos, 0, wd, G)}


def intended(ls=LS, pos=POS, wd=WD):
    xy = {0.0, P / 2, R, G / 2, ls, pos - wd / 2, pos + wd / 2, wd / 2}
    xy |= {-v for v in xy}
    z = {0.0, T, T / 2, -H, -H - GT, LID}
    return dict(xy=np.array(sorted(round(v, 6) for v in xy)), z=np.array(sorted(round(v, 6) for v in z)), dir=np.array([-1.0, 0.0, 1.0]))


def export_and_check():
    OUT.mkdir(parents=True, exist_ok=True)
    parts = build()
    comp = cq.Compound.makeCompound([parts["METAL_PEC"], parts["SUBSTRATE"], parts["AIR"], parts["DX_SHEET"], parts["DY_SHEET"]])
    pre = CK.check(comp, pos=POS, ls=LS, g=G, wd=WD)
    asm = cq.Assembly()
    cols = {"METAL_PEC": (0.9, 0.7, 0.2), "SUBSTRATE": (0.25, 0.55, 0.25), "AIR": (0.7, 0.8, 1.0), "DX_SHEET": (0.9, 0.1, 0.1), "DY_SHEET": (0.1, 0.1, 0.9)}
    for nm, shp in parts.items():
        asm.add(shp, name=nm, color=cq.Color(*cols[nm], 1))
    path = OUT / "cell_geometry.step"
    asm.export(str(path), exportType="STEP")
    raw = path.read_text(encoding="utf-8")
    snapped, stats = snap_step_text(raw, intended(), description="sidequest RT5880 coded cell geometry, mm, 0.001 mm grid")
    path.write_text(snapped, encoding="utf-8", newline="\n")
    post = CK.check_file(path, pos=POS, ls=LS, g=G, wd=WD)
    text = path.read_text(encoding="utf-8")
    bad = inexact_literals(text, intended())
    vs = verify_snap(raw, snapped)
    re = cq.importers.importStep(str(path))
    bb = re.val().BoundingBox() if hasattr(re.val(), "BoundingBox") else None
    dims = {"bbox_x": [round(bb.xmin, 6), round(bb.xmax, 6)], "bbox_z": [round(bb.zmin, 6), round(bb.zmax, 6)]}
    exact = {"ok": (len(bad) == 0 and vs.get("collapsed_edges", 1) == 0 and not vs.get("distinct_merge_violations")), "inexact_literals": len(bad), "n_snapped": stats["n_snapped"], "max_snap_mm": stats["max_snap_mm"], **dims}
    pre["checks"]["EXACT_DIMS"] = exact
    post["checks"]["EXACT_DIMS"] = exact
    return pre, post


# ---------------------------------------------------------------- predictions
def y_air_lid(beta, om, k0, L=LID * 1e-3):
    kz = np.sqrt(k0**2 - beta**2 + 0j)
    if kz.imag > 0:
        kz = -kz
    return -1j * om * EPS0 / (kz * np.tan(kz * L))


def stack_b(beta, f, lid=True):
    om, k0 = 2 * np.pi * f, 2 * np.pi * f / C0
    ya = y_air_lid(beta, om, k0) if lid else sm._y_air(beta, om, k0)
    return (ya + sm._y_slab(beta, om, gc.SUB.eps_complex, k0, gc.SUB.thickness)).imag


def f_bare(lid=True):
    return brentq(lambda f: stack_b(BETA, f, lid), 9e9, 11.7e9) / 1e9


def _job(a):
    import warnings as w

    w.simplefilter("ignore")
    ls, sx, sy, f = a
    c = C.cell("rt5880", sx, sy, diam=6.0, ls=ls, pos=POS, wd=WD)
    return C.bloch_loaded(c, f, 64, BETA, 0.0)


def resonance(yxx_grid, lid=True):
    cs = CubicSpline(FG, np.array(yxx_grid).imag)
    try:
        fr = brentq(lambda f: stack_b(BETA, f, lid) + float(cs(f)), FG[0], FG[-1], xtol=1e3)
    except ValueError:
        return None
    return fr


def predictions():
    jobs = []
    variants = {"A_rim": LS, "B_norim": 3.1}
    states = list(itertools.product(("OFF", "ON"), repeat=2))
    for var, ls in variants.items():
        for sx, sy in states:
            for f in FG:
                jobs.append((ls, sx, sy, f))
    for f in FG:
        jobs.append((LS, None, None, f))  # slits open, no diode at all
    with Pool(4) as pool:
        res = pool.map(_job, jobs)
    it = iter(res)
    tens = {}
    for var in variants:
        for sx, sy in states:
            tens[(var, sx, sy)] = np.array([next(it) for _ in FG])
    tens[("open", None, None)] = np.array([next(it) for _ in FG])
    out = {"bare_substrate_GHz_with_lid": f_bare(True), "bare_substrate_GHz_open_air": f_bare(False), "beta_rad_per_m": BETA, "rows": [], "label": "MODEL_ONLY (loaded Bloch tensor, mesh 64, transverse resonance with the 30 mm PEC lid, +-0.15 GHz)"}
    for var in variants:
        for sx, sy in states:
            for axis in ("x", "y"):
                if axis == "x":
                    yt = tens[(var, sx, sy)]
                else:
                    yt = tens[(var, sy, sx)][:, ::-1, ::-1]  # y propagation of (sx, sy) = x propagation of the mirrored state (xx <-> yy): the propagation-axis entry is its Y_xx
                fr = resonance(yt[:, 1, 1] if axis == "y" else yt[:, 0, 0])
                row = {"variant": var, "state": f"{sx}/{sy}", "axis": axis, "f_GHz": None if fr is None else round(fr / 1e9, 4)}
                if fr is not None:
                    cs = {k: CubicSpline(FG, yt[:, k[0], k[1]]) for k in ((0, 0), (1, 1), (0, 1))}
                    k0 = 2 * np.pi * fr / C0
                    row.update({"beta_over_k0": round(BETA / k0, 4), "Y_xx_mS": [complex(cs[(0, 0)](fr)).real * 1e3, complex(cs[(0, 0)](fr)).imag * 1e3], "Y_yy_mS": [complex(cs[(1, 1)](fr)).real * 1e3, complex(cs[(1, 1)](fr)).imag * 1e3],
                                "Y_xy_mS": [complex(cs[(0, 1)](fr)).real * 1e3, complex(cs[(0, 1)](fr)).imag * 1e3]})
                out["rows"].append(row)
    fo = resonance(tens[("open", None, None)][:, 0, 0])
    out["open_slits_no_diode_f_GHz"] = None if fo is None else round(fo / 1e9, 4)
    out["grid"] = {f"{k[0]}_{k[1]}_{k[2]}": [[[v.real * 1e3, v.imag * 1e3] for v in row] for row in t.reshape(len(FG), 4)] for k, t in tens.items()}
    return out


if __name__ == "__main__":
    pre, post = export_and_check()
    print("PRE ", {k: v["ok"] for k, v in pre["checks"].items()}, pre["PASS"], flush=True)
    print("POST", {k: v["ok"] for k, v in post["checks"].items()}, post["PASS"], flush=True)
    (OUT / "geometry_checks.json").write_text(json.dumps({"before_export": pre, "re_imported_step": post}, indent=1, default=str), encoding="utf-8", newline="\n")
    if not (pre["PASS"] and post["PASS"]):
        sys.exit("geometry checks failed: no predictions, no package")
    with open(OUT / "states.csv", "w", encoding="utf-8", newline="\n") as fh:
        fh.write("state,DX_SHEET_element,DX_value,DY_SHEET_element,DY_value,note\n")
        for sx, sy in itertools.product(("OFF", "ON"), repeat=2):
            el = lambda s: ("C", "28.16 fF") if s == "OFF" else ("R", "2.44 ohm")  # noqa: E731
            fh.write(f"x{sx}_y{sy},{el(sx)[0]},{el(sx)[1]},{el(sy)[0]},{el(sy)[1]},R only = ON / C only = OFF (MA4AGP907 measured 10.54 GHz); current line across the 0.4 mm gap\n")
    d = predictions()
    (OUT / "predictions.json").write_text(json.dumps(d, indent=1, default=float), encoding="utf-8", newline="\n")
    print("bare substrate (lid):", d["bare_substrate_GHz_with_lid"], "open:", d["open_slits_no_diode_f_GHz"])
    for r in d["rows"]:
        print(r["variant"], r["state"], r["axis"], r["f_GHz"], r.get("beta_over_k0"))
