"""Leaf cell with the exact Table I of Zhao, Pan, Shen, IEEE TAP 74(1) 2026: a 13.3, g 0.5, l 0.6, r1 6.4, r2 6.0; two lens leaves side by side along x (long axis y), pads l x l, diode bridge g x l on the x axis.
Lens = overlap of two equal discs (left r1, right r2) of width W along x (W = (a - g - s)/2, s = intercell clearance, not in Table I: DEFAULT 0.8, assumption). Y_sheet -> X_s (open and PEC lid lambda/4). MODEL_ONLY.
    python sidequest_leaf_cell2.py <tag> '<json list of cases>'   keys: f er tand h a g l r1 r2 s mesh state diode"""
import json, math, sys, warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent)); warnings.simplefilter("ignore")
import numpy as np
import sidequest_leaf_cell as L1
from src.composite_geometry import AxisRect
from src.geometry import UnitCell
from src.hologram.coded_cell_eval import sheet_tensor
from src.loads import LumpedLoad
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice
from src.polygon_geometry import Poly, PolyComposite
MM = 1e-3


def leaf(xin, W, r, sgn):
    d = 2 * r - W; xm = xin + sgn * W / 2
    yt = math.sqrt(r * r - d * d / 4); t = math.asin(yt / r)
    a1 = [(xm + sgn * (d / 2 - r * math.cos(a)) * -1 if False else xm + sgn * (r * math.cos(a) - d / 2) * 1, r * math.sin(a)) for a in np.linspace(-t, t, 41)]
    return a1, d, yt


def lens_poly(xin, W, r, sgn):
    # outer arc: disc centred at xm + sgn*(-d/2)... build both arcs explicitly; outer boundary away from the gap, inner boundary towards the gap
    d = 2 * r - W; xm = xin + sgn * W / 2; yt = math.sqrt(r * r - d * d / 4); t = math.asin(yt / r)
    th = np.linspace(-t, t, 41)
    outer = [(xm + sgn * (d / 2 - r * math.cos(a)) * -1 * -1 if False else xm - sgn * (d / 2) + sgn * r * math.cos(a) * 0 + 0, 0) for a in []]
    # disc A centred at (xm - sgn*d/2): its far side (sgn direction) boundary is the outer arc;  disc B centred at (xm + sgn*d/2): its near side (-sgn direction) is the inner arc
    outer = [(xm - sgn * d / 2 + sgn * r * math.cos(a), r * math.sin(a)) for a in th]
    inner = [(xm + sgn * d / 2 - sgn * r * math.cos(a), r * math.sin(a)) for a in th[::-1]]
    pts = outer + inner
    return Poly(tuple((x * MM, y * MM) for x, y in pts)), yt


def cell(c):
    a, g, l, r1, r2 = c["a"], c["g"], c["l"], c["r1"], c["r2"]
    W = (a - g - c["s"]) / 2
    pl, yt1 = lens_poly(-g / 2, W, r1, -1); pr, yt2 = lens_poly(+g / 2, W, r2, +1)
    metal = [pl, pr, AxisRect(-(g / 2 + l / 2) * MM, 0, l * MM, l * MM), AxisRect((g / 2 + l / 2) * MM, 0, l * MM, l * MM), AxisRect(0, 0, g * MM, l * MM)]
    f = c["f"]
    if c["diode"] == "smp":
        zf = (lambda ff: complex(2.0, 2 * math.pi * ff * 0.45e-9)) if c["state"] == "ON" else (lambda ff: complex(0.0, 2 * math.pi * ff * 0.45e-9 - 1 / (2 * math.pi * ff * 0.23e-12)))
    else:
        zf = lambda ff: complex(2.44, 0.0)  # noqa: E731
    ld = LumpedLoad(0.0, 0.0, g * MM, l * MM, "x", z_of_f=zf)
    sub = Substrate("s", c["er"], c["tand"], c["h"] * MM)
    return UnitCell(RectangularLattice(a * MM, a * MM), PolyComposite(metal=tuple(metal)), sub, Metal(), backing="pec", loads=(ld,)), sub, (W, yt1, yt2)


def job(c):
    warnings.simplefilter("ignore")
    uc, sub, dims = cell(c)
    y = sheet_tensor(uc, c["f"], c["mesh"]); H = 3e8 / c["f"] / 4
    rx = L1.xs_from_y(complex(y[0, 0]), sub, c["f"], H); ry = L1.xs_from_y(complex(y[1, 1]), sub, c["f"], H)
    return {**c, "W": dims[0], "tips_mm": [2 * dims[1], 2 * dims[2]], "Y_mS": [[v.real * 1e3, v.imag * 1e3] for v in y.ravel()], "x": rx, "y": ry}


if __name__ == "__main__":
    cases = json.loads(sys.argv[2])
    with ProcessPoolExecutor(min(len(cases), 32)) as ex:
        res = list(ex.map(job, cases))
    for r in res:
        print(r["state"], "mesh", r["mesh"], "W %.2f tips %.2f/%.2f" % (r["W"], *r["tips_mm"]), "Yxx %.2f%+.2fj Yyy %.2f%+.2fj" % (*r["Y_mS"][0], *r["Y_mS"][3]), "Xx_lid %.1f Xy_lid %.1f ohm (open %.1f/%.1f) alpha_x %.2f" % (r["x"]["X_lid_ohm"], r["y"]["X_lid_ohm"], r["x"]["X_open_ohm"], r["y"]["X_open_ohm"], r["x"]["alpha_np_per_m"]), flush=True)
    Path(f"results/sidequest_leaf/{sys.argv[1]}.json").write_text(json.dumps(res, indent=1), encoding="utf-8", newline="\n")
