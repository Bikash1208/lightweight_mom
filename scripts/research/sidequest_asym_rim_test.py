"""Hypothesis test: HFSS shows diodes that switch only together (C/C 10.670, R/C 10.665, C/R 10.690, R/R 11.006 GHz at 90 deg, x run) and the user saw the slit-end rims missing at the TOP and RIGHT ends of the HFSS model.
MoM with the rim removed at the top end of the vertical slit and at the right end of the horizontal slit (rims at the bottom and left kept), diodes as exported (x diode at y = +0.8, y diode at x = +0.8, footprint 1.0):
eigenfrequency (lid-aware transverse resonance with the loaded Bloch tensor) for the four R/C states; OFF = 28.16 fF load or OFF = open (calibrated). Output results/sidequest_hybrid/asym_rim_test.json"""
import itertools
import json
import sys
import warnings
from multiprocessing import Pool
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

import sidequest_cell_search as C
import sidequest_gate_cell as gc
import sidequest_rt5880_cell_package as K
from src.composite_geometry import AxisRect
from src.geometry import UnitCell
from src.materials import Metal
from src.periodic_bc import RectangularLattice
from src.polygon_geometry import PolyComposite, regular_polygon
from src.ris import diodes_measured as dm

MM = 1e-3
G, POS, WD = 0.4, 0.8, 1.0
FG = np.array([9.6, 10.0, 10.4, 10.8, 11.2, 11.6]) * 1e9


def build(sx, sy, vx=(-2.8, 3.1), hy=(-2.8, 3.1), off="C"):
    disc = regular_polygon(0.0, 0.0, 3.0 * MM, 48, 0.0)
    ap, loads = [], []
    for axis, (lo, hi), state in (("x", vx, sx), ("y", hy, sy)):
        win = (POS - WD / 2, POS + WD / 2)
        segs = [(lo, win[0]), (win[1], hi)]
        if state == "OFF" and off == "open":
            segs = [(lo, hi)]  # open footprint: nothing in the slit (calibrated OFF)
        for a, b in segs:
            ap.append(AxisRect(0.0, (a + b) / 2 * MM, G * MM, (b - a) * MM) if axis == "x" else AxisRect((a + b) / 2 * MM, 0.0, (b - a) * MM, G * MM))
        if not (state == "OFF" and off == "open"):
            c = POS
            loads.append(dm.load(gc.STATES[state], 0.0, c * MM, G * MM, WD * MM, "x") if axis == "x" else dm.load(gc.STATES[state], c * MM, 0.0, WD * MM, G * MM, "y"))
    return UnitCell(RectangularLattice(6.4e-3, 6.4e-3), PolyComposite(metal=(disc,), apertures=tuple(ap)), gc.SUB, Metal(), backing="pec", loads=tuple(loads))


def job(a):
    warnings.simplefilter("ignore")
    geom, off, sx, sy, f = a
    vx, hy = ((-2.8, 3.1), (-2.8, 3.1)) if geom == "asym" else ((-2.8, 2.8), (-2.8, 2.8))
    return complex(C.bloch_loaded(build(sx, sy, vx, hy, off), f, 64, K.BETA, 0.0)[0, 0])


if __name__ == "__main__":
    states = list(itertools.product(("OFF", "ON"), repeat=2))
    jobs = [(g, off, sx, sy, f) for g in ("asym",) for off in ("C", "open") for sx, sy in states for f in FG]
    with Pool(4) as pool:
        res = pool.map(job, jobs)
    out, it = {}, iter(res)
    print("HFSS: C/C 10.670  R/C 10.665  C/R 10.690  R/R 11.006")
    for g in ("asym",):
        for off in ("C", "open"):
            row = {}
            for sx, sy in states:
                ys = np.array([next(it) for _ in FG])
                cs = CubicSpline(FG, ys.imag)
                try:
                    fr = brentq(lambda f: K.stack_b(K.BETA, f) + float(cs(f)), FG[0], FG[-1], xtol=1e3) / 1e9
                except ValueError:
                    fr = None
                row[f"{sx}/{sy}"] = fr
            out[f"{g}_OFF={off}"] = row
            print(g, "OFF=" + off, {k: (None if v is None else round(v, 3)) for k, v in row.items()}, flush=True)
    Path("results/sidequest_hybrid/asym_rim_test.json").write_text(json.dumps(out, indent=1), encoding="utf-8", newline="\n")
