"""Redesign candidates for independent per-axis switching (x diodes on one line, y diodes on another): both halves of each slit bridged (2 diodes per axis, 4 per cell), with the rim (C1) and without (C4: four islands, diodes the only bridges).
Per state (x diodes OFF/ON, y diodes OFF/ON): x-run eigenfrequency at 90 deg (lid-aware transverse resonance with the loaded Bloch tensor, mesh 64), Y tensor at f_res; y run by x<->y symmetry. OFF = 28.16 fF load or OFF = open (HFSS: OFF ~ open).
Footprint 0.4 x 1.0 at 0.3-1.3 mm from the slit crossing on both sides (sheets keep >= 0.14 mm clearance).  Output results/sidequest_hybrid/independent_axes.json"""
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

import sidequest_cell_independence as CI
import sidequest_cell_search as C
import sidequest_rt5880_cell_package as K

FG = np.array([9.6, 10.0, 10.4, 10.8, 11.2, 11.6]) * 1e9
TOP, BOT = (0.3, 1.3), (-1.3, -0.3)


def cell(name, sx, sy, off):
    ls = 2.8 if name == "C1_rim_4diodes" else 3.1
    def wins(s):
        if s == "OFF" and off == "open":
            return [("OPEN", *TOP), ("OPEN", *BOT)]
        return [(s, *TOP), (s, *BOT)]
    return CI.build(ls, wins(sx), wins(sy))


def job(a):
    warnings.simplefilter("ignore")
    name, off, sx, sy, f = a
    return C.bloch_loaded(cell(name, sx, sy, off), f, 64, K.BETA, 0.0)


if __name__ == "__main__":
    names = ("C1_rim_4diodes", "C4_norim_4diodes")
    states = list(itertools.product(("OFF", "ON"), repeat=2))
    jobs = [(n, off, sx, sy, f) for n in names for off in ("open", "C") for sx, sy in states for f in FG]
    with Pool(4) as pool:
        res = pool.map(job, jobs)
    it = iter(res)
    out = {"HFSS_reference_one_cell_two_diodes_per_cell": {"C/C": 10.670, "R/C": 10.665, "C/R": 10.690, "R/R": 11.006}, "rows": []}
    for n in names:
        for off in ("open", "C"):
            tens = {}
            for sx, sy in states:
                tens[(sx, sy)] = np.array([next(it) for _ in FG])
            fr, yr = {}, {}
            for (sx, sy), t in tens.items():
                cs = CubicSpline(FG, t[:, 0, 0].imag)
                try:
                    f0 = brentq(lambda f: K.stack_b(K.BETA, f) + float(cs(f)), FG[0], FG[-1], xtol=1e3)
                except ValueError:
                    f0 = None
                fr[(sx, sy)] = f0
                if f0:
                    yr[(sx, sy)] = {k: [complex(CubicSpline(FG, t[:, i, j])(f0)).real * 1e3, complex(CubicSpline(FG, t[:, i, j])(f0)).imag * 1e3] for k, (i, j) in {"xx": (0, 0), "yy": (1, 1), "xy": (0, 1)}.items()}
            # per-axis dY at the common frequency 11.0 GHz (grid point), cross coupling
            def y11(sx, sy, i, j):
                return complex(CubicSpline(FG, tens[(sx, sy)][:, i, j])(10.8e9)).imag * 1e3
            base = y11("OFF", "OFF", 0, 0)
            row = {"candidate": n, "OFF": off, "f_x_run_GHz": {f"{a}/{b}": (None if v is None else round(v / 1e9, 4)) for (a, b), v in fr.items()},
                   "dYxx_at_10p8GHz_mS": {"x_ON_only": y11("ON", "OFF", 0, 0) - base, "y_ON_only": y11("OFF", "ON", 0, 0) - base, "both_ON": y11("ON", "ON", 0, 0) - base},
                   "dYyy_at_10p8GHz_mS": {"x_ON_only": y11("ON", "OFF", 1, 1) - y11("OFF", "OFF", 1, 1), "y_ON_only": y11("OFF", "ON", 1, 1) - y11("OFF", "OFF", 1, 1), "both_ON": y11("ON", "ON", 1, 1) - y11("OFF", "OFF", 1, 1)},
                   "Y_at_f_res": {f"{a}/{b}": v for (a, b), v in yr.items()}}
            out["rows"].append(row)
            print(n, off, row["f_x_run_GHz"], "dYxx", {k: round(v, 2) for k, v in row["dYxx_at_10p8GHz_mS"].items()}, "dYyy", {k: round(v, 2) for k, v in row["dYyy_at_10p8GHz_mS"].items()}, flush=True)
    Path("results/sidequest_hybrid/independent_axes.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8", newline="\n")
