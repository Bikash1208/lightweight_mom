"""Crossed two-diode cell: split disc (diameter 6.0 mm) cut by a vertical slit (x diode bridges it) and a horizontal slit (y diode bridges it), both 0/90 deg, measured MA4AGP907.
Topologies: T1 full-length slits (4 islands, x diode at y = +ydiode, y diode at x = +xdiode, one island floating), T2 slits of half-length `ls` (rim connects the islands).
Sheet tensor from the MoM (normal incidence), 4 states OFF/OFF, ON/OFF, OFF/ON, ON/ON. MODEL_ONLY."""
import itertools
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np

import sidequest_coded_cell as cc
import sidequest_gate_cell as gc
from src.composite_geometry import AxisRect
from src.geometry import UnitCell
from src.hologram.coded_cell_eval import sheet_tensor
from src.materials import Metal
from src.periodic_bc import RectangularLattice
from src.polygon_geometry import PolyComposite, regular_polygon
from src.ris import diodes_measured as dm

MM = 1e-3


def cross_cell(sx, sy, diam=6.0, g=0.4, wd=1.2, ls=None, pos=1.6, four=False):
    """sx, sy in {OFF, ON}. ls: half-length of the slits (None = through the disc). Diodes at (0, +pos) [x] and (+pos, 0) [y]; four=True adds the mirror diodes (0, -pos), (-pos, 0)."""
    r = diam / 2
    half = (r + 0.1) if ls is None else ls
    disc = regular_polygon(0.0, 0.0, r * MM, 48, 0.0)
    ps = (pos, -pos) if four else (pos,)
    cuts = sorted([(-half, 1e9)] + [(c - wd / 2, c + wd / 2) for c in ps], key=lambda t: t[0])  # footprint intervals along the slit
    segs, a = [], -half
    for lo, hi in sorted((c - wd / 2, c + wd / 2) for c in ps):
        if lo > a:
            segs.append((a, lo))
        a = hi
    if a < half:
        segs.append((a, half))
    ap = []
    for a0, a1 in segs:
        ap.append(AxisRect(0.0, (a0 + a1) / 2 * MM, g * MM, (a1 - a0) * MM))  # vertical slit pieces
        ap.append(AxisRect((a0 + a1) / 2 * MM, 0.0, (a1 - a0) * MM, g * MM))  # horizontal slit pieces
    lds = tuple(dm.load(gc.STATES[sx], 0.0, c * MM, g * MM, wd * MM, "x") for c in ps) + tuple(dm.load(gc.STATES[sy], c * MM, 0.0, g * MM, wd * MM, "y") for c in ps)
    return UnitCell(RectangularLattice(gc.P, gc.P), PolyComposite(metal=(disc,), apertures=tuple(ap)), gc.SUB, Metal(), backing="pec", loads=lds)


def job(p):
    warnings.simplefilter("ignore")
    ls, pos, mesh, four = p
    out = {"ls": ls, "pos": pos, "mesh": mesh, "four": four, "states": {}}
    for sx, sy in itertools.product(("OFF", "ON"), repeat=2):
        y = sheet_tensor(cross_cell(sx, sy, ls=ls, pos=pos, four=four), gc.F0, mesh)
        out["states"][f"{sx}/{sy}"] = [[complex(v).real * 1e3, complex(v).imag * 1e3] for v in y.ravel()]
    return out


def analyse(o):
    Y = {k: np.array([complex(a, b) for a, b in v]).reshape(2, 2) * 1e-3 for k, v in o["states"].items()}
    dx = 0.5 * ((Y["ON/OFF"][0, 0] - Y["OFF/OFF"][0, 0]) + (Y["ON/ON"][0, 0] - Y["OFF/ON"][0, 0]))  # mean over the other diode's state
    dy = 0.5 * ((Y["OFF/ON"][1, 1] - Y["OFF/OFF"][1, 1]) + (Y["ON/ON"][1, 1] - Y["ON/OFF"][1, 1]))
    iso = np.mean([0.5 * (y[0, 0] + y[1, 1]) for y in Y.values()])
    res = {"dXX_mS": [dx.real * 1e3, dx.imag * 1e3], "dYY_mS": [dy.real * 1e3, dy.imag * 1e3]}
    for nm, d in (("x", dx), ("y", dy)):
        try:
            bk, ab, aa = cc.alpha_pair(iso - d / 2, iso + d / 2)
        except Exception:  # noqa: BLE001
            bk, ab, aa = None, None, None
        res[f"beta_k0_{nm}"], res[f"alpha_beam_{nm}"], res[f"alpha_abs_{nm}"] = bk, ab, aa
    res["Y_iso_mS"] = [iso.real * 1e3, iso.imag * 1e3]
    return res


if __name__ == "__main__":
    grid = [(None, 1.6, 64, False), (None, 1.6, 32, False), (None, 1.6, 32, True), (None, 1.6, 64, True)]
    with ProcessPoolExecutor(4) as ex:
        res = list(ex.map(job, grid))
    for o in res:
        o["analysis"] = analyse(o)
        a = o["analysis"]
        print(o["four"], o["mesh"], "dXX %.2f%+.2fj dYY %.2f%+.2fj" % (*a["dXX_mS"], *a["dYY_mS"]), "alpha x/y", a["alpha_beam_x"], a["alpha_beam_y"], "beta", a["beta_k0_x"], flush=True)
        print("   ", {k: np.round(complex(v[0][0], v[0][1]), 2) for k, v in o["states"].items()})
    Path("results/sidequest_hybrid/coded_cross_scan.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline="\n")
