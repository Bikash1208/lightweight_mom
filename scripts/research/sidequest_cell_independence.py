"""Why did the MoM predict a single diode to switch (x diode alone +0.46 GHz) when HFSS shows none (R/C 10.665, C/R 10.690, C/C 10.670, R/R 11.006 GHz at 90 deg)? Experiments on the exported geometry
(RT5880 disc, slits 0.4 x 5.6, rim 0.2, MoM mesh 64, loaded Bloch tensor at kt = beta = 245.4 rad/m, f = 10.8 GHz): Y_xx for ideal elements (open / metal strap / diode loads) in each slit,
strap length scan, both halves of a slit bridged, other-axis (cross) coupling, mesh check.  Output results/sidequest_hybrid/cell_independence.json
Element codes per slit: list of (state, y_or_x_start, length) windows in mm measured from the slit crossing outwards on the +side (and mirrored for 'both'); state OFF/ON = MA4AGP907 loads, SHORT = metal strap, OPEN = nothing."""
import json
import sys
import warnings
from multiprocessing import Pool
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np

import sidequest_cell_search as C
import sidequest_gate_cell as gc
from src.composite_geometry import AxisRect
from src.geometry import UnitCell
from src.materials import Metal
from src.periodic_bc import RectangularLattice
from src.polygon_geometry import PolyComposite, regular_polygon
from src.ris import diodes_measured as dm

MM = 1e-3
G = 0.4
BETA = np.radians(90.0) / 6.4e-3
F = 10.8e9


def build(ls, xwin, ywin):
    """xwin: windows on the vertical slit [(state, y0, y1), ...] (y in mm, any sign); ywin: windows on the horizontal slit [(state, x0, x1), ...]."""
    disc = regular_polygon(0.0, 0.0, 3.0 * MM, 48, 0.0)
    ap, loads = [], []
    for axis, wins in (("x", xwin), ("y", ywin)):
        cuts = sorted([(a, b) for (s, a, b) in wins if s != "OPEN"])
        segs, cur = [], -ls
        for a, b in cuts:
            if a > cur:
                segs.append((cur, a))
            cur = max(cur, b)
        if cur < ls:
            segs.append((cur, ls))
        for a, b in segs:
            ap.append(AxisRect(0.0, (a + b) / 2 * MM, G * MM, (b - a) * MM) if axis == "x" else AxisRect((a + b) / 2 * MM, 0.0, (b - a) * MM, G * MM))
        for s, a, b in wins:
            if s in ("OFF", "ON"):
                c, w = (a + b) / 2, b - a
                loads.append(dm.load(gc.STATES[s], 0.0, c * MM, G * MM, w * MM, "x") if axis == "x" else dm.load(gc.STATES[s], c * MM, 0.0, w * MM, G * MM, "y"))
    return UnitCell(RectangularLattice(6.4e-3, 6.4e-3), PolyComposite(metal=(disc,), apertures=tuple(ap)), gc.SUB, Metal(), backing="pec", loads=tuple(loads))


def job(a):
    warnings.simplefilter("ignore")
    name, ls, xwin, ywin, mesh = a
    y = C.bloch_loaded(build(ls, xwin, ywin), F, mesh, BETA, 0.0)
    return name, [[complex(v).real * 1e3, complex(v).imag * 1e3] for v in y.ravel()]


def cases():
    L = 2.8
    top = lambda s, a=0.3, b=1.3: (s, a, b)  # noqa: E731
    bot = lambda s, a=0.3, b=1.3: (s, -b, -a)  # noqa: E731
    cs = []
    cs.append(("open/open", L, [], [], 64))
    for s in ("SHORT", "ON", "OFF"):
        cs.append((f"x_top_{s}", L, [top(s)], [], 64))
        cs.append((f"y_right_{s}", L, [], [top(s)], 64))
        cs.append((f"x_top_{s}+y_right_{s}", L, [top(s)], [top(s)], 64))
        cs.append((f"x_both_{s}", L, [top(s), bot(s)], [], 64))
        cs.append((f"x_both_{s}+y_both_{s}", L, [top(s), bot(s)], [top(s), bot(s)], 64))
    for w in (0.5, 1.5, 2.0, 2.5):
        cs.append((f"x_top_SHORT_len{w}", L, [("SHORT", 0.3, 0.3 + w)], [], 64))
    cs.append(("x_full_SHORT_0.2..2.8", L, [("SHORT", 0.2, 2.8), ("SHORT", -2.8, -0.2)], [], 64))
    for mesh in (32, 64, 96):
        cs.append((f"mesh{mesh}_open", L, [], [], mesh))
        cs.append((f"mesh{mesh}_x_top_SHORT_0.4..1.2", L, [("SHORT", 0.4, 1.2)], [], mesh))
    cs.append(("norim_open", 3.1, [], [], 64))
    cs.append(("norim_x_top_SHORT", 3.1, [top("SHORT")], [], 64))
    cs.append(("norim_x_both_SHORT", 3.1, [top("SHORT"), bot("SHORT")], [], 64))
    cs.append(("norim_x_both+y_both_SHORT", 3.1, [top("SHORT"), bot("SHORT")], [top("SHORT"), bot("SHORT")], 64))
    return cs


if __name__ == "__main__":
    cs = cases()
    with Pool(4) as pool:
        res = dict(pool.map(job, cs))
    base = res["open/open"][0][1]
    out = {"f_GHz": F / 1e9, "Y": res}
    print("Y_xx (mS imag) open/open:", round(base, 3))
    for k, v in res.items():
        print(f"{k:34s} Y_xx {v[0][1]:8.3f}  dY_xx {v[0][1] - res['open/open'][0][1] if not k.startswith(('mesh','norim')) else float('nan'):8.3f}  Y_yy {v[3][1]:8.3f}  Y_xy {v[1][1]:7.3f}", flush=True)
    Path("results/sidequest_hybrid/cell_independence.json").write_text(json.dumps(out, indent=1), encoding="utf-8", newline="\n")
