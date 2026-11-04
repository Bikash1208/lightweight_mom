"""Rotated 1-diode cell for the polarisation-switched coded hologram (uses the Phase 28 rotated-load fix, `loads.default_transverse_conducts`).
Split patch a x h with the bridge/diode footprint (gap g, width wd) turned by chi; the patch must fit the 6.4 mm cell at any chi, so it is smaller than the
x-only pair of CODED_HOLOGRAM.md. Sheet tensor from the MoM (normal incidence, mesh 32), measured MA4AGP907 ON 10 mA / OFF 0 V. MODEL_ONLY."""
import json
import math
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np

import sidequest_gate_cell as gc
from src.geometry import UnitCell
from src.hologram.coded_cell_eval import sheet_tensor
from src.materials import Metal
from src.periodic_bc import RectangularLattice
from src.polygon_geometry import Poly, PolyComposite, oriented_rect
from src.ris import diodes_measured as dm


def half_disc(chi_rad, side, diam, g, n=15):
    """Half of a disc of diameter `diam` cut by a gap `g` through the centre, perpendicular to the diode axis (angle chi); side = +-1."""
    r = diam / 2
    th0 = math.acos(g / 2 / r)
    pts = [(r * math.cos(t), r * math.sin(t)) for t in np.linspace(-th0, th0, n)]
    c, s_ = math.cos(chi_rad), math.sin(chi_rad)
    return Poly(tuple((side * (c * px - s_ * py), side * (s_ * px + c * py)) for px, py in pts))


def rot_cell(state, chi_deg, diam=6.0e-3, g=0.4e-3, wd=1.2e-3):
    t = math.radians(chi_deg)
    halves = tuple(half_disc(t, s, diam, g) for s in (-1, 1))
    geo = PolyComposite(metal=halves + (oriented_rect(0.0, 0.0, g, wd, t),))
    d = "x" if chi_deg % 180 == 0 else "y" if chi_deg % 180 == 90 else float(chi_deg)
    ld = dm.load(gc.STATES[state], 0.0, 0.0, g, wd, d)
    return UnitCell(RectangularLattice(gc.P, gc.P), geo, gc.SUB, Metal(), backing="pec", loads=(ld,))


def uv(y, chi_deg):
    t = math.radians(chi_deg)
    u, v = np.array([math.cos(t), math.sin(t)]), np.array([-math.sin(t), math.cos(t)])
    return u @ y @ u, v @ y @ v, u @ y @ v


if __name__ == "__main__":
    out = {}
    for mesh in (32,):
        for chi in (0, 22.5, 45):
            row = {}
            for st in ("OFF", "ON"):
                y = sheet_tensor(rot_cell(st, chi), gc.F0, mesh) * 1e3
                yu, yv, yuv = uv(y, chi)
                row[st] = {"y_u_mS": [yu.real, yu.imag], "y_v_mS": [yv.real, yv.imag], "y_uv_mS": [yuv.real, yuv.imag]}
            out[f"m{mesh}_chi{chi}"] = row
            print(mesh, chi, {s: (np.round(complex(*r['y_u_mS']), 3), np.round(complex(*r['y_v_mS']), 3), np.round(complex(*r['y_uv_mS']), 3)) for s, r in row.items()}, flush=True)
    Path("results/sidequest_hybrid/coded_pol_cell.json").write_text(json.dumps(out, indent=1), encoding="utf-8", newline="\n")
