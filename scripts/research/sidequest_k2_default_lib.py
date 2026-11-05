"""Scalar cell library for the K2-default board (one control line, both diodes of a cell on the line), MODEL_ONLY. MoM sheet tensor (normal incidence, mesh 64) at 10.5 GHz, RT5880 1.575 mm, 6.4 mm lattice, measured MA4AGP907 (OFF_0V / ON_10mA).
Cell = disc r 3.0 mm, cross slit (width 0.4) of half-length ls, footprint 0.4 x 1.0 mm at 0.8 mm on each axis (x load at (0, +pos) across the vertical slit, y load at (+pos, 0)), both loads in the SAME state (scalar cell).
Modes: 'open' = diode OFF = full slit, NO load (HFSS C/C 10.670 GHz = MoM open slit 10.662, the lumped 28 fF is ineffective in HFSS; the MoM with a 28 fF load would give 10.525), ON = R-only loads, 'strap' (footprint metal, no load). Family A = rim cells (ls <= 2.8): ON lowers Y.
Family B = through-slit cells (ls > 3.0, four islands, diode the only bridge): ON raises Y.     python sidequest_k2_default_lib.py [workers]  ->  results/sidequest_hybrid/k2_default_lib.json"""
import json
import sys
from dataclasses import replace
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import sidequest_gate_cell as gc
from src.composite_geometry import AxisRect
from src.geometry import UnitCell
from src.hologram.coded_cell_eval import sheet_tensor
from src.materials import Metal
from src.periodic_bc import RectangularLattice
from src.polygon_geometry import PolyComposite, regular_polygon
from src.ris import diodes_measured as dm

MM = 1e-3
F = 10.5e9
G, WD, POS, R = 0.4, 1.0, 0.8, 3.0
R_ON = 2.44  # ohm: ON = this resistance only (HFSS Impedance boundary, Rs = R_ON * WD / G = 6.1 ohm/sq); the measured MA4AGP907 ON impedance 2.44 - j3.84 ohm is NOT used (not what HFSS gets)


def load(mode, x, y, direction):
    l = dm.load(gc.STATES[mode], x, y, G * MM, WD * MM, direction)
    return replace(l, z_of_f=(lambda f: complex(R_ON))) if mode == "ON" else l  # OFF = measured open diode, Z = -538j ohm at 10.5 GHz (28.2 fF)


def cell(ls, pos, mode):
    """mode: OFF, ON (loads on both axes), strap (footprint metal, no load), open (full slit, no load), none (no slit)."""
    disc = regular_polygon(0.0, 0.0, R * MM, 48, 0.0)
    ap, loads = [], []
    if mode != "none":
        foot = mode in ("OFF", "ON", "strap")
        segs = [(-ls, ls)]
        if foot:
            segs = []
            if pos - WD / 2 > -ls:
                segs.append((-ls, min(pos - WD / 2, ls)))
            if pos + WD / 2 < ls:
                segs.append((max(pos + WD / 2, -ls), ls))
        for a0, a1 in segs:
            ap.append(AxisRect(0.0, (a0 + a1) / 2 * MM, G * MM, (a1 - a0) * MM))
            ap.append(AxisRect((a0 + a1) / 2 * MM, 0.0, (a1 - a0) * MM, G * MM))
        if mode in ("OFF", "ON"):
            loads.append(load(mode, 0.0, pos * MM, "x"))
            loads.append(load(mode, pos * MM, 0.0, "y"))
    return UnitCell(RectangularLattice(gc.P, gc.P), PolyComposite(metal=(disc,), apertures=tuple(ap)), gc.SUB, Metal(), backing="pec", loads=tuple(loads))


VARIANTS = ([("A", 0.0, POS, "none"), ("A", 1.0, POS, "open")]  # solid disc / plain slit
            + [("A", ls, POS, "open") for ls in (1.8, 2.2, 2.6, 2.8)]
            + [("A", 2.8, POS, "ON")] + [("A", ls, POS, "strap") for ls in (2.2, 2.8)]
            + [("B", ls, pos, md) for ls in (2.9, 3.1) for pos in (0.8, 1.4, 2.0) for md in ("open", "ON")])


def job(v):
    warnings.simplefilter("ignore")
    fam, ls, pos, mode = v
    y = sheet_tensor(cell(ls, pos, mode), F, 64)
    return {"family": fam, "ls": ls, "pos": pos, "mode": mode, "Y_mS": [[complex(a).real * 1e3, complex(a).imag * 1e3] for a in y.ravel()]}


if __name__ == "__main__":
    nw = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    with ProcessPoolExecutor(nw) as ex:
        res = list(ex.map(job, VARIANTS))
    for r in res:
        print(r["family"], r["ls"], r["pos"], r["mode"], "Yxx %.2f%+.2fj Yyy %.2f%+.2fj" % (*r["Y_mS"][0], *r["Y_mS"][3]), flush=True)
    Path("results/sidequest_hybrid/k2_default_lib.json").write_text(json.dumps({"f_hz": F, "mesh": 64, "footprint_mm": [G, WD], "rows": res}, indent=1), encoding="utf-8", newline="\n")
