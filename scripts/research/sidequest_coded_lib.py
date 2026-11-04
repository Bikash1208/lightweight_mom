"""Cell sub-library for the crossed rim cell (x element varied; y element fixed: slit 5.6 mm, diode OFF). MoM sheet tensor, mesh 64, 10.54 GHz, measured MA4AGP907.
Switched variants (diode OFF / ON): slit half-length ls and diode position pos. Static variants (no diode): open slit, shorted slit (metal in the footprint), or no slit.
    python sidequest_coded_lib.py   -> results/sidequest_hybrid/coded_lib.json (Y_xx, Y_yy per variant)"""
import json
import sys
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
G, WD = 0.4, 1.2


def slit_pieces(ls, pos, footprint_metal):
    """Vertical slit pieces (x in +-G/2, y in [-ls, ls]) minus the footprint [pos-WD/2, pos+WD/2] when footprint_metal."""
    if not footprint_metal:
        return [(-ls, ls)]
    segs = []
    if pos - WD / 2 > -ls:
        segs.append((-ls, pos - WD / 2))
    if pos + WD / 2 < ls:
        segs.append((pos + WD / 2, ls))
    return segs


def lib_cell(ls, pos, mode, diam=6.0):
    disc = regular_polygon(0.0, 0.0, diam / 2 * MM, 48, 0.0)
    ap, loads = [], []
    if mode != "none":
        for a0, a1 in slit_pieces(ls, pos, mode in ("OFF", "ON", "short")):
            ap.append(AxisRect(0.0, (a0 + a1) / 2 * MM, G * MM, (a1 - a0) * MM))
        if mode in ("OFF", "ON"):
            loads.append(dm.load(gc.STATES[mode], 0.0, pos * MM, G * MM, WD * MM, "x"))
    # fixed y element: slit 5.6 mm, diode OFF at x = +1.6
    for a0, a1 in slit_pieces(2.8, 1.6, True):
        ap.append(AxisRect((a0 + a1) / 2 * MM, 0.0, (a1 - a0) * MM, G * MM))
    loads.append(dm.load(gc.STATES["OFF"], 1.6 * MM, 0.0, G * MM, WD * MM, "y"))
    return UnitCell(RectangularLattice(gc.P, gc.P), PolyComposite(metal=(disc,), apertures=tuple(ap)), gc.SUB, Metal(), backing="pec", loads=tuple(loads))


VARIANTS = ([("sw", ls, pos, st) for ls in (2.6, 2.8) for pos in (0.8, 1.2, 1.6, 2.0) if pos + WD / 2 <= ls for st in ("OFF", "ON")]
            + [("st", ls, 0.0, md) for ls in (1.0, 2.0, 2.8) for md in ("open", "short")] + [("st", 0.0, 0.0, "none")])


def job(v):
    warnings.simplefilter("ignore")
    kind, ls, pos, mode = v
    y = sheet_tensor(lib_cell(ls, pos, mode), gc.F0, 64)
    return {"kind": kind, "ls": ls, "pos": pos, "mode": mode, "Y_mS": [[complex(a).real * 1e3, complex(a).imag * 1e3] for a in y.ravel()]}


if __name__ == "__main__":
    with ProcessPoolExecutor(3) as ex:
        res = list(ex.map(job, VARIANTS))
    for r in res:
        print(r["kind"], r["ls"], r["pos"], r["mode"], "Yxx %.2f%+.2fj Yyy %.2f%+.2fj" % (*r["Y_mS"][0], *r["Y_mS"][3]), flush=True)
    Path("results/sidequest_hybrid/coded_lib.json").write_text(json.dumps(res, indent=1), encoding="utf-8", newline="\n")
