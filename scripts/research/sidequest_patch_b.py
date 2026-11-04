"""Square-patch (reference / transition) sheet susceptance at mesh 64 (0.1 mm grid), 10.54 GHz."""
import json, sys, warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sidequest_gate_cell as gc
from src.composite_geometry import AxisRect, CompositeGeometry
from src.geometry import RectangularLattice, UnitCell
from src.hologram.coded_cell_eval import sheet_tensor
from src.materials import Metal


def job(a):
    warnings.simplefilter("ignore")
    uc = UnitCell(RectangularLattice(gc.P, gc.P), CompositeGeometry(metal=(AxisRect(0, 0, a * 1e-3, a * 1e-3),)), gc.SUB, Metal(), backing="pec")
    y = sheet_tensor(uc, gc.F0, 64)
    return {"a_mm": a, "y_ms": [complex(y[0, 0]).real * 1e3, complex(y[0, 0]).imag * 1e3], "yy_ms": complex(y[1, 1]).imag * 1e3}


if __name__ == "__main__":
    with ProcessPoolExecutor(4) as ex:
        res = list(ex.map(job, [5.3, 5.4, 5.5, 5.7, 5.8, 5.9]))
    Path("results/sidequest_hybrid/patch_b_mesh64.json").write_text(json.dumps(res, indent=1), encoding="utf-8", newline="\n")
    for r in res:
        print(r)
