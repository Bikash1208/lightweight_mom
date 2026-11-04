"""Rim robustness / mesh convergence of the crossed two-diode cell: rim = metal left at each slit end (ls = 3.0 - rim). Writes results/sidequest_hybrid/coded_rim/<rim>_<mesh>.json.
    python sidequest_coded_rim.py <mesh> <workers> <states: 2|4> rim1 rim2 ..."""
import itertools
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import sidequest_coded_cross_cell as X
import sidequest_gate_cell as gc
from src.hologram.coded_cell_eval import sheet_tensor

OUT = Path("results/sidequest_hybrid/coded_rim")


def job(p):
    warnings.simplefilter("ignore")
    rim, mesh, sx, sy = p
    y = sheet_tensor(X.cross_cell(sx, sy, ls=round(3.0 - rim, 4)), gc.F0, mesh)
    return p, [[complex(v).real * 1e3, complex(v).imag * 1e3] for v in y.ravel()]


if __name__ == "__main__":
    mesh, workers, nst = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    rims = [float(r) for r in sys.argv[4:]]
    OUT.mkdir(parents=True, exist_ok=True)
    states = list(itertools.product(("OFF", "ON"), repeat=2)) if nst == 4 else [("OFF", "OFF"), ("ON", "OFF")]
    jobs = [(r, mesh, a, b) for r in rims for a, b in states]
    with ProcessPoolExecutor(workers) as ex:
        for (rim, m, a, b), y in ex.map(job, jobs):
            f = OUT / f"rim{rim}_m{m}.json"
            d = json.loads(f.read_text()) if f.exists() else {"rim": rim, "mesh": m, "states": {}}
            d["states"][f"{a}/{b}"] = y
            f.write_text(json.dumps(d, indent=1), encoding="utf-8", newline="\n")
            print(rim, m, a, b, flush=True)
