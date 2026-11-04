"""Side quest step 2: characterisation of the selected gate cell (dog-bone, one MA4AGP907 along x). MODEL_ONLY, mesh 32 and 64, F = n.
    python scripts/research/sidequest_gate_final.py > log; results/sidequest_hybrid/gate_cell_final.json
"""
import json, sys, time, warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import sidequest_gate_cell as gc
from src.hologram.coded_cell_eval import sheet_tensor

CAND = dict(e=1.2, s=0.8, h=6.0, g=0.4, lt=5.6, wd=0.8)  # mm; aligned to 0.2 mm (mesh 32) and 0.1 mm (mesh 64)
B_REF_10P54 = 7.82e-3


def cell(state, scale=None):
    c = CAND
    return gc.make_dogbone(c["e"] * 1e-3, c["s"] * 1e-3, c["h"] * 1e-3, c["g"] * 1e-3, state, scale=scale, lt=c["lt"] * 1e-3, wd=c["wd"] * 1e-3)


def solve(job):
    warnings.simplefilter("ignore")
    state, mesh, f, scale = job
    t = time.time()
    y = sheet_tensor(cell(state, scale), f, mesh)
    return {"state": state, "mesh": mesh, "f_ghz": f / 1e9, "scale": scale, "y_ms": [[[complex(v).real * 1e3, complex(v).imag * 1e3] for v in row] for row in y], "sec": time.time() - t}


def line(yxx, f, rows):
    old = (gc.F0, gc.K0, gc.W0)
    gc.F0, gc.K0, gc.W0 = f, gc.free_space_wavenumber(f), gc.angular_frequency(f)
    try:
        yref = 1j * B_REF_10P54 * f / 10.54e9  # assumption: the hologram mean sheet is capacitive-like, b ~ f
        return gc.strip_transmission(yref, yxx, rows)
    finally:
        gc.F0, gc.K0, gc.W0 = old


if __name__ == "__main__":
    jobs = [(st, 32, f * 1e9, None) for st in ("OFF", "ON") for f in np.round(np.arange(10.0, 11.01, 0.1), 2)]
    jobs += [(st, m, 10.54e9, None) for st in ("OFF", "ON") for m in (64,)]
    jobs += [("ON", 64, 10.54e9, {"R": 0.8}), ("ON", 64, 10.54e9, {"R": 1.2}), ("OFF", 64, 10.54e9, {"C": 0.8}), ("OFF", 64, 10.54e9, {"C": 1.2})]
    jobs += [("ON", 32, 10.54e9, {"R": 0.8}), ("ON", 32, 10.54e9, {"R": 1.2}), ("OFF", 32, 10.54e9, {"C": 0.8}), ("OFF", 32, 10.54e9, {"C": 1.2})]
    jobs += [(st, 64, f * 1e9, None) for st in ("OFF", "ON") for f in (10.2, 10.9)]
    t0 = time.time()
    with ProcessPoolExecutor(4) as ex:
        rows = list(ex.map(solve, jobs))
    Path("results/sidequest_hybrid/gate_cell_final_raw.json").write_text(json.dumps({"cand_mm": CAND, "rows": rows}, indent=1), encoding="utf-8", newline="\n")
    print("done", round(time.time() - t0), "s")
