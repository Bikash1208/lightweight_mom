"""Retune the dog-bone ON state so that the no-bound-mode band is centred on 10.54 GHz (mesh 32 sweep, then mesh-64 check)."""
import json, sys, warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import sidequest_gate_cell as gc
from src.hologram.coded_cell_eval import sheet_tensor

FREQS = np.round(np.arange(10.0, 11.21, 0.05), 3)


def job(p):
    warnings.simplefilter("ignore")
    e, s, h, lt, g, wd, mesh, fl = p
    cell = gc.make_dogbone(e * 1e-3, s * 1e-3, h * 1e-3, g * 1e-3, "ON", lt=lt * 1e-3, wd=wd * 1e-3)
    ys = [complex(sheet_tensor(cell, f * 1e9, mesh)[0, 0]) for f in fl]
    flags = gc.no_mode_band([f * 1e9 for f in fl], ys)
    return {"p": p[:6], "mesh": mesh, "freqs": list(map(float, fl)), "y_ms": [[y.real * 1e3, y.imag * 1e3] for y in ys], "nomode": [bool(x) for x in flags]}


def band(r):
    fl, nm = np.array(r["freqs"]), np.array(r["nomode"])
    if not nm.any():
        return None
    idx = np.flatnonzero(nm)
    return float(fl[idx[0]]), float(fl[idx[-1]]), int(len(idx))


if __name__ == "__main__":
    pts = [(e, s, h, lt, 0.4, 0.8, 32, FREQS) for e in (0.8, 1.2, 1.6) for s in (0.8, 1.2) for h in (5.2, 5.6, 6.0) for lt in (5.6, 6.0)]
    with ProcessPoolExecutor(4) as ex:
        res = list(ex.map(job, pts))
    Path("results/sidequest_hybrid/gate_retune_mesh32.json").write_text(json.dumps(res), encoding="utf-8", newline="\n")
    for r in res:
        b = band(r)
        print(r["p"], b, None if b is None else round((b[0] + b[1]) / 2, 3))
