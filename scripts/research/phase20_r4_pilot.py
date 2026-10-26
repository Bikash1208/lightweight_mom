"""Phase 20 R4 item 2: pilot for the mesh-correction map. At 10 GHz, theta = 30 (phi every 45 deg) and theta = 45 (phi 0, 90, 180, 270), all four states, meshes 60 / 90 / 120:
how large is the correction r(120) - r(60) in r_RR (phase, dB) and is it smooth in phi (a low-order Fourier fit in phi, per state)? Also records wall time and peak memory per
point. MODEL_ONLY.

    python scripts/research/phase20_r4_pilot.py <out.json>
"""

import json
import math
import resource
import sys
import time
import warnings
from pathlib import Path

import numpy as np
from phase20_r4_symmetry import PL, ROOT, cell, wrap

from src.basis_functions import MeshSettings
from src.ris.cells import madp_kit
from src.ris.polarization import circular_reflection
from src.sparameters import compute_s_parameters_oblique_over_loads

POINTS = [(30.0, p) for p in range(0, 360, 45)] + [(45.0, p) for p in (0, 90, 180, 270)]
MESHES = (60, 90, 120)


def main(out):
    c = cell()
    kit = madp_kit()
    sets = [c.unit_cell(st, kit).loads for st in c.states()]
    base = c.unit_cell(c.states()[0], kit)
    rows = []
    for th, phi in POINTS:
        for n in MESHES:
            t = time.time()
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res = compute_s_parameters_oblique_over_loads(
                    base,
                    10e9,
                    MeshSettings(n, n, n // 2 + 1),
                    PL,
                    math.radians(th),
                    math.radians(phi),
                    sets,
                    warn_grating=False,
                )
            dt = time.time() - t
            for k, r in enumerate(res):
                cr = circular_reflection(r.s11)
                rows.append(
                    {
                        "theta": th,
                        "phi": phi,
                        "mesh": n,
                        "state": k,
                        "rr_re": float(cr[0, 0].real),
                        "rr_im": float(cr[0, 0].imag),
                        "lr_abs": float(abs(cr[1, 0])),
                        "seconds_all_states": dt,
                        "maxrss_gb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9,
                    }
                )
        print(th, phi, "done", flush=True)
    Path(out).write_text(json.dumps(rows, indent=1))


def analyse(path):
    rows = json.loads(Path(path).read_text())
    key = {(r["theta"], r["phi"], r["mesh"], r["state"]): complex(r["rr_re"], r["rr_im"]) for r in rows}
    out = {}
    for th in (30.0, 45.0):
        for k in range(4):
            phis = sorted({r["phi"] for r in rows if r["theta"] == th})
            d_ph = np.array([wrap(math.degrees(np.angle(key[(th, p, 120, k)] / key[(th, p, 60, k)]))) for p in phis])
            d_db = np.array([20 * math.log10(abs(key[(th, p, 120, k)]) / abs(key[(th, p, 60, k)])) for p in phis])
            d9 = np.array([wrap(math.degrees(np.angle(key[(th, p, 120, k)] / key[(th, p, 90, k)]))) for p in phis])
            out[f"theta{th:g}_state{k}"] = {
                "phi": phis,
                "dphase_120_60": d_ph.tolist(),
                "ddb_120_60": d_db.tolist(),
                "dphase_120_90": d9.tolist(),
            }
    return out


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ROOT / "results/phase20_r4/pilot.json")
