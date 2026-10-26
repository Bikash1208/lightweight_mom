"""Phase 20 R5-1 (SERVER JOB `r5-tolerance-mom`): MoM prediction of the diode-tolerance cases the user runs in HFSS: baseline CP-2, MADP diode OFF capacitance C_off = 0.0225 / 0.025 / 0.030 pF
(= C_T x 0.9 / 1 / 1.2), four one-hot states, theta = 0 and 30 deg, phi = 0, f = 9.5-11.5 GHz in 0.05 GHz steps, meshes 60 and 90 at F = 2n. Predictions are written BEFORE the HFSS comparison. MODEL_ONLY.

    python scripts/research/phase20_r5_tolerance_mom.py <out.npz> [--workers N] [--smoke]
Run with OPENBLAS_NUM_THREADS=1 per worker.
"""

import json
import math
import os
import sys
import warnings
from multiprocessing import Pool
from pathlib import Path

for _k in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "1")

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

from phase20_r3_finalist import build  # noqa: E402
from provenance import finish, provenance  # noqa: E402

from src.basis_functions import MeshSettings  # noqa: E402
from src.excitation import ReferencePlanes  # noqa: E402
from src.ris.cells import DiodeKit, madp_kit  # noqa: E402
from src.ris.polarization import circular_reflection, to_structure_frame  # noqa: E402
from src.sparameters import compute_s_parameters_oblique_over_loads  # noqa: E402

COFFS = (0.0225, 0.025, 0.030)  # pF
THETAS = (0.0, 30.0)
PL = ReferencePlanes(0.0, None, 1.575e-3)


def _cell():
    d = json.loads((ROOT / "results/phase20_r3/opt/cp2_madp_final.json").read_text())
    return build("cp2", d["params_mm"], d["period_mm"])


def solve(task):
    f, th, n = task
    cell = _cell()
    nom = madp_kit()
    kits = [DiodeKit(nom.name, nom.on, nom.off.scaled(C_T=c / 0.025)) for c in COFFS]
    sets = [cell.unit_cell(st, k).loads for k in kits for st in cell.states()]
    base = cell.unit_cell(cell.states()[0], nom)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = compute_s_parameters_oblique_over_loads(
            base, f * 1e9, MeshSettings(n, n, 2 * n), PL, math.radians(th), 0.0, sets, warn_grating=False
        )
    c = np.array([to_structure_frame(circular_reflection(r.s11), 0.0) for r in res])  # (3*4, 2, 2)
    return c.reshape(len(COFFS), 4, 2, 2)


def main(out, workers, smoke):
    freqs = np.round(np.arange(10.5, 10.51, 0.05), 3) if smoke else np.round(np.arange(9.5, 11.5001, 0.05), 3)
    meshes = (24, 36) if smoke else (60, 90)
    rec = provenance({"job": "r5-tolerance-mom", "coffs_pF": COFFS, "meshes": meshes, "smoke": smoke}, workers=workers)
    tasks = [(float(f), th, n) for n in sorted(meshes, reverse=True) for th in THETAS for f in freqs]
    with Pool(workers) as pool:
        res = pool.map(solve, tasks, chunksize=1)
    arr = np.zeros((len(meshes), len(COFFS), len(freqs), len(THETAS), 4, 2, 2), dtype=complex)
    for (f, th, n), r in zip(tasks, res, strict=True):
        arr[meshes.index(n), :, int(np.argmin(abs(freqs - f))), THETAS.index(th)] = r
    np.savez_compressed(
        out,
        c_structure=arr,
        coffs_pf=np.array(COFFS),
        f_ghz=freqs,
        theta_deg=np.array(THETAS),
        meshes=np.array(meshes),
        provenance=np.array([json.dumps(finish(rec))]),
    )
    print("wrote", out, arr.shape)


if __name__ == "__main__":
    wk = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 4
    main(sys.argv[1], wk, "--smoke" in sys.argv)
