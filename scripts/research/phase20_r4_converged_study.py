"""Phase 20 R4 (finding): mesh convergence of the CP-2 cell at 10 GHz with a CONVERGED Floquet cutoff (F = 2n: meshes 60 / 90 / 120 -> F 120 / 180 / 240), all four states,
theta 0 and 30 deg (phi 0). Compares with the minimal-cutoff values used in R3 and in the first R4 grid. MODEL_ONLY.

    python scripts/research/phase20_r4_converged_study.py <out.json>
"""

import json
import math
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import phase20_r4_build_library as b

from src.basis_functions import MeshSettings
from src.ris.polarization import circular_reflection
from src.sparameters import compute_s_parameters_oblique_over_loads


def run(f, th, ph, n, big_f):
    st = b._setup()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = compute_s_parameters_oblique_over_loads(
            st["base"],
            f * 1e9,
            MeshSettings(n, n, big_f),
            st["planes"],
            math.radians(th),
            math.radians(ph),
            st["sets"],
            warn_grating=False,
        )
    return np.array([circular_reflection(r.s11) for r in res])


def main(out):
    rows = []
    for th in (0.0, 30.0):
        for n in (60, 90, 120):
            for label, big_f in (("minimal", n // 2 + 1), ("2n", 2 * n)):
                t = time.monotonic()
                c = run(10.0, th, 0.0, n, big_f)
                rows.append(
                    {
                        "theta": th,
                        "mesh": n,
                        "F": big_f,
                        "label": label,
                        "c": [[[[z.real, z.imag] for z in r] for r in s] for s in c],
                        "seconds": time.monotonic() - t,
                    }
                )
                ph = np.degrees(np.angle(c[:, 0, 0]))
                print(
                    f"theta {th:g} mesh {n} F {big_f:3d} ({label}): r_RR phase {np.round(ph, 1)}  dB {np.round(20 * np.log10(abs(c[:, 0, 0])), 2)}  opp {np.round(20 * np.log10(np.maximum(abs(c[:, 1, 0]), abs(c[:, 0, 1]))), 1)}  [{time.monotonic() - t:.0f} s]",
                    flush=True,
                )
    Path(out).write_text(json.dumps(rows))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "../../results/phase20_r4/converged_study.json")
