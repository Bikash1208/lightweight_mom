"""Phase 20 R4 item 1: is state k at azimuth phi the same as state 0 at phi - 45 k deg (in the incidence-plane TE/TM frame)? Rotating the whole structure and
the incident direction by 45 k deg maps state-0's stub set onto a different set (the three OFF stubs move) and the square lattice onto a rotated one,
so the relation can only hold approximately; this script measures how well, at mesh 60 (the R3 search mesh), for the registered CP-2 cell.

    python scripts/research/phase20_r4_symmetry.py <out.json>
"""

import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

from phase20_r3_finalist import build  # noqa: E402

from src.basis_functions import MeshSettings  # noqa: E402
from src.excitation import ReferencePlanes  # noqa: E402
from src.ris.cells import madp_kit  # noqa: E402
from src.ris.polarization import circular_reflection  # noqa: E402
from src.sparameters import compute_s_parameters_oblique_over_loads  # noqa: E402

PL = ReferencePlanes(0.0, None, 1.575e-3)
N = 60
FREQS = (9.0e9, 10.0e9, 11.0e9)
THETAS = (0.0, 30.0, 45.0)
PHIS = tuple(range(0, 360, 15))


def wrap(x):
    return (x + 180.0) % 360.0 - 180.0


def cell():
    d = json.loads((ROOT / "results/phase20_r3/opt/cp2_madp_final.json").read_text())
    return build("cp2", d["params_mm"], d["period_mm"])


def solve_all_states(c, f, th, phi):
    kit = madp_kit()
    sets = [c.unit_cell(st, kit).loads for st in c.states()]
    base = c.unit_cell(c.states()[0], kit)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = compute_s_parameters_oblique_over_loads(
            base, f, MeshSettings(N, N, N // 2 + 1), PL, math.radians(th), math.radians(phi), sets, warn_grating=False
        )
    return [r.s11 for r in res]


def main(out):
    c = cell()
    table = {}
    for f in FREQS:
        for th in THETAS:
            phis = PHIS
            for phi in phis:
                table[(f, th, phi)] = solve_all_states(c, f, th, phi)
    rows = []
    for (f, th, phi), s_k in table.items():
        for k in (1, 2, 3):
            ref = table[(f, th, (phi - 45 * k) % 360)][0]
            s = s_k[k]
            cr, rr = circular_reflection(s), circular_reflection(ref)
            rows.append(
                {
                    "f_ghz": f / 1e9,
                    "theta": th,
                    "phi": phi,
                    "k": k,
                    "co_tm_phase_err": wrap(math.degrees(np.angle(s[1, 1] / ref[1, 1]))),
                    "co_tm_db_err": 20 * math.log10(abs(s[1, 1]) / abs(ref[1, 1])),
                    "co_te_phase_err": wrap(math.degrees(np.angle(s[0, 0] / ref[0, 0]))),
                    "rr_phase_err": wrap(math.degrees(np.angle(cr[0, 0] / rr[0, 0]))),
                    "rr_db_err": 20 * math.log10(abs(cr[0, 0]) / abs(rr[0, 0])),
                    "max_abs_dS": float(np.abs(s - ref).max()),
                }
            )
    Path(out).write_text(json.dumps(rows, indent=1))
    for key in ("rr_phase_err", "rr_db_err", "co_tm_phase_err", "co_te_phase_err", "max_abs_dS"):
        v = np.array([abs(r[key]) for r in rows])
        print(f"{key} median {np.median(v):.3f} 95% {np.percentile(v, 95):.3f} max {v.max():.3f}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ROOT / "results/phase20_r4/symmetry.json")
