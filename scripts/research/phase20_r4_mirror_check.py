"""Phase 20 R4 full-grid, mirror check (server): at mesh 60, F = 2n, f = 10.5 GHz, state 0, is S(theta, -phi) the mirror image of S(theta, phi)?

The structure's state 0 has its ON stub on the x axis, so mirroring y -> -y maps phi -> -phi; in the incidence-plane (TE, TM) frame the co-pol terms S_TE,TE and S_TM,TM
stay and the cross terms S_TE,TM, S_TM,TE change sign. The stub set {0, 45, 90, 135} has no partner for the OFF stubs, so this is approximate (R4: SYMMETRY.md).

    python scripts/research/phase20_r4_mirror_check.py <out.json>
"""

import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

import phase20_r4_build_library as lib  # noqa: E402

F_GHZ = 10.5
MESH = 60
PAIRS = [(15.0, 30.0), (15.0, 60.0), (30.0, 30.0), (30.0, 60.0)]  # (theta, phi); the mirror partner is 360 - phi


def main(out):
    jobs = []
    for th, ph in PAIRS:
        jobs += [(F_GHZ, th, ph, MESH), (F_GHZ, th, 360.0 - ph, MESH)]
    with Pool(len(jobs)) as pool:
        res = pool.map(lib.solve_point, jobs, chunksize=1)
    sol = {(j[1], j[2]): r[0] for j, r in zip(jobs, res, strict=True)}  # state 0, 2x2 (TE, TM)
    rows = []
    for th, ph in PAIRS:
        a, b = sol[(th, ph)], sol[(th, 360.0 - ph)].copy()
        b[0, 1] *= -1
        b[1, 0] *= -1
        rows.append(
            {
                "theta": th,
                "phi": ph,
                "partner_phi": 360.0 - ph,
                "max_abs_dS": float(np.abs(a - b).max()),
                "max_rel_dS": float(np.abs(a - b).max() / np.abs(a).max()),
                "S": [[[z.real, z.imag] for z in row] for row in a],
                "S_partner_mirrored": [[[z.real, z.imag] for z in row] for row in b],
            }
        )
    worst = max(r["max_rel_dS"] for r in rows)
    Path(out).write_text(json.dumps({"f_ghz": F_GHZ, "mesh": MESH, "floquet": lib.floquet_order(MESH), "state": 0, "max_rel_dS": worst, "rows": rows}, indent=1))
    for r in rows:
        print(f"theta {r['theta']:g} phi {r['phi']:g}/{r['partner_phi']:g}: max|dS| {r['max_abs_dS']:.4g} (rel {r['max_rel_dS']:.4g})")
    print(f"max relative |dS| = {worst:.4g}  ->  {'MIRROR_USABLE' if worst <= 1e-6 else 'MIRROR_NOT_USABLE (keep phi 0-360)'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ROOT / "results/infra/r4_full_grid_mirror_check.json")
