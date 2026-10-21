"""Phase 20 R4 item 1b: does state k have a mirror symmetry about its own stub axis (angle 45 k deg)? Mirroring maps incidence azimuth phi -> 2*45k - phi and, in
the incidence-plane (TE, TM) frame, S_TE,TE and S_TM,TM are unchanged while S_TE,TM and S_TM,TE change sign. The stub set {0,45,90,135} is NOT mirror symmetric
about any of its members (the OFF stubs have no partner), so this can only be approximate; measured here at mesh 60.

    python scripts/research/phase20_r4_mirror.py <out.json>
"""

import json
import sys
from pathlib import Path

import numpy as np
from phase20_r4_symmetry import FREQS, ROOT, N, cell, solve_all_states, wrap  # noqa: F401

THETAS = (0.0, 30.0, 45.0)
PHIS = tuple(range(0, 360, 15))


def main(out):
    c = cell()
    table = {}
    for f in (10.0e9,):
        for th in THETAS:
            for phi in PHIS:
                table[(f, th, phi)] = solve_all_states(c, f, th, phi)
    rows = []
    for (f, th, phi), s_k in table.items():
        for k in range(4):
            phi_m = (90 * k - phi) % 360
            ref = table[(f, th, phi_m)][k].copy()
            ref[0, 1] *= -1
            ref[1, 0] *= -1
            s = s_k[k]
            rows.append(
                {
                    "theta": th,
                    "phi": phi,
                    "k": k,
                    "tm_phase_err": wrap(float(np.degrees(np.angle(s[1, 1] / ref[1, 1])))),
                    "tm_db_err": float(20 * np.log10(abs(s[1, 1]) / abs(ref[1, 1]))),
                    "te_phase_err": wrap(float(np.degrees(np.angle(s[0, 0] / ref[0, 0])))),
                    "max_abs_dS": float(np.abs(s - ref).max()),
                }
            )
    Path(out).write_text(json.dumps(rows, indent=1))
    for k in range(4):
        for th in THETAS:
            sel = [r for r in rows if r["k"] == k and r["theta"] == th]
            v = np.array([abs(r["tm_phase_err"]) for r in sel])
            d = np.array([r["max_abs_dS"] for r in sel])
            print(
                f"k {k} theta {th:g}: TM phase err median {np.median(v):.1f} max {v.max():.1f} deg; max|dS| median {np.median(d):.3f} max {d.max():.3f}"
            )


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ROOT / "results/phase20_r4/mirror.json")
