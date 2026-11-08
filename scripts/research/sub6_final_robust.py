"""SUB6 FINAL BOARD: robustness of the chosen design (200 mm, N 240, calibrated, switched-cell table of SUB6_MESH) against the island-cell values: ALL switched-cell alpha (OFF and ON) and |Delta| (OFF and ON, phase kept)
scaled by s_alpha / s_delta, NO re-optimisation (sites, design phases, slot rotations and directions fixed; the marches / phases / weights are recomputed with the scaled values).  MODEL_ONLY.
    python sub6_final_robust.py        (SUB6_MESH=96 for the mesh-96 table)
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sub6_final_cal as FC
import sub6_final_design as FD
import sub6_board_island2 as X
import sub6_board_model4 as R4

GATES = dict(bal=("<=", 3.0), off_ARA3=("<=", 3.0), off_ARB3=("<=", 3.0), on_ARA3=("<=", 3.0), supp=(">=", 12.0), sll_off=("<=", -13.0), sll_on=("<=", -13.0))


def gate_ok(m):
    return {k: bool(m[k] <= v if op == "<=" else m[k] >= v) for k, (op, v) in GATES.items()}


def main():
    fd = FD.final_design("cal")
    D, sa = fd["D"], fd["sa"]
    mask = D["mask"]
    rows = []
    for sa_s, sd_s, so, sn in ((1.0, 1.0, 1.0, 1.0), (0.85, 0.85, 1.0, 1.0), (1.15, 1.15, 1.0, 1.0), (1.0, 1.0, 0.85, 1.15), (1.0, 1.0, 1.15, 0.85), (1.0, 1.0, 0.85, 1.0), (1.0, 1.0, 1.0, 0.85), (1.0, 1.0, 1.15, 1.0), (1.0, 1.0, 1.0, 1.15)):
        s2 = dict(sa)
        for k in ("alpha_off", "alpha_on"):
            s2[k] = np.where(mask, sa[k] * sa_s * (so if k == 'alpha_off' else sn), sa[k])
        for k in ("D_off", "D_on"):
            s2[k] = np.where(mask, sa[k] * sd_s, sa[k])
        ph = R4.physics(D, s2)
        m = X.metrics_state(D, s2, ph, fd["direction"], fd["flips"], FC.LAUNCH)
        g = gate_ok(m)
        rows.append(dict(s_alpha=sa_s, s_delta=sd_s, s_alpha_off_extra=so, s_alpha_on_extra=sn, all_gates_pass=all(g.values()), failed=[k for k, v in g.items() if not v], **{k: (float(v) if not isinstance(v, str) else v) for k, v in m.items()}))
        r = rows[-1]
        print(f"alpha x{sa_s} Delta x{sd_s} offx{so} onx{sn}: GA {r['off_GA']:.2f} GB {r['off_GB']:.2f} bal {r['bal']:.2f} ARA/B {r['off_ARA3']:.2f}/{r['off_ARB3']:.2f} ON GA {r['on_GA']:.2f} supp {r['supp']:.1f} AR_on {r['on_ARA3']:.2f} SLL {r['sll_off']:.1f}/{r['sll_on']:.1f} eta {r['eta_off']:.3f}/{r['eta_on']:.3f} gate {r['gate']} margin {r['margin']:.2f} pass={r['all_gates_pass']} {r['failed']}", flush=True)
    req = [r for r in rows if r["s_alpha"] == r["s_delta"] and r["s_alpha"] != 1.0]  # joint x0.85 / x1.15 (required)
    out = dict(label="MODEL_ONLY / NOT_HFSS_VALIDATED", note="|Delta| does not enter the board metrics of this model (only alpha, through the amplitude ratio sqrt(alpha_ON/alpha_OFF)); Delta x s is therefore inert, the extra rows scale alpha_OFF and alpha_ON separately (ratio change) as the real stress test.", mesh=FC.MESH, config=FD.CONFIG, rows=rows, required_pass=bool(all(r["all_gates_pass"] for r in req)), corners_pass=bool(all(r["all_gates_pass"] for r in rows)))
    (FC.OUTD / "robustness.json").write_text(json.dumps(out, indent=1), encoding="utf-8", newline="\n")
    print("REQUIRED (joint x0.85, x1.15) pass:", out["required_pass"], " all 9 pass:", out["corners_pass"])


if __name__ == "__main__":
    main()
