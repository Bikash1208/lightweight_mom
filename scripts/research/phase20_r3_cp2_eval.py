"""CP-2 4-state evaluation (PB ladder, handedness-preserving r_RR, tolerance corners) of a given geometry.

python scripts/research/phase20_r3_cp2_eval.py <period_mm> <madp|smp> R line width   (mm)
"""

import sys
from multiprocessing import Pool

from phase20_r3_opt import corner_kits, lp_responses, lp_score
from phase20_r3_opt_cp2 import make_cell

from src.ris.cells import madp_kit, smp1345_kit


def evaluate(cell, kit, pool, thetas=(0.0, 30.0)):
    resp = lp_responses(cell, corner_kits(kit), thetas, pool, basis="circular")
    return resp, lp_score(resp)


if __name__ == "__main__":
    period = float(sys.argv[1]) * 1e-3
    kit = madp_kit() if sys.argv[2] == "madp" else smp1345_kit()
    cell = make_cell([float(a) for a in sys.argv[3:6]], period)
    print(cell.line, cell.stub_extent())
    with Pool(6) as pool:
        resp, (score, metrics) = evaluate(cell, kit, pool)
    for k, (p, a, x) in resp.items():
        print(
            k,
            "phase",
            p.round(0),
            "|rRR|dB",
            (20 * __import__("numpy").log10(a)).round(1),
            "opp-hand dB",
            (20 * __import__("numpy").log10(x)).round(1),
        )
    print(score, metrics)
