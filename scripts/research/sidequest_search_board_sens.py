"""Board sensitivity to the switched-cell alpha (eps_r +-0.3 and mesh move it): RT5880 (round-11 option sets, A0 in a grid) and RO3010 (E1 static cells, alpha_E1 13.1 Np/m, A0 in a grid). Output cell_search_board_sens.json."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_search_board as S
import sidequest_coded_improve as I
import sidequest_coded_improve_run as R

if __name__ == "__main__":
    VS_BIT, PAIRS_BIT = I.bit_sets()
    kw = dict(R.BASE)
    for it in ("I1", "I2", "I5"):
        kw = R.ITEMS[it](kw)
    out = {}
    S.GS = (2.0,)
    for a0 in (0.39, 0.7, 1.0, 1.84, 3.43):
        S.configure(10.54e9, a0, 0.094 * a0 / 1.84 + 0.0, 0.15, 1.079, 0.627)
        S.go(f"RT5880_A0_{a0}", R.bd, kw["vs"], kw["pairs"], mu=kw["mu"], gs=(2.0,))
        out[f"RT5880_A0_{a0}"] = S.LOG[f"RT5880_A0_{a0}"]
    for a0 in (0.85, 1.15, 1.6, 3.0, 7.125):
        S.configure(10.05e9, a0, 2.517, 0.6, 1.14, 0.476)
        bd = I.Board()
        tens = I.e1_options(float(np.sqrt(13.1 / a0)))
        S.go(f"RO3010_A0_{a0}", bd, VS_BIT, PAIRS_BIT, tens=tens, gs=(1.0,))
        out[f"RO3010_A0_{a0}"] = S.LOG[f"RO3010_A0_{a0}"]
    Path("results/sidequest_hybrid/cell_search_board_sens.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8", newline="\n")
