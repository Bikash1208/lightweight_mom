"""Phase 20 R3 item 3: full target table (theta 0/30 nominal, tolerance corners at theta 0, band) of the BEST optimiser result of an LP family, at the search mesh.
The LP cells miss the phase target by a wide margin (see the stage-1 files), so no finalist mesh correction is run for them: this table only quantifies the miss.

    python scripts/research/phase20_r3_lp_report.py <lp1|lp2|lp3> <madp|smp> <stage1.json> <out.json>
"""

import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from phase20_r3_opt import corner_kits, lp_responses, lp_score

from src.ris.cells import madp_kit, smp1345_kit


def build(family, x, names, codes=None):
    kw = {k: float(v) * 1e-3 for k, v in zip(names, x, strict=True)}
    if family == "lp3":
        from phase20_r3_opt_lp3 import SubsetCell

        from src.ris.cells import Lp3Cell

        return SubsetCell(Lp3Cell(**kw), tuple(tuple(c) for c in codes))
    from src.ris.cells import Lp1Cell, Lp2Cell

    return {"lp1": Lp1Cell, "lp2": Lp2Cell}[family](**kw)


def main(family, kit_name, src, out):
    d = json.loads(Path(src).read_text())
    best = min(d["runs"], key=lambda r: r["J"])
    kit = madp_kit() if kit_name == "madp" else smp1345_kit()
    cell = build(family, best["x_mm"], d["names"], best.get("codes"))
    with Pool(8) as pool:
        resp = lp_responses(cell, corner_kits(kit), (0.0, 30.0), pool)
    score, metrics = lp_score(resp)
    table = {
        f"{k}@{th:g}": {
            "phase_deg": p.tolist(),
            "amp_db": (20 * np.log10(a)).tolist(),
            "xpol_db": float(20 * np.log10(max(x.max(), 1e-9))),
        }
        for (k, th), (p, a, x) in resp.items()
    }
    rec = {
        "family": family,
        "kit": kit.name,
        "x_mm": best["x_mm"],
        "names": d["names"],
        "codes": best.get("codes"),
        "J_stage1": best["J"],
        "J_full": score,
        "metrics": metrics,
        "table": table,
        "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "PROVISIONAL_DATASHEET", "VIA_OMITTED_APPROXIMATION"],
    }
    Path(out).write_text(json.dumps(rec, indent=1))
    print(json.dumps({"x_mm": rec["x_mm"], "metrics": metrics}))


if __name__ == "__main__":
    main(*sys.argv[1:5])
