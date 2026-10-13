"""Phase 25 follow-up 2b: the 5x5 / 4x4 plane-wave locality check (same method as
`phase25_forward_validation.py` Part II) for the candidate pairs of `phase25_gap_pair_study.py`.

    python scripts/research/phase25_pair_locality.py <out.json>

* `P6.4_gap1.0`: 6.4 mm cell, 5.33 mm square (mesh-aligned stand-in for the 1.0 mm-gap pair), N = 5 (D = 1.12 lambda).
* `P8.0_gap1.14`: 8.0 mm cell, 6.857 mm square (24 cells of the 0.2857 mm mesh; the 6.85 mm recommended square),
  N = 4 (D = 32 mm, the same 1.12 lambda, orders (+-1,0),(0,+-1) propagate).
The patterns are crops of the hex4 1-bit designs (blocks of ports 1 and 3) and a seeded random map; they are
patterns, not designs for these pairs. The empty cell and the square are re-solved at each mesh.
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

import phase25_forward_validation as fv  # noqa: E402

JOBS = (
    ("P6.4_gap1.0", 6.4, 12, 10, 5),
    ("P6.4_gap1.0", 6.4, 24, 20, 5),
    ("P8.0_gap1.14", 8.0, 14, 12, 4),
    ("P8.0_gap1.14", 8.0, 28, 24, 4),
)

if __name__ == "__main__":
    res = fv._synthesize()
    out = {"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED"]}
    for name, period, mesh, cells, n in JOBS:
        out.setdefault(name, {})[f"mesh{mesh}_N{n}"] = fv.part_two(
            res, mesh, n, period_mm=period, patch_cells=cells, tag=f"{name} "
        )
        Path(sys.argv[1]).write_text(json.dumps(out, indent=1))
