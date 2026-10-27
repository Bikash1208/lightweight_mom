"""Phase 28 item 3: summary JSON of the screening searches (what was tried, the best numbers, counts), from the raw search result files.

    python scripts/research/phase28_search_summary.py <out.json> <family:period_mm:file> [...]

`family` is `octagon_stub` (ris.cells.Cp2Cell scaled), `sector_screen`, `sector_balanced` or `sector_refined` (SectorHubCell). Only samples with a computed tensor count as valid.
"""

import json
import sys
from pathlib import Path


def _score(r, key):
    s = r["states"]
    vals = [(s[k].get(key) or 0.0) for k in s]
    return min(vals)


def summarise(family, period, path):
    rows = json.loads(Path(path).read_text(encoding="utf-8"))
    ok = [r for r in rows if "states" in r]
    best_dx = max((_score(r, "dx_aniso_ohm") for r in ok), default=0.0)
    best_al = max((_score(r, "alpha_beam_np_per_m") for r in ok), default=0.0)
    out = {
        "family": family,
        "period_mm": period,
        "n_samples": len(rows),
        "n_valid": len(ok),
        "best_dx_own_ohm": best_dx,
        "best_alpha_np_per_m": best_al,
        "source": Path(path).name,
    }
    bal = [r for r in ok if r.get("balance") is not None and r["balance"] < 0.15]
    out["n_balanced"] = len(bal)
    out["best_balanced_alpha_np_per_m"] = max((r["score"] for r in bal), default=0.0)
    return out


if __name__ == "__main__":
    out = Path(sys.argv[1])
    fams = []
    for spec in sys.argv[2:]:
        fam, per, path = spec.split(":", 2)
        fams.append(summarise(fam, float(per), path))
    out.write_text(json.dumps({"families": fams, "notes": []}, indent=1), encoding="utf-8")
    print("wrote", out)
