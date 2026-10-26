"""Phase 25 follow-up 3b: the larger-block half of `hfss_exports/cell_pair_locality_check/`.

    python scripts/research/phase25_locality_blocks.py [package_dir]

The 1x2 / 2x2 supercells in the package are sub-wavelength (period < lambda): only the specular order
propagates, so energy conservation forces |S11| ~ 1 and only its PHASE is informative. The blocks here are the
ones with propagating Floquet orders (D = 1.12 lambda at 10.54 GHz) whose per-order power fractions the MoM
check (`results/phase25/forward_validation.json`, `pair_locality.json`) compares with the isolated-cell
model: 5x5 blocks of the current pair and of the locality-best candidate, 4x4 of the separation-best (8 mm) candidate; STEP for each pattern + one reference CSV.
"""

import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.export.array_export import write_hfss_step  # noqa: E402

THICK_MM = 1.575


def main(pkg):
    pkg = Path(pkg)
    fv = json.loads((ROOT / "results/phase25/forward_validation.json").read_text())
    pl = json.loads((ROOT / "results/phase25/pair_locality.json").read_text())
    sources = {
        "current_P6.4_L5.85": (6.4, 5.85, fv["part_two_plane_wave_supercell_mom"]["mesh24"]),
        "candidate_locality_best_P6.4_L5.33": (6.4, 5.33, pl["P6.4_gap1.0"]["mesh24_N5"]),
        "candidate_separation_best_P8.0_L6.85": (8.0, 6.85, pl["P8.0_gap1.14"]["mesh28_N4"]),
    }
    rows = []
    for pair, (period, side, m) in sources.items():
        for name, c in m["cases"].items():
            if name == "uniform_patch":
                continue
            lev = np.array(c["states"])
            n = lev.shape[0]
            grid = [[(side * 1e-3, side * 1e-3) if lev[i, j] else (0.0, 0.0) for j in range(n)] for i in range(n)]
            fname = f"block{n}x{n}_{name}.step"
            write_hfss_step(grid, period * 1e-3, period * 1e-3, THICK_MM * 1e-3, pkg / pair / fname)
            for order, r in c["propagating_orders"].items():
                rows.append(
                    {
                        "pair": pair,
                        "block_step": fname,
                        "states_row_major_1_is_square": "".join(str(int(v)) for v in lev.ravel()),
                        "order_m_n": order,
                        "frequency_GHz": 10.54,
                        "mom_mesh_cells_per_period": m.get("mesh_cells_per_6p4mm_cell"),
                        "mom_square_side_mm": m["patch_side_mm_represented"],
                        "mom_power_fraction": r["mom_power_fraction"],
                        "isolated_cell_model_power_fraction_renormalised": r["forward_power_fraction_renormalised"],
                        "error_percentage_points": r["error_percentage_points"],
                    }
                )
    with open(pkg / "mom_order_power_fractions.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.5f}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"wrote {len(rows)} rows, STEP blocks in {pkg}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else str(ROOT / "hfss_exports" / "cell_pair_locality_check"))
