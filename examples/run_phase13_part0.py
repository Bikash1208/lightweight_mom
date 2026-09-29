#!/usr/bin/env python3
"""Phase 13 Part 0: carry-overs from PHASE12D_ADDENDUM.md. Writes
``results/phase13_part0/part0_carry_overs.json``. Expensive (several real 7-point Richardson
mesh-convergence families plus real HFSS Touchstone ingestion) -- run once; see
``docs/phases/PHASE13_WEBAPP.md`` for the numbers already extracted from a run of this script."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.phase13.part0_run import run_part0, write_part0_json

if __name__ == "__main__":
    result = run_part0()
    out_dir = Path(__file__).resolve().parent.parent / "results" / "phase13_part0"
    path = write_part0_json(result, out_dir)
    print(f"Phase 13 Part 0 written to {path}")
    print(
        f"  G4 phase_inf: {result.convergence_g4['phase']['value_inf']:.2f} deg "
        f"(order {result.convergence_g4['phase']['order']:.2f})"
    )
    for r in result.resonance_restatement:
        print(
            f"  {r['label']}: {r['cells']} cells, dx={r['dx_mm']:.3f}mm, "
            f"scaled correction {r['scaled_correction_deg']:.2f} deg"
        )
    for v in result.order_power_verdicts:
        print(
            f"  N={v['n_cells']}: predicted intended {v['predicted_intended_fraction']:.3f} "
            f"vs measured {v['measured_intended_fraction']:.3f} "
            f"(diff {v['absolute_difference']:.3f})"
        )
