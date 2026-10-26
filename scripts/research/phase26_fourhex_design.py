"""Phase 26: designs and model predictions of the two 4-hexagon tensor-CP boards (shared 6.4 mm lattice).

    python scripts/research/phase26_fourhex_design.py <out_dir> [--quick] [--scale-s 1.0,1.7]

Board A (static): beams (theta 30; phi 45/135/225/315) on ports 1..4, each hexagon with its own phase reference.
Board B (Butler-ready): all four hexagons aimed at (30, 225), ONE global phase reference (the board centre), so that
equal-phase drive of the four ports is a single coherent aperture and the Butler phase steps steer it.

Per port the model gives: main pattern (cell currents of the design, grounded-slab Green function) and the
`pickup` pattern of the unconverted surface wave crossing the shared edges into the neighbours' cells (cell current
`A exp(-j phi_in) (cos chi, sin chi)`, chi = 2 a_cell - psi_in, alpha of the cell's own design: ASSUMPTION), with a
sensitivity scale `s` on the unconverted edge power (s = 1: the design march; s = 1.7: upper reading of the HFSS
|E| line of the single hexagon, `hfss_results/hex150_tensor_cp_single/efield_line_magE_phi225_1mm.csv`).
Gain anchoring: E_j is normalised to the fraction `LAUNCH eta_j (1 + r_j)` of the accepted power radiated by the cells,
+ the HFSS offset of the single hexagon (-0.68 dB). Labels: MODEL_ONLY, NOT_HFSS_VALIDATED, SOLVER_SUPERCELL.
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from src.hologram.tensor_cp_design import *  # noqa: E402,F401,F403  (the design flow lives in the module; this script is a wrapper)
from src.hologram.tensor_cp_design import FEEDS, design_fourhex  # noqa: E402


def main(out_dir, quick=False, scales=(1.0, 1.7)):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    res, cells, fields = design_fourhex(quick=quick, scales=scales)
    for board, c in cells.items():
        (out_dir / f"fourhex_{board}_cells.json").write_text(json.dumps({"board": board, "cells": c}))
    (out_dir / "fourhex_design.json").write_text(json.dumps(res, indent=1, default=float))
    th = np.arange(0.0, 90.0 + 1e-9, 3.0 if quick else 1.0)
    ph = np.arange(0.0, 360.0, 6.0 if quick else 1.0)
    np.savez_compressed(out_dir / "_fields_scratch.npz", theta=th, phi=ph,
                        **{f"{b}_{s}_{j}_{c}": fields[(b, s)][j][i] for (b, s) in [(k[0], k[1]) for k in fields if len(k) == 2] for j in FEEDS for i, c in enumerate(("et", "ep"))})


if __name__ == "__main__":
    a = sys.argv[1:]
    sc_arg = tuple(float(v) for v in a[a.index("--scale-s") + 1].split(",")) if "--scale-s" in a else (1.0, 1.7)
    main(a[0], quick="--quick" in a, scales=sc_arg)
