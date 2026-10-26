"""Phase 20 R3 item 2: how the three STARTING cells (LP-1, LP-2, CP-2; `src/ris/cells.py`) respond at 10 GHz. NOT designs; labelled MODEL_ONLY,
PROVISIONAL_DATASHEET. Writes `results/phase20_r3/cells_start.json`.

    python scripts/research/phase20_r3_cells_start.py [out_dir]
"""

import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.basis_functions import MeshSettings  # noqa: E402
from src.excitation import ReferencePlanes  # noqa: E402
from src.ris.cells import Cp2Cell, Lp1Cell, Lp2Cell, ideal_kit, madp_kit, smp1345_kit  # noqa: E402
from src.ris.polarization import L, R, circular_reflection  # noqa: E402
from src.sparameters import compute_s_parameters_oblique  # noqa: E402

PLANES = ReferencePlanes(0.0, None, 1.575e-3)
MESH = {"lp1": 48, "lp2": 48, "cp2": 54}


def solve(cell, state, kit, theta, f=10e9):
    n = MESH[cell.kind]
    uc = cell.unit_cell(state, kit)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return compute_s_parameters_oblique(uc, f, MeshSettings(n, n, n // 2 + 1), PLANES, math.radians(theta), 0.0)


def _c(z):
    return [float(abs(z)), float(math.degrees(np.angle(z)))]


def run(out_dir):
    out = {
        "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "PROVISIONAL_DATASHEET", "VIA_OMITTED_APPROXIMATION"],
        "f_ghz": 10.0,
        "cells": {},
    }
    for cell in (Lp1Cell(), Lp2Cell(), Cp2Cell()):
        kits = (smp1345_kit(), madp_kit()) + ((ideal_kit(),) if cell.kind == "cp2" else ())
        for kit in kits:
            for th in (0.0, 30.0):
                for st in cell.states():
                    s = solve(cell, st, kit, th)
                    row = {
                        "te": _c(s.s11[0, 0]),
                        "tm": _c(s.s11[1, 1]),
                        "te_tm_cross": float(abs(s.s11[0, 1])),
                        "power_out": [float(v) for v in s.power_out_per_incident()[0]],
                    }
                    c = circular_reflection(s.s11)
                    row["circular"] = {
                        "r_RR": _c(c[R, R]),
                        "r_LL": _c(c[L, L]),
                        "r_LR": _c(c[L, R]),
                        "r_RL": _c(c[R, L]),
                    }
                    out["cells"].setdefault(cell.kind, {}).setdefault(kit.name, {}).setdefault(f"theta{th:g}", {})[
                        "".join("1" if x == "ON" else "0" for x in st)
                    ] = row
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    (Path(out_dir) / "cells_start.json").write_text(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else ROOT / "results" / "phase20_r3")
