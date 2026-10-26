"""Phase 26 A: in-context contrast of the 1-bit pair versus the stripe period (supercell MoM, plane wave at
normal incidence), the evidence behind `context_contrast_ratio` of `phase26_hex150_leakage.py`.

    python scripts/research/phase26_context_contrast.py <out.json>

Straight stripes of the design pair (square | empty, 5.867 mm square, mesh 24, 10.54 GHz, grounded RT/duroid
5880 1.575 mm) with 1, 2 and 3 cells per half-period and duty 1/3, 1/2, 2/3 (every real period of the
hex150 hologram is 18.6-54 mm = 2.9-8.4 cells; the Phase 25 stripes were the shortest, 2 cells). The
homogenised sheet admittance of each supercell comes from its specular order; the in-context contrast
ratio is `kappa = (Y_MoM - Y_empty) / (duty (Y_square - Y_empty))` with Y_square from the isolated cell.
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

import phase25_block_locality as bl  # noqa: E402
from src.dispersion.sheet_model import sheet_admittance_from_gamma  # noqa: E402
from src.materials import Substrate  # noqa: E402

CASES = {
    "stripes_1of2_period2": np.array([[1, 0]] * 2),
    "stripes_1of3_period3": np.array([[1, 0, 0]] * 3),
    "stripes_2of3_period3": np.array([[1, 1, 0]] * 3),
    "stripes_2of4_period4": np.array([[1, 1, 0, 0]] * 4),
}


def main(out):
    res = bl.run(24, cases=CASES, n=None, tag="ctx ") if False else None
    # `run` assumes one n for every case: drive it per case
    result = {"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED"], "cases": {}}
    sub = Substrate("user", bl.fv.EPS_R, bl.fv.TAND, bl.fv.H_MM * 1e-3)
    for name, lev in CASES.items():
        n = lev.shape[0]
        r = bl.run(24, cases={name: lev}, n=n, tag="ctx ")
        g0 = r["isolated_empty"]
        g1 = r["isolated_square"]
        gm = r["cases"][name]["orders"]["0,0"]["mom"]
        c = lambda m, p: m * np.exp(1j * np.radians(p))  # noqa: E731
        f = bl.F0
        y0, y1, ym = (sheet_admittance_from_gamma(c(*g), sub, f) for g in (g0, g1, gm))
        duty = float(lev[0].mean())
        kappa = float(np.real((ym - y0) / (duty * (y1 - y0))))
        result["cases"][name] = {
            "period_cells": n,
            "duty": duty,
            "specular_mom_mag_phase": gm,
            "kappa": kappa,
            "n_unknowns": r["cases"][name]["n_unknowns"],
            "seconds": r["cases"][name]["seconds"],
        }
        print(name, "kappa", round(kappa, 3), flush=True)
        Path(out).write_text(json.dumps(result, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
