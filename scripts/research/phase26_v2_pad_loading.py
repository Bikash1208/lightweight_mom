"""Phase 26 V2 launch pad: mean loading Y(D) of unslotted circles and of the hologram's slotted cells.

    python scripts/research/phase26_v2_pad_loading.py <out.json>

Mean sheet admittance `Ybar = Im tr(Y_s)/2` [mS] of the Bloch sheet tensor (periodic MoM, mesh 24, pitch 6.4 mm,
RT5880 1.575 mm) of (a) unslotted circles D = 3.5 .. 6.2 mm at kt = 1.10 k0 (the pad's own bound wave is a little
tighter than the hologram's 1.06 k0) and (b) each slot size of the hologram at 8 slot angles, kt = 1.06 k0.
SOLVER_INTERNAL, MODEL_ONLY. Used to grade the pad diameter so that Ybar falls from the high value at the pin to the
local mean of the hologram.
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import phase26_ramp_locality as rl  # noqa: E402
from src.basis_functions import MeshSettings  # noqa: E402
from src.dispersion import tensor_sheet as ts  # noqa: E402
from src.dispersion.slotted_cells import SlottedCircularPatch  # noqa: E402
from src.geometry import UnitCell  # noqa: E402
from src.materials import Metal  # noqa: E402
from src.periodic_bc import RectangularLattice  # noqa: E402

K0, F, SUB, P = rl.K0, rl.F, rl.SUB, rl.P


def ybar(geo, kt_over_k0, mesh_n=24, angle_deg=0.0):
    uc = UnitCell(RectangularLattice(P, P), geo, SUB, Metal(), backing="pec")
    y = ts.bloch_sheet_tensor(uc, F, MeshSettings(mesh_n, mesh_n, mesh_n), kt_over_k0 * K0, 0.0)
    return float(0.5 * (y[0, 0] + y[1, 1]).imag * 1e3)


def main(out):
    res = {"labels": ["SOLVER_INTERNAL", "MODEL_ONLY"], "pitch_mm": P * 1e3, "circles": [], "slotted": {}}
    for d in np.arange(3.5, 6.21, 0.1):
        v = ybar(SlottedCircularPatch(d * 1e-3, 0.0, 0.0), 1.10)
        res["circles"].append({"D_mm": round(float(d), 2), "ybar_mS": v})
        print("D", round(d, 2), round(v, 2), flush=True)
    for L, W in ((4.4, 1.2), (4.6, 1.6), (4.8, 1.6), (4.8, 2.0)):
        vals = [ybar(SlottedCircularPatch(6.0e-3, L * 1e-3, W * 1e-3, np.radians(a)), 1.06) for a in np.arange(0, 180, 22.5)]
        res["slotted"][f"{L}x{W}"] = {"ybar_mS_by_angle": vals, "mean_mS": float(np.mean(vals))}
        print("slot", L, W, round(float(np.mean(vals)), 2), flush=True)
    Path(out).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
