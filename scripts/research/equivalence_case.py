"""Cross-machine equivalence case (infra brief rule 5): the baseline CP-2 cell (`results/phase20_r3/opt/cp2_madp_final.json`), MADP nominal diodes, mesh 60, Floquet cutoff F = 2n = 120,
the four one-hot states, f = 10.5 GHz, (theta, phi) = (0, 0) and (30, 0). Writes the complex TE/TM S11 of every state with a provenance record.

    python scripts/research/equivalence_case.py <out.json>
"""

import json
import math
import sys
import time
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

from phase20_r3_finalist import build  # noqa: E402
from provenance import finish, provenance  # noqa: E402

from src.basis_functions import MeshSettings  # noqa: E402
from src.excitation import ReferencePlanes  # noqa: E402
from src.ris.cells import madp_kit  # noqa: E402
from src.sparameters import compute_s_parameters_oblique_over_loads  # noqa: E402

MESH, FREQ_GHZ, ANGLES = 60, 10.5, ((0.0, 0.0), (30.0, 0.0))


def main(out):
    d = json.loads((ROOT / "results/phase20_r3/opt/cp2_madp_final.json").read_text())
    cell = build("cp2", d["params_mm"], d["period_mm"])
    kit = madp_kit()
    planes = ReferencePlanes(0.0, None, 1.575e-3)
    sets = [cell.unit_cell(st, kit).loads for st in cell.states()]
    base = cell.unit_cell(cell.states()[0], kit)
    rec = provenance(
        {
            "case": "cp2_baseline_mesh60_F2n_10p5GHz",
            "params": {"mesh": MESH, "floquet": 2 * MESH, "f_ghz": FREQ_GHZ, "angles_deg": ANGLES},
        },
        workers=1,
    )
    cases = []
    for th, ph in ANGLES:
        t0 = time.monotonic()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = compute_s_parameters_oblique_over_loads(
                base,
                FREQ_GHZ * 1e9,
                MeshSettings(MESH, MESH, 2 * MESH),
                planes,
                math.radians(th),
                math.radians(ph),
                sets,
                warn_grating=False,
            )
        for k, r in enumerate(res):
            cases.append(
                {
                    "state": k + 1,
                    "theta_deg": th,
                    "phi_deg": ph,
                    "s11": [[[float(z.real), float(z.imag)] for z in row] for row in r.s11],
                    "seconds": time.monotonic() - t0,
                }
            )
    Path(out).write_text(json.dumps({"provenance": finish(rec), "cases": cases}, indent=1))
    print("wrote", out, "cases", len(cases))


if __name__ == "__main__":
    main(sys.argv[1])
