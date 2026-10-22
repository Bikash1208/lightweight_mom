"""Phase 20 R4 item 5: diode tolerance corners of the library on a reduced grid: f = 9.5, 10, 10.5 GHz; theta = 0, 30 deg (phi = 0 at theta 0; phi = 0, 45, 90 at 30 deg);
MADP nominal (L_s = 0), C_T -20 % / +20 %, L_s = 0.1 and 0.2 nH (ASSUMED_SENSITIVITY_NOT_DEVICE_DATA). Each corner is solved at the grid mesh (60) and corrected with the
library's NOMINAL mesh correction at the same (f, theta) (an approximation: the correction map is that of the nominal diode). MODEL_ONLY.

    python scripts/research/phase20_r4_corners.py        # writes results/phase20_r4/tolerance_corners.json
"""

import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

from phase20_r3_finalist import build  # noqa: E402
from phase20_r4_metrics import point_metrics  # noqa: E402

from src.basis_functions import MeshSettings  # noqa: E402
from src.excitation import ReferencePlanes  # noqa: E402
from src.ris.cells import DiodeKit, madp_kit  # noqa: E402
from src.ris.polarization import circular_reflection, to_structure_frame  # noqa: E402
from src.ris.state_library import StateLibrary  # noqa: E402
from src.sparameters import compute_s_parameters_oblique_over_loads  # noqa: E402

FREQS = (9.5, 10.0, 10.5)
POINTS = [(0.0, 0.0)] + [(30.0, p) for p in (0.0, 45.0, 90.0)]
PL = ReferencePlanes(0.0, None, 1.575e-3)


def kits():
    nom = madp_kit()
    out = {"nominal (L_s = 0)": nom}
    out["C_T -20 %"] = DiodeKit(nom.name, nom.on, nom.off.scaled(C_T=0.8))
    out["C_T +20 %"] = DiodeKit(nom.name, nom.on, nom.off.scaled(C_T=1.2))
    out["L_s = 0.1 nH (assumed)"] = madp_kit(0.1e-9)
    out["L_s = 0.2 nH (assumed)"] = madp_kit(0.2e-9)
    return out


def main():
    d = json.loads((ROOT / "results/phase20_r3/opt/cp2_madp_final.json").read_text())
    cell = build("cp2", d["params_mm"], d["period_mm"])
    lib = StateLibrary.load(ROOT / "results/phase20_r4/cp2_library_v1.npz")
    ks = kits()
    sets = [cell.unit_cell(st, kit).loads for kit in ks.values() for st in cell.states()]
    base = cell.unit_cell(cell.states()[0], madp_kit())
    rows = []
    for f in FREQS:
        for th, ph in POINTS:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res = compute_s_parameters_oblique_over_loads(
                    base,
                    f * 1e9,
                    MeshSettings(60, 60, 120),
                    PL,
                    math.radians(th),
                    math.radians(ph),
                    sets,
                    warn_grating=False,
                )
            # nominal-diode mesh correction at this (f, theta): library c_corr - c_raw at the nearest phi index (it is phi independent)
            ip = int(np.argmin(np.abs(lib.phi - ph)))
            i_f, i_t = int(np.argmin(np.abs(lib.f - f))), int(np.argmin(np.abs(lib.theta - th)))
            corr = lib.a["correction"][i_f, i_t, ip]  # (state, 2, 2)
            nom_c = None
            for ik, name in enumerate(ks):
                c = np.array(
                    [to_structure_frame(circular_reflection(res[ik * 4 + k].s11), np.radians(ph)) for k in range(4)]
                )
                c = c + corr
                m = point_metrics(c)
                rr_phase = np.degrees(np.angle(c[:, 0, 0]))
                if nom_c is None:
                    nom_c = rr_phase
                rows.append(
                    {
                        "f_ghz": f,
                        "theta": th,
                        "phi": ph,
                        "corner": name,
                        **m,
                        "rr_phase_deg": rr_phase.tolist(),
                        "rr_phase_shift_vs_nominal_deg": ((rr_phase - nom_c + 180.0) % 360.0 - 180.0).tolist(),
                        "rr_db": (20 * np.log10(np.abs(c[:, 0, 0]))).tolist(),
                    }
                )
    Path(ROOT / "results/phase20_r4/tolerance_corners.json").write_text(json.dumps(rows, indent=1))
    for r in rows:
        print(
            f"f {r['f_ghz']:g} th {r['theta']:g} phi {r['phi']:g} {r['corner']:<24} dev {r['dev']:5.1f} step {r['step']:5.1f} "
            f"amp_min {r['amp_min']:5.2f} spread {r['spread']:4.2f} opp {r['opp']:6.1f}"
        )


if __name__ == "__main__":
    main()
