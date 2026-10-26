"""Phase 20 R5-1: compare the user's HFSS diode-tolerance run (CSV of `read_hfss_reim` + a C_off column) with the MoM prediction `phase20_r5_tolerance_mom.py` (mesh 90, library nominal correction added).
COMPARISON ONLY: nothing is calibrated. Reports per C_off, state and theta: HFSS - MoM r_RR phase and level at 10.5 GHz and over 10.0-11.0 GHz, and the sensitivity d r_RR / d C_off of both.

    python scripts/research/phase20_r5_tolerance_compare.py <hfss_tol.csv> <mom.npz> [--plane-offset-mm D] [--out report.json]
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.ris.calibrated_library import read_hfss_reim_tolerance  # noqa: E402
from src.ris.polarization import circular_reflection, to_structure_frame  # noqa: E402
from src.ris.state_library import StateLibrary  # noqa: E402


def db(z):
    return 20 * np.log10(max(abs(z), 1e-12))


def main(hfss_csv, mom_npz, plane_mm=0.0, out=None):
    rows = read_hfss_reim_tolerance(hfss_csv, plane_mm)
    z = np.load(mom_npz)
    lib = StateLibrary.load(ROOT / "results/phase20_r4/cp2_library_v1.npz")
    meshes = list(z["meshes"])
    im = len(meshes) - 1  # finest mesh
    fz, thz, cz = z["f_ghz"], z["theta_deg"], z["coffs_pf"]
    report = {"labels": ["MODEL_ONLY", "HFSS_COMPARISON_ONLY_NOT_CALIBRATED"], "mom_mesh": int(meshes[im]), "rows": []}
    for (coff, st, f, th), m in sorted(rows.items()):
        if not (10.0 - 1e-9 <= f <= 11.0 + 1e-9) or th not in thz:
            continue
        ic, jf, jt = int(np.argmin(abs(cz - coff))), int(np.argmin(abs(fz - f))), int(np.argmin(abs(thz - th)))
        if abs(cz[ic] - coff) > 1e-6 or abs(fz[jf] - f) > 1e-6:
            continue
        corr = (
            lib.structure_circular(f, th, 0.0, st, key="c_corr")[0]
            - lib.structure_circular(f, th, 0.0, st, key="c_raw")[0]
        )
        mom = z["c_structure"][im, ic, jf, jt, st - 1] + corr
        hf = to_structure_frame(circular_reflection(m), 0.0)
        report["rows"].append(
            {
                "coff_pF": coff, "state": st, "f_ghz": f, "theta_deg": th,
                "rr_phase_diff_deg": float(np.degrees(np.angle(hf[0, 0] / mom[0, 0]))),
                "rr_level_diff_db": db(hf[0, 0]) - db(mom[0, 0]),
                "opp_hfss_db": max(db(hf[0, 1]), db(hf[1, 0])), "opp_mom_db": max(db(mom[0, 1]), db(mom[1, 0])),
            }
        )  # fmt: skip
    r = [x for x in report["rows"] if abs(x["f_ghz"] - 10.5) < 1e-9]
    if r:
        ph = np.array([abs(x["rr_phase_diff_deg"]) for x in r])
        report["at_10p5_ghz"] = {
            "n": len(r),
            "median_abs_phase_diff_deg": float(np.median(ph)),
            "max_abs_phase_diff_deg": float(ph.max()),
        }
    if out:
        Path(out).write_text(json.dumps(report, indent=1))
    print(json.dumps({k: v for k, v in report.items() if k != "rows"}, indent=1))
    return report


if __name__ == "__main__":
    a = sys.argv
    main(
        a[1],
        a[2],
        float(a[a.index("--plane-offset-mm") + 1]) if "--plane-offset-mm" in a else 0.0,
        a[a.index("--out") + 1] if "--out" in a else None,
    )
