"""Phase 20 R5-1: does the converged Floquet cutoff (F = 2n) close the gap between the MoM and HFSS diode-tolerance sensitivity?

Compares, at 10.5 GHz (and over 10.25-10.75 GHz), the CHANGE of the CP-2 state phases and ladder steps with the OFF capacitance (C_off 0.0225 / 0.030 pF versus the nominal 0.025 pF)
in the HFSS tolerance run and in the server's MoM prediction (`phase20_r5_tolerance_mom.npz`, meshes 60 and 90, F = 2n). Differences are taken within each data set, so the reference-plane
correction and the library's mesh correction cancel. The earlier R4 corner table (C_T +-20 %, minimum-cutoff-era statement) is quoted from `results/phase20_r4/tolerance_corners.json`.
COMPARISON ONLY.

    python scripts/research/phase20_r5_tolerance_gap.py <hfss_tol.csv> <mom.npz> <out.json>
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.ris.calibrated_library import read_hfss_reim_tolerance  # noqa: E402
from src.ris.polarization import circular_reflection, to_structure_frame  # noqa: E402

NOM = 0.025


def wrap(x):
    return (x + 180.0) % 360.0 - 180.0


def steps(rr):
    """Ladder steps 1->2, 2->3, 3->4, 4->1 [deg] from the four co-handed responses (wrapped to 0..360)."""
    ph = np.degrees(np.angle(rr))
    return np.array([(ph[(i + 1) % 4] - ph[i]) % 360.0 for i in range(4)])


def main(hfss_csv, mom_npz, out):
    h = read_hfss_reim_tolerance(hfss_csv, 0.0)
    z = np.load(mom_npz)
    f_arr, th_arr, cz, meshes = z["f_ghz"], z["theta_deg"], z["coffs_pf"], list(z["meshes"])

    def hf_rr(coff, f, th):
        return np.array([to_structure_frame(circular_reflection(h[(round(coff, 6), s, round(f, 6), th)]), 0.0)[0, 0] for s in range(1, 5)])

    def mom_rr(im, coff, f, th):
        ic, jf, jt = int(np.argmin(abs(cz - coff))), int(np.argmin(abs(f_arr - f))), int(np.argmin(abs(th_arr - th)))
        return z["c_structure"][im, ic, jf, jt, :, 0, 0]

    res = {"labels": ["MODEL_ONLY", "COMPARISON_ONLY"], "rows": []}
    for th in (0.0, 30.0):
        for f in (10.25, 10.5, 10.75):
            for coff in (0.0225, 0.030):
                rec = {"theta": th, "f_ghz": f, "coff_pF": coff}
                s_h0, s_h = steps(hf_rr(NOM, f, th)), steps(hf_rr(coff, f, th))
                rec["hfss_step_change_deg"] = [float(wrap(a - b)) for a, b in zip(s_h, s_h0, strict=True)]
                rec["hfss_state_phase_change_deg"] = [float(wrap(a)) for a in np.degrees(np.angle(hf_rr(coff, f, th) / hf_rr(NOM, f, th)))]
                for im, mesh in enumerate(meshes):
                    s_m0, s_m = steps(mom_rr(im, NOM, f, th)), steps(mom_rr(im, coff, f, th))
                    rec[f"mom_mesh{int(mesh)}_step_change_deg"] = [float(wrap(a - b)) for a, b in zip(s_m, s_m0, strict=True)]
                    rec[f"mom_mesh{int(mesh)}_state_phase_change_deg"] = [float(wrap(a)) for a in np.degrees(np.angle(mom_rr(im, coff, f, th) / mom_rr(im, NOM, f, th)))]
                    rec[f"mom_mesh{int(mesh)}_max_abs_step_change_deg"] = float(np.abs(rec[f"mom_mesh{int(mesh)}_step_change_deg"]).max())
                rec["hfss_max_abs_step_change_deg"] = float(np.abs(rec["hfss_step_change_deg"]).max())
                res["rows"].append(rec)
    res["r4_corner_table_note"] = "results/phase20_r4/tolerance_corners.json (MoM, C_T +-20 %, R4 mesh/cutoff): the 17-28 deg step errors quoted in R3/R4"
    Path(out).write_text(json.dumps(res, indent=1))
    for r in res["rows"]:
        if r["f_ghz"] == 10.5:
            print(r["theta"], r["coff_pF"], "HFSS", round(r["hfss_max_abs_step_change_deg"], 1), "| MoM", {int(m): round(r[f"mom_mesh{int(m)}_max_abs_step_change_deg"], 1) for m in meshes})


if __name__ == "__main__":
    main(*sys.argv[1:4])
