"""Phase 20 R5-3: beam table, gain bandwidth and diode-tolerance sensitivity of the plane-wave CP-2 RIS (MODEL_ONLY).

    python scripts/research/phase20_r5_beam_table.py results/phase20_r5/beam_table.json

Normal incidence, 20 x 20 cells at 15 mm, design frequency 10.5 GHz, steering theta_s = 0..50 deg (+ 60 deg marked) in phi_s = 0 and 90 deg. Variants: plain,
dither (seed 0 in the table, spread over seeds 0-4 reported), tilt, refined; ideal continuous / ideal 2-bit references. Bandwidth: fixed plain and refined
codings swept over 10.0-11.0 GHz (step 0.05). Sensitivity: the same codings with C_off = 0.0225 / 0.025 / 0.030 pF from the user's HFSS tolerance run.
All numbers are MODEL_ONLY with the library + run-2 HFSS delta; nothing is tuned to the C_off cases.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from scripts.research.provenance import finish, provenance  # noqa: E402
from src.ris import coding as cd  # noqa: E402
from src.ris.calibrated_library import CalibratedLibrary, ToleranceBand, read_hfss_reim_tolerance  # noqa: E402
from src.ris.state_library import OutOfDomain  # noqa: E402

F0 = 10.5
ANGLES = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60]
PLANES = [0.0, 90.0]
VARIANTS = ["plain", "dither", "tilt", "refined"]
FREQS = [round(10.0 + 0.05 * i, 2) for i in range(21)]
COFFS = [0.0225, 0.025, 0.030]
TOL_CSV = ROOT / "hfss_results" / "phase20_r3_cp_cell" / "hfss_cp2_tolerance_Coff_0225_025_030.csv"
TOL_PLANE_MM = -0.374  # the run-2 residual plane correction used for the delta table (same HFSS project)
STEP = 0.005


def row(m, g_t):
    return {
        "gain_target_dbi": round(g_t, 2),
        "peak_gain_dbi": round(m.peak_gain_dbi, 2),
        "peak_theta_deg": round(m.peak_theta_deg, 2),
        "peak_phi_deg": round(m.peak_phi_deg, 1),
        "pointing_error_deg": round(m.pointing_error_deg, 2),
        "directivity_dbi": round(m.directivity_dbi, 2),
        "radiated_over_incident": round(m.radiated_over_incident, 3),
        "sll_db": round(m.sll_db, 1),
        "sll_at_deg": [round(m.sll_theta_deg, 1), round(m.sll_phi_deg, 1)],
        "cross_peak_db": round(m.cross_peak_db, 1),
        "cross_peak_at_deg": [round(m.cross_peak_theta_deg, 1), round(m.cross_peak_phi_deg, 1)],
        "cross_at_specular_db": round(m.cross_at_specular_db, 1),
        "axial_ratio_db": round(m.axial_ratio_db, 2),
    }


def _num(v):
    return None if not np.isfinite(v) else round(float(v), 2)


def main(out_path):
    t0 = time.time()
    cal = CalibratedLibrary.load()
    cal.tolerance = ToleranceBand.from_hfss(read_hfss_reim_tolerance(TOL_CSV, plane_offset_mm=TOL_PLANE_MM))
    result = {"provenance": provenance(), "design_f_ghz": F0, "labels": list(cd.ASSUMPTIONS), "table": [], "bandwidth": [], "sensitivity": []}
    base_gain, codings = {}, {}
    for phi in PLANES:
        for th in ANGLES:
            note = ""
            if th == 60:
                note = f"near band centre only: first grating lobe enters at {cd.grating_onset_deg(F0):.1f} deg at 10.5 GHz; no band margin above 10.54 GHz"
            entry = {"theta_s": th, "phi_s": phi, "note": note, "variants": {}}
            for k, v in cd.ideal_references(th, phi, F0, step=STEP).items():
                entry["variants"][k] = row(v, v.peak_gain_dbi)
            for var in VARIANTS:
                try:
                    c = cd.design(cal, var, F0, th, phi)
                    m, gt = cd.evaluate(cal, c, step=STEP)
                except OutOfDomain as e:
                    entry["variants"][var] = {"status": "NOT_AVAILABLE_FROM_CURRENT_DATA", "reason": str(e)}
                    continue
                entry["variants"][var] = row(m, gt)
                codings[(var, th, phi)] = c
                if var == "dither":
                    gs = [cd.evaluate(cal, cd.design(cal, "dither", F0, th, phi, seed=s), step=STEP)[1] for s in range(5)]
                    entry["variants"][var]["gain_target_over_seeds_dbi"] = [round(min(gs), 2), round(float(np.mean(gs)), 2), round(max(gs), 2)]
                if var == "refined":
                    entry["variants"][var]["refine_log"] = c.log
            if th == 0:
                base_gain[phi] = entry["variants"]["plain"]["gain_target_dbi"]
            for var in entry["variants"].values():
                if "gain_target_dbi" in var:
                    var["scan_loss_db"] = round(base_gain[phi] - var["gain_target_dbi"], 2)
            result["table"].append(entry)
            print(f"table phi={phi} theta={th} done {time.time() - t0:.0f}s", flush=True)

    for th in (0, 30, 50):
        for var in ("plain", "refined"):
            c = codings[(var, th, 0.0)]
            gt, gp, tp, ok = [], [], [], []
            for f in FREQS:
                try:
                    m, g = cd.evaluate(cal, c, f_ghz=f, step=STEP)
                    gt.append(g), gp.append(m.peak_gain_dbi), tp.append(m.peak_theta_deg), ok.append(True)
                except OutOfDomain:
                    gt.append(np.nan), gp.append(np.nan), tp.append(np.nan), ok.append(False)
            good = np.array(ok)
            fr = np.array(FREQS)[good]
            result["bandwidth"].append({
                "theta_s": th, "phi_s": 0.0, "variant": var, "freqs": FREQS,
                "gain_at_design_direction_dbi": [_num(g) for g in gt],
                "peak_gain_dbi": [_num(g) for g in gp],
                "peak_theta_deg": [_num(g) for g in tp],
                "grating_onset_deg": [round(cd.grating_onset_deg(f), 1) for f in FREQS],
                "bandwidth_design_direction": cd.gain_bandwidth(fr, np.array(gt)[good]) if good.any() else {},
                "bandwidth_peak": cd.gain_bandwidth(fr, np.array(gp)[good]) if good.any() else {},
                "out_of_domain_freqs": [f for f, o in zip(FREQS, ok) if not o],
            })
    print(f"bandwidth done {time.time() - t0:.0f}s", flush=True)

    for phi in PLANES:
        for th in (0, 30, 50):
            for var in ("plain", "refined"):
                c = codings[(var, th, phi)]
                rows_ = {str(coff): row(*(lambda r: (r[0], r[1]))(cd.evaluate(cal, c, coff=coff, step=STEP))) for coff in COFFS}
                g_all = [r["gain_target_dbi"] for r in rows_.values()]
                result["sensitivity"].append({
                    "theta_s": th, "phi_s": phi, "variant": var, "by_coff_pf": rows_,
                    "gain_target_spread_db": round(max(g_all) - min(g_all), 2),
                    "label": "MODEL_ONLY with HFSS C_off differences added to the calibrated states (HFSS tolerance run, theta 0 data)",
                })
    result["runtime_s"] = round(time.time() - t0, 1)
    finish(result["provenance"])
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    Path(out_path).write_text(json.dumps(result, indent=1, default=float))
    print("wrote", out_path, result["runtime_s"], "s")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "results/phase20_r5/beam_table.json")
