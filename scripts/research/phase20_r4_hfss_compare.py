"""Phase 20 R4 item 6: compare the user's CP-2 HFSS export with the state library. COMPARISON ONLY: nothing is calibrated here (calibration needs the user's approval and
`src.ris.hfss_calibration.apply_calibration(..., approved_by_user=True)`).

    python scripts/research/phase20_r4_hfss_compare.py <hfss_export.csv> [--tm-shift-180] [--library path.npz] [--out report.json]

The CSV is in the long format of `hfss_exports/phase20_r3_cp_cell/mom_reference.csv` with the columns `state, f_ghz, theta_deg, phi_deg, quantity, hfss_mag_dB, hfss_phase_deg`
(quantity: S_TE_TE, S_TM_TM, S_TE_TM, S_TM_TE). Run the bare-slab convention check of the package README first and pass --tm-shift-180 if it shows the 180 degree TM shift.
"""

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.ris.hfss_calibration import compare, proposed_offsets, read_hfss_rows  # noqa: E402
from src.ris.state_library import StateLibrary  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--tm-shift-180", action="store_true")
    ap.add_argument("--library", default=str(ROOT / "results/phase20_r4/cp2_library_v1.npz"))
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    lib = StateLibrary.load(a.library)
    comps, summary = compare(lib, read_hfss_rows(a.csv, tm_shift_180=a.tm_shift_180))
    report = {
        "labels": ["MODEL_ONLY", "HFSS_COMPARISON_ONLY_NOT_CALIBRATED"],
        "summary": summary,
        "candidate_offsets_not_applied": proposed_offsets(comps),
        "points": [asdict(c) for c in comps],
    }
    if a.out:
        Path(a.out).write_text(json.dumps(report, indent=1))
    print(
        json.dumps(
            {"summary": summary, "candidate_offsets_not_applied": report["candidate_offsets_not_applied"]}, indent=1
        )
    )


if __name__ == "__main__":
    main()
