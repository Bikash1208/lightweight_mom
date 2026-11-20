"""Score the three 2026 R1 full-board gate patterns from exported HFSS FFD files."""
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from full_board_gate_v31_extract import ffd
from sidequest_sll_def import analyze


def main():
    work = Path(sys.argv[1])
    out = Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    states = {}
    for state in ("BOTH", "A_ONLY", "B_ONLY"):
        metrics = json.loads((work / f"extract_{state}_001" / "metrics.json").read_text())
        theta, phi, et, ep = ffd(Path(metrics["ffd_path"]))
        plus = np.abs(et + 1j * ep) ** 2
        minus = np.abs(et - 1j * ep) ** 2
        rh, lh = (plus, minus) if metrics["ffd_rhcp_convention"] == "et_plus_j_ep" else (minus, plus)
        for hand, power in (("RHCP", rh), ("LHCP", lh)):
            peak = (metrics[hand]["peak_theta_deg"], metrics[hand]["peak_phi_deg"])
            pattern = analyze(theta, phi, power, peak=peak)
            metrics[hand]["sll_first_null_db"] = pattern["SLL_db"]
            metrics[hand]["sll_direction_deg"] = pattern["SLL_dir"]
            metrics[hand]["hpbw_mean_deg"] = pattern["HPBW_mean_deg"]
            metrics[hand]["wide_angle_relative_db"] = (
                metrics[hand]["wide_angle_max_gain_db_theta60_to90"] - metrics[hand]["peak_gain_db"])
        states[state] = metrics
    both = states["BOTH"]
    for state in ("A_ONLY", "B_ONLY"):
        row = states[state]
        for hand in ("RHCP", "LHCP"):
            row[hand]["peak_change_vs_both_db"] = row[hand]["peak_gain_db"] - both[hand]["peak_gain_db"]
        suppressed, survivor = ("LHCP", "RHCP") if state == "A_ONLY" else ("RHCP", "LHCP")
        row["suppressed_hand"] = suppressed
        row["survivor_hand"] = survivor
        row["suppressed_drop_db"] = -row[suppressed]["peak_change_vs_both_db"]
        row["survivor_change_db"] = row[survivor]["peak_change_vs_both_db"]
        row["gate_pass"] = row["suppressed_drop_db"] >= 10 and row["survivor_change_db"] >= 0
    result = {"label": "HFSS_2026_R1_V31", "frequency_ghz": 5.68, "states": states,
              "overall_gate_pass": states["A_ONLY"]["gate_pass"] and states["B_ONLY"]["gate_pass"],
              "sll_definition": "Highest same-hand lobe outside direction-dependent first null; if no null, exclude 2x mean HPBW. Above horizon only.",
              "wide_angle_definition": "Maximum same-hand gain for 60 <= theta <= 90 deg."}
    (out / "metrics_summary.json").write_text(json.dumps(result, indent=2))
    print("A_ONLY", states["A_ONLY"]["suppressed_drop_db"], states["A_ONLY"]["survivor_change_db"])
    print("B_ONLY", states["B_ONLY"]["suppressed_drop_db"], states["B_ONLY"]["survivor_change_db"])


if __name__ == "__main__":
    main()
