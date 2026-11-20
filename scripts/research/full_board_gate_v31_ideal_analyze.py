"""Compare the one v3.1 A-only-IDEAL HFSS solve with the saved BOTH reference."""
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from full_board_gate_v31_extract import ffd
from sidequest_sll_def import analyze


def main():
    both_path, ideal_path, profile_path, manifest_path, output_path = map(Path, sys.argv[1:6])
    both = json.loads(both_path.read_text())["states"]["BOTH"]
    ideal = json.loads(ideal_path.read_text())
    manifest = json.loads(manifest_path.read_text())
    profile = profile_path.read_text(errors="replace")
    theta, phi, et, ep = ffd(Path(ideal["ffd_path"]))
    plus = np.abs(et + 1j * ep) ** 2
    minus = np.abs(et - 1j * ep) ** 2
    rh, lh = (plus, minus) if ideal["ffd_rhcp_convention"] == "et_plus_j_ep" else (minus, plus)
    for hand, power in (("RHCP", rh), ("LHCP", lh)):
        peak = (ideal[hand]["peak_theta_deg"], ideal[hand]["peak_phi_deg"])
        pattern = analyze(theta, phi, power, peak=peak)
        ideal[hand]["sll_first_null_db"] = pattern["SLL_db"]
        ideal[hand]["sll_direction_deg"] = pattern["SLL_dir"]
        ideal[hand]["peak_change_vs_both_db"] = ideal[hand]["peak_gain_db"] - both[hand]["peak_gain_db"]
    first = profile.split("Adaptive Pass 1", 1)[1].split("Adaptive Pass 2", 1)[0]
    memory = re.search(r"Matrix Solve\s+[^\n]*?\s([\d.]+)\s+G\s+Type:", first)
    assert memory, "FIRST_PASS_MEMORY_NOT_FOUND"
    assert manifest["ideal_plain_replacements"] == 234 and manifest["blocked_vias"] == 0
    assert len([r for r in manifest["rows"] if r["ideal_plain_replacement"]]) == 234
    result = {"label": "HFSS_2026_R1_V31_A_ONLY_IDEAL", "status": "COMPLETE", "frequency_ghz": 5.68,
              "changed_cells": 234, "gate_discs": 16, "vias": 0,
              "first_pass_hfss_memory_gb": float(memory.group(1)),
              "ideal": ideal,
              "both_reference": {"RHCP": both["RHCP"], "LHCP": both["LHCP"],
                                 "s11_db": both["s11_db"],
                                 "accepted_power_fraction": both["accepted_power_fraction"]},
              "B_peak_drop_db": -ideal["LHCP"]["peak_change_vs_both_db"],
              "A_peak_drop_db": -ideal["RHCP"]["peak_change_vs_both_db"],
              "B_design_direction_drop_db": both["LHCP"]["at_30_phi180_db"] - ideal["LHCP"]["at_30_phi180_db"],
              "A_design_direction_drop_db": both["RHCP"]["at_30_phi0_db"] - ideal["RHCP"]["at_30_phi0_db"]}
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2))
    print("B_drop", result["B_peak_drop_db"], "A_drop", result["A_peak_drop_db"],
          "first_memory", result["first_pass_hfss_memory_gb"])


if __name__ == "__main__":
    main()
