"""Score the one B-ALL-SHORT HFSS state against saved v3.1 references."""
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from full_board_gate_v31_extract import ffd
from sidequest_sll_def import analyze


def gain_at(theta, phi, intensity, peak, target_theta, target_phi):
    i = int(np.flatnonzero(np.isclose(theta, target_theta))[0])
    j = int(np.flatnonzero(np.isclose(phi, target_phi))[0])
    ip = int(np.flatnonzero(np.isclose(theta, peak["peak_theta_deg"]))[0])
    jp = int(np.flatnonzero(np.isclose(phi, peak["peak_phi_deg"]))[0])
    return float(peak["peak_gain_db"] + 10 * np.log10(intensity[i, j] / intensity[ip, jp]))


def main():
    both_path, off_path, short_path, ffd_path, geometry_path, output_path = map(Path, sys.argv[1:7])
    both = json.loads(both_path.read_text())["states"]["BOTH"]
    off = json.loads(off_path.read_text())["b_all_off_ideal_258"]
    short = json.loads(short_path.read_text())
    geometry = json.loads(geometry_path.read_text())
    assert geometry["centre_short_strip_cells"] == 250
    assert geometry["centre_short_strip_width_mm"] == 0.5
    assert geometry["resulting_slot_lengths_mm"] == [4.11, 4.11]
    assert geometry["blocked_vias"] == 0 and geometry["gate_sites"] == 16
    theta, phi, et, ep = ffd(ffd_path)
    plus = np.abs(et + 1j * ep) ** 2
    minus = np.abs(et - 1j * ep) ** 2
    rh, lh = (plus, minus) if short["ffd_rhcp_convention"] == "et_plus_j_ep" else (minus, plus)
    for hand, power in (("RHCP", rh), ("LHCP", lh)):
        peak = short[hand]
        pattern = analyze(theta, phi, power, peak=(peak["peak_theta_deg"], peak["peak_phi_deg"]))
        peak["sll_first_null_db"] = pattern["SLL_db"]
        peak["sll_direction_deg"] = pattern["SLL_dir"]
        peak["at_32_phi0_db"] = gain_at(theta, phi, power, peak, 32, 0)
        peak["at_36_phi180_db"] = gain_at(theta, phi, power, peak, 36, 180)
    gates = {
        "global_lhcp_le_minus1_dbi": short["LHCP"]["peak_gain_db"] <= -1,
        "lhcp_36_180_le_minus1_dbi": short["LHCP"]["at_36_phi180_db"] <= -1,
        "rhcp_32_0_ge_8p32_dbi": short["RHCP"]["at_32_phi0_db"] >= 8.32,
    }
    result = {"label": "HFSS_2026_R1_V31_B_ALL_SHORT", "frequency_ghz": 5.68,
              "changed_cells": 250, "strip_width_mm": 0.5, "slot_lengths_mm": [4.11, 4.11],
              "all_gates_pass": all(gates.values()), "gates": gates,
              "states": {"BOTH": both, "B_ALL_OFF_IDEAL": off, "B_ALL_SHORT": short}}
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2))
    print("B_ALL_SHORT", short["LHCP"]["at_36_phi180_db"], short["LHCP"]["peak_gain_db"],
          short["RHCP"]["at_32_phi0_db"], gates)


if __name__ == "__main__":
    main()
