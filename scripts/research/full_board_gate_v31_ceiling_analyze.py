"""Score the v3.1 B-all-off ceiling and conditional all-plain HFSS floor."""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from full_board_gate_v31_extract import ffd
from sidequest_sll_def import analyze


def add_sll(metrics):
    theta, phi, et, ep = ffd(Path(metrics["ffd_path"]))
    plus = np.abs(et + 1j * ep) ** 2
    minus = np.abs(et - 1j * ep) ** 2
    rh, lh = (plus, minus) if metrics["ffd_rhcp_convention"] == "et_plus_j_ep" else (minus, plus)
    for hand, power in (("RHCP", rh), ("LHCP", lh)):
        peak = (metrics[hand]["peak_theta_deg"], metrics[hand]["peak_phi_deg"])
        pattern = analyze(theta, phi, power, peak=peak)
        metrics[hand]["sll_first_null_db"] = pattern["SLL_db"]
        metrics[hand]["sll_direction_deg"] = pattern["SLL_dir"]
    return metrics


def main():
    both_path, prior_path, b_path, plain_path, points_path, output_path = map(Path, sys.argv[1:7])
    both = json.loads(both_path.read_text())["states"]["BOTH"]
    prior = json.loads(prior_path.read_text())["ideal"]
    b = add_sll(json.loads(b_path.read_text()))
    plain = add_sll(json.loads(plain_path.read_text()))
    points = json.loads(points_path.read_text())
    assert points["status"] == "COMPLETE" and points["solves_started"] == 0
    for state in (prior, b, plain):
        state["B_peak_drop_vs_both_db"] = both["LHCP"]["peak_gain_db"] - state["LHCP"]["peak_gain_db"]
        state["A_peak_drop_vs_both_db"] = both["RHCP"]["peak_gain_db"] - state["RHCP"]["peak_gain_db"]
    assert b["B_peak_drop_vs_both_db"] < 13, "CONDITIONAL_ALL_PLAIN_NOT_TRIGGERED"
    result = {"label": "HFSS_2026_R1_V31_CEILING", "frequency_ghz": 5.68,
              "all_plain_condition_met": True, "both": both, "a_only_ideal_234": prior,
              "b_all_off_ideal_258": b, "all_plain_516": plain,
              "gain_at_reference_directions": points["samples"],
              "sll_definition": "Highest same-hand lobe outside first null, fallback 2x mean HPBW, above horizon."}
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2))
    print("B_ALL drop", b["B_peak_drop_vs_both_db"], "ALL_PLAIN LHCP at 36/180",
          points["samples"]["ALL_PLAIN"]["th36_ph180"]["LHCP_gain_dbi"])


if __name__ == "__main__":
    main()
