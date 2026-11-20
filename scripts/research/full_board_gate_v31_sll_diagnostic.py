"""Re-score stored v3.1 HFSS far fields; no AEDT process or solve is started."""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from full_board_gate_v31_extract import direction_separation, ffd
from sidequest_sll_def import analyze


STATE_KEYS = (
    ("BOTH", "both", "results/sub6_meander_model/full_board_gate/v31/BOTH.ffd"),
    ("A_ONLY_IDEAL_234", "a_only_ideal_234", "results/sub6_meander_model/full_board_gate/v31_ideal/A_ONLY_IDEAL.ffd"),
    ("B_ALL_OFF_IDEAL_258", "b_all_off_ideal_258", "results/sub6_meander_model/full_board_gate/v31_ceiling/B_ALL_OFF_IDEAL.ffd"),
    ("ALL_PLAIN", "all_plain_516", "results/sub6_meander_model/full_board_gate/v31_ceiling/ALL_PLAIN.ffd"),
)
TARGETS = ((44, 144), (44, 176), (32, 0), (36, 180))
EXCLUSION_CENTRES = ((44, 144), (44, 176))


def db(value):
    return 10 * math.log10(max(float(value), 1e-300))


def sample_gain(theta, phi, gain, target_theta, target_phi):
    i = int(np.flatnonzero(np.isclose(theta, target_theta, rtol=0, atol=1e-9))[0])
    j = int(np.flatnonzero(np.isclose(phi, target_phi, rtol=0, atol=1e-9))[0])
    return float(gain[i, j])


def candidates(theta, phi, gain, main_mask, peak_gain, n=5):
    allowed = (~main_mask) & (theta[:, None] <= 90)
    allowed[:, -1] = False  # Duplicate 360° phi column.
    result = []
    for rank in range(1, n + 1):
        work = np.where(allowed, gain, -np.inf)
        i, j = np.unravel_index(int(np.argmax(work)), work.shape)
        if not np.isfinite(work[i, j]):
            break
        th, ph, level = float(theta[i]), float(phi[j]), float(gain[i, j])
        result.append({"rank": rank, "theta_deg": th, "phi_deg": ph,
                       "gain_dbi": level, "db_below_peak": float(peak_gain - level)})
        allowed &= direction_separation(theta, phi, th, ph) >= 15
    return result


def main():
    root, metrics_path, output_dir = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    output_dir.mkdir(parents=True, exist_ok=True)
    metrics = json.loads(metrics_path.read_text())
    output = {"label": "HFSS_2026_R1_V31_STORED_FFD_ONLY", "frequency_ghz": 5.68,
              "sll_definition": "Highest same-hand gain above horizon outside the first-null main-lobe mask; fallback 2x mean HPBW.",
              "second_definition": "Next highest grid candidate outside the same main-lobe mask and at least 15 deg from the first SLL.",
              "states": {}, "both_sll_region_test": {}}
    csv_rows = []
    for label, key, ffd_rel in STATE_KEYS:
        state = metrics[key]
        theta, phi, et, ep = ffd(root / ffd_rel)
        plus = np.abs(et + 1j * ep) ** 2
        minus = np.abs(et - 1j * ep) ** 2
        rh, lh = (plus, minus) if state["ffd_rhcp_convention"] == "et_plus_j_ep" else (minus, plus)
        output["states"][label] = {}
        for hand, power in (("RHCP", rh), ("LHCP", lh)):
            peak = state[hand]
            peak_th, peak_ph, peak_gain = (float(peak["peak_theta_deg"]),
                                           float(peak["peak_phi_deg"]), float(peak["peak_gain_db"]))
            i = int(np.flatnonzero(np.isclose(theta, peak_th))[0])
            j = int(np.flatnonzero(np.isclose(phi, peak_ph))[0])
            gain = db(power[i, j])
            gain_grid = 10 * np.log10(np.maximum(power, 1e-300)) + peak_gain - gain
            main_mask, _ = analyze(theta, phi, power, peak=(peak_th, peak_ph), return_mask=True)
            lobes = candidates(theta, phi, gain_grid, main_mask, peak_gain)
            assert len(lobes) >= 2
            prior_sll = peak.get("sll_first_null_db")
            if prior_sll is not None:
                assert abs(-lobes[0]["db_below_peak"] - float(prior_sll)) < 0.02, (label, hand, prior_sll, lobes[0])
            gains = {f"th{th}_ph{ph}": sample_gain(theta, phi, gain_grid, th, ph) for th, ph in TARGETS}
            record = {"peak": {"theta_deg": peak_th, "phi_deg": peak_ph, "gain_dbi": peak_gain},
                      "lobes": lobes, "gain_at_targets_dbi": gains}
            if label == "BOTH":
                exclude = np.zeros_like(main_mask)
                for th, ph in EXCLUSION_CENTRES:
                    exclude |= direction_separation(theta, phi, th, ph) <= 20
                masked = np.where((~main_mask) & (~exclude) & (theta[:, None] <= 90), gain_grid, -np.inf)
                masked[:, -1] = -np.inf
                k = np.unravel_index(int(np.argmax(masked)), masked.shape)
                record["sll_excluding_residual_regions"] = {
                    "theta_deg": float(theta[k[0]]), "phi_deg": float(phi[k[1]]),
                    "gain_dbi": float(masked[k]), "db_below_peak": float(peak_gain - masked[k])}
            output["states"][label][hand] = record
            for kind, lobe in [("peak", record["peak"])] + [(f"sll_{v['rank']}", v) for v in lobes]:
                csv_rows.append({"state": label, "hand": hand, "kind": kind,
                                 "theta_deg": lobe["theta_deg"], "phi_deg": lobe["phi_deg"],
                                 "gain_dbi": lobe["gain_dbi"],
                                 "db_below_peak": lobe.get("db_below_peak", 0.0)})
            if label == "BOTH":
                v = record["sll_excluding_residual_regions"]
                csv_rows.append({"state": label, "hand": hand, "kind": "sll_excluding_44deg_regions",
                                 "theta_deg": v["theta_deg"], "phi_deg": v["phi_deg"],
                                 "gain_dbi": v["gain_dbi"], "db_below_peak": v["db_below_peak"]})
    both = output["states"]["BOTH"]
    for hand in ("RHCP", "LHCP"):
        first = both[hand]["lobes"][0]
        other = both["LHCP" if hand == "RHCP" else "RHCP"]["lobes"][0]
        sep = float(direction_separation(np.array([first["theta_deg"]]),
                                         np.array([first["phi_deg"]]),
                                         other["theta_deg"], other["phi_deg"])[0, 0])
        distances = {f"th{th}_ph{ph}": float(direction_separation(
            np.array([first["theta_deg"]]), np.array([first["phi_deg"]]), th, ph)[0, 0])
            for th, ph in EXCLUSION_CENTRES}
        output["both_sll_region_test"][hand] = {
            "distance_to_other_hand_sll_deg": sep, "same_region_as_other_hand_within_20deg": sep <= 20,
            "distance_to_residual_centres_deg": distances,
            "within_20deg_of_either_residual_centre": min(distances.values()) <= 20}
    (output_dir / "sll_diagnostic.json").write_text(json.dumps(output, indent=2))
    with (output_dir / "lobes.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(csv_rows[0]))
        writer.writeheader()
        writer.writerows(csv_rows)
    print("BOTH", output["both_sll_region_test"])


if __name__ == "__main__":
    main()
