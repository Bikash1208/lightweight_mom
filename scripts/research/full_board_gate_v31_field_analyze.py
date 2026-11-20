"""Compare stored complex-E radial samples from the three unchanged v3.1 HFSS projects."""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

PHIS = (0, 45, 90, 135, 180, 225, 270, 315)
STATES = ("BOTH", "A_ONLY", "B_ONLY")
GROUPS = {"x_lt_0": (135, 180, 225), "x_gt_0": (315, 0, 45),
          "axis_y_positive": (90,), "axis_y_negative": (270,)}


def mean_power(rows, state, phis, r_min, r_max, column):
    vals = [float(row[column]) for row in rows if row["state"] == state
            and int(row["phi_deg"]) in phis and r_min <= int(row["r_mm"]) <= r_max]
    assert vals and np.isfinite(vals).all()
    return float(np.mean(vals))


def ratio_db(rows, state, phis, r_min, r_max, column):
    both = mean_power(rows, "BOTH", phis, r_min, r_max, column)
    case = mean_power(rows, state, phis, r_min, r_max, column)
    return 10 * math.log10(case / both)


def main():
    work, out, layout = [Path(v) for v in sys.argv[1:4]]
    out.mkdir(parents=True, exist_ok=True)
    all_rows = []
    for state in STATES:
        with (work / state / f"{state}_radial.csv").open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        assert len(rows) == 8 * 151 and rows[0]["state"] == state
        all_rows.extend(rows)
    with (out / "radial_profiles.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(all_rows[0]))
        writer.writeheader()
        writer.writerows(all_rows)
    with layout.open(newline="", encoding="utf-8") as stream:
        sites = list(csv.DictReader(stream))
    result = {"label": "HFSS_2026_R1_V31_POSTPROCESS_ONLY", "frequency_ghz": 5.68,
              "sample_height_mm_above_top_cu": 2.0, "sample_z_mm": 2.035,
              "radial_step_mm": 1, "radii_mm": [0, 150], "averaging_r_mm": [40, 140],
              "states": {}, "groups": GROUPS,
              "nan_note": "r=0 lies inside the PEC feed pin and is NaN in every field line."}
    for state in ("A_ONLY", "B_ONLY"):
        blocked_beam = "B" if state == "A_ONLY" else "A"
        via_sites = [{"i": int(site["i"]), "j": int(site["j"]),
                      "x_mm": float(site["x_mm"]), "y_mm": float(site["y_mm"]),
                      "gate_for": site["gate_for"]}
                     for site in sites if site["gate_for"] == blocked_beam]
        assert len(via_sites) == 8
        entry = {"blocked_beam": blocked_beam, "via_sites": via_sites, "per_phi": {}, "groups": {}}
        for phi in PHIS:
            entry["per_phi"][str(phi)] = {
                "total_e2_ratio_db_r40_140": ratio_db(all_rows, state, (phi,), 40, 140, "e_total_sq"),
                "ez2_ratio_db_r40_140": ratio_db(all_rows, state, (phi,), 40, 140, "e_z_sq"),
                "total_e2_ratio_db_r6_29": ratio_db(all_rows, state, (phi,), 6, 29, "e_total_sq")}
        for name, phis in GROUPS.items():
            entry["groups"][name] = {
                "total_e2_ratio_db_r40_140": ratio_db(all_rows, state, phis, 40, 140, "e_total_sq"),
                "ez2_ratio_db_r40_140": ratio_db(all_rows, state, phis, 40, 140, "e_z_sq"),
                "total_e2_ratio_db_r6_29": ratio_db(all_rows, state, phis, 6, 29, "e_total_sq")}
        result["states"][state] = entry
    (out / "field_ratios.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    for state in ("A_ONLY", "B_ONLY"):
        print(state, result["states"][state]["groups"])


if __name__ == "__main__":
    main()
