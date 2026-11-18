"""No-HFSS anchor frequency-slope and first-order validity audit."""

import csv
import json
import math
from pathlib import Path

from scripts.research.sub6_match_supercell import load_modes
from scripts.research.sub6_supercell_model_check import R, beta_pair, estimate

C = 299792458.0
P_SUPER = 0.0928
P_CELL = 0.0116
PHASES = (15, 30, 45)
MODES = (4, 4, 3)


def main():
    baseline = json.loads((R / "folded_qbase.json").read_text(encoding="utf-8"))["rows"]
    picked = {p: load_modes("ANCHOR", p)[m - 1] for p, m in zip(PHASES, MODES)}
    slope_vg = {}
    for a, b in ((15, 30), (30, 45)):
        slope_vg[f"{a}-{b}"] = abs((picked[b]["f_ghz"] - picked[a]["f_ghz"]) * 1e9
                                   / (b - a)) * 360 * P_SUPER
    source_85 = {}
    for direction in ("across", "along"):
        data = json.loads((R / f"ANCHOR_PML_{direction}_85_air40.035_modes.json")
                          .read_text(encoding="utf-8"))["modes"]
        source_85[direction] = next(m["Q"] for m in data if m["Q"] < 100000)
    qbase_85 = (source_85["across"] + source_85["along"]) / 2
    rows = []
    for p, mode_index in zip(PHASES, MODES):
        mode = picked[p]
        f, q = mode["f_ghz"], mode["Q"]
        point = next(row for row in baseline if row["cell"] == "ANCHOR"
                     and row["super_phase_deg"] == p and row["fold_sign"] == -1
                     and row["branch_index"] == 1)
        qb = point["Q_base"]
        phase_cell = point["per_cell_phase_deg"]
        beta_ratio = phase_cell * C / (360 * P_CELL * f * 1e9)
        harmonic_angle = math.degrees(math.asin(beta_ratio - C / (P_SUPER * f * 1e9)))
        pair = beta_pair("ANCHOR", f)
        chain = estimate(f, pair, harmonic_angle)
        local = (["15-30"] if p == 15 else ["30-45"] if p == 45 else
                 ["15-30", "30-45"])
        for interval in local:
            vg = slope_vg[interval]
            alpha = math.pi * f * 1e9 * (1 / q - 1 / qb) / vg
            qb_clamp_alpha = math.pi * f * 1e9 * (1 / q - 1 / qbase_85) / vg
            rows.append({"phase_deg": p, "mode": mode_index, "f_ghz": f,
                         "Q": q, "Q_base": qb, "Q_base_85_clamp": qbase_85,
                         "Q_base_extrapolated": point["phase_extrapolated"],
                         "per_cell_phase_deg": phase_cell, "beta_over_k0": beta_ratio,
                         "harmonic_angle_deg": harmonic_angle,
                         "local_slope_interval": interval, "v_g_m_s": vg,
                         "hfss_alpha_local_np_m": alpha,
                         "hfss_alpha_Qbase85_np_m": qb_clamp_alpha,
                         "model_alpha_own_angle_np_m": chain["alpha_np_m"],
                         "model_over_hfss": chain["alpha_np_m"] / alpha,
                         "B_x_mS": chain["B_x_mS"], "B_y_mS": chain["B_y_mS"],
                         "B_mean_mS": chain["B_mean_mS"],
                         "delta_mS": chain["delta_mS"],
                         "delta_over_abs_B": chain["delta_mS"] / abs(chain["B_mean_mS"]),
                         "x_beta_bracketed": pair["x"]["bracketed"],
                         "y_beta_bracketed": pair["y"]["bracketed"]})
    (R / "anchor_resume4_diag.json").write_text(json.dumps({
        "label": "HFSS modal alpha with local slopes; beta-B-alpha MODEL_ONLY",
        "qbase_85_by_direction": source_85,
        "local_vg_m_s": slope_vg, "rows": rows}, indent=2) + "\n", encoding="utf-8")
    with (R / "anchor_resume4_diag.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(row["phase_deg"], row["local_slope_interval"],
              f"f={row['f_ghz']:.6f}", f"vg={row['v_g_m_s']/1e6:.3f}",
              f"alpha_HFSS={row['hfss_alpha_local_np_m']:.4f}",
              f"alpha_model={row['model_alpha_own_angle_np_m']:.4f}",
              f"delta/B={row['delta_over_abs_B']:.3f}",
              f"Qbase={row['Q_base']:.1f}")


if __name__ == "__main__":
    main()
