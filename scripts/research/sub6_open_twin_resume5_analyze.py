"""OFF-twin raw-mode inventory and conservative gate audit from saved HFSS data."""

import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment


ROOT = Path(__file__).resolve().parents[2]
R = ROOT / "results" / "sub6_meander_model" / "supercell_hfss"
C0, P_SUPER, P_CELL = 299792458.0, 0.0928, 0.0116
PHASES = (15, 30)
TETS = {15: 91679, 30: 131335}


def modes(phase):
    path = R / f"SUPER_D3p6_OPEN_ph{phase}_PML40_resume5_modes.json"
    return json.loads(path.read_text(encoding="utf-8"))["modes"]


def profile(cell, phase, tag, mode):
    path = R / f"SUPER_{cell}_ph{phase}_PML40_{tag}_mode{mode}_MagE.fld"
    if not path.exists():
        return None
    samples = np.loadtxt(path, skiprows=2)[:, 3].reshape(8, 2)
    value = np.sqrt(np.mean(samples * samples, axis=1))
    norm = np.linalg.norm(value)
    return value / norm if norm > 0 else None


def cosine(a, b):
    return float(a @ b) if a is not None and b is not None else None


def baseline(phase, branch):
    point = (720 - phase) / 8
    q = []
    for direction in ("across", "along"):
        values = []
        for p in (85, 95):
            path = R / f"D3p6_OPEN_PML_{direction}_{p}_air40.035_resume5_modes.json"
            values.append(json.loads(path.read_text(encoding="utf-8"))["modes"][branch-1]["Q"])
        q.append(values[0] + (point - 85) / 10 * (values[1] - values[0]))
    return sum(q) / 2


def write_csv(name, rows):
    with (R / name).open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def main():
    by_phase = {p: modes(p) for p in PHASES}
    candidates = {p: [m for m in by_phase[p]
                      if 5.4 <= m["f_ghz"] <= 6.5 and m["Q"] < 1e5]
                  for p in PHASES}
    a, b = candidates[15], candidates[30]
    profiles = {(p, m["mode"]): profile("D3p6_OPEN", p, "resume5", m["mode"])
                for p, group in candidates.items() for m in group}
    cost = np.array([[1 - cosine(profiles[15, x["mode"]], profiles[30, y["mode"]])
                      + 0.3 * ((x["f_ghz"] - y["f_ghz"]) / .15) ** 2
                      for y in b] for x in a])
    ia, ib = linear_sum_assignment(cost)
    pairs = {a[i]["mode"]: b[j]["mode"] for i, j in zip(ia, ib)}
    inverse = {v: k for k, v in pairs.items()}
    all_rows = []
    for phase in PHASES:
        for m in by_phase[phase]:
            mode, f = m["mode"], m["f_ghz"]
            partner = (pairs if phase == 15 else inverse).get(mode)
            p_other = 30 if phase == 15 else 15
            other = by_phase[p_other][partner-1] if partner else None
            corr = cosine(profiles.get((phase, mode)), profiles.get((p_other, partner)))
            slope = ((other["f_ghz"] - f) / (p_other - phase)) if other else None
            vg = abs(slope * 1e9) * 360 * P_SUPER if slope is not None else None
            folded = {}
            for sign, label in ((-1, "minus"), (+1, "plus")):
                phase_cell = (720 + sign * phase) / 8
                beta = phase_cell * C0 / (360 * P_CELL * f * 1e9)
                harmonic = beta - C0 / (P_SUPER * f * 1e9)
                folded[f"beta_{label}_over_k0"] = beta
                folded[f"harmonic_{label}_angle_deg"] = (
                    math.degrees(math.asin(harmonic)) if abs(harmonic) <= 1 else None)
            all_rows.append({"phase_deg": phase, "mode": mode, "f_ghz": f,
                             "Q": m["Q"], "max_solved_tets": TETS[phase],
                             **folded, "partner_phase_deg": p_other if partner else None,
                             "field_matched_partner_mode": partner,
                             "field_cosine": corr, "slope_ghz_per_deg": slope,
                             "vg_abs_m_s": vg,
                             "window": "IN_5P4_6P5" if 5.4 <= f <= 6.5 else "ABOVE_6P5",
                             "classification": ("BOX_Q" if m["Q"] >= 1e5 else
                                                "LOW_Q" if m["Q"] <= 60 else
                                                "CANDIDATE")})
    write_csv("off_all_modes_resume5.csv", all_rows)
    gate_tracks = ((1, 1), (2, 2), (3, 3), (4, 5), (5, 4))
    on = next(x for x in json.loads((R / "d3p6_resume4_tracks.json").read_text(encoding="utf-8"))
              if x["track"] == "x_low_mean_candidate" and x["phase_deg"] == 15)
    fop, gate = on["f_ghz"], on["alpha_np_m"] / 12
    gate_rows = []
    for ma, mb in gate_tracks:
        x, y = by_phase[15][ma-1], by_phase[30][mb-1]
        vg = abs((y["f_ghz"] - x["f_ghz"]) * 1e9 / 15) * 360 * P_SUPER
        for phase, m in ((15, x), (30, y)):
            q1, q2 = baseline(phase, 1), baseline(phase, 2)
            qbase_conservative = min(q1, q2)
            alpha = math.pi * m["f_ghz"] * 1e9 * (
                1 / m["Q"] - 1 / qbase_conservative) / vg
            on_field = profile("B2p2_D3p6", phase, "resume4", 3)
            off_field = profiles[phase, m["mode"]]
            gate_rows.append({"off_track": f"{ma}->{mb}", "phase_deg": phase,
                              "mode": m["mode"], "f_ghz": m["f_ghz"],
                              "Q": m["Q"], "v_g_abs_m_s": vg,
                              "Q_base_branch1": q1, "Q_base_branch2": q2,
                              "Q_base_conservative_min": qbase_conservative,
                              "alpha_if_SW_np_m": alpha,
                              "ON_xlow_field_cosine": cosine(on_field, off_field),
                              "brackets_5p664": min(x["f_ghz"], y["f_ghz"]) <= fop <=
                                                  max(x["f_ghz"], y["f_ghz"]),
                              "ON_alpha_np_m": on["alpha_np_m"],
                              "OFF_gate_np_m": gate})
    write_csv("off_gate_candidates_resume5.csv", gate_rows)
    summary = {"method": "HFSS modes; 16-point |E| field-cosine plus frequency Hungarian pairing",
               "fold_warning": "Both +/- unfolded beta values are kinematic for non-SW modes",
               "fop_ghz": fop, "on_alpha_np_m": on["alpha_np_m"],
               "off_gate_np_m": gate,
               "qbase": {str(p): {str(branch): baseline(p, branch)
                                   for branch in (1, 2)} for p in PHASES},
               "pairs": pairs,
               "candidate_alpha_note": "Uses smaller of branch-1/2 PML Q to minimize extracted alpha; not a clean-SW claim"}
    spur = []
    for direction in ("across", "along"):
        pec = json.loads((R / f"D3p6_OPEN_PEC_{direction}_85_air40.035_resume5_modes.json")
                         .read_text(encoding="utf-8"))["modes"]
        p85 = json.loads((R / f"D3p6_OPEN_PML_{direction}_85_air40.035_resume5_modes.json")
                         .read_text(encoding="utf-8"))["modes"]
        p95 = json.loads((R / f"D3p6_OPEN_PML_{direction}_95_air40.035_resume5_modes.json")
                         .read_text(encoding="utf-8"))["modes"]
        for index in (0, 1):
            vg_unit = abs((p95[index]["f_ghz"] - p85[index]["f_ghz"]) * 1e9 / 10) * 360 * P_CELL
            alpha_spur = math.pi * p85[index]["f_ghz"] * 1e9 * (
                1 / p85[index]["Q"] - 1 / pec[index]["Q"]) / vg_unit
            spur.append({"direction": direction, "branch_index": index+1,
                         "f_PEC_ghz": pec[index]["f_ghz"],
                         "f_PML_ghz": p85[index]["f_ghz"],
                         "Q_PEC": pec[index]["Q"], "Q_PML": p85[index]["Q"],
                         "vg_unit_m_s": vg_unit, "alpha_spur_np_m": alpha_spur,
                         "pass_le_0p1": alpha_spur <= .1})
    summary["pml_spur"] = spur
    (R / "off_resume5_summary.json").write_text(json.dumps(summary, indent=2) + "\n",
                                                encoding="utf-8")
    print("gate", fop, gate, "pairs", pairs)
    for row in gate_rows:
        print(row["off_track"], row["phase_deg"], row["mode"],
              f"f={row['f_ghz']:.6f}", f"Q={row['Q']:.1f}",
              f"alpha_if_SW={row['alpha_if_SW_np_m']:.3f}",
              f"cos={row['ON_xlow_field_cosine']:.3f}")


if __name__ == "__main__":
    main()
