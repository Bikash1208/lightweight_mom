import csv, json, math, re
from pathlib import Path
R = Path(__file__).resolve().parents[2] / "results" / "sub6_gate" / "board"
names = ("BOARD_REF", "BOARD_OFF", "BOARD_ON")
rows = {k: json.loads((R / "rows" / f"{k}_6p0.json").read_text()) for k in names}
out = []
for k, r in rows.items():
    t = r["mesh_raw"]["convergence_text"]
    passes = re.search(r"Completed : (\d+)", t).group(1); conv = re.search(r"Converged : (\w+)", t).group(1)
    tab = re.findall(r"^\s*\d+\|\s*(\d+)\|\s*([\w.]+)\|", t, re.M)
    fm, ff = r["field_metrics"], r["far_field"]
    s60 = next(x for x in r["S11"] if abs(x["f_ghz"] - 6.0) < 1e-6)
    L = fm["lines"]
    rec = {"run": k, "passes": passes, "converged": conv, "last_dS": tab[-1][1], "tets": tab[-1][0], "seconds": round(r["solve_seconds"]),
           "counts_ok": r["build"]["object_counts"] == r["build"]["expected_counts"],
           "S11_6p0_abs": round(s60["abs_S11"], 4), "S11_6p0_dB": round(20 * math.log10(s60["abs_S11"]), 2),
           "phasorE_z3_behind": round(L["PROBE_P3"]["CMag_E"]["mean"], 2), "phasorE_z3_mirror": round(L["PROBE_M3"]["CMag_E"]["mean"], 2),
           "ratio_dB_phasorE_z3": round(fm["behind_over_mirror_dB_CMag_E_z3"], 3),
           "phasorE_z1_behind": round(L["PROBE_P1"]["CMag_E"]["mean"], 2), "phasorE_z1_mirror": round(L["PROBE_M1"]["CMag_E"]["mean"], 2),
           "ratio_dB_phasorE_z1": round(fm["behind_over_mirror_dB_CMag_E_z1"], 3),
           "Ez_z1_behind": round(L["PROBE_P1"]["Mag_Ez"]["mean"], 2), "Ez_z1_mirror": round(L["PROBE_M1"]["Mag_Ez"]["mean"], 2),
           "ratio_dB_Ez_z1": round(20 * math.log10(L["PROBE_P1"]["Mag_Ez"]["mean"] / L["PROBE_M1"]["Mag_Ez"]["mean"]), 3),
           "MagE_ph0_z3_ratio_dB": round(fm["behind_over_mirror_dB_Mag_E_z3"], 3), "MagE_ph0_z1_ratio_dB": round(fm["behind_over_mirror_dB_Mag_E_z1"], 3),
           "frac_phi_m90_90": round(ff["upper_phi_m90_90_fraction"], 4), "frac_phi_90_270": round(ff["upper_phi_90_270_fraction"], 4),
           "P_rad_W": round(ff["radiated_power_W_per_1W_incident"], 4), "P_acc_W": round(ff["accepted_power_W_per_1W_incident"], 4),
           "rad_over_acc": round(ff["radiated_over_accepted"], 4)}
    out.append(rec)
with (R / "board_summary_6p0.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
# radial profiles CSV
with (R / "radial_E_6p0.csv").open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["x_mm"] + [f"{k}_phasorE_V_per_m" for k in names] + [f"{k}_MagE_phase0_V_per_m" for k in names])
    for side in ("RAD_M", "RAD_P"):
        xs = rows["BOARD_REF"]["field_metrics"]["radial"][side]["x_mm"]
        order = range(len(xs)) if side == "RAD_P" else range(len(xs) - 1, -1, -1)
        for i in order:
            w.writerow([xs[i]] + [round(rows[k]["field_metrics"]["radial"][side]["CMag_E"][i], 3) for k in names] +
                       [round(rows[k]["field_metrics"]["radial"][side]["Mag_E"][i], 3) for k in names])
for o in out: print(o)
sw = lambda k: "; ".join(f"{x['f_ghz']}:{x['abs_S11']:.3f}" for x in rows[k]["S11"])
for k in names: print(k, "S11 sweep", sw(k))
