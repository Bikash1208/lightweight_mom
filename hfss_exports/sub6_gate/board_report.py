import csv, json, math, re
from pathlib import Path
R = Path(__file__).resolve().parents[2] / "results" / "sub6_gate" / "board"
rows = {k: json.loads((R / "rows" / f"{k}.json").read_text()) for k in ("BOARD_REF", "BOARD_OFF", "BOARD_ON")}
out, lines = [], []
for k, r in rows.items():
    t = r["mesh_raw"]["convergence_text"]
    tab = re.findall(r"^\s*(\d+)\|\s*(\d+)\|\s*([\w.]+)\|", t, re.M)
    conv = re.search(r"Converged : (\w+)", t).group(1)
    s58 = next(x for x in r["S11"] if abs(x["f_ghz"] - 5.8) < 1e-6)
    fl, ff, fp, fm = r["field_lines"], r["far_field"], r["field_lines"]["PROBE_P"], r["field_lines"]["PROBE_M"]
    rec = {"run": k, "passes": len(tab), "converged": conv, "last_dS": tab[-1][2], "tets": tab[-1][1], "cores": r["cores"], "seconds": round(r["solve_seconds"]),
           "counts": r["build"]["object_counts"], "expected": r["build"]["expected_counts"],
           "S11_5p8_abs": round(s58["abs_S11"], 4), "S11_5p8_dB": round(20 * math.log10(s58["abs_S11"]), 2),
           "E_behind_V_m": round(fp["mean_mag_E_V_per_m"], 3), "E_mirror_V_m": round(fm["mean_mag_E_V_per_m"], 3),
           "behind_over_mirror_dB": round(r["behind_over_mirror_dB"], 3),
           "frac_phi_m90_90_upper": round(ff["upper_phi_m90_90_fraction"], 4), "frac_phi_90_270_upper": round(ff["upper_phi_90_270_fraction"], 4),
           "P_rad_W": round(ff["radiated_power_W_per_1W_incident"], 4), "P_acc_W": round(ff["accepted_power_W_per_1W_incident"], 4),
           "rad_over_acc": round(ff["radiated_over_accepted"], 4), "upper_frac_of_rad": round(ff["upper_fraction_of_radiated"], 4)}
    out.append(rec)
    lines.append(k + " S11 sweep: " + "; ".join(f"{x['f_ghz']}:{x['abs_S11']:.3f}" for x in r["S11"]))
with (R / "board_summary.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
for o in out: print(o)
print("\n".join(lines))
ref = rows["BOARD_REF"]["field_lines"]
for k in ("BOARD_OFF", "BOARD_ON"):
    fl = rows[k]["field_lines"]
    print(k, "behind vs REF behind dB", round(20*math.log10(fl["PROBE_P"]["mean_mag_E_V_per_m"]/ref["PROBE_P"]["mean_mag_E_V_per_m"]),2),
          "mirror vs REF mirror dB", round(20*math.log10(fl["PROBE_M"]["mean_mag_E_V_per_m"]/ref["PROBE_M"]["mean_mag_E_V_per_m"]),2))
