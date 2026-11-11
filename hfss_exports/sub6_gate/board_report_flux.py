import csv, json, math, re
from pathlib import Path
R = Path(__file__).resolve().parents[2] / "results" / "sub6_gate" / "board"
runs = {}
for k in ("BOARD_REF", "BOARD_OFF", "BOARD_ON"):
    runs[k] = {"flux": json.loads((R / "refield" / f"{k}_5.8GHz_flux.json").read_text()), "row": json.loads((R / "rows" / f"{k}.json").read_text()), "src": "5.8 saved project"}
for k in ("BOARD_C1", "BOARD_C2", "BOARD_C3"):
    r = json.loads((R / "rows" / f"{k}.json").read_text()); runs[k] = {"flux": r["flux"], "row": r, "src": "new solve"}
cols = ["run", "source", "passes", "converged", "dS", "tets", "S11_5p8", "counts_ok"]
planes = [f"{s}{x}_{b}" for x in ("75", "45", "25") for b in ("lo", "hi") for s in ("P", "M")]
out = []
for k, v in runs.items():
    t = v["row"]["mesh_raw"]["convergence_text"]
    tab = re.findall(r"^\s*\d+\|\s*(\d+)\|\s*([\w.]+)\|", t, re.M)
    s58 = next(x for x in v["row"]["S11"] if abs(x["f_ghz"] - 5.8) < 1e-6)
    rec = {"run": k, "source": v["src"], "passes": re.search(r"Completed : (\d+)", t).group(1), "converged": re.search(r"Converged : (\w+)", t).group(1),
           "dS": tab[-1][1], "tets": tab[-1][0], "S11_5p8": round(s58["abs_S11"], 4),
           "counts_ok": v["row"]["build"]["object_counts"] == v["row"]["build"]["expected_counts"]}
    f = v["flux"]["planes"]
    for p in planes:
        rec["flux_W_" + p] = round(f[p], 5)
    for x in ("75", "45", "25"):
        for b in ("lo", "hi"):
            pp, mm = f[f"P{x}_{b}"], f[f"M{x}_{b}"]
            rec[f"ratio_dB_{x}_{b}"] = round(10 * math.log10(pp / mm), 2) if pp > 0 and mm > 0 else "n/a(<=0)"
    ff = v["row"].get("far_field")
    if ff:
        rec["P_rad_W"] = round(ff["radiated_power_W_per_1W_incident"], 4); rec["P_acc_W"] = round(ff["accepted_power_W_per_1W_incident"], 4)
    out.append(rec)
fields = []
for r in out:
    for k in r:
        if k not in fields: fields.append(k)
with (R / "flux_summary_5p8.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(out)
for r in out:
    print(r["run"], r["source"], "passes", r["passes"], r["converged"], "dS", r["dS"], "tets", r["tets"], "S11", r["S11_5p8"], r["counts_ok"])
    print("   x75 lo P/M %.4f/%.4f (%s dB)  hi %.4f/%.4f (%s) | x45 lo %.4f/%.4f (%s) hi %.4f/%.4f (%s) | x25 lo %.4f/%.4f (%s) hi %.4f/%.4f (%s)" % (
        r["flux_W_P75_lo"], r["flux_W_M75_lo"], r["ratio_dB_75_lo"], r["flux_W_P75_hi"], r["flux_W_M75_hi"], r["ratio_dB_75_hi"],
        r["flux_W_P45_lo"], r["flux_W_M45_lo"], r["ratio_dB_45_lo"], r["flux_W_P45_hi"], r["flux_W_M45_hi"], r["ratio_dB_45_hi"],
        r["flux_W_P25_lo"], r["flux_W_M25_lo"], r["ratio_dB_25_lo"], r["flux_W_P25_hi"], r["flux_W_M25_hi"], r["ratio_dB_25_hi"]))
