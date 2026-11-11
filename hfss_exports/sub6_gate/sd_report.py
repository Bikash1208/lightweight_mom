import csv, json, math, re
from pathlib import Path
R = Path(__file__).resolve().parents[2] / "results" / "sub6_gate"
rows = {p.stem: json.loads(p.read_text()) for p in sorted((R / "rows").glob("*.json")) if p.stem.startswith(("NR_", "SDW2P9"))}
def conv(r):
    t = r["mesh_raw"]["convergence_text"] or ""
    return (re.search(r"Completed : (\d+)", t).group(1), re.search(r"Converged : (\w+)", t).group(1),
            float(re.search(r"Current   : ([\d.eE+-]+)", t).group(1)), re.findall(r"^\s*\d+\|\s*(\d+)\|", t, re.M)[-1])
def at(r, f): return next(d for d in r["per_frequency"] if abs(d["f_ghz"] - f) < 1e-6)
db = lambda z: 20 * math.log10(max(abs(complex(*z)), 1e-12))
out, lines = [], []
for rid, r in rows.items():
    p, c, ds, tets = conv(r)
    d = at(r, 5.8); zo = d["port_zo_ohm"]
    row = r["row"]
    rec = {"run_id": rid, "kind": row["kind"], "variant": row.get("variant", ""), "n_gate": row["n_gate"] if row["kind"] == "SDGATE" else "",
           "state": row["state"] or "", "n_static": row.get("n_static", ""), "passes": p, "converged": c, "dS": round(ds, 4), "tets": tets,
           "readback_ok": r["readback_ok"], "S21_dB_5p8": round(db(d["S21"]), 3), "abs_S11_5p8": round(abs(complex(*d["S11"])), 4),
           "absorbed_5p8": round(d["absorbed_fraction"], 4), "Zo_P1_ohm": f"{zo['Zo(P1)'][0]:.2f}{zo['Zo(P1)'][1]:+.2f}j"}
    out.append(rec)
# reference-normalised SD rows
ref = rows["NR_REF_H20"]
for rec in out:
    if rec["kind"] == "SDGATE":
        rec["S21_over_REF_dB_5p8"] = round(rec["S21_dB_5p8"] - db(at(ref, 5.8)["S21"]), 3)
        r = rows[rec["run_id"]]
        win = [x for x in r["per_frequency"] if 5.6 - 1e-6 <= x["f_ghz"] <= 6.0 + 1e-6]
        rec["min_S21_dB_5p6_6p0"] = round(min(db(x["S21"]) for x in win), 3)
        rec["max_S11_5p6_6p0"] = round(max(abs(complex(*x["S11"])) for x in win), 4)
fields = sorted({k for r in out for k in r}, key=lambda k: list(out[0]).index(k) if k in out[0] else 99)
with (R / "summary_5p8GHz_SD_NR.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(out)
# per-frequency CSV for NR rows
with (R / "nr_bench_freq.csv").open("w", newline="") as f:
    w = csv.writer(f); w.writerow(["run_id", "f_ghz", "S21_dB", "abs_S11", "absorbed"])
    for rid, r in rows.items():
        if rid.startswith("NR_"):
            for x in r["per_frequency"]:
                w.writerow([rid, x["f_ghz"], round(db(x["S21"]), 4), round(abs(complex(*x["S11"])), 5), round(x["absorbed_fraction"], 5)])
for rec in out: print(rec)
