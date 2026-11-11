import csv, json, math, re
from pathlib import Path
rows = Path(__file__).resolve().parents[2] / "results" / "sub6_gate" / "rows"
order = ["BENCH_D1_slab_H20", "BENCH_D2_n2_H20", "BENCH_D2_n4_H20", "REF_H20", "BENCH_D2_n8_H20",
         "BENCH_D3_n6_H10", "REF_H15", "BENCH_D3_n6_H30", "BENCH_D3_n6_H40",
         "BENCH_D4_n6_slab5_H20", "BENCH_D4_n6_slab15_H20"]
C = 299792458.0
out, lines = [], []
for rid in order:
    r = json.loads((rows / f"{rid}.json").read_text())
    t = r["mesh_raw"]["convergence_text"] or ""
    passes = re.search(r"Completed : (\d+)", t).group(1); conv = re.search(r"Converged : (\w+)", t).group(1)
    ds = re.search(r"Current   : ([\d.eE+-]+)", t).group(1)
    tets = re.findall(r"^\s*\d+\|\s*(\d+)\|", t, re.M)[-1]
    win = [d for d in r["per_frequency"] if 5.6 - 1e-6 <= d["f_ghz"] <= 6.0 + 1e-6]
    worst_s21 = min(20 * math.log10(abs(complex(*d["S21"]))) for d in win)
    worst_s11 = max(abs(complex(*d["S11"])) for d in win)
    lines.append(f"{rid}: passes {passes} conv {conv} dS {float(ds):.4f} tets {tets} | 5.6-6.0: min S21 {worst_s21:.2f} dB, max|S11| {worst_s11:.3f}"
                 f" PASS={worst_s21 > -1.5 and worst_s11 < 0.3}")
    for d in r["per_frequency"]:
        f = d["f_ghz"]; k0 = 2 * math.pi * f * 1e9 / C
        g = d.get("port_gamma_1_per_m")
        rec = {"run_id": rid, "f_ghz": f, "S21_dB": round(20 * math.log10(abs(complex(*d["S21"]))), 3),
               "abs_S11": round(abs(complex(*d["S11"])), 4), "absorbed": round(d["absorbed_fraction"], 4),
               "alpha_per_m": "" if not g else round(g[0], 5), "beta_over_k0": "" if not g else round(g[1] / k0, 4)}
        out.append(rec)
    if "gamma_all_1_per_m" in r["per_frequency"][0]:
        for g in r["per_frequency"][0]["gamma_all_1_per_m"]:
            seq = [(d["f_ghz"], d["gamma_all_1_per_m"][g]) for d in r["per_frequency"]]
            lines.append("  " + g + " (alpha/m, beta/k0): " + "; ".join(
                f"{f:.2f}:{a:.2f}/{b/(2*math.pi*f*1e9/C):.3f}" for f, (a, b) in seq if round(f * 100) % 10 in (0, 5) and f in (5.5, 5.7, 5.8, 6.0, 6.1)))
p = rows.parents[0] / "bench_diagnosis.csv"
with p.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
(rows.parents[0] / "bench_summary.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
