"""Compare two island_rows.json (mesh 64 vs mesh 96) for the ids present in both: alpha, beta/k0, Delta per state, relative change |x96 - x64| / |x64|.
Prints one line per row and state, the max relative change per row, and HARD_STOP if any > 10 %.   python sub6_island_meshcompare.py <rows64.json> <rows96.json> [ids,...]"""
import json
import sys

a = {r["row"]["id"]: r for r in json.load(open(sys.argv[1]))["rows"]}
b = {r["row"]["id"]: r for r in json.load(open(sys.argv[2]))["rows"]}
ids = sys.argv[3].split(",") if len(sys.argv) > 3 else [k for k in b if k in a]
c = lambda z: complex(z[0], z[1]) if isinstance(z, list) else complex(z)
rel = lambda x, y: abs(y - x) / abs(x) if abs(x) > 0 else float("nan")
worst = 0.0
print("id,state,aU64,aU96,d_alpha,b64,b96,d_beta,|D|64,|D|96,d_Delta")
for k in ids:
    if k not in a or k not in b:
        print(k, "MISSING in", "64" if k not in a else "96")
        continue
    m = 0.0
    for st in ("OFF_BAR50", "ON_BAR50"):
        s4, s9 = a[k]["states"][st], b[k]["states"][st]
        da, db, dd = rel(s4["alpha_np_per_m"], s9["alpha_np_per_m"]), rel(s4["beta_over_k0"], s9["beta_over_k0"]), rel(c(s4["Delta_mS"]), c(s9["Delta_mS"]))
        m = max(m, da, db, dd)
        print(f"{k},{st},{s4['alpha_np_per_m']:.3g},{s9['alpha_np_per_m']:.3g},{da:.1%},{s4['beta_over_k0']:.4g},{s9['beta_over_k0']:.4g},{db:.1%},{abs(c(s4['Delta_mS'])):.3g},{abs(c(s9['Delta_mS'])):.3g},{dd:.1%}")
    print(f"{k} MAX_REL_CHANGE {m:.1%} leak64={a[k]['leakage']} leak96={b[k]['leakage']}")
    worst = max(worst, m)
print("WORST", f"{worst:.1%}", "HARD_STOP" if worst > 0.10 else "WITHIN_10PCT")
