"""Compare two island_rows.json (mesh 64 vs mesh 96) for the ids present in both.  Rule 3 AMENDED (user decision 2026-11-08):
  * OFF_BAR50 alpha, beta/k0, Delta: relative change |x96 - x64| / |x64| must be <= 10 %.
  * ON_BAR50 alpha: PASS if |alpha| < 0.3 Np/m at BOTH meshes (absolute criterion for near-zero values); otherwise the 10 % relative criterion applies.
  * ON_BAR50 beta/k0 and Delta: reported, not gated.
Prints one line per row and state, a verdict per row (PASS / FAIL) and HARD_STOP if any row fails.
    python sub6_island_meshcompare.py <rows64.json> <rows96.json> [ids,...]"""
import json
import sys

LIM, A_ABS = 0.10, 0.3
a = {r["row"]["id"]: r for r in json.load(open(sys.argv[1]))["rows"]}
b = {r["row"]["id"]: r for r in json.load(open(sys.argv[2]))["rows"]}
ids = sys.argv[3].split(",") if len(sys.argv) > 3 else [k for k in b if k in a]
c = lambda z: complex(z[0], z[1]) if isinstance(z, list) else complex(z)
rel = lambda x, y: abs(y - x) / abs(x) if abs(x) > 0 else float("nan")
fails = []
print("id,state,alpha64,alpha96,d_alpha,beta64,beta96,d_beta,|Delta|64,|Delta|96,d_Delta,gated")
for k in ids:
    if k not in a or k not in b:
        print(k, "MISSING in", "64" if k not in a else "96")
        fails.append(k)
        continue
    worst, why = 0.0, ""
    for st in ("OFF_BAR50", "ON_BAR50"):
        s4, s9 = a[k]["states"][st], b[k]["states"][st]
        a4, a9 = s4["alpha_np_per_m"], s9["alpha_np_per_m"]
        da, db, dd = rel(a4, a9), rel(s4["beta_over_k0"], s9["beta_over_k0"]), rel(c(s4["Delta_mS"]), c(s9["Delta_mS"]))
        if st == "OFF_BAR50":
            gated = max(da, db, dd)
            tag = "OFF: alpha,beta,Delta"
        else:
            on_abs = abs(a4) < A_ABS and abs(a9) < A_ABS
            gated = 0.0 if on_abs else da
            tag = "ON: alpha ABS<0.3 both meshes" if on_abs else "ON: alpha relative (>=0.3 at one mesh)"
        if gated > worst:
            worst, why = gated, f"{st} {tag}"
        print(f"{k},{st},{a4:.3g},{a9:.3g},{da:.1%},{s4['beta_over_k0']:.4g},{s9['beta_over_k0']:.4g},{db:.1%},{abs(c(s4['Delta_mS'])):.3g},{abs(c(s9['Delta_mS'])):.3g},{dd:.1%},{tag}")
    ok = worst <= LIM
    if not ok:
        fails.append(k)
    print(f"{k} GATED_MAX {worst:.1%} ({why or 'n/a'}) -> {'PASS' if ok else 'FAIL'}  leak64={a[k]['leakage']} leak96={b[k]['leakage']}")
print("VERDICT", "HARD_STOP failing rows: " + ",".join(fails) if fails else "ALL_PASS")
