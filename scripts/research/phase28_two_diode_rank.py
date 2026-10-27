"""Ranks the two-diode screen: for each (base, roles) the four states (A, B) give tensor vectors v = delta exp(j 2 axis); a usable 2-bit set needs the four vectors at 90 deg spacing (axes 45 deg apart) with equal delta and equal b_bar.
Quantities: ang_err (largest deviation of an axis gap from 45 deg), d_spread, b_spread (max-min over mean), alpha_min (own loading, mesh 32). python scripts/research/phase28_two_diode_rank.py <screen.jsonl> <out.json>"""
import json, sys, math, collections
rows = [json.loads(l) for l in open(sys.argv[1])]
g = collections.defaultdict(dict)
for r in rows:
    g[(r["base"], r["roles"])][(r["a"], r["b"])] = r
out = []
for (base, roles), st in g.items():
    if len(st) != 4 or any(s.get("alpha_beam_np_per_m") is None or s.get("delta_ms") is None for s in st.values()):
        continue
    ss = [st[k] for k in sorted(st)]
    ang = sorted((2 * math.radians(s["axis_deg"])) % (2 * math.pi) for s in ss)
    gaps = [(ang[(i + 1) % 4] - ang[i]) % (2 * math.pi) for i in range(4)]
    ang_err = max(abs(math.degrees(x) / 2 - 45.0) for x in gaps)
    d = [s["delta_ms"] for s in ss]; b = [s["b_bar_ms"] for s in ss]; al = [s["alpha_beam_np_per_m"] for s in ss]
    out.append({"base": base, "roles": roles, "ang_err_axis_deg": ang_err, "delta_spread": (max(d) - min(d)) / (sum(d) / 4), "b_spread": (max(b) - min(b)) / (sum(b) / 4),
                "b_mean_ms": sum(b) / 4, "delta_min_ms": min(d), "alpha_min": min(al), "alpha_max": max(al),
                "states": {f"{k[0]}{k[1]}": {"axis_deg": st[k]["axis_deg"], "delta_ms": st[k]["delta_ms"], "b_bar_ms": st[k]["b_bar_ms"], "alpha": st[k]["alpha_beam_np_per_m"]} for k in sorted(st)}})
json.dump(out, open(sys.argv[2], "w"))
print(len(out), "configs complete")
def show(title, key, flt):
    c = sorted([o for o in out if flt(o)], key=key)[:8]
    print(title, len(c))
    for o in c: print(" ", o["base"], o["roles"], "ang_err %.1f dspr %.2f bspr %.2f bmean %.2f dmin %.2f alpha_min %.3f alpha_max %.3f" % (o["ang_err_axis_deg"], o["delta_spread"], o["b_spread"], o["b_mean_ms"], o["delta_min_ms"], o["alpha_min"], o["alpha_max"]))
show("balanced (ang<=12, dspread<=.3, bspread<=.2): best alpha_min", lambda o: -o["alpha_min"], lambda o: o["ang_err_axis_deg"] <= 12 and o["delta_spread"] <= .3 and o["b_spread"] <= .2)
show("loose (ang<=20, dspread<=.5, bspread<=.4): best alpha_min", lambda o: -o["alpha_min"], lambda o: o["ang_err_axis_deg"] <= 20 and o["delta_spread"] <= .5 and o["b_spread"] <= .4)
show("best angular quality, any", lambda o: o["ang_err_axis_deg"], lambda o: True)
show("best by ang_err with alpha_min>=0.5", lambda o: o["ang_err_axis_deg"], lambda o: o["alpha_min"] >= 0.5)
