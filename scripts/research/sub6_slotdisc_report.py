"""CSV + summary of rows_*.jsonl.  MODEL_ONLY."""
import csv, json, sys, collections, cmath, math
from pathlib import Path
R = Path(__file__).resolve().parents[2] / "results/sub6_slotdisc"
rows = []
for f in sorted(R.glob("rows_*.jsonl")):
    rows += [json.loads(l) for l in f.read_text().splitlines() if l.strip()]
C = lambda z: complex(z[0], z[1])
out = []
for d in rows:
    r, s = d["row"], d["states"]
    rho = C(d["rho_nom"])
    out.append(dict(id=r["id"], sub=r["sub_name"], p=r["p"], gap=r["gap"], L_over_2r=r["ratio"], r=round(r["r"], 3), L=round(r["L"], 3), W=r["W"],
        beta_NOCUT=s["NOCUT"]["beta_over_k0"], beta_OFF=s["OFF_BAR50"]["beta_over_k0"], beta_ON=s["ON_BAR50"]["beta_over_k0"],
        Dimag_NOCUT=C(s["NOCUT"]["Delta_mS"]).imag, Dimag_OPEN=C(s["OFF_OPEN"]["Delta_mS"]).imag, Dimag_OFF=C(s["OFF_BAR50"]["Delta_mS"]).imag, Dimag_ON=C(s["ON_BAR50"]["Delta_mS"]).imag,
        mean_OFF=C(s["OFF_BAR50"]["mean_mS"]).imag, mean_ON=C(s["ON_BAR50"]["mean_mS"]).imag,
        A_perp_OFF=C(s["OFF_BAR50"]["A_perp_mS"]).real, A_perp_ON=C(s["ON_BAR50"]["A_perp_mS"]).real, dA_perp=C(d["dA_perp"]).real if isinstance(d["dA_perp"], list) else d["dA_perp"],
        rho_re=rho.real, rho_im=rho.imag, rho_abs=abs(rho), rho_arg_deg=math.degrees(cmath.phase(rho)),
        rho_open_re=C(d["rho_open"]).real, rho_ON_over_NOCUT_re=C(d["rho_vs_nocut"]).real,
        alpha_NOCUT=s["NOCUT"]["alpha_np_per_m"], alpha_OPEN=s["OFF_OPEN"]["alpha_np_per_m"], alpha_OFF=d["alpha_off"], alpha_ON=d["alpha_on"], alpha_ratio=d["alpha_ratio_on_over_off"],
        LEAKAGE=d["leakage"], resonance=("PROXY_FLIP:" + ",".join(d["resonance_flip"])) if d["resonance_flag"] else "", nodes=len(d["nodes"]), seconds=round(d["seconds"])))
out.sort(key=lambda x: -max(x["alpha_ON"], x["alpha_OFF"]))
with open(R / "rows_summary.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0])); w.writeheader()
    for x in out: w.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in x.items()})
print(len(out), collections.Counter(x["LEAKAGE"] for x in out))
print("resonance proxy flags:", sum(bool(x["resonance"]) for x in out))
for x in out[:12]:
    print(x["id"], "aON %.2f aOFF %.3f ratio %.0f" % (x["alpha_ON"], x["alpha_OFF"], x["alpha_ratio"]), "DimOFF %.2f DimON %.2f rho %.1f arg %.0f" % (x["Dimag_OFF"], x["Dimag_ON"], x["rho_abs"], x["rho_arg_deg"]), "rhoO %.2f" % x["rho_open_re"], "dA %.2f" % x["dA_perp"], x["LEAKAGE"])
# rho in the cell-scan range -1.4..-3.5 (real part), arg 160-180
sel = [x for x in out if -3.6 <= x["rho_re"] <= -1.3 and abs(x["rho_im"]) < 0.3 * abs(x["rho_re"])]
print("rows with rho_nom real in [-3.5,-1.3]:", len(sel), [(x["id"], round(x["rho_re"], 2), round(x["alpha_ON"], 2), round(x["alpha_OFF"], 3)) for x in sel][:10])
sel = [x for x in out if -3.6 <= x["rho_open_re"] <= -1.3]
print("rows with rho_ON/OPEN in [-3.5,-1.3]:", len(sel), [(x["id"], round(x["rho_open_re"], 2), round(x["alpha_ON"], 2), round(x["alpha_OFF"], 3)) for x in sel][:10])
import statistics as st
print("Dimag_OFF range %.2f..%.2f, Dimag_ON %.2f..%.2f, alpha_ON median %.2f max %.2f, alpha_OFF median %.3f max %.3f" % (min(x["Dimag_OFF"] for x in out), max(x["Dimag_OFF"] for x in out), min(x["Dimag_ON"] for x in out), max(x["Dimag_ON"] for x in out), st.median(x["alpha_ON"] for x in out), max(x["alpha_ON"] for x in out), st.median(x["alpha_OFF"] for x in out), max(x["alpha_OFF"] for x in out)))
print("alpha_NOCUT range %.2f..%.2f; alpha_OPEN max %.2f" % (min(x["alpha_NOCUT"] for x in out), max(x["alpha_NOCUT"] for x in out), max(x["alpha_OPEN"] for x in out)))
