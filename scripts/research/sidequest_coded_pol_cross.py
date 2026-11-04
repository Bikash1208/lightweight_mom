"""K = 2 RHCP/LHCP design with the SOLVED crossed two-diode cell (coded_cross_scan.json, mesh 64, x and y diode 0/90 deg). Same model and optimiser as sidequest_coded_pol.py (route R2)."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_arch_study as st
import sidequest_coded_cell as cc
import sidequest_coded_pol as P

VARIANT = sys.argv[1] if len(sys.argv) > 1 else "rim"  # "rim": slits 5.6 mm long (0.2 mm rim), "through": slits cut through the disc
fn = "coded_cross_rim_scan.json" if VARIANT == "rim" else "coded_cross_scan.json"
scan = json.loads(Path("results/sidequest_hybrid/" + fn).read_text(encoding="utf-8"))
a = [o for o in scan if o["mesh"] == 64 and not o["four"]][0]["analysis"]
delta = complex(*a["dXX_mS"]) * 1e-3
iso = complex(*a["Y_iso_mS"]) * 1e-3
P.BETA_K0, P.A0, P.AABS0 = cc.alpha_pair(iso - delta / 2, iso + delta / 2)
P.BETA = P.BETA_K0 * st.K0
P.DELTA = delta

if __name__ == "__main__":
    print("cell: beta/k0 %.3f alpha0 %.3f alpha_abs0 %.3f |dY| %.2f mS" % (P.BETA_K0, P.A0, P.AABS0, abs(delta) * 1e3))
    o, base, rows = P.run("R2")
    out = []
    for lam in (0.0, 0.3, 0.5, 0.8):
        rr = [r for r in rows if r["lam"] == lam and r["loss_db"] < 1.0]
        r = max(rr, key=lambda r: r["static_frac"])
        f = P.full("R2", lam, r["eps_B_deg"], f"lam{lam}")
        out.append(f)
        s = f["states"]
        print("R2x", lam, "eps_B %.1f static %.3f diodes %d lines %d/%d" % (r["eps_B_deg"], f["static_frac"], f["diodes"], f["lines_with_inverters"], f["lines_no_inverters"]),
              "gain", [round(x["gain_dbi"], 1) for x in s], "ideal", round(s[0]["G_ideal_dbi"], 1), "SLL", [round(x["SLL_db"], 1) for x in s], "AR", [round(x["AR_db"], 1) for x in s], "XPD", [round(x["XPD_db"], 1) for x in s], "eta", round(s[0]["converted_beam_power"], 4), flush=True)
    Path("results/sidequest_hybrid/coded_pol_cross_%s_results.json" % VARIANT).write_text(json.dumps({"cell": {"beta_k0": P.BETA_K0, "alpha0": P.A0, "alpha_abs0": P.AABS0}, "designs": out}, indent=1, default=float), encoding="utf-8", newline="\n")
