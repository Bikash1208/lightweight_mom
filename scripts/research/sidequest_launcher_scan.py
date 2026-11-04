"""Scan of the parasitic-pin launcher (2 quadrants: ON = pins in the back half; 4 quadrants: ON = pins outside the chosen 90-deg window). MODEL_ONLY."""
import json
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np

import sidequest_launcher_model as L


def run2(N, top, rho, z_on=None, z_off=None, rho_scale=1.0):
    on = [np.cos(np.radians((m + 0.5) * 360 / N)) < 0 for m in range(N)]
    r = L.solve(N, rho * 1e-3, top, on, z_on=z_on, z_off=z_off)
    fb, frac = L.fb_2q(r)
    return r, fb, L.sectors(r, 0, 90) / r["p_acc"]


def run4(N, top, rho):
    on = [abs((np.degrees(np.radians((m + 0.5) * 360 / N)) + 180) % 360 - 180) > 45 for m in range(N)]
    r = L.solve(N, rho * 1e-3, top, on)
    f = L.sectors(r, 0, 45)
    oth = max(L.sectors(r, c, 45) for c in (90, 180, 270))
    return r, float(10 * np.log10(f / max(oth, 1e-30))), f / r["p_acc"], f / (L.sectors(r, 0, 180) + 1e-30)


if __name__ == "__main__":
    rows = []
    for N in (4, 6, 8, 12):
        for top in np.arange(3.0, 9.01, 1.0):
            for rho in np.arange(4.0, 36.0, 1.0):
                r, fb, eta = run2(N, top, rho)
                rows.append({"N": N, "top": float(top), "rho": float(rho), "fb_db": fb, "s11_db": r["s11_db"], "eta_front": eta, "sw_over_acc": r["p_sw"] / r["p_acc"], "diode_loss": r["p_diode"] / r["p_acc"], "space_loss": r["p_space"] / r["p_acc"]})
    Path("results/sidequest_hybrid/launcher_scan2q.json").write_text(json.dumps(rows), encoding="utf-8", newline="\n")
    ok = [r for r in rows if r["fb_db"] >= 10 and r["s11_db"] <= -10]
    print("n feasible (FB>=10, S11<=-10):", len(ok), "of", len(rows))
    for N in (4, 6, 8, 12):
        c = [r for r in rows if r["N"] == N]
        bf = max(c, key=lambda r: r["fb_db"])
        cf = [r for r in c if r["fb_db"] >= 10 and r["s11_db"] <= -10]
        be = max(cf, key=lambda r: r["eta_front"]) if cf else None
        print("N", N, "best FB:", {k: round(v, 2) for k, v in bf.items()}, "| best eta with FB>=10,S11<=-10:", be and {k: round(v, 2) for k, v in be.items()})
