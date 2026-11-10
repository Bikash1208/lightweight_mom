"""FOLLOW-UP DIAGNOSTIC item 2: 5.8 GHz model pattern from the STATIC cells only and from the SWITCHED cells only (identity and rot180 convention) vs HFSS OFF / ON.  MODEL_ONLY.
    SUB6_MESH=96 python diag_split_static_switched.py <off.ffd> <on.ffd> <out.json>
Complex least squares on the RHCP field over theta <= 60 (all phi): HFSS_state ~ a * S_static + b * S_switched(state)  (a, b complex), residual fraction = incoherent part; also with the switched response of the OTHER state.
"""
import json
import math
import os
import sys
from pathlib import Path

os.environ.setdefault("SUB6_MESH", "96")
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import numpy as np

import diag_hfss_vs_model as G
import sub6_final_design as FD
import sub6_board_model4 as R4


def main():
    ffd_off, ffd_on, outp = sys.argv[1:4]
    H = {"off": G.load_ffd(ffd_off), "on": G.load_ffd(ffd_on)}
    fd = FD.final_design("cal")
    D, sa, ph = fd["D"], fd["sa"], fd["ph"]
    S = D["S"]
    fl = fd["flips"].astype(float)
    mask = D["mask"] & D["live"]
    live = D["live"]
    cur = {st: R4.currents(D, sa, ph, fd["direction"], fl, st)[:2] for st in ("off", "on")}
    comp = {}
    for st in ("off", "on"):
        jx, jy = cur[st]
        comp[("static", st)] = G.model_field(S, np.where(mask, 0, jx), np.where(mask, 0, jy))
        comp[("switched", st)] = G.model_field(S, np.where(mask, jx, 0), np.where(mask, jy, 0))
    comp = {k: (G.pad_to_180(v[0]), G.pad_to_180(v[1])) for k, v in comp.items()}
    n = 31          # theta 0 ... 60 (2 deg rows)

    def rh(e, kind):
        et, ep = e
        if kind == "rot180":
            et, ep = G.transform(et, ep, "rot", 180)
        return G.circ(et, ep)[0][:n, :180]

    def power_share(e):
        et, ep = e
        return float(np.sum((np.abs(et[:n, :180]) ** 2 + np.abs(ep[:n, :180]) ** 2) * np.sin(np.radians(G.TH_DEG[:n]))[:, None]))

    res = dict(label="MODEL_ONLY split static / switched", power_static_vs_switched_OFF=dict(static=power_share(comp[("static", "off")]), switched=power_share(comp[("switched", "off")])),
               power_static_vs_switched_ON=dict(static=power_share(comp[("static", "on")]), switched=power_share(comp[("switched", "on")])), rows=[])
    for kind in ("id", "rot180"):
        for st in ("off", "on"):
            h = G.circ(*H[st])[0][:n, :180]
            hv = h.ravel()
            # region of the HFSS beam
            tt, pp = np.meshgrid(G.TH_DEG[:n], G.PH_DEG[:180], indexing="ij")
            near = (np.hypot(tt - 24, ((pp - 178 + 180) % 360) - 180) <= 12).ravel()
            for swmode in ("own", "other"):
                sw_state = st if swmode == "own" else {"off": "on", "on": "off"}[st]
                a = rh(comp[("static", st)], kind).ravel()
                b = rh(comp[("switched", sw_state)], kind).ravel()
                A = np.stack([a, b], 1)
                coef, *_ = np.linalg.lstsq(A, hv, rcond=None)
                fit = A @ coef
                fa = coef[0] * a
                fb = coef[1] * b
                row = dict(convention=kind, hfss_state=st, switched_response=swmode, coef_static=complex(coef[0]), coef_switched=complex(coef[1]),
                           ratio_abs_static_over_switched=float(abs(coef[0]) / max(abs(coef[1]), 1e-30)), residual_fraction_power=float(np.sum(np.abs(hv - fit) ** 2) / np.sum(np.abs(hv) ** 2)),
                           static_only_ccorr=float(abs(np.vdot(a, hv)) / math.sqrt(np.vdot(a, a).real * np.vdot(hv, hv).real)), switched_only_ccorr=float(abs(np.vdot(b, hv)) / math.sqrt(np.vdot(b, b).real * np.vdot(hv, hv).real)),
                           frac_beam_region_power_explained_by_static=float(np.sum(np.abs(fa[near]) ** 2) / np.sum(np.abs(hv[near]) ** 2)), frac_beam_region_power_explained_by_switched=float(np.sum(np.abs(fb[near]) ** 2) / np.sum(np.abs(hv[near]) ** 2)),
                           frac_beam_region_residual=float(np.sum(np.abs((hv - fit)[near]) ** 2) / np.sum(np.abs(hv[near]) ** 2)),
                           frac_total_power_in_beam_region_hfss=float(np.sum(np.abs(hv[near]) ** 2) / np.sum(np.abs(hv) ** 2)))
                # peaks of the model components (after the convention transform), RHCP power, theta 6..60
                for nm, e in (("static", rh(comp[("static", st)], kind)), ("switched", rh(comp[("switched", sw_state)], kind))):
                    p = np.abs(e) ** 2
                    p[:3] = 0
                    i = np.unravel_index(np.argmax(p), p.shape)
                    row["peak_" + nm] = (int(G.TH_DEG[i[0]]), int(G.PH_DEG[i[1]]))
                res["rows"].append(row)
    Path(outp).write_text(json.dumps(res, indent=1, default=str), encoding="utf-8", newline="\n")
    print(json.dumps(res["power_static_vs_switched_OFF"]), json.dumps(res["power_static_vs_switched_ON"]))
    for r in res["rows"]:
        print("%-7s %-3s sw=%-5s |a|/|b| %.2f resid %.2f | static-only ccorr %.2f switched-only ccorr %.2f | beam region: static %.2f switched %.2f resid %.2f (beam region = %.2f of HFSS power) | peaks static %s switched %s" % (
            r["convention"], r["hfss_state"], r["switched_response"], r["ratio_abs_static_over_switched"], r["residual_fraction_power"], r["static_only_ccorr"], r["switched_only_ccorr"],
            r["frac_beam_region_power_explained_by_static"], r["frac_beam_region_power_explained_by_switched"], r["frac_beam_region_residual"], r["frac_total_power_in_beam_region_hfss"], r["peak_static"], r["peak_switched"]))


if __name__ == "__main__":
    main()
