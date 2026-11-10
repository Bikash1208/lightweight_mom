"""ITERATION 1 study v2 (model, relative trends; MODEL_ONLY / NOT_HFSS_VALIDATED), user decision 2026-11-10: reliability first.
    ROT180=1 SUB6_MESH=96 python jobs1_iter2_study.py <out.json> <alpha_off> <alpha_on>
k = 1 (single island cell, proven) primary, k = 2 optional; N = 80 / 100 / 120 / 150; circular board 220 / 280 / 360 mm (R 110 / 140 / 180); partitions inner / outer / sector_axis / sector_perp / field_weighted;
rot180 frame (A designed at phi 180, B at phi 0 in the model, ROT180=1).  Targets (from the full patterns): |G_A - G_B| <= 3 dB, B suppression G_B(OFF) - G_B(ON) >= 10 dB (12 desired), |G_A(ON) - G_A(OFF)| <= 1 dB, then max min(G_A, G_B).
Robustness of the best candidates: switched alpha x0.85 / x1.15 (ON ratio kept) and static alpha x0.85 / x1.15 (amplitude x sqrt); only designs meeting the targets in ALL five cases survive.
Everything as in jobs1_iter1_study.py (same engine, assumptions and gain model).
"""
import json
import math
import os
import sys
from pathlib import Path

os.environ.setdefault("ROT180", "1")
os.environ.setdefault("SUB6_MESH", "96")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1]))
import numpy as np

import jobs1_iter1_study as S

M, R4, B = S.M, S.R4, S.B


def full(bd, e):
    cur = e["_cur"]
    (dAo, dBo), _ = S.pattern_D(S.green_ref, bd.x, bd.y, *cur["off"], [S.BEAM_A, S.BEAM_B])
    (dAn, dBn), _ = S.pattern_D(S.green_ref, bd.x, bd.y, *cur["on"], [S.BEAM_A, S.BEAM_B])
    eo, en = S.LAUNCH * e["p_conv_off"], S.LAUNCH * e["p_conv_on"]
    g = lambda d, eta: d + 10 * math.log10(max(eta, 1e-9))
    r = dict(G_A_off=g(dAo, eo), G_B_off=g(dBo, eo), G_A_on=g(dAn, en), G_B_on=g(dBn, en), D_A_off=dAo, D_B_off=dBo, eta_off=eo, eta_on=en)
    r["balance"] = r["G_A_off"] - r["G_B_off"]
    r["supp"] = r["G_B_off"] - r["G_B_on"]
    r["dA"] = r["G_A_on"] - r["G_A_off"]
    r["meets"] = bool(abs(r["balance"]) <= 3.0 and r["supp"] >= 10.0 and abs(r["dA"]) <= 1.0)
    r["meets_12dB"] = bool(r["supp"] >= 12.0 and r["meets"])
    r["score"] = min(r["G_A_off"], r["G_B_off"])
    return r


def main():
    outp = sys.argv[1]
    a_off, a_on = float(sys.argv[2]), float(sys.argv[3])
    r_on = math.sqrt(a_on / a_off)
    T = R4.tables(); tab = M.static_table(T); bf = float(np.mean(T["st"]["b"]))
    S.green_ref = B.get_side(0.2).green
    radii = [0.110, 0.140, 0.180]
    ns = [int(v) for v in os.environ.get("JOB_NS", "80,100,120,150").split(",")]
    ks = [int(v) for v in os.environ.get("JOB_KS", "1,2").split(",")]
    places = ["inner", "outer", "sector_axis", "sector_perp", "field_weighted"]
    boards = {}
    cand = []
    for R in radii:
        bd = boards[R] = S.Board(R, S.green_ref, T, tab, bf)
        for n in ns:
            for k in ks:
                for place in places:
                    best = None
                    for cap in S.CAPS:
                        e = bd.evaluate(place, n, k, a_off, r_on, cap, need_fields=True)
                        if e is None:
                            continue
                        # relaxed level pre-filter score: violation of the three targets (level based)
                        viol = max(0.0, abs(e["balance_A_minus_B"]) - 3.0) + max(0.0, 10.0 - e["suppression_dB"]) + max(0.0, abs(e["dA_off_to_on"]) - 1.0)
                        key = (-viol, min(e["A_off"], e["B_off"]))
                        if best is None or key > best[0]:
                            best = (key, e)
                    if best is None:
                        continue
                    e = best[1]
                    e["R_mm"] = R * 1e3
                    cand.append(e)
            print("R", R, "N", n, "done", flush=True)
    # pre-filter by level violation, then full patterns for the 40 best-by-(violation, level)
    def viol(e):
        return max(0.0, abs(e["balance_A_minus_B"]) - 3.0) + max(0.0, 10.0 - e["suppression_dB"]) + max(0.0, abs(e["dA_off_to_on"]) - 1.0)
    cand.sort(key=lambda e: (viol(e) > 1.5, viol(e) if viol(e) > 1.5 else -min(e["A_off"], e["B_off"])))
    top = cand[:40]
    rows = []
    for e in top:
        bd = boards[e["R_mm"] / 1e3]
        f = full(bd, e)
        rows.append(dict(R_mm=e["R_mm"], diameter_mm=2 * e["R_mm"], N=e["N"], k=e["k"], place=e["place"], cap=e["cap"], n_switched_cells=e["n_switched_cells"], vias=e["vias"], level_balance=e["balance_A_minus_B"], **f))
    rows.sort(key=lambda r: (not r["meets"], -r["score"]))
    # robustness for the best ten that meet the targets (else the ten closest)
    rob = []
    for r in rows[:10]:
        bd = boards[r["R_mm"] / 1e3]
        base = bd.evaluate(r["place"], r["N"], r["k"], a_off, r_on, r["cap"], need_fields=True)
        els = base["_els"]
        scen = {}
        for name, sw_s, st_s, ron_s in (("sw x0.85", 0.85, 1.0, 1.0), ("sw x1.15", 1.15, 1.0, 1.0), ("static x0.85", 1.0, 0.85, 1.0), ("static x1.15", 1.0, 1.15, 1.0), ("ON amplitude x1.2 (extra, not in the 5)", 1.0, 1.0, 1.2)):
            e = bd.evaluate(r["place"], r["N"], r["k"], a_off * sw_s, r_on * ron_s, r["cap"], need_fields=True, static_scale=st_s, elements=els)
            scen[name] = full(bd, e)
        r["robust"] = {k_: dict(balance=v["balance"], supp=v["supp"], dA=v["dA"], G_A=v["G_A_off"], G_B=v["G_B_off"], meets=v["meets"]) for k_, v in scen.items()}
        r["robust_all_meet"] = bool(r["meets"] and all(v["meets"] for kk, v in scen.items() if not kk.startswith("ON amplitude")))
        r["robust_incl_ON_x1.2"] = bool(r["meets"] and all(v["meets"] for v in scen.values()))
        r["worst_balance"] = max([abs(r["balance"])] + [abs(v["balance"]) for v in scen.values()])
        r["worst_supp"] = min([r["supp"]] + [v["supp"] for v in scen.values()])
        r["worst_dA"] = max([abs(r["dA"])] + [abs(v["dA"]) for v in scen.values()])
        r["worst_score"] = min([r["score"]] + [min(v["G_A_off"], v["G_B_off"]) for v in scen.values()])
        rob.append(r)
        print("%d mm k%d N%d %-14s cap %.1f | G %.1f/%.1f bal %.1f supp %.1f dA %+.1f | meets %s robust %s | worst bal %.1f supp %.1f dA %.1f" % (r["diameter_mm"], r["k"], r["N"], r["place"], r["cap"], r["G_A_off"], r["G_B_off"], r["balance"], r["supp"], r["dA"], r["meets"], r["robust_all_meet"], r["worst_balance"], r["worst_supp"], r["worst_dA"]), flush=True)
    Path(outp).write_text(json.dumps(dict(label="MODEL_ONLY / NOT_HFSS_VALIDATED iteration 1 v2", alpha_off=a_off, alpha_on=a_on, r_on=r_on, rot180_frame=os.environ.get("ROT180") == "1", ranked=rows, robust=rob,
                                          all_candidates=[{k_: v for k_, v in e.items() if not k_.startswith("_")} for e in cand]), indent=1, default=float), encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
