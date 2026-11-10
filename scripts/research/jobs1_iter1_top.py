"""top-3 combinations per N from iter1_study.json, re-evaluated with full patterns (D, gains).  python jobs1_iter1_top.py <study.json> <out.json> [r_on]"""
import json, math, os, sys
from pathlib import Path
os.environ.setdefault("SUB6_MESH", "96")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1]))
import numpy as np
import jobs1_iter1_study as S
M, R4, B = S.M, S.R4, S.B
st = json.load(open(sys.argv[1])); r_on = float(sys.argv[3]) if len(sys.argv) > 3 else st["r_on_amplitude"]
T = R4.tables(); tab = M.static_table(T); bf = float(np.mean(T["st"]["b"])); green = B.get_side(0.2).green
out = {}
boards = {}
for N in (40, 50):
    rows = [r for r in st["rows"] if r["N"] == N]
    bal = [r for r in rows if r["balanced_within_3dB"]]
    pool = sorted(bal, key=lambda r: -min(r["A_off"], r["B_off"]))
    seen, top = set(), []
    for r in pool:
        key = (r["R_mm"], r["k"], r["place"])
        if key in seen: continue
        seen.add(key); top.append(r)
        if len(top) == 3: break
    near = sorted([r for r in rows if not r["balanced_within_3dB"]], key=lambda r: abs(r["balance_A_minus_B"]))[:2]
    res = []
    for r in top + near:
        R = r["R_mm"] / 1e3
        if R not in boards: boards[R] = S.Board(R, green, T, tab, bf)
        bd = boards[R]
        e = bd.evaluate(r["place"], N, r["k"], r["alpha_cell"], r_on, r["cap"], need_fields=True)
        (dAoff, dBoff), _ = S.pattern_D(green, bd.x, bd.y, *e["_cur"]["off"], [S.BEAM_A, S.BEAM_B])
        (dAon, dBon), _ = S.pattern_D(green, bd.x, bd.y, *e["_cur"]["on"], [S.BEAM_A, S.BEAM_B])
        eta_off, eta_on = S.LAUNCH * e["p_conv_off"], S.LAUNCH * e["p_conv_on"]
        g = lambda D, eta: D + 10 * math.log10(max(eta, 1e-9))
        e.pop("_cur")
        e.update(R_mm=r["R_mm"], D_A_off_dBi=dAoff, D_B_off_dBi=dBoff, D_A_on_dBi=dAon, D_B_on_dBi=dBon, eta_rad_off=eta_off, eta_rad_on=eta_on,
                 G_A_off_dBi=g(dAoff, eta_off), G_B_off_dBi=g(dBoff, eta_off), G_A_on_dBi=g(dAon, eta_on), G_B_on_dBi=g(dBon, eta_on), balanced=r["balanced_within_3dB"])
        res.append(e)
    out[str(N)] = res
    print("N", N)
    for e in res:
        print(" R%3d k%d %-14s cap %.1f | A/B off %.1f/%.1f bal %.1f supp %.1f dA(off->on) %+.1f | G_A %.1f G_B %.1f dBi (on: A %.1f B %.1f) eta %.2f/%.2f | cells %d diodes %d vias %d %s" % (e["R_mm"], e["k"], e["place"], e["cap"], e["A_off"], e["B_off"], e["balance_A_minus_B"], e["suppression_dB"], e["dA_off_to_on"], e["G_A_off_dBi"], e["G_B_off_dBi"], e["G_A_on_dBi"], e["G_B_on_dBi"], e["eta_rad_off"], e["eta_rad_on"], e["n_switched_cells"], e["diodes"], e["vias"], "" if e["balanced"] else "(NOT balanced)"))
Path(sys.argv[2]).write_text(json.dumps(out, indent=1, default=float), encoding="utf-8", newline="\n")
