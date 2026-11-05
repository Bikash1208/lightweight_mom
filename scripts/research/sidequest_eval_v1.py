"""Re-evaluate the delivered r 2.30 package design (v1, designed with a single beta/k0 1.404) with (a) the SLL metric of the earlier design runs (exclusion 1.2 HPBW of the actual beam, about 45 deg) and (b) the stricter metric of the v2 runs (exclusion 21 deg around each wanted beam), at the design beta and at the real per-state beta (1.33 / 1.26).  MODEL_ONLY."""
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import sidequest_ro3006_design as RD
import sidequest_ro3006_v2 as V

cells = json.loads(Path("results/sidequest_hybrid/ro3006_board_design_r2.3.json").read_text())["cells"]
lib = RD.load_lib(96, True)
np_, K, D, bd = RD.engine(lib, 0.71)
lv = lib["levels"]
pairs = lib["pairs"]
n = bd.n
m = np.zeros((2, 2, n))
nd = np.zeros(n)
for i, (ix, iy, kd, tk, fr) in enumerate(cells):
    if kd == "SW":
        m[0, :, i], m[1, :, i] = pairs[fr][0], pairs[fr][1]
        nd[i] = 2
    else:
        v = 1.0 if (tk == "open" and fr == 1.0) else (lv["solid"] if tk == "none" else (lv["strap_%s" % ("0.79" if fr < 1 else "1.0")] if tk == "strap" else lv["open_%s" % fr]))
        m[:, :, i] = v
wanted = {0: [0], 1: [1]}
k0 = D.st.K0
for label, betas in (("design beta 1.404 / 1.404", (1.404, 1.404)), ("real beta 1.33 / 1.26", (1.33, 1.26))):
    K.A0 = lib["A0"] if label.startswith("design") else sum(V.lib_for(True, betas)["A0_states"]) / 2
    ev_old = RD.evaluate(np_, K, D, bd, m, nd, 0.71) if label.startswith("design") else None
    ev = V.evaluate_gen(np_, K, D, bd, m, nd, 0.71, [b * k0 for b in betas], wanted)
    s = ev["states"]
    print(label, "| gain own S1 %.1f S2 %.1f | wrong-beam S1 %.1f S2 %.1f | SLL(21 deg excl.) S1 %.1f S2 %.1f" % (s[0]["beam_gain_dbi"]["b1_(25,0)"], s[1]["beam_gain_dbi"]["b2_(25,180)"], s[0]["beam_gain_dbi"]["b2_(25,180)"], s[1]["beam_gain_dbi"]["b1_(25,0)"], s[0]["SLL_db"], s[1]["SLL_db"]), "| edge %.0f%%" % (100 * ev["edge"]))
    if ev_old:
        so = ev_old["states"]
        print("   old metric (1.2 HPBW exclusion): SLL_ap S1 %.1f S2 %.1f, gain S1 %.1f S2 %.1f" % (so[0]["SLL_aperture_db"], so[1]["SLL_aperture_db"], so[0]["gain_design_dir"], so[1]["gain_design_dir"]))
