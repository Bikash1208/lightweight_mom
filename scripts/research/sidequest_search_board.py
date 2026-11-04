"""K = 2 board for the two cell-search routes at <= 766 diodes (1 line + inverter per axis), same engine as sidequest_coded_improve (RHCP (25, 0) / LHCP (25, 180), hexagon side 150 mm, centre probe, 6.4 mm superpixel grid).
Primary RO3010 0.635 mm (eps_r 11.2, 10.05 GHz): switched cell = notch superpixel (alpha_beam from the MoM cell, first-order formula about the cell's own mean), static cells = E1-type tensor cells (any axis angle, |m| <= m_E1 =
sqrt(alpha_E1 / A0)), alpha_E1 = 13.1 Np/m (E1 library first-order maximum; 9.2 Floquet-ratio sensitivity).  Secondary RT5880: rim disc cell (diode position 0.8 mm) with the round-11 static/pair option sets.
MODEL_ONLY. Output results/sidequest_hybrid/cell_search_board.json."""
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_coded_improve_run as R  # builds the RT5880 board / library at the default 10.54 GHz
import sidequest_coded_improve as I

C0 = 299792458.0
LOG = {}
GS = (1.0, 1.5, 2.0, 2.7)
LIMIT_SLL = -12.7


def configure(f0, a0, aa0, aa_bg, beta_k0, launch):
    I.st.K0 = 2 * np.pi * f0 / C0
    I.A0, I.AA0, I.AA_BG, I.BETA_K0, I.LAUNCH0 = a0, aa0, aa_bg, beta_k0, launch
    I.BETA = beta_k0 * I.st.K0
    I.LAM0 = C0 / f0
    I.DMAX_DBI = 10 * np.log10(4 * np.pi * I.AREA / I.LAM0**2)


def go(label, bd, vs, pairs, tens=None, mu=0.0, launch=None, e1_scale=None, gs=GS):
    t = time.time()
    best = None
    for g in gs:
        r, aux = I.full_run(vs, pairs, [0.0, 0.0], mu, g, launch or I.LAUNCH0, (1.0, 1.0), None, 766, bd, None, 1.0, tens, e1_scale)
        gm = np.mean([x["gain_dbi"] for x in r["states"]])
        sll = max(x["SLL_db"] for x in r["states"])
        ok = sll <= LIMIT_SLL
        if best is None or (ok, gm) > (best[0], best[1]):
            best = (ok, gm, r, aux, g)
    _, _, r, aux, g = best
    s = r["states"]
    LOG[label] = {"gain": [round(x["gain_dbi"], 2) for x in s], "SLL": [round(x["SLL_db"], 1) for x in s], "AR": [round(x["AR_db"], 2) for x in s], "XPD": [round(x["XPD_db"], 1) for x in s],
                  "ap_eff_pct": round(100 * r["ap_eff"][0], 2), "edge_pct": round(100 * r["edge"], 1), "SW_to_beam_pct": round(100 * r["SW_to_beam"], 1), "static_pct": round(100 * r["static_frac"], 1),
                  "diodes": r["diodes"], "lines": r["lines_with_inverters"], "g": g, "e1_cell_pct": round(100 * r["e1_cell_frac"], 1), "ON_loss": round(r["ON_loss"], 3)}
    print(label, LOG[label], round(time.time() - t), flush=True)
    return r, aux


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    VS_BIT, PAIRS_BIT = I.bit_sets()
    if which in ("all", "rt"):
        # secondary: RT5880 rim disc, Bloch cell, diode position 0.8 mm (alpha 1.84), round-11 option sets (I1 + I2 + I5)
        configure(10.54e9, 1.836, 0.094, 0.15, 1.079, 0.627)
        kw = dict(R.BASE)
        go("RT5880_bits", R.bd, VS_BIT, PAIRS_BIT)
        kw = R.ITEMS["I1"](kw)
        kw = R.ITEMS["I2"](kw)
        kw = R.ITEMS["I5"](kw)
        go("RT5880_I1I2I5", R.bd, kw["vs"], kw["pairs"], mu=kw["mu"])
    if which in ("all", "ro"):
        configure(10.05e9, 7.125, 2.517, 0.6, 1.14, 0.476)
        bd = I.Board()
        mmax = float(np.sqrt(13.1 / 7.125))
        tens = I.e1_options(mmax)
        go("RO3010_bits", bd, VS_BIT, PAIRS_BIT)
        r1, aux1 = go("RO3010_E1", bd, VS_BIT, PAIRS_BIT, tens=tens)
        go("RO3010_E1_loss", bd, VS_BIT, PAIRS_BIT, tens=tens, mu=0.05)
        go("RO3010_E1_launch0627", bd, VS_BIT, PAIRS_BIT, tens=tens, launch=0.627)
        go("RO3010_E1_floquet", bd, VS_BIT, PAIRS_BIT, tens=I.e1_options(float(np.sqrt(9.2 / 7.125))))
        _, _, m, nd, loss, alpha, aa = aux1
        scale, info = I.taper_scale(bd, alpha, aa)
        LOG["taper_info"] = [{k: (float(v) if not isinstance(v, bool) else v) for k, v in d.items()} for d in info]
        go("RO3010_E1_taper", bd, VS_BIT, PAIRS_BIT, tens=tens, e1_scale=scale)
    Path(f"results/sidequest_hybrid/cell_search_board_{which}.json").write_text(json.dumps(LOG, indent=1, default=float), encoding="utf-8", newline="\n")
