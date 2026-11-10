"""Smallest board at which the model's static beam A dominates the feed's direct radiation (model, MODEL_ONLY / NOT_HFSS_VALIDATED).
All-static circular board (ring exclusion 12 mm), best static cap per radius, A designed at phi 0 (the HFSS frame is rot180).  G = D x 0.656 x converted fraction; feed reference = 7.6 dBi (HFSS 9-cell total peak, feed + rim).
    python jobs1_min_board.py <out.json>"""
import json, math, os, sys
from pathlib import Path
os.environ.setdefault("SUB6_MESH", "96")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1]))
import numpy as np
import jobs1_iter1_study as S
S.R_EXCL = 0.012
M, R4, B = S.M, S.R4, S.B
T = R4.tables(); tab = M.static_table(T); bf = float(np.mean(T["st"]["b"])); green = B.get_side(0.2).green
rows = []
for R in (0.040, 0.050, 0.060, 0.070, 0.080, 0.090, 0.100, 0.110, 0.140):
    bd = S.Board(R, green, T, tab, bf)
    best = None
    for cap in S.CAPS:
        d = bd.d[cap]
        (DA,), P = S.pattern_D(green, bd.x, bd.y, d["jxA"], d["jyA"], [S.BEAM_A])
        pconv = float(np.sum(d["ampA"][d["live"]] ** 2) * S.PITCH**2 / (2 * math.pi))
        eta = S.LAUNCH * pconv
        G = DA + 10 * math.log10(max(eta, 1e-9))
        r = dict(diameter_mm=2 * R * 1e3, cells=int(bd.x.size), live=int(d["live"].sum()), cap=cap, D_A_dBi=DA, eta_rad=eta, G_A_dBi=G, beam_power_fraction=eta, feed_direct_power_fraction=1 - S.LAUNCH)
        if best is None or G > best["G_A_dBi"]:
            best = r
    # peak position of the all-static pattern
    d = bd.d[best["cap"]]
    th = np.radians(np.arange(0, 61, 2.0)); ph = np.radians(np.arange(0, 360, 4.0))
    et, ep = S.far_grid(green, bd.x, bd.y, d["jxA"], d["jyA"], th, ph)
    rh = np.abs(S.copol_component(et, ep, "rhcp")); lh = np.abs(S.copol_component(et, ep, "lhcp"))
    i = np.unravel_index(np.argmax(rh), rh.shape)
    best.update(peak_theta=float(np.degrees(th[i[0]])), peak_phi=float(np.degrees(ph[i[1]])), XPD_dB=float(20 * math.log10(rh[i] / max(lh[i], 1e-12))))
    rows.append(best)
    print("%3d mm cells %3d cap %.1f  D %.1f  eta %.2f  G %.1f dBi  peak (%.0f, %.0f) XPD %.0f dB" % (best["diameter_mm"], best["cells"], best["cap"], best["D_A_dBi"], best["eta_rad"], best["G_A_dBi"], best["peak_theta"], best["peak_phi"], best["XPD_dB"]), flush=True)
Path(sys.argv[1]).write_text(json.dumps(rows, indent=1), encoding="utf-8")
