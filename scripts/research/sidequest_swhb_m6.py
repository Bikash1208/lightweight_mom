"""M6 check of SWHB against the old HFSS boards (r 2.30 RO3006 v2 and add-a-beam (b)).  MODEL_ONLY.   python sidequest_swhb_m6.py"""
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_sll_def as SLD
import sidequest_swhb as S

RES = ROOT / "results/sidequest_hybrid"
P = 6.4e-3


def old_cells(name, cal=1.0):
    d = json.loads((RES / name).read_text())
    lib = json.loads((RES / "ro3006_r2.3_mesh96_lib.json").read_text())["rows"]

    def Y(kind, fr):
        r = next(x for x in lib if x["kind"] == kind and abs(x["fr"] - fr) < 1e-9)
        return complex(*r["Y_mS"]) * 1e-3

    yon_ref = Y("ON", 1.0)
    xy, ys = [], [[], []]
    for ix, iy, typ, kind, f in d["cells"]:
        xy.append((ix * P, iy * P))
        if typ == "ST":
            y0 = y1 = Y(kind, f)
        else:
            y0, y1 = Y("open", f), Y("ON", f)
        for s, y in ((0, y0), (1, y1)):
            ys[s].append(yon_ref + cal * (y - yon_ref))
    return np.array(xy), [np.array(a) for a in ys], d


def summarize(env, xy, y, eta=1.0):
    r = S.run_state(env, xy, y, lambda kr: 0j, (0j, 0j), P, eta)
    if not r["ok"]:
        return r
    u = r["u"]
    th, ph = r["th"], r["ph"]
    prad = r["P_rad_unit"]
    TH, _ = np.meshgrid(th, ph, indexing="ij")
    d_peak = 4 * math.pi * u.max() / prad
    hi = (np.degrees(TH) >= 60)
    frac_hi = float(np.sum(u * np.sin(TH) * hi) / np.sum(u * np.sin(TH)))
    return {"beta_k0": r["beta_k0"], "alpha_loss": r["alpha_loss"], "kappa": r["kappa"], "P_rad": prad, "P_edge": r["P_edge"], "P_abs": 1 - prad - r["P_edge"], "D_peak_dBi": 10 * math.log10(d_peak), "G_peak_dBi": 10 * math.log10(4 * math.pi * u.max()), "frac_60_90": frac_hi,
            "edge_dB_vs_cyl": 10 * math.log10(max(r["P_edge"], 1e-12)), "peak_dir": [float(np.degrees(th[np.unravel_index(u.argmax(), u.shape)[0]])), float(np.degrees(ph[np.unravel_index(u.argmax(), u.shape)[1]]))]}


if __name__ == "__main__":
    env = S.Env(6.15, 0.002, 1.27e-3)
    for name in ("ro3006_board_design_r2.3_v2.json", "ro3006_board_design_r2.3_addbeam_b.json"):
        for cal in (1.0, 0.75):
            xy, ys, d = old_cells(name, cal)
            for s in (0, 1):
                r = summarize(env, xy, ys[s])
                print(name[-20:], "cal", cal, "state", s + 1, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
