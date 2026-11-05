"""Azimuthal harmonic content (e^{j n phi}) of the HFSS far field of board_K2_default on the rings theta 12 / 24 deg: common mode (E1+E2)/2, difference mode (E1-E2)/2, and each state.  Not a model.  -> results/sidequest_hybrid/hfss_ffd_harmonics.json"""
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import sidequest_ffd_analysis as A

d = {k: A.load(p) for k, p in A.FILES.items()}
(e1t, e1p), (e2t, e2p) = d[1], d[2]
modes = {"state1": (e1t, e1p), "state2": (e2t, e2p), "common": ((e1t + e2t) / 2, (e1p + e2p) / 2), "difference": ((e1t - e2t) / 2, (e1p - e2p) / 2)}
out = {}
for th in (12, 24):
    i = th // 2
    for nm, (et, ep) in modes.items():
        a, b = et[i, :-1], ep[i, :-1]
        n = a.size
        ph = np.radians(np.arange(-180, 180, 2.0))
        p = {}
        for k in range(-8, 9):
            e = np.exp(-1j * k * ph)
            p[k] = float((abs(np.sum(a * e) / n) ** 2 + abs(np.sum(b * e) / n) ** 2))
        tot = sum(p.values())
        tot_all = float(np.mean(abs(a) ** 2 + abs(b) ** 2))
        frac = {k: round(v / tot_all, 3) for k, v in p.items()}
        top = sorted(frac.items(), key=lambda t: -t[1])[:6]
        out[f"theta{th}_{nm}"] = {"fraction_by_harmonic_n": frac, "top": top, "captured_by_|n|<=8": round(tot / tot_all, 3)}
        print(f"theta {th} {nm:10s} top harmonics (n, fraction):", top, "captured", round(tot / tot_all, 3))
Path("results/sidequest_hybrid/hfss_ffd_harmonics.json").write_text(json.dumps(out, indent=1), encoding="utf-8", newline=chr(10))
