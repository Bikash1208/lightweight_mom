"""K = 2 RHCP/LHCP board result for each (rim, mesh) of the crossed rim cell (data results/sidequest_hybrid/coded_rim/rim*_m*.json from sidequest_coded_rim.py).
Per-axis dY = Y_xx(ON/OFF) - Y_xx(OFF/OFF) (other diode OFF); sheet mean = mean of (Yxx+Yyy)/2 over OFF/OFF and ON/OFF; alpha/beta from the repo first-order leakage (sidequest_coded_cell.alpha_pair).
Board: sidequest_coded_pol (route R2, lambda 0 and 0.5, eps_B scanned). Prints/saves alpha, dY, beta, static %, diodes, gain, SLL, AR, XPD."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_arch_study as st
import sidequest_coded_cell as cc
import sidequest_coded_pol as P

D = Path("results/sidequest_hybrid/coded_rim")


def cell_numbers(d):
    y = {k: np.array([complex(a, b) for a, b in v]).reshape(2, 2) * 1e-3 for k, v in d["states"].items()}
    dx = y["ON/OFF"][0, 0] - y["OFF/OFF"][0, 0]
    iso = np.mean([0.5 * (v[0, 0] + v[1, 1]) for k, v in y.items() if k in ("OFF/OFF", "ON/OFF")])
    bk, ab, aa = cc.alpha_pair(iso - dx / 2, iso + dx / 2)
    return dx, iso, bk, ab, aa


def board(bk, ab, aa, dx):
    P.BETA_K0, P.A0, P.AABS0, P.DELTA = bk, ab, aa, dx
    P.BETA = bk * st.K0
    out = {}
    for lam in (0.0, 0.5):
        best = None
        grid = np.radians(np.arange(0, 360, 22.5)) if lam else [0.0]
        base_amp = None
        for eb in grid:
            B, dd, fc = P.converge("R2", [0.0, eb], lam)
            amp = np.array(P.amplitude(B, P.realise(B, dd)))
            if base_amp is None and not lam:
                base_amp = amp
            sm = P.summarize(dd)
            if lam == 0:
                best = (0.0, np.degrees(eb))
            else:
                B0, d0, _ = P.converge("R2", [0.0, 0.0], 0.0)
                a0 = np.array(P.amplitude(B0, P.realise(B0, d0)))
                loss = float(-20 * np.log10(amp.min() / a0.min()))
                if loss < 1.0 and (best is None or sm["static_frac"] > best[0]):
                    best = (sm["static_frac"], np.degrees(eb))
        f = P.full("R2", lam, float(best[1]), f"lam{lam}")
        out[str(lam)] = {k: f[k] for k in ("static_frac", "diodes", "lines_with_inverters", "lines_no_inverters")} | {"eps_B_deg": float(best[1]), "states": f["states"]}
    return out


if __name__ == "__main__":
    res = []
    for fn in sorted(D.glob("rim*_m*.json"), key=lambda p: (json.loads(p.read_text())["rim"], json.loads(p.read_text())["mesh"])):
        d = json.loads(fn.read_text())
        if "ON/OFF" not in d["states"] or "OFF/OFF" not in d["states"]:
            continue
        if len(sys.argv) > 1 and int(sys.argv[1]) != d["mesh"]:
            continue
        dx, iso, bk, ab, aa = cell_numbers(d)
        b = board(bk, ab, aa, dx)
        s = b["0.5"]["states"]
        row = {"rim": d["rim"], "mesh": d["mesh"], "dY_mS": [dx.real * 1e3, dx.imag * 1e3], "Y_mean_mS": [iso.real * 1e3, iso.imag * 1e3], "beta_k0": bk, "alpha_beam": ab, "alpha_abs": aa, "board": b}
        res.append(row)
        print(d["rim"], d["mesh"], "dY %.2f%+.2fj alpha %.3f abs %.3f beta %.3f | lam.5 static %.3f diodes %d gain %.2f/%.2f SLL %.1f/%.1f AR %.2f/%.2f XPD %.1f/%.1f" % (
            *row["dY_mS"], ab, aa, bk, b["0.5"]["static_frac"], b["0.5"]["diodes"], s[0]["gain_dbi"], s[1]["gain_dbi"], s[0]["SLL_db"], s[1]["SLL_db"], s[0]["AR_db"], s[1]["AR_db"], s[0]["XPD_db"], s[1]["XPD_db"]), flush=True)
    Path("results/sidequest_hybrid/coded_rim_board.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline="\n")
