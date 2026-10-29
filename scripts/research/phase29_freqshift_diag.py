"""Phase 29: frequency-shift check of the two-feed forward model (MODEL_ONLY, analysis only, no re-synthesis). The two-feed 150 mm design (psi2 = 180) is fixed (cell slots and orientations); at frequency f
the model is re-evaluated with k0 -> k0 f / 10.54 GHz in the radiation (Green function, far-field phases) and in the surface-wave march, with the cell table (alpha, beta/k0, kappa) frozen at its 10.54 GHz values
(non-dispersive-cell approximation; the table exists at 10.54 GHz only). Mode A: beta/k0 fixed (beta*k0 scales with f). Mode B: the phase constant is shifted so beta*k0 is unchanged (the design-frequency surface wave) while
the radiation is at f. Also reports, for mode A, the gain-peak frequency. Writes results/phase29/freqshift_diag.json."""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts/research"))
RES = ROOT / "results/phase29"
F0 = 10.54e9


def run(freqs_ghz, psi2=180.0):
    import phase29_tensor_ladder as L
    import phase29_weakmod as W
    from phase29_pointB_sweep import port_metrics

    from src.frequency import free_space_wavenumber
    from src.green_function import LayeredGreenFunction
    from src.hologram import shared_aperture_multiport as sa
    from src.hologram import shared_aperture_tensor as st
    from src.hologram import tensor_cp_design as t

    L.init()
    lut, table = L._S["lut"], L._S["table"]
    d0 = W.design("n2", [psi2], {"mode": "table", "M": None}, 240)
    slot, ia = d0["slot"], np.rint(d0["axis_rad"] / lut.step).astype(int)
    poly = sa.hexagon_m(d0["r_mm"])
    cxy = np.c_[d0["x_m"], d0["y_m"]]
    th, ph = np.arange(0.0, 90.0 + 1e-9, 1.0), np.arange(0.0, 360.0, 1.0)
    out = {}
    K0_0, F_0 = t.K0, t.F
    for mode in ("A", "B"):
        for fg in freqs_ghz:
            s = fg * 1e9 / F0
            kfreq = free_space_wavenumber(fg * 1e9)
            k_march = K0_0 * s if mode == "A" else K0_0
            ana = [
                st.march_analysis_tensor(table, lut, poly, d0["feeds_m"][j], k_march, d0["beams"][j][0], d0["beams"][j][1], cxy, slot, ia, 240)
                for j in range(len(d0["beams"]))
            ]
            cur = [
                st.cell_current(ana[j], lut, table, d0["feeds_m"][j], d0["x_m"], d0["y_m"], slot, ia, d0["beams"][j][0], d0["beams"][j][1])
                for j in range(len(ana))
            ]
            d = dict(d0, J=cur, eta_conv=[a["eta_conv"] for a in ana])
            t.K0, t.F = kfreq, fg * 1e9
            try:
                F = sa.port_fields(d, LayeredGreenFunction(t.SUB, fg * 1e9, backing="pec"), th, ph)
                pm = port_metrics(F, th, ph, d0["beams"])
            finally:
                t.K0, t.F = K0_0, F_0
            out[f"{mode}_{fg:.2f}"] = {
                "gain": [round(p["peak_gain_dbi"], 1) for p in pm],
                "sll": [round(p["sll_db"], 1) for p in pm],
                "peak_theta_phi": [[round(x, 1) for x in p["peak_theta_phi"]] for p in pm],
                "eta_conv": [round(e, 3) for e in d["eta_conv"]],
            }
            print(mode, fg, json.dumps(out[f"{mode}_{fg:.2f}"]), flush=True)
    return out


if __name__ == "__main__":
    S = run([10.0, 10.3, 10.54, 10.7, 10.83, 11.0])
    (RES / "freqshift_diag.json").write_text(json.dumps(S, indent=1))
