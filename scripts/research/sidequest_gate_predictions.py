"""Predictions for the gate test board (MODEL_ONLY) -> hfss_exports/sidequest_gate_test/predictions.json"""
import json
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np

import sidequest_gate_cell as gc
import sidequest_gate_match as gm
import sidequest_gate_reflection as gr
import sidequest_gate_ring as R

R_G = 0.016  # gate front face of the test board along the axis (cell edge at 16.0 mm)
YREF_T, YT_T = 1j * 7.49e-3, 1j * 5.44e-3  # actual test-board cells (mesh 64)
R.YREF = YREF_T


def stack(y_gate_xx, y_gate_yy, psi, trans=True):
    c, s = np.cos(np.radians(psi)), np.sin(np.radians(psi))
    y = y_gate_xx * c * c + y_gate_yy * s * s
    return gm.cascade(YREF_T, ([(YT_T, 1)] if trans else []) + [(y, 3)] + ([(YT_T, 1)] if trans else []))


def on_best(psi):
    c, s = np.cos(np.radians(psi)), np.sin(np.radians(psi))
    y = R.YXX_ON * c * c + R.YYY * s * s
    zr = gc.zm(gc.beta_complex(YREF_T))
    best = None
    for bg in gc.roots_all(y):
        zg = gc.zm(bg); th = bg * 3 * gc.P
        mg = np.array([[np.cos(th), 1j * zg * np.sin(th)], [1j * np.sin(th) / zg, np.cos(th)]])
        mt, _ = gm.abcd(YT_T, 1)
        m = mt @ mg @ mt
        den = m[0, 0] + m[0, 1] / zr + m[1, 0] * zr + m[1, 1]
        s21, s11 = 2 / den, (m[0, 0] + m[0, 1] / zr - m[1, 0] * zr - m[1, 1]) / den
        if best is None or abs(s21) > abs(best[1]):
            best = (s11, s21, bg)
    return best


def board_model(state, gam_b_phase=None):
    """Harmonic model of the test board: gate on the -x half (psi in (90, 270) deg). state 'OFF' (pass stack), 'ON_real' (leaky stack), 'ON_ideal' (|Gamma| = 1, T = 0, phase swept)."""
    br = gc.beta_complex(YREF_T)
    pr = R.psi_rel_of(gr.PSI)
    P = R.PSI_REL
    if state == "OFF":
        tab = [stack(R.YXX, R.YYY, p) for p in P]
        s11 = np.array([a for a, _ in tab]); s21 = np.array([b for _, b in tab])
    elif state == "ON_real":
        tab = [on_best(p) for p in P]
        s11 = np.array([a for a, _, _ in tab]); s21 = np.array([b for _, b, _ in tab])
    if state in ("OFF", "ON_real"):
        g11 = (np.interp(pr, P, s11.real) + 1j * np.interp(pr, P, s11.imag)) * np.exp(2j * br.real * gc.P)
        g21 = np.interp(pr, P, np.abs(s21)) * np.exp(1j * (np.interp(pr, P, np.unwrap(np.angle(s21))) + br.real * 5 * gc.P))
    else:
        g11 = np.full(gr.PSI.shape, np.exp(1j * gam_b_phase)); g21 = np.zeros(gr.PSI.shape, complex)
    gate = gr.OPEN == 0
    gam = np.where(gate, g11, 0.0).astype(complex)
    tr = np.where(gate, g21, 1.0).astype(complex)
    s11f, c = R.solve2(R_G, gam, tr)
    f = sum(cn * (1j) ** int(n) * np.exp(1j * n * gr.PSI) for n, cn in zip(gr.NS, c, strict=True))
    pf, pb = float(np.mean(gr.OPEN * np.abs(f) ** 2)), float(np.mean((1 - gr.OPEN) * np.abs(f) ** 2))
    return {"S11_db": float(20 * np.log10(abs(s11f) + 1e-12)), "P_front": pf, "P_back": pb, "front_over_back_db": float(10 * np.log10(pf / max(pb, 1e-12)))}


if __name__ == "__main__":
    pred = {"label": "MODEL_ONLY (single-mode line model + cylindrical-harmonic feed model, infinite plane: board edges, pin direct radiation and the finite 83 x 70 mm size are NOT modelled)",
            "cells_mm": {"ref_square": 5.8, "trans_square": 5.4, "gate": "dog-bone e1.6 s1.2 h6.0 g0.4 lt6.0 wd0.8"},
            "sheets_ms_10p54": {"ref_b": 7.49, "trans_b": 5.44, "gate_OFF_xx": "0.001+2.5j", "gate_OFF_yy": "6.6j", "gate_ON_xx": "100.1+68.0j", "gate_ON_yy": "6.6j"}}
    pred["stack_OFF"] = [{"psi_rel_deg": p, "RL_db": float(-20 * np.log10(abs(stack(R.YXX, R.YYY, p)[0]))), "IL_db": float(-20 * np.log10(abs(stack(R.YXX, R.YYY, p)[1])))} for p in (0, 15, 30, 45)]
    pred["gate_OFF_without_transition"] = [{"psi_rel_deg": p, "RL_db": float(-20 * np.log10(abs(stack(R.YXX, R.YYY, p, False)[0]))), "IL_db": float(-20 * np.log10(abs(stack(R.YXX, R.YYY, p, False)[1])))} for p in (0, 45)]
    pred["stack_ON_real"] = [{"psi_rel_deg": p, "S21_db": float(20 * np.log10(abs(on_best(p)[1]))), "S11_db": float(20 * np.log10(abs(on_best(p)[0]))), "mode_beta_k0": float(on_best(p)[2].real / gc.K0), "alpha_np_m": float(-on_best(p)[2].imag)} for p in (0, 15, 30, 45)]
    pred["board"] = {"OFF": board_model("OFF"), "ON_real": board_model("ON_real"),
                     "ON_ideal_block_phase_sweep": [{"phase_deg": float(d), **board_model("ON_ideal", np.radians(d))} for d in range(0, 360, 45)]}
    pred["falsifiers"] = [
        "ON-state surface |E| 1 mm above the board at x = -(25..38) mm vs OFF: model says only 3-5 dB lower (leaky block). A drop of >= 15 dB would falsify the model (the cell then blocks as originally claimed).",
        "OFF-state |E| along -x equal to +x within 1 dB (transition rows matched). A step of > 2 dB at the gate would falsify the pass-state match.",
        "S11 ON vs OFF: model -21 dB (leaky block) / -5 to -13 dB (ideal block, phase unknown); S11 better than -10 dB in both states is expected only for the leaky block.",
    ]
    Path("hfss_exports/sidequest_gate_test/predictions.json").write_text(json.dumps(pred, indent=1, default=float), encoding="utf-8", newline="\n")
    print(json.dumps({k: pred[k] for k in ("board", "stack_ON_real", "stack_OFF")}, indent=1, default=float)[:3500])
