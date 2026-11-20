"""Shared loaders for the S4 baseline track (Phase 1: stored data only, no AEDT). MODEL/STORED_HFSS_2026R1 labels in REPORT.md."""
import json, math, sys, os
from pathlib import Path
import numpy as np
os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")
sys.dont_write_bytecode = True

REF = Path(r"C:\work\lwmom_s4_ref")
GATE = REF / "results/sub6_meander_model/full_board_gate"
sys.path.insert(0, str(REF / "scripts/research"))
from full_board_gate_v31_extract import ffd, direction_separation  # reuse reader (imports pyaedt at top; fine)

_ms = json.loads((GATE / "v31/metrics_summary.json").read_text())["states"]
_ce = json.loads((GATE / "v31_ceiling/metrics.json").read_text())
_sh = json.loads((GATE / "v31_short/metrics.json").read_text())["states"]
STATES = {
    "BOTH": (GATE / "v31/BOTH.ffd", _ms["BOTH"]),
    "A_ONLY": (GATE / "v31/A_ONLY.ffd", _ms["A_ONLY"]),
    "B_ONLY": (GATE / "v31/B_ONLY.ffd", _ms["B_ONLY"]),
    "A_ONLY_IDEAL": (GATE / "v31_ideal/A_ONLY_IDEAL.ffd", _ce["a_only_ideal_234"]),
    "B_ALL_OFF_IDEAL": (GATE / "v31_ceiling/B_ALL_OFF_IDEAL.ffd", _ce["b_all_off_ideal_258"]),
    "ALL_PLAIN": (GATE / "v31_ceiling/ALL_PLAIN.ffd", _ce["all_plain_516"]),
    "B_ALL_SHORT": (GATE / "v31_short/B_ALL_SHORT.ffd", _sh["B_ALL_SHORT"]),
}
_cache = {}

def load(state):
    """Return dict: theta, phi, Et, Ep (complex, gain-scaled so |E|^2*... see below), gain_R, gain_L (dBi grids), meta.
    Scale: HFSS peak gain of RHCP anchors absolute level; same scale applied to both hands (consistency checked)."""
    if state in _cache: return _cache[state]
    path, m = STATES[state]
    th, ph, et, ep = ffd(path)
    plus, minus = np.abs(et + 1j*ep)**2, np.abs(et - 1j*ep)**2
    rh, lh = (plus, minus) if m["ffd_rhcp_convention"] == "et_plus_j_ep" else (minus, plus)
    sc = {}
    for hand, p in (("RHCP", rh), ("LHCP", lh)):
        i = int(np.argmin(abs(th - m[hand]["peak_theta_deg"]))); j = int(np.argmin(abs(ph - m[hand]["peak_phi_deg"])))
        sc[hand] = m[hand]["peak_gain_db"] - 10*math.log10(p[i, j])
    assert abs(sc["RHCP"] - sc["LHCP"]) < 0.02, (state, sc)   # same scale for both hands => convention consistent
    s = sc["RHCP"]
    # amplitude scale so that |E|^2 -> gain linear for each hand: G_hand = 10^(s/10) * power_hand
    amp = 10**(s/20)
    d = dict(theta=th, phi=ph, Et=et*amp, Ep=ep*amp, rh=rh*10**(s/10), lh=lh*10**(s/10), m=m, scale_db=s,
             rhcp_conv=m["ffd_rhcp_convention"])
    _cache[state] = d
    return d

def sphere_weights(theta, phi):
    """dOmega weights (trapezoid) on the (theta,phi) grid; phi includes duplicate 360 column."""
    t = np.deg2rad(theta); p = np.deg2rad(phi)
    wt = np.zeros_like(t); wt[1:-1] = (t[2:] - t[:-2])/2; wt[0] = (t[1]-t[0])/2; wt[-1] = (t[-1]-t[-2])/2
    wp = np.zeros_like(p); wp[1:-1] = (p[2:] - p[:-2])/2; wp[0] = (p[1]-p[0])/2; wp[-1] = (p[-1]-p[-2])/2
    return (np.sin(t)*wt)[:, None] * wp[None, :]

def db(x): return 10*np.log10(np.maximum(x, 1e-300))
