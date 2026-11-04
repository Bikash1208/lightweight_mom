"""Method 1 gain estimate: the SWITCHED LAUNCHER (probe + 12 pins, top 6 mm, ring radius 10.5 mm, back-half pins ON) feeds the 2-quadrant board from the centre; the launched SW
angular power P(psi) of the network model weights every ray. Same hologram/cap assumptions as `sidequest_twoprobe_model.py` (alpha cap 3.4 Np/m, ideal cells). MODEL_ONLY.
The launched wave is treated as cylindrical about the probe (the ring's phase-centre shift is ignored)."""
import json
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np

import sidequest_arch_study as st
import sidequest_launcher_model as L
import sidequest_twoprobe_model as tp

N, TOP, RHO = 12, 6.0, 10.5


def run(active=1, n=N, top=TOP, rho=RHO, on_back_half=True):
    poly, x, y, _ = st.aperture("hex")
    on = [(np.cos(np.radians((m + 0.5) * 360 / n)) < 0) == (active == 1) for m in range(n)]
    r = L.solve(n, rho * 1e-3, top, on)
    psi_g, p = r["psi"], r["p_psi"]
    p_rel = p / np.mean(p)
    beams = {1: (np.radians(25), 0.0), -1: (np.radians(25), np.pi)}
    pols = {1: "RHCP", -1: "LHCP"}
    from src.hologram.leaky_forward import polygon_ray_distance

    jx = np.zeros(x.shape, complex)
    jy = np.zeros(x.shape, complex)
    w_tot = 0.0
    for i in range(x.size):
        rho_i = np.hypot(x[i], y[i])
        if rho_i < st.RHO_KEEP:
            continue
        psi = np.arctan2(y[i], x[i])
        side_active = (x[i] * active) > 0  # this cell belongs to the active half
        beam = beams[1] if x[i] > 0 else beams[-1]
        pol = pols[1] if x[i] > 0 else pols[-1]
        ell = float(polygon_ray_distance(poly, (0.0, 0.0), np.array([psi]))[0])
        # own-half ray length to the boundary (x = 0 is the half boundary; rays from the centre never cross it inside the half)
        a = min(-np.log(1 - tp.ETA_D) / (2 * max(ell - st.RHO_KEEP, 1e-3)), tp.ALPHA_CAP)
        pw = float(np.interp(psi, psi_g, p_rel))
        amp = np.sqrt(2 * a / rho_i * pw) * np.exp(-a * (rho_i - st.RHO_KEEP))
        j = st.jones_to_current(st.JONES[pol], *beam)
        ph = -st.K0 * np.sin(beam[0]) * (x[i] * np.cos(beam[1]) + y[i] * np.sin(beam[1]))
        cur = amp * np.exp(1j * ph)
        jx[i], jy[i] = j[0] * cur, j[1] * cur
        w_tot += 2 * a * pw * np.exp(-2 * a * (rho_i - st.RHO_KEEP)) / rho_i * st.PITCH**2 / (2 * np.pi)
    et, ep = st.far_field(x, y, jx, jy)
    th0, ph0 = beams[active]
    tho, pho = beams[-active]
    mt = st.metrics(et, ep, th0, ph0, pols[active])
    lk = st.leak_db(et, ep, tho, pho, mt["_inten_peak"])
    gain = mt["D_dbi"] + 10 * np.log10((r["p_sw"] / r["p_acc"]) * w_tot)
    fb, _ = L.fb_2q(r) if active == 1 else (None, None)
    return {"gain_dbi": float(gain), "D_dbi": mt["D_dbi"], "SLL_db": mt["SLL_db"], "other_beam_level_db": lk, "isolation_db": -lk, "s11_db": r["s11_db"], "sw_over_acc": r["p_sw"] / r["p_acc"],
            "converted_fraction": w_tot, "fb_launcher_db": fb}


if __name__ == "__main__":
    out = {}
    for act in (1,):
        out["state_A"] = run(act)
        print(act, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in out["state_A"].items()}, flush=True)
    # reference: bare probe (all pins OFF) feeding the same board: both beams on
    out["bare_probe_reference"] = None
    Path("results/sidequest_hybrid/launcher_board_gain.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8", newline="\n")
