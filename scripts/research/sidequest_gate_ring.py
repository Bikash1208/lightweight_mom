"""Step 4: feed S11 and open-quadrant gain of A-2q with the REAL pass-state (matched by 1 transition row each side) and a closed block gate. MODEL_ONLY.
Same cylindrical-harmonic model as sidequest_gate_reflection.py, but the open half is a thin interface with the complex S11(psi), S21(psi) of the
5-row stack [transition | 3 gate rows OFF | transition] (single-mode line model, x-cells for |psi| < 45 deg, y-cells beyond), the blocked half an
ideal reflector (Gamma = +1 at the gate front, T = 0, step-1 convention; its phase is NOT known from the cell model and is absorbed in the radius).
Interface at the gate front face r_g; the inner transition row is the row just inside r_g, hologram starts r_g + 25.6 mm."""
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
from sidequest_arch_study import LAUNCH, ETA_D

YREF = 1j * 7.82e-3
YXX, YYY = 0.0009e-3 + 2.5e-3j, 6.6e-3j  # OFF state, mesh 64, 10.54 GHz (xx along the diode axis)
BT = 5.64e-3  # transition-row susceptance [S] (optimised for psi_rel 0..45 deg)
PSI_REL = np.linspace(0, 45, 10)


def pass_stack(psi_rel, bt=BT):
    c, s = np.cos(np.radians(psi_rel)), np.sin(np.radians(psi_rel))
    y = YXX * c * c + YYY * s * s
    tr = [(1j * bt, 1)]
    return gm.cascade(YREF, tr + [(y, 3)] + tr)


def tables(bt=BT):
    s = [pass_stack(p, bt) for p in PSI_REL]
    return np.array([a for a, _ in s]), np.array([b for _, b in s])


def psi_rel_of(psi):
    a = np.abs(np.arctan2(np.sin(psi), np.cos(psi)))  # |psi| in [0, pi]
    a = np.where(a > np.pi / 2, np.pi - a, a)  # fold to [0, 90 deg] (mirror about the y axis)
    return np.degrees(np.minimum(a, np.pi / 2 - a))


def interface_functions(r_g, ideal=False, gam_b=1.0, bt=BT):
    br = gc.beta_complex(YREF)
    d_front = 1 * gc.P  # inner transition row
    l_tot = 5 * gc.P
    s11, s21 = tables(bt)
    pr = psi_rel_of(gr.PSI)
    g11 = np.interp(pr, PSI_REL, s11.real) + 1j * np.interp(pr, PSI_REL, s11.imag)
    g21 = np.interp(pr, PSI_REL, np.abs(s21)) * np.exp(1j * (np.interp(pr, PSI_REL, np.unwrap(np.angle(s21))) + br.real * l_tot))
    g11 = g11 * np.exp(2j * br.real * d_front)  # S11 referenced to r_g instead of the stack front plane
    gam = np.where(gr.OPEN > 0, 0.0 if ideal else g11, gam_b)
    tr = np.where(gr.OPEN > 0, 1.0 if ideal else g21, 0.0)
    return gam.astype(complex), tr.astype(complex)


def solve2(r, gam_psi, t_psi):
    h1, h2 = gr.hankel1(gr.NS, gr.BETA * r), gr.hankel2(gr.NS, gr.BETA * r)
    Tg = gr.TOEP(gr.FOUR(gam_psi))
    Tt = gr.TOEP(gr.FOUR(t_psi))
    e0 = (gr.NS == 0).astype(float)
    amat = np.eye(len(e0)) - np.outer(e0, e0)  # matched probe: A_0 = 1 (s_p = 0), A_n = B_n
    lhs = np.diag(h1) - Tg @ np.diag(h2) @ amat
    b = np.linalg.solve(lhs, Tg @ (h2 * e0))
    a = amat @ b + e0
    c = np.linalg.solve(np.diag(h2), Tt @ (h2 * a))
    return LAUNCH * b[gr.N], c


def run(r_g, ideal=False, gain=True, gam_b=1.0, bt=BT):
    gam, tr = interface_functions(r_g, ideal, gam_b, bt)
    s11, c = solve2(r_g, gam, tr)
    out = {"r_g_mm": r_g * 1e3, "S11_db": float(20 * np.log10(abs(s11) + 1e-12)), "p_open": gr.p_open(c)}
    if gain:
        mt = gr.pattern_gain(c, r_g + gc.P, rows=3)  # hologram starts at r_g + 6.4 + 19.2 mm
        out.update({"D_dbi": mt["D_dbi"], "SLL_db": mt["SLL_db"], "gain_dbi": mt["D_dbi"] + 10 * np.log10(LAUNCH * out["p_open"] * ETA_D)})
    return out


if __name__ == "__main__":
    s11, s21 = tables()
    print("pass stack RL/IL vs psi_rel:", [(round(float(p)), round(float(-20 * np.log10(abs(a))), 1), round(float(-20 * np.log10(abs(b))), 2)) for p, a, b in zip(PSI_REL, s11, s21)])
    res = {"bt_ms": BT * 1e3, "pass_RL_db": [float(-20 * np.log10(abs(a))) for a in s11], "pass_IL_db": [float(-20 * np.log10(abs(b))) for b in s21], "sweep": []}
    lam_g = 2 * np.pi / gr.BETA
    r0 = 0.01505
    for dr in np.linspace(-lam_g / 4, lam_g / 4, 17):
        r = r0 + dr
        if r < 0.0065:
            continue
        real, ideal = run(r), run(r, ideal=True)
        res["sweep"].append({"dr_mm": float(dr * 1e3), "real": real, "ideal": ideal})
        print(round(dr * 1e3, 1), "real", {k: round(v, 2) for k, v in real.items()}, "| ideal S11 %.1f gain %.1f" % (ideal["S11_db"], ideal["gain_dbi"]), flush=True)
    res["centre"] = {"real": run(r0), "ideal": run(r0, ideal=True)}
    Path("results/sidequest_hybrid/gate_ring_real_pass.json").write_text(json.dumps(res, indent=1), encoding="utf-8", newline="\n")
    print("centre", res["centre"])


# ------------------------------------------------------------------ REAL block state (all modal roots of the ON sheet)
YXX_ON = (100.1 + 68.0j) * 1e-3  # retuned cell, ON, mesh 64, 10.54 GHz


def on_stack(psi_rel, bt=BT):
    """S11, S21 (complex) of [transition | 3 ON gate rows | transition]; the gate section uses the root with the LARGEST transmission."""
    c, s = np.cos(np.radians(psi_rel)), np.sin(np.radians(psi_rel))
    y = YXX_ON * c * c + YYY * s * s
    best = None
    for bg in gc.roots_all(y):
        zg = gc.zm(bg)
        th = bg * 3 * gc.P
        mg = np.array([[np.cos(th), 1j * zg * np.sin(th)], [1j * np.sin(th) / zg, np.cos(th)]])
        mt, _ = gm.abcd(1j * bt, 1)
        m = mt @ mg @ mt
        zr = gc.zm(gc.beta_complex(YREF))
        a, b, cc, d = m[0, 0], m[0, 1], m[1, 0], m[1, 1]
        den = a + b / zr + cc * zr + d
        s21, s11 = 2 / den, (a + b / zr - cc * zr - d) / den
        if best is None or abs(s21) > abs(best[1]):
            best = (s11, s21)
    return best


def run_real_block(r_g):
    br = gc.beta_complex(YREF)
    tab = [on_stack(p) for p in PSI_REL]
    s11 = np.array([a for a, _ in tab])
    s21 = np.array([b for _, b in tab])
    pr = psi_rel_of(gr.PSI)
    g11 = (np.interp(pr, PSI_REL, s11.real) + 1j * np.interp(pr, PSI_REL, s11.imag)) * np.exp(2j * br.real * gc.P)
    g21 = np.interp(pr, PSI_REL, np.abs(s21)) * np.exp(1j * (np.interp(pr, PSI_REL, np.unwrap(np.angle(s21))) + br.real * 5 * gc.P))
    gam_p, tr_p = interface_functions(r_g, False)
    gam = np.where(gr.OPEN > 0, gam_p, g11)
    tr = np.where(gr.OPEN > 0, tr_p, g21)
    s11f, c = solve2(r_g, gam, tr)
    f = sum(cn * (1j) ** int(n) * np.exp(1j * n * gr.PSI) for n, cn in zip(gr.NS, c, strict=True))
    p_open = float(np.mean(gr.OPEN * np.abs(f) ** 2))
    p_blk = float(np.mean((1 - gr.OPEN) * np.abs(f) ** 2))
    return {"S11_db": float(20 * np.log10(abs(s11f))), "p_open": p_open, "p_blocked_side": p_blk, "isolation_db": float(10 * np.log10(p_open / max(p_blk, 1e-12))),
            "on_stack_S21_db": [float(20 * np.log10(abs(x))) for x in s21], "on_stack_S11_db": [float(20 * np.log10(abs(x))) for x in s11]}
