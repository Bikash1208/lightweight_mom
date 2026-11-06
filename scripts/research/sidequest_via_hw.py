"""Via-switched GB cell with the REAL bias hardware (RT5880 1.575, 3 Cu layers): cascade model.  MODEL_ONLY.
Cascade from the patch to RF ground (all per square cell):
  patch --L_via--> node V (C_V) --[C_s chip + ESL]--> node A (C_A, choke Y_ch) --diode (ON R / OFF C_j)--> node B (cathode pad: C_B to L2, RF ground).
Hardware: the anti-pad is a slot in L2 under via land + C_s + island A (their capacitance to L2 is only fringing, ~15 fF each); the cathode pad B (RF ground, about 2.3 pF to L2 over the 0.1 mm bond ply);
choke = microstrip lambda/4-class line (0.15 mm, over the bond ply) from node A to the bypass pad C (RF ground, about 2.2 pF) - the string bias path.   python sidequest_via_hw.py"""
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from scipy.optimize import brentq

import sidequest_cellB_model as M
import sidequest_via_cell as V

RES = ROOT / "results/sidequest_hybrid"
MU0, EPS0, C0 = M.MU0, M.EPS0, M.C0
H = V.H
PLY = dict(er=3.52, tand=0.004, d=0.101e-3)  # Rogers RO4450F (datasheet mirrors; bonding to PTFE RT5880 to be confirmed)
HW = dict(L_via_extra=0.05e-9, C_V=15e-15, C_A=15e-15, L_cs=0.25e-9, R_cs=0.2, C_B=2.3e-12, L_B=0.05e-9, C_C=2.2e-12, L_C=0.05e-9, w_ch=0.15e-3, t_cu=35e-6)
DIODE = dict(R=2.44, C=28.16e-15)


def microstrip(w, h, er, t=35e-6):
    u = w / h
    ee = (er + 1) / 2 + (er - 1) / 2 * (1 + 12 / u) ** -0.5
    z0 = 60 / math.sqrt(ee) * math.log(8 / u + u / 4) if u <= 1 else 120 * math.pi / (math.sqrt(ee) * (u + 1.393 + 0.667 * math.log(u + 1.444)))
    return z0, ee


def choke_Zin(f, length, hw=HW, ply=PLY, w=None, d=None, er=None):
    w = hw["w_ch"] if w is None else w
    d = ply["d"] if d is None else d
    er = ply["er"] if er is None else er
    z0, ee = microstrip(w, d, er)
    om = 2 * math.pi * f
    k = om / C0 * math.sqrt(ee) * (1 - 0.5j * ply["tand"] * (er / ee) * (ee - 1) / (er - 1))  # complex propagation constant (dielectric loss)
    zl = 1.0 / (1j * om * hw["C_C"]) + 1j * om * hw["L_C"]
    t_ = np.tan(k * length)
    return z0 * (zl + 1j * z0 * t_) / (z0 + 1j * zl * t_)


def branch(state, f, C_s, length, hw=HW, kappa=1.0, dC=1.0, dR=1.0, dL=1.0, dCs=1.0, dPad=1.0, dChoke=None, ply=PLY):
    om = 2 * math.pi * f
    zB = 1.0 / (1j * om * hw["C_B"] * dPad) + 1j * om * hw["L_B"]
    zd = DIODE["R"] * dR if state == "ON" else 1.0 / (1j * om * DIODE["C"] * dC)
    z_ch = choke_Zin(f, length * (1.0 if dChoke is None else dChoke[0]), hw, ply, w=None if dChoke is None else hw["w_ch"] * dChoke[1], d=None if dChoke is None else ply["d"] * dChoke[2], er=None if dChoke is None else ply["er"] * dChoke[3])
    yA = 1j * om * hw["C_A"] * dPad + 1.0 / z_ch
    zA = 1.0 / (yA + 1.0 / (zd + zB))
    zcs = hw["R_cs"] + 1j * om * hw["L_cs"] - 1j / (om * C_s * dCs)
    yV = 1j * om * hw["C_V"] * dPad
    zV = 1.0 / (yV + 1.0 / (zcs + zA))
    zin = 1j * om * (kappa * MU0 * H * dL + hw["L_via_extra"]) + zV
    return 1.0 / zin


def patch():
    rows = json.loads((RES / "cellB_rt5880_a.json").read_text())
    return M.Patch([r for r in rows if r["lam"] == "RT5880" and r["mesh"] == 64 and r["psi"] == 0.0 and abs(r["gap"] - 0.4) < 1e-6])


def beta_mean(P, C_s, length, **kw):
    yon, yoff = branch("ON", V.F0, C_s, length, **kw), branch("OFF", V.F0, C_s, length, **kw)
    r = V.root(P, 0.5 * (yon + yoff))
    return None if r is None else r[0]


def knob(P, length, target=1.062, **kw):
    f = lambda c: (lambda b: 1.0 if b is None else b)(beta_mean(P, c * 1e-12, length, **kw)) - target
    cs = np.geomspace(0.15, 6.0, 120)
    v = np.array([f(c) for c in cs])
    ix = np.where(np.diff(np.sign(v)) != 0)[0]
    return [brentq(f, cs[i], cs[i + 1], xtol=1e-7) * 1e-12 for i in ix]


def margins(C_s, length, window=1.5, **kw):
    """series / parallel resonances of the branch (|Y| extrema) within +-window GHz of 10.5 GHz, for ON and OFF."""
    out = {}
    fs = np.linspace(4e9, 17e9, 1300)
    for st in ("ON", "OFF"):
        y = np.array([branch(st, f, C_s, length, **kw) for f in fs])
        a = np.abs(y)
        pk = [fs[i] / 1e9 for i in range(1, len(fs) - 1) if (a[i] > a[i - 1] and a[i] > a[i + 1]) or (a[i] < a[i - 1] and a[i] < a[i + 1])]
        out[st] = pk
    return out


def main():
    P = patch()
    # (0) the old ideal-branch design (no hardware): C_s 0.614 pF, step 15.0 mS
    import sidequest_via_design as D
    yon0, yoff0 = D.br("ON", "MA4AGP907 (measured)", 0.614e-12, 0.0), D.br("OFF", "MA4AGP907 (measured)", 0.614e-12, 0.0)
    print("ideal branch (no hardware): Y_ON %.2f%+.2fj Y_OFF %.2f mS step %.1f beta_mean %.3f" % (yon0.real * 1e3, yon0.imag * 1e3, yoff0.imag * 1e3, (yoff0 - yon0).imag * 1e3, 1.062))
    # (1) hardware added, C_s unchanged
    for length_mm in (3.0, 4.0, 4.9, 6.0):
        L = length_mm * 1e-3
        yon, yoff = branch("ON", V.F0, 0.614e-12, L), branch("OFF", V.F0, 0.614e-12, L)
        b = beta_mean(P, 0.614e-12, L)
        print("hardware, C_s 0.614 unchanged, choke %.1f mm: Y_ON %.2f%+.2fj Y_OFF %.2f%+.2fj step %.1f beta_mean %s Y_choke %s" % (length_mm, yon.real * 1e3, yon.imag * 1e3, yoff.real * 1e3, yoff.imag * 1e3, (yoff - yon).imag * 1e3, None if b is None else round(b, 3),
              np.round(1 / choke_Zin(V.F0, L) * 1e3, 2)))
    # (2) re-tuned C_s per choke length
    for length_mm in (3.0, 3.5, 4.0, 4.5, 4.9, 5.4, 6.0):
        L = length_mm * 1e-3
        for C in knob(P, L):
            yon, yoff = branch("ON", V.F0, C, L), branch("OFF", V.F0, C, L)
            mg = margins(C, L)
            near = {k: [round(x, 2) for x in v if abs(x - 10.5) < 1.5] for k, v in mg.items()}
            print("choke %.1f mm retuned C_s %.3f pF: Y_ON %.2f%+.2fj Y_OFF %.2f%+.2fj step %.1f | extrema within 1.5 GHz: %s" % (length_mm, C * 1e12, yon.real * 1e3, yon.imag * 1e3, yoff.real * 1e3, yoff.imag * 1e3, (yoff - yon).imag * 1e3, near))


if __name__ == "__main__" and (len(sys.argv) < 2 or sys.argv[1] == "tl"):
    main()


# ---------------------------------------------------------------- resistive bias feed (recommended): ballast resistor from node A to the bias bus, cathode to ground by a blind via
RB = dict(R_b=680.0, C_p=25e-15, Z_bus=45.0, L_bv=0.15e-9, C_A=35e-15, C_V=15e-15, L_cs=0.25e-9, R_cs=0.2, L_via_extra=0.05e-9)


def branch2(state, f, C_s, rb=RB, kappa=1.0, dC=1.0, dR=1.0, dL=1.0, dCs=1.0, dRb=1.0, dCp=1.0, dCA=1.0, dCV=1.0):
    om = 2 * math.pi * f
    zd = DIODE["R"] * dR if state == "ON" else 1.0 / (1j * om * DIODE["C"] * dC)
    zB = (1j * om * rb["L_bv"] + 0.05) if rb.get("cathode", "via") == "via" else (1.0 / (1j * om * rb.get("C_B", 2.3e-12) * dCV) + 1j * om * 0.05e-9)
    zR = 1.0 / (1.0 / (rb["R_b"] * dRb) + 1j * om * rb["C_p"] * dCp)
    y_ch = 1.0 / (zR + rb["Z_bus"])
    yA = 1j * om * rb["C_A"] * dCA + y_ch
    zA = 1.0 / (yA + 1.0 / (zd + zB))
    zcs = rb["R_cs"] + 1j * om * rb["L_cs"] - 1j / (om * C_s * dCs)
    yV = 1j * om * rb["C_V"] * dCV
    zV = 1.0 / (yV + 1.0 / (zcs + zA))
    return 1.0 / (1j * om * (kappa * MU0 * H * dL + rb["L_via_extra"]) + zV)


def beta_mean2(P, C_s, **kw):
    yon, yoff = branch2("ON", V.F0, C_s, **kw), branch2("OFF", V.F0, C_s, **kw)
    r = V.root(P, 0.5 * (yon + yoff))
    return None if r is None else r[0]


def knob2(P, target=1.062, **kw):
    f = lambda c: (lambda b: 1.0 if b is None else b)(beta_mean2(P, c * 1e-12, **kw)) - target
    cs = np.geomspace(0.15, 8.0, 160)
    v = np.array([f(c) for c in cs])
    ix = np.where(np.diff(np.sign(v)) != 0)[0]
    return [brentq(f, cs[i], cs[i + 1], xtol=1e-7) * 1e-12 for i in ix]


def extrema(C_s, st, rb=RB, f0=4e9, f1=17e9, **kw):
    fs = np.linspace(f0, f1, 1300)
    a = np.array([abs(branch2(st, f, C_s, rb, **kw)) for f in fs])
    return [round(fs[i] / 1e9, 2) for i in range(1, len(fs) - 1) if (a[i] > a[i - 1] and a[i] > a[i + 1]) or (a[i] < a[i - 1] and a[i] < a[i + 1])]


def main2():
    P = patch()
    print("== resistive bias feed ==")
    yon0, yoff0 = branch2("ON", V.F0, 0.614e-12, dict(RB, R_b=1e12, C_p=0.0, C_A=15e-15, L_bv=0.0)), branch2("OFF", V.F0, 0.614e-12, dict(RB, R_b=1e12, C_p=0.0, C_A=15e-15, L_bv=0.0))
    print("pads / anti-pad / blind via only (no feed): Y_ON %.2f%+.2fj Y_OFF %.2f%+.2fj step %.1f" % (yon0.real * 1e3, yon0.imag * 1e3, yoff0.real * 1e3, yoff0.imag * 1e3, (yoff0 - yon0).imag * 1e3))
    out = {}
    for Rb in (330.0, 680.0, 1000.0, 2200.0, 4700.0):
        rb = dict(RB, R_b=Rb)
        yon, yoff = branch2("ON", V.F0, 0.614e-12, rb), branch2("OFF", V.F0, 0.614e-12, rb)
        b0 = beta_mean2(P, 0.614e-12, rb=rb)
        cs = knob2(P, rb=rb)
        line = "R_b %5.0f: raw (C_s 0.614): Y_ON %.2f%+.2fj Y_OFF %.2f%+.2fj step %.1f beta_mean %s" % (Rb, yon.real * 1e3, yon.imag * 1e3, yoff.real * 1e3, yoff.imag * 1e3, (yoff - yon).imag * 1e3, None if b0 is None else round(b0, 3))
        for C in cs[:1]:
            yon, yoff = branch2("ON", V.F0, C, rb), branch2("OFF", V.F0, C, rb)
            ex = {st: [x for x in extrema(C, st, rb) if abs(x - 10.5) < 1.5] for st in ("ON", "OFF")}
            line += " | retuned C_s %.3f pF: Y_ON %.2f%+.2fj Y_OFF %.2f%+.2fj step %.1f extrema<1.5GHz %s" % (C * 1e12, yon.real * 1e3, yon.imag * 1e3, yoff.real * 1e3, yoff.imag * 1e3, (yoff - yon).imag * 1e3, ex)
            out[Rb] = dict(C_s_pF=C * 1e12, step=(yoff - yon).imag * 1e3, Y_ON=[yon.real * 1e3, yon.imag * 1e3], Y_OFF=[yoff.real * 1e3, yoff.imag * 1e3])
        print(line)
    return out


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "rb":
    main2()
