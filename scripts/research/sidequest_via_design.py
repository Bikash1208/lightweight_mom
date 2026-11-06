"""Via-switched GB cell on RT5880 1.575 with packaged PIN diodes: knob search (C_s, L_ser, patch gap) for a 15-20 mS step, mean sheet = tensor mean (beta/k0 1.062), both branch resonances >= 1.5 GHz from 10.5 GHz, minimum tolerance sensitivity.  MODEL_ONLY.
Diode lumped models (R_on, C_j(OFF), L_package) are from distributor / datasheet snippets, NOT vendor S-parameters (not retrievable here): see DIODES.   python sidequest_via_design.py"""
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_cellB_model as M
import sidequest_via_cell as V

RES = ROOT / "results/sidequest_hybrid"
W = V.W0
MU0 = M.MU0
H = V.H
# name: (R_on ohm, C_j OFF F, L_package H, source note)
DIODES = {
    "MA4AGP907 (measured)": (2.44, 28.16e-15, 0.0, "user measured model (package-free chip)"),
    "SMP1345-079LF": (1.5, 0.19e-12, 0.7e-9, "repo src/ris/diodes.py from Skyworks SMP1345 Series datasheet 200046U (2018-06-19): R_s 1.5 ohm typ (10 mA, 100 MHz), C_T 0.19 pF typ (1 V, 1 MHz), L_s 0.7 nH (SC-79)"),
    "MADP-000907-14020 (chip)": (5.2, 0.025e-12, 0.0, "repo src/ris/diodes.py from MACOM MADP-000907-14020 datasheet Rev V8: R_s 5.2 ohm typ (10 mA, 1 GHz), C_T 0.025 pF typ (-10 V, 1 MHz), L_s NOT_AVAILABLE (0 used)"),
    "BAR63-02L": (1.0, 0.25e-12, 0.6e-9, "Infineon BAR63-02L: Ls 0.6 nH (preliminary datasheet), rf <= 1 ohm and CT 0.21-0.3 pF from third-party listings (UNVERIFIED)"),
    "BAR90-02LRH": (1.3, 0.19e-12, 0.4e-9, "Infineon BAR90 family: Ls 0.4 nH, rf 1.3 ohm at 3 mA (sibling), CT 0.19 pF (BAR90-02LS) (UNVERIFIED for the -02LRH)"),
    "BAR50-02V": (4.5, 0.30e-12, 0.7e-9, "Infineon BAR50-02V: rf 4.5 ohm at 10 mA 100 MHz, CT 0.24-0.4 pF (listings disagree; mid value used), Ls NOT found (SC-79 value of the SMP1345 assumed)"),
}


def br(state, d, C_s, L_ser, kappa=1.0, C_pad=25e-15, dR=1.0, dC=1.0, dL=1.0, dCs=1.0, dLs=1.0, f=V.F0):
    R, Cj, Lp, _ = DIODES[d]
    w = 2 * math.pi * f
    L = kappa * MU0 * H * dL + Lp + L_ser * dLs
    x = w * L - 1.0 / (w * C_s * dCs)
    if state == "ON":
        return 1.0 / (R * dR + 1j * x)
    return 1.0 / (1j * (x - 1.0 / (w * (Cj * dC + C_pad))))


def resonances(d, C_s, L_ser, C_pad=25e-15):
    R, Cj, Lp, _ = DIODES[d]
    L = MU0 * H + Lp + L_ser
    f_on = 1.0 / (2 * math.pi * math.sqrt(L * C_s))
    ceq = 1.0 / (1.0 / C_s + 1.0 / (Cj + C_pad))
    f_off = 1.0 / (2 * math.pi * math.sqrt(L * ceq))
    return f_on / 1e9, f_off / 1e9


def patch_table():
    rows = json.loads((RES / "cellB_rt5880_explore.json").read_text())
    G = {}
    for r in rows:
        G.setdefault(round(r["gap"], 3), []).append(r)
    gaps = sorted(G)
    yp = {g: M.Patch(G[g]).y_sw(1.062).imag for g in gaps}
    return gaps, yp


def search(M_target=7.55e-3, step_lo=14.0, step_hi=21.0, margin=1.5, log=print):
    gaps, yp = patch_table()
    out = {}
    cs_grid = np.geomspace(0.05e-12, 8e-12, 90)
    ls_grid = np.linspace(0.0, 3.0e-9, 31)
    dyp = {g: (yp[gaps[min(i + 1, len(gaps) - 1)]] - yp[gaps[max(i - 1, 0)]]) / (gaps[min(i + 1, len(gaps) - 1)] - gaps[max(i - 1, 0)]) for i, g in enumerate(gaps)}
    for d in DIODES:
        best = None
        n_feas = 0
        for g in gaps[:4]:
            for C_s in cs_grid:
                for L_ser in ls_grid:
                    yon, yoff = br("ON", d, C_s, L_ser), br("OFF", d, C_s, L_ser)
                    step_s = (yoff - yon).imag * 1e3
                    step = abs(step_s)
                    mean = yp[g] + 0.5 * (yon + yoff).imag
                    if abs(mean - M_target) > 0.25e-3 or not (step_lo <= step <= step_hi):
                        continue
                    fo, ff = resonances(d, C_s, L_ser)
                    if abs(fo - 10.5) < margin or abs(ff - 10.5) < margin:
                        continue
                    n_feas += 1
                    # sensitivities: change of mean sheet and step [mS]
                    dm, ds = [], []
                    for kw in (dict(dC=0.8), dict(dC=1.2), dict(dR=0.8), dict(dR=1.2), dict(dL=0.92), dict(dL=1.08), dict(dCs=0.9), dict(dCs=1.1), dict(dLs=0.85), dict(dLs=1.15), dict(C_pad=20e-15), dict(C_pad=30e-15)):
                        a, b = br("ON", d, C_s, L_ser, **kw), br("OFF", d, C_s, L_ser, **kw)
                        dm.append(abs(0.5 * ((a + b).imag - (yon + yoff).imag)) * 1e3)
                        ds.append(abs(abs((b - a).imag) - abs((yoff - yon).imag)) * 1e3)
                    for sg in (-0.06, 0.06):
                        dm.append(abs(dyp[g] * sg) * 1e3)
                        ds.append(0.0)
                    rss_m, rss_s = math.sqrt(sum(x * x for x in dm)), math.sqrt(sum(x * x for x in ds))
                    merit = rss_m / 1.7 + rss_s / (0.12 * step)  # 1.7 mS ~ 0.02 in beta/k0; 12 % step ~ 25 % in alpha
                    if best is None or merit < best["merit"]:
                        best = dict(diode=d, polarity=('A: ON lowers Y' if step_s > 0 else 'B: ON raises Y'), gap=g, C_s_pF=C_s * 1e12, L_ser_nH=L_ser * 1e9, step_mS=step, mean_mS=mean * 1e3, Y_on_mS=[yon.real * 1e3, yon.imag * 1e3], Y_off_mS=[yoff.real * 1e3, yoff.imag * 1e3], f_res_on_GHz=fo, f_res_off_GHz=ff,
                                    rss_dmean_mS=rss_m, rss_dstep_mS=rss_s, merit=merit, worst_dmean_mS=max(dm), worst_dstep_mS=max(ds))
        out[d] = {"n_feasible": n_feas, "best": best}
        log(d, "feasible:", n_feas, None if best is None else {k: (round(v, 3) if isinstance(v, float) else v) for k, v in best.items()})
    return out


if __name__ == "__main__":
    res = search()
    (RES / "via_design_search.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline="\n")
