"""Writes results/sub6_board_island_hfsscal/HFSSCAL.md from points_hfsscal.csv and calibration.json (MODEL_ONLY; dry run while the ON data are missing)."""
import csv
import json
from pathlib import Path

import numpy as np

D = Path(__file__).resolve().parents[2] / "results" / "sub6_board_island_hfsscal"
R = list(csv.DictReader(open(D / "points_hfsscal.csv")))
cal = json.load(open(D / "calibration.json"))
g = lambda r, k: float(r[k])  # noqa: E731
f = lambda r, k, d=2: f"{g(r, k):.{d}f}"  # noqa: E731
V3 = [r for r in R if r["variant"] == "V3"]
by = {(r["mode"], r["N_sw"], r["placement"], r["taper"], r["excl"]): r for r in V3}
tc = "| direction | HFSS points (phase, f, Q) | MoM f at that phase | shift a (GHz) | beta/k0 at 5.8 GHz: MoM -> calibrated | d phi / d f (deg/GHz) | alpha_x = omega/(2 Q v_g) (Np/m) |\n|---|---|---|---|---|---|---|\n"
for d in ("x", "y"):
    c = cal["OFF"][d]
    pts = "; ".join(f"{p['phase_deg']} deg, {p['f_GHz']} GHz, Q {p['Q']}" for p in c["points"])
    tc += f"| {c['direction']} | {pts} | {', '.join(f'{v:.3f}' for v in c['f_MoM_at_phase_GHz'])} | {c['shift_a_GHz']:+.3f} | {c['beta_over_k0_MoM_5p8']:.4f} -> {c['beta_over_k0_5p8']:.4f} | {c['dphi_df_deg_per_GHz']:.1f} | {c['alpha_x_Np_per_m']:.3f} |\n"
t1 = "| N_sw | mode | MEETS_ALL of 12 | best design (pl / taper / excl) | gate (margin) | OFF G_A / G_B / bal | OFF AR_A / AR_B | ON G_A / B supp / AR_A | SLL off / on | eta off / on | absorbed OFF | switched beta OFF / ON (mean) | alpha_x (Np/m) |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|\n"
for n in ("120", "135", "240"):
    for mode in ("base", "cal"):
        rs = [r for r in V3 if r["mode"] == mode and r["N_sw"] == n]
        b = max(rs, key=lambda r: g(r, "margin"))
        t1 += (f"| {n} | {mode}{' (ON = MoM placeholder)' if mode == 'cal' else ' (MoM, same single row)'} | {sum(r['gate'] == 'MEETS_ALL' for r in rs)} | {b['placement']} / {float(b['taper']):g} / {int(float(b['excl']))} | {b['gate']} ({f(b, 'margin')}) | "
               f"{f(b, 'off_GA', 1)} / {f(b, 'off_GB', 1)} / {f(b, 'bal')} | {f(b, 'off_ARA3')} / {f(b, 'off_ARB3')} | {f(b, 'on_GA', 1)} / {f(b, 'supp', 1)} / {f(b, 'on_ARA3')} | {f(b, 'sll_off', 1)} / {f(b, 'sll_on', 1)} | "
               f"{f(b, 'eta_off', 3)} / {f(b, 'eta_on', 3)} | {f(b, 'absorbed_frac_off', 3)} | {f(b, 'sw_beta_off_mean', 3)} / {f(b, 'sw_beta_on_mean', 3)} | {f(b, 'sw_alpha_x_off_mean')} |\n")
t2 = "| N_sw | paired designs | mean change OFF G_A (dB) | mean change ON G_A (dB) | mean change of the gate margin (dB) | mean eta_off ratio cal/base | designs MEETS_ALL base -> cal |\n|---|---|---|---|---|---|---|\n"
for n in ("120", "135", "240"):
    d = [(g(by[("cal", n, pl, t, ex)], "off_GA") - g(by[("base", n, pl, t, ex)], "off_GA"), g(by[("cal", n, pl, t, ex)], "on_GA") - g(by[("base", n, pl, t, ex)], "on_GA"),
          g(by[("cal", n, pl, t, ex)], "margin") - g(by[("base", n, pl, t, ex)], "margin"), g(by[("cal", n, pl, t, ex)], "eta_off") / g(by[("base", n, pl, t, ex)], "eta_off"))
         for pl in ("a", "b") for t in ("0.5000", "0.8000", "1.0000") for ex in ("0", "1") if ("cal", n, pl, t, ex) in by]
    a = np.array(d).mean(0)
    mb = sum(r["gate"] == "MEETS_ALL" for r in V3 if r["mode"] == "base" and r["N_sw"] == n)
    mc = sum(r["gate"] == "MEETS_ALL" for r in V3 if r["mode"] == "cal" and r["N_sw"] == n)
    t2 += f"| {n} | {len(d)} | {a[0]:+.2f} | {a[1]:+.2f} | {a[2]:+.2f} | {a[3]:.3f} | {mb} -> {mc} |\n"
md = f"""# SUB6 board (island cells) with HFSS-calibrated switched cells: DRY RUN (stage: sub6-board-island-hfsscal)

**MODEL_ONLY / DRY RUN. ON = MoM placeholder (no HFSS ON data yet), OFF calibrated from ONE phase per direction. Do not use for decisions.** Code `scripts/research/sub6_board_hfsscal.py` (imports `sub6_board_island2.py` and `sub6_board_model4.py` unchanged); data `hfss_data.json` (the user's HFSS numbers, quoted from the chat, not a pushed commit), `calibration.json`, `points_hfsscal.csv` (288 rows = 72 designs x V0-V3).

## Method
* Cell: the single island row `ISL_W2.9_q0.8_Li0.5_c0.6` (the row HFSS was run for; `base` and `cal` both use only this row, so they differ only by the calibration; the 9-row board-island numbers are not directly comparable).
* HFSS (ISL_vias2: vias + bias stubs, FR4 0.8, OFF, 85 deg): along 5.963 GHz Q 250, across 5.566 GHz Q 189 (surface-wave modes); not used: via resonance 4.54 / 4.70 GHz (Q about 125), air mode 6.023 GHz.
* Per state and direction the MoM dispersion f_M(phase) is corrected, f_cal = f_M + a (+ b (phi - phi_ref) once two phases exist); beta/k0 at 5.8 GHz from f_cal; v_g from its slope; alpha_x = omega / (2 Q v_g) = f (d beta/d f) / (2 Q). The MoM leakage alpha of the row is kept (MoM-only); alpha_x is absorbed, not radiated.
* On the board: beta(theta') = beta_along cos^2 theta' + beta_across sin^2 theta' with theta' = ray angle - slot axis of the hologram (3 iterations, axes before the 90-degree flips; the flips are not fed back); alpha_x likewise. ON: MoM beta and alpha, alpha_x = 0.
{tc}
## Result (V3, 200 mm, 12 designs per N: a/b x taper 0.5/0.8/1.0 x ring exclusion off/on)
{t1}
Paired change (cal - base):
{t2}
## Reading
The OFF calibration lowers the switched-cell beta (mean 1.104 -> 1.085) and adds 0.58 Np/m of absorbed loss (2.7-4.6 % of the launched power); it changes the gains by <= 0.05 dB and the gate margin by -0.15 ... +0.05 dB on average, but N 120 loses all 4 passing designs (SLL misses by 0.1-0.2 dB), N 135 keeps 3, N 240 goes 11 -> 12 of 12. The ON-state numbers are the MoM ones and the ON calibration (beta shift ON-OFF, ON loss) will matter more: the OFF-ON beta change of the switched cells moves from 0.10 (MoM) to 0.08 now but is untested.

## Limits
One HFSS point per direction (shift only, MoM slope assumed for v_g and for the phase-to-frequency map); Q includes every loss of the HFSS cell (diode R, FR4, copper, the bottom radiation boundary through the anti-pads) and is treated as extra absorption of the switched cells only; the static cells are not calibrated; the theta' weighting (cos^2) is a model choice; the flips are not fed back into beta; HFSS v2 has the 2.0 x 0.3 mm bias stubs. Waiting for: OFF 95 deg, ON 80 / 85 / 95 deg.
"""
(D / "HFSSCAL.md").write_text(md, encoding="utf-8", newline="\n")
print("ok")
