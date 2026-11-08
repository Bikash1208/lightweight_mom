"""Writes results/sub6_board_island/followup/FOLLOWUP.md from budget_best4.json and followup_points.csv (MODEL_ONLY)."""
import csv
import json
from pathlib import Path

D = Path(__file__).resolve().parents[2] / "results" / "sub6_board_island" / "followup"
bud = json.load(open(D / "budget_best4.json"))
R = list(csv.DictReader(open(D / "followup_points.csv")))
g = lambda r, k: float(r[k])  # noqa: E731
t1 = "| point (side N pl taper excl) | G_A OFF / ON | dG (ON-OFF) = dD + d10log(eta) | eta OFF / ON | A cone | B cone | other radiated | unconverted | est. sheet absorption (switched part) | G_A ON: beta shift off / on | switched dbeta, static dphi rms |\n|---|---|---|---|---|---|---|---|---|---|---|\n"
for b in bud:
    p, o, n = b["point"], b["off"], b["on"]
    t1 += (f"| {int(p['side'] * 1e3)} {p['N']} {p['pl']} {p['taper']} {p['excl']} | {o['G_A']:.2f} / {n['G_A']:.2f} | {b['delta_G_ON_minus_OFF']:+.2f} = {b['delta_D']:+.2f} {b['delta_eta_db']:+.2f} | {o['eta']:.3f} / {n['eta']:.3f} | "
           f"{o['frac_A_cone']:.3f} / {n['frac_A_cone']:.3f} | {o['frac_B_cone']:.3f} / {n['frac_B_cone']:.4f} | {o['frac_other_radiated']:.3f} / {n['frac_other_radiated']:.3f} | {o['unconverted']:.3f} / {n['unconverted']:.3f} | "
           f"{o['absorbed_est']:.4f} ({o['absorbed_switched_est']:.4f}) / {n['absorbed_est']:.4f} ({n['absorbed_switched_est']:.4f}) | {b['G_A_ON_no_beta_shift']:.2f} / {b['G_A_ON_with_beta_shift']:.2f} | {b['sw_beta_shift_mean']:+.3f}, {b['dphi_rms_deg']:.1f} deg |\n")
mean = lambda rs, k: sum(g(r, k) for r in rs) / len(rs)  # noqa: E731
t2 = "| variant | points | MEETS_ALL (-3 dB AR) | mean OFF AR_A / AR_B (-3 dB beam) | mean ON AR_A | mean min(G_A) | mean B supp ON | mean SLL off / on | mean rotated |\n|---|---|---|---|---|---|---|---|---|\n"
names = {"V0": "V0 round 'island' (w = 0, single-direction rotations)", "V1": "V1 w = 0.5 (averaged beta)", "V2": "V2 w* (best min G_A)", "V3": "V3 = best of V2 and multi-direction rotations"}
for v in ("V0", "V1", "V2", "V3"):
    rs = [r for r in R if r["variant"] == v]
    t2 += (f"| {names[v]} | {len(rs)} | {sum(r['gate'] == 'MEETS_ALL' for r in rs)} | {mean(rs, 'off_ARA3'):.2f} / {mean(rs, 'off_ARB3'):.2f} | {mean(rs, 'on_ARA3'):.2f} | {mean(rs, 'minGA'):.2f} | {mean(rs, 'supp'):.1f} | "
           f"{mean(rs, 'sll_off'):.1f} / {mean(rs, 'sll_on'):.1f} | {mean(rs, 'rotated'):.0f} |\n")
t3 = "| side | N_sw | best V3 config (pl / taper / excl / w) | OFF G_A / G_B / bal | OFF AR_A / AR_B (-3 dB) | ON G_A / B supp / AR_A | SLL off / on | eta off / on | rotated | V3 gate (margin) | V0 best margin | MEETS_ALL V0 -> V3 (of 12) |\n|---|---|---|---|---|---|---|---|---|---|---|---|\n"
for key in sorted({(r["side_mm"], r["N_sw"]) for r in R}, key=lambda k: (-float(k[0]), int(k[1]))):
    rs3 = [r for r in R if r["variant"] == "V3" and (r["side_mm"], r["N_sw"]) == key]
    rs0 = [r for r in R if r["variant"] == "V0" and (r["side_mm"], r["N_sw"]) == key]
    b = max(rs3, key=lambda r: g(r, "margin"))
    v0 = max(rs0, key=lambda r: g(r, "margin"))
    t3 += (f"| {int(float(key[0]))} | {key[1]} | {b['placement']} / {float(b['taper']):g} / {'yes' if int(float(b['excl'])) else 'no'} / {float(b['w']):g} | {g(b, 'off_GA'):.1f} / {g(b, 'off_GB'):.1f} / {g(b, 'bal'):.2f} | {g(b, 'off_ARA3'):.2f} / {g(b, 'off_ARB3'):.2f} | "
           f"{g(b, 'on_GA'):.1f} / {g(b, 'supp'):.1f} / {g(b, 'on_ARA3'):.2f} | {g(b, 'sll_off'):.1f} / {g(b, 'sll_on'):.1f} | {g(b, 'eta_off'):.2f} / {g(b, 'eta_on'):.2f} | {b['rotated']} | {b['gate']} ({g(b, 'margin'):.2f}) | {g(v0, 'margin'):.2f} | "
           f"{sum(r['gate'] == 'MEETS_ALL' for r in rs0)} -> {sum(r['gate'] == 'MEETS_ALL' for r in rs3)} |\n")
v2 = [r for r in R if r["variant"] == "V2"]
v0 = [r for r in R if r["variant"] == "V0"]
dm = [g(a, "minGA") - g(b, "minGA") for a, b in zip(v2, v0)]
md = f"""# SUB6 board-island follow-up (stage: sub6-board-island, follow-up)

**MODEL_ONLY / PENDING_HFSS_CELL_CHECK / NOT_HFSS_VALIDATED.** Code `scripts/research/sub6_board_island2.py` (+ `sub6_board_model4.py` unchanged); data `budget_best4.json`, `followup_points.csv` (72 designs x 4 variants), this file.

## 1. Why ON G_A < OFF G_A (best 200 / 220 mm points, N 135 and N 240)
{t1}
G = D_A + 10 log10(launch x eta): the ON deficit is the **conversion efficiency**, not the beta shift. In ON the island cells stop radiating, the static cells keep their (alpha-fixed) amplitudes, so eta falls from 0.16-0.39 to 0.06-0.24 and 76-94 % of the accepted power reaches the aperture edge unconverted; the directivity rises by +2.7 to +3.2 dB (one beam instead of two), which only partly compensates. At N 240 the two terms nearly cancel (dG +0.6 / -0.1 dB). **The switched-cell beta shift (-0.07 to -0.09; static-cell phase error rms 10-21 deg) costs only 0.03-0.12 dB of G_A,ON** (columns "beta shift off / on"). Radiated shares come from the far-field power split (cone <= 15 deg around each beam) times the bookkeeping eta, so they are model-level (eta from the alpha bookkeeping, shares from the pattern integral). Residual B in ON: 0.0006-0.0033 of the accepted power (0.2-5 % of the radiated power). Absorption is an ESTIMATE from Re(Y_mean) of each row (alpha_abs = Re Y / (d Im Y_tot/d beta), integrated along the rays with the leakage-only ray power); it is 1-3.6 % of the accepted power (ON diodes: about 1.4-3.2 %), is not subtracted from the gains, and the diode share in ON is about 2x the OFF value.

## 2. Beam-A design phase (averaged / optimised beta)
Since the beta shift is not the cause, the gain change is small: designing the static axes with the phase P_d = P_OFF + w (P_ON - P_OFF) and choosing w* from {{0, .25, .5, .75, 1}} to maximise min(G_A,OFF, G_A,ON) changes min(G_A) by {sum(dm) / len(dm):+.2f} dB on average (max {max(dm):+.2f} dB over the 72 designs); w = 0.5 alone gives {mean([r for r in R if r['variant'] == 'V1'], 'minGA') - mean(v0, 'minGA'):+.2f} dB and loses MEETS_ALL points (SLL). Not worth adopting. The real lever for G_A,ON is eta_on (more static conversion in ON, or fewer unused island cells), not the phase.

## 3. AR in OFF (rotation optimisation at several directions)
Rotation optimisation of ALL sites (90-degree flips, round-2 method, same costs) with the cross-pol constraint (LHCP/RHCP <= -20 dB) imposed at seven directions per beam (centre + 6 on a 5-degree ring) for A (OFF, ON) and B (OFF); candidates w in {{0, w*}}, ring weight in {{0.1, 0.3, 1}}, plus V2 itself; the design with the best gate margin is kept (V3).
{t2}
{t3}
The round-island rotation fixed the cross-pol only at the beam centre, leaving the -3 dB-beam AR of beam B at 2.2-3.1 dB; the multi-direction constraint brings OFF AR_A / AR_B to 1.7-2.5 dB at the cost of SLL (the optimiser does not see SLL, so the selection by gate margin keeps only designs that still meet -13 dB) and of 2x more rotated slots (mean {mean([r for r in R if r['variant'] == 'V3'], 'rotated'):.0f}, max {max(g(r, 'rotated') for r in R if r['variant'] == 'V3'):.0f}; the cap is 200).

## Limits
MODEL_ONLY; mesh-64 island cells (switched table 9 rows); the selection of the best of several optimiser runs by the gate margin is a design search on the evaluation metric; rotated slots near each other with no coupling model; ON residual-B cancellation relies on the mean rho_d (-0.33); G_A numbers are referred to accepted power (launch 0.656 with the standard SMA pin of sub6-launcher). Stopped here.
"""
(D / "FOLLOWUP.md").write_text(md, encoding="utf-8", newline="\n")
print("ok")
