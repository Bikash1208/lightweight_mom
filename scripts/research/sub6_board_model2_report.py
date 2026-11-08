"""Writes results/sub6_board_model/round2/BOARD_MODEL_R2.md from points_r2.csv, tolerance_C.csv, provenance_run2.json (MODEL_ONLY)."""
import csv
import json
from collections import defaultdict
from pathlib import Path

D = Path(__file__).resolve().parents[2] / "results" / "sub6_board_model" / "round2"
R = list(csv.DictReader(open(D / "points_r2.csv")))
TOL = list(csv.DictReader(open(D / "tolerance_C.csv")))
prov = json.load(open(D / "provenance_run2.json"))
tb = json.load(open(D / "tolerance_best_point.json"))
g = lambda r, k: float(r[k])  # noqa: E731
f = lambda r, k, d=2: f"{float(r[k]):.{d}f}"  # noqa: E731
n = len(R)
m_peak = sum(r["gate"] == "MEETS_ALL" for r in R)
m_3 = sum(r["gate_3dB_beam_AR"] == "MEETS_ALL" for r in R)
groups = defaultdict(list)
for r in R:
    groups[(int(g(r, "side_mm")), int(g(r, "N_sw")))].append(r)
t1 = "| side | N_sw (f) | best config pl / taper / ring-excl | OFF G_A / G_B / bal | OFF AR A / B (peak; -3 dB beam) | ON G_A / B supp | ON AR_A (peak; -3 dB beam) | SLL off / on | eta off / on | rotated slots | gate peak-AR / -3dB-AR (margin) |\n|---|---|---|---|---|---|---|---|---|---|---|\n"
for k in sorted(groups, key=lambda k: (-k[0], k[1])):
    r = max(groups[k], key=lambda r: g(r, "worst_margin_3dB_beam_db"))
    t1 += (f"| {k[0]} | {k[1]} ({f(r, 'f_sw')}) | {r['placement']} / {float(r['taper']):g} / {'yes' if int(g(r, 'excl')) else 'no'} | {f(r, 'off_GA', 1)} / {f(r, 'off_GB', 1)} / {f(r, 'off_balance')} | "
           f"{f(r, 'off_AR_A')} / {f(r, 'off_AR_B')}; {f(r, 'off_AR_A_3dB')} / {f(r, 'off_AR_B_3dB')} | {f(r, 'on_GA', 1)} / {f(r, 'on_suppression', 1)} | {f(r, 'on_AR_A')}; {f(r, 'on_AR_A_3dB')} | "
           f"{f(r, 'sll_off', 1)} / {f(r, 'sll_on', 1)} | {f(r, 'eta_off')} / {f(r, 'eta_on')} | {r['flip_count']} ({r['flip_switched']} sw) | {r['gate']} / {r['gate_3dB_beam_AR']} ({f(r, 'worst_margin_3dB_beam_db')}) |\n")


def mean(rs, k):
    return sum(g(r, k) for r in rs) / len(rs)


t2 = "| subset | points | ON AR_A: no A term / A term, unrotated / after joint rotation | SLL_on: same three | SLL_off after | B supp ON: no A / unrotated / after | MEETS_ALL peak-AR / -3 dB-AR |\n|---|---|---|---|---|---|---|\n"
subs = [("all", lambda r: True)]
subs += [(f"taper {t}", (lambda t: lambda r: abs(g(r, "taper") - t) < 1e-9)(t)) for t in (0.5, 0.8, 1.0)]
subs += [("ring exclusion off", lambda r: int(g(r, "excl")) == 0), ("ring exclusion on", lambda r: int(g(r, "excl")) == 1)]
subs += [("200 mm", lambda r: int(g(r, "side_mm")) == 200), ("220 mm", lambda r: int(g(r, "side_mm")) == 220), ("placement a", lambda r: r["placement"] == "a"), ("placement b", lambda r: r["placement"] == "b")]
for name, fn in subs:
    rs = [r for r in R if fn(r)]
    t2 += (f"| {name} | {len(rs)} | {mean(rs, 'noA_on_AR_A'):.2f} / {mean(rs, 'preflip_on_AR_A'):.2f} / {mean(rs, 'on_AR_A'):.2f} | {mean(rs, 'noA_sll_on'):.1f} / {mean(rs, 'preflip_sll_on'):.1f} / {mean(rs, 'sll_on'):.1f} | {mean(rs, 'sll_off'):.1f} | "
           f"{mean(rs, 'noA_on_supp'):.1f} / {mean(rs, 'preflip_on_supp'):.1f} / {mean(rs, 'on_suppression'):.1f} | {sum(r['gate'] == 'MEETS_ALL' for r in rs)} / {sum(r['gate_3dB_beam_AR'] == 'MEETS_ALL' for r in rs)} |\n")
t3 = "| diode C (pF) | switched OFF alpha mean (Np/m) | OFF G_A / G_B | balance | B contrast OFF-ON (dB) | ON B suppression | OFF AR A / B (peak; -3 dB) | ON AR_A (peak; -3 dB) | SLL off / on | gate |\n|---|---|---|---|---|---|---|---|---|---|\n"
for r in TOL:
    t3 += (f"| {g(r, 'C_pF'):.2f} | {f(r, 'sw_alpha_off_mean', 4)} | {f(r, 'off_GA', 1)} / {f(r, 'off_GB', 1)} | {f(r, 'off_balance')} | {g(r, 'off_GB') - g(r, 'on_GB'):.1f} | {f(r, 'on_suppression', 1)} | "
           f"{f(r, 'off_AR_A')} / {f(r, 'off_AR_B')}; {f(r, 'off_AR_A_3dB')} / {f(r, 'off_AR_B_3dB')} | {f(r, 'on_AR_A')}; {f(r, 'on_AR_A_3dB')} | {f(r, 'sll_off', 1)} / {f(r, 'sll_on', 1)} | {r['gate']} |\n")
bp = tb["point"]
rot = [g(r, "flip_count") / g(r, "N_cells") for r in R]
md = f"""# SUB6 board model, round 2: polarisation-aware synthesis (stage: sub6-board-model-r2)

**MODEL_ONLY / PENDING_HFSS_CELL_CHECK / NOT_HFSS_VALIDATED.** Same cells (`results/sub6_slotdisc/rows_RT5880.jsonl`, mesh-64 MoM, no HFSS check), same architecture and gates as `BOARD_MODEL.md` (round 1, untouched). Code `scripts/research/sub6_board_model2.py` (+ `sub6_board_model.py`, `sub6_board_req.py` unchanged); data `points_r2.csv` ({n} points), `tolerance_C.csv`, `layout_<side>mm_N<n>_<a|b>_t0.8.csv` (taper 0.8 with ring exclusion; columns include `rotated90`), `provenance_run2.json`. Launch fraction {prov['launch_fraction']:.4f}. {prov['seconds']:.0f} s for {n} points on {prov['workers']} workers.

## What changed
1. **Full tensor in the far field.** Cell current J = |E| [ (Delta/2) e(2a - phi_R) + A_eff e(phi_R - 90 deg) ], A_eff = A_perp sin^2(theta') + A_par cos^2(theta'), theta' = phi_R - a. The Delta term is the engine's (amplitude sqrt(2 alpha W / rho)); the A term is added with the complex ratio A_eff/(Delta/2) of the cell state (OFF_BAR50 / ON_BAR50; NOCUT rows have A = 0, so only switched cells carry it; mean |ratio| OFF about 5, ON about 0.05, but the OFF amplitude itself is tiny). A 90 deg rotation of a slot flips the Delta term and swaps sin^2 <-> cos^2. **Assumption (mine): the table alpha is the Delta-term leakage and the A term is added on top.** The A term alone has wavevector beta (bound), so it matters only through the discrete sampling of the sites.
2. **Joint axis choice (joint_opt2 style).** Greedy single 90-degree rotations (psi + 90 deg, same in OFF and ON) over ALL live sites, static and switched, with the corrected RHCP index (`copol_component(...,'rhcp')`). Cost: maximise |R_A|^2 in OFF and ON and |R_B|^2 in OFF, penalise |R_B|^2 in ON; constraints: LHCP/RHCP <= -{prov['xpd_db']:.0f} dB at beam A (both states) and beam B (OFF), co-pol floors {prov['floor_a_db']} dB (A) / {prov['floor_b_db']} dB (B) below the unrotated sums, B-in-ON no higher than max(unrotated, {prov['supp_db']:.0f} dB below B OFF); at most {prov['max_flips']} rotations. Rotated slots per point: {min(rot) * 100:.0f}-{max(rot) * 100:.0f} % of the cells (CP_STEP3 risk 3: mutual coupling of rotated cells not modelled).
3. **ON SLL.** Taper c in {{0.5, 0.8, 1.0}} (g = 1 - c (rho/rho_edge)^2 on the march target; round 1 used 0.5); switched cells excluded from lattice rings 1-2 around the feed (ring = max(|i|,|j|) <= 2) and from the outer ring (cells with a missing 8-neighbour) when "ring exclusion" is on. Both placements are restricted to the allowed sites (placement b = the stage-1b optimiser with the same swap rule).
4. **Sweep.** N_sw 180/200/220/240/260 at 200 mm ({int(g(R[0], 'N_cells'))} cells) and 220 mm, a/b, taper x3, exclusion x2 = {n} points. 220 mm is not scaled by the cell count here (N_sw as given).

## Gates (balance <= 3 dB, AR <= 3 dB, B suppression >= 12 dB, SLL <= -13 dB; SLL = worst of both states)
**Peak-AR gate (the job's gate, `metrics`: AR at the RHCP peak within 6 deg of the beam): MEETS_ALL {m_peak} of {n}. With AR taken as the worst over the -3 dB RHCP contour of the beam (stricter, my addition): MEETS_ALL {m_3} of {n}.** The joint rotation puts the cross-pol at the beam-A direction at its floor (-{prov['xpd_db']:.0f} dB gives AR 1.74 dB), so peak-AR margins are all about 1.26 dB; the -3 dB-beam AR is 0.2-1.2 dB worse. Best point by the stricter margin: {bp['side_mm'][:3]} mm, N_sw {bp['N_sw']}, placement {bp['placement']}, taper {float(bp['taper']):.2g}, ring exclusion {'yes' if int(float(bp['excl'])) else 'no'}.

Best configuration per (side, N_sw), chosen by the stricter margin:

{t1}
Gain is referred to accepted power x model conversion eta, so absolute G_A are model numbers; ON G_A is above OFF because the switched cells add conversion in ON (eta_on > eta_off).

## What each lever did (means over the subset)
{t2}
"no A term" = the round-1 current model; "A term, unrotated" = tensor with A but the designed axes; "after" = after the joint rotations. B suppression after rotation is held near the 17 dB cap (unrotated 19-29 dB): the rotations trade suppression for AR and SLL, the gate (12 dB) is still met.

## Tolerance of the diode capacitance at the best point (switched table restricted to the 5 rows of `results/sub6_slotdisc_ls_tol`)
{t3}
Source data `input_ls_tol.csv` (copied from github/compute/sub6-slotdisc-ls-tol 6d9c0ba). Design: the best point rebuilt with only the 5 ls_tol rows as switched cells (C = 0.15 pF nominal), then OFF alpha, beta and Delta of those rows replaced by the C = 0.12 / 0.20 values, sites, slot axes and rotations unchanged; A_perp kept at its nominal value (not in ls_tol). C acts on the OFF state only, so ON B suppression and ON AR do not change; the C effect is the OFF leakage of the switched cells (OFF alpha mean {f(TOL[0], 'sw_alpha_off_mean', 3)} at 0.12 pF and {f(TOL[2], 'sw_alpha_off_mean', 3)} at 0.20 pF against {f(TOL[1], 'sw_alpha_off_mean', 3)} nominal; at 0.20 pF one row exceeds the 0.1 Np/m filter), seen as OFF balance / OFF B and the B contrast OFF-ON.

## Limits
Cells are mesh-64 MoM, not HFSS-checked (2 of 80 rows have an HFSS package). The A-term model (see 1), the AR/XPD target at a single direction, and rotated slots next to unrotated ones without coupling are my modelling choices. Mixed gap / W neighbours: coupling not modelled. Stopped here: no cell choice, no board design, no HFSS package.
"""
(D / "BOARD_MODEL_R2.md").write_text(md, encoding="utf-8", newline="\n")
print(m_peak, m_3, n)
