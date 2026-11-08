"""Writes results/sub6_board_edge/BOARD_EDGE.md from points_edge.csv (MODEL_ONLY)."""
import csv
import json
from collections import Counter
from pathlib import Path

D = Path(__file__).resolve().parents[2] / "results" / "sub6_board_edge"
R = list(csv.DictReader(open(D / "points_edge.csv")))
prov = json.load(open(D / "provenance_run_edge.json"))
g = lambda r, k: float(r[k])  # noqa: E731
f = lambda r, k, d=2: f"{g(r, k):.{d}f}"  # noqa: E731
n = len(R) // 4
t1 = "| side | N_sw | place | best V3 config (taper / excl / w) | OFF G_A / G_B / bal | OFF AR_A / AR_B (-3 dB) | ON G_A / B supp / AR_A | SLL off / on | eta off / on | rotated | lines (x) len min-max | string V min-max | diodes | gate (margin) |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n"
for key in sorted({(r["side_mm"], r["N_sw"]) for r in R}, key=lambda k: (-float(k[0]), int(float(k[1])))):
    for pl in ("a", "b"):
        rs = [r for r in R if r["variant"] == "V3" and (r["side_mm"], r["N_sw"]) == key and r["placement"] == pl]
        b = max(rs, key=lambda r: g(r, "margin"))
        t1 += (f"| {int(float(key[0]))} | {int(float(key[1]))} | {pl} | {float(b['taper']):g} / {'yes' if int(float(b['excl'])) else 'no'} / {float(b['w']):g} | {f(b, 'off_GA', 1)} / {f(b, 'off_GB', 1)} / {f(b, 'bal')} | {f(b, 'off_ARA3')} / {f(b, 'off_ARB3')} | "
               f"{f(b, 'on_GA', 1)} / {f(b, 'supp', 1)} / {f(b, 'on_ARA3')} | {f(b, 'sll_off', 1)} / {f(b, 'sll_on', 1)} | {f(b, 'eta_off')} / {f(b, 'eta_on')} | {b['rotated']} | {b['n_lines']} ({b['n_lines_x']}) {b['line_len_min']}-{b['line_len_max']} | "
               f"{f(b, 'string_V_min', 1)}-{f(b, 'string_V_max', 1)} V | {b['n_diodes_in_strings']} | {b['gate']} ({f(b, 'margin')}) |\n")
t2 = "| variant | designs | MEETS_ALL | mean margin | best margin | gate counts |\n|---|---|---|---|---|---|\n"
for v in ("V0", "V1", "V2", "V3"):
    rs = [r for r in R if r["variant"] == v]
    t2 += f"| {v} | {len(rs)} | {sum(r['gate'] == 'MEETS_ALL' for r in rs)} | {sum(g(r, 'margin') for r in rs) / len(rs):.2f} | {max(g(r, 'margin') for r in rs):.2f} | {dict(Counter(r['gate'] for r in rs).most_common(4))} |\n"
t3 = "| side | N_sw | place | G_A ON, beta shift OFF (OFF phase) / ON | cost (dB) | sw d_beta (ON-OFF) | static dphi rms / max (deg) | theta(slot, line) mean (deg), fraction >= 67.5 | theta-aware amplitude factor (mean, fraction < 0.5) |\n|---|---|---|---|---|---|---|---|---|\n"
for key in sorted({(r["side_mm"], r["N_sw"]) for r in R}, key=lambda k: (-float(k[0]), int(float(k[1])))):
    for pl in ("a", "b"):
        rs = [r for r in R if r["variant"] == "V3" and (r["side_mm"], r["N_sw"]) == key and r["placement"] == pl]
        b = max(rs, key=lambda r: g(r, "margin"))
        t3 += (f"| {int(float(key[0]))} | {int(float(key[1]))} | {pl} | {f(b, 'on_GA_beta_shift_off', 1)} / {f(b, 'on_GA_beta_shift_on', 1)} | {f(b, 'beta_shift_cost_dB')} | {f(b, 'sw_dbeta_on_minus_off', 3)} | {f(b, 'static_dphi_rms_deg', 1)} / {f(b, 'static_dphi_max_deg', 0)} | "
               f"{f(b, 'theta_rel_mean_deg', 0)}, {f(b, 'theta_rel_frac_ge67', 2)} | {f(b, 'theta_aware_amp_factor_mean')}, {f(b, 'theta_aware_frac_below_half')} |\n")
meets = [r for r in R if r["gate"] == "MEETS_ALL"]
md = f"""# SUB6 board model with EDGE-DIODE switched strings (stage: sub6-board-edge)

**MODEL_ONLY / PENDING_HFSS_CELL_CHECK / NOT_HFSS_VALIDATED.** Code `scripts/research/sub6_board_edge.py` (imports `sub6_board_island2.py` and `sub6_board_model4.py` unchanged; V3 method as the island follow-up); data `points_edge.csv` ({len(R)} rows = {n} designs x V0-V3), `provenance_run_edge.json`, input `input_edge_diode.json/csv` (compute/sub6-edge-diode-cell 34c4a37). Launch fraction {prov['launch_fraction']}. {prov['seconds']:.0f} s on {prov['workers']} workers.

## Setup
STATIC = NOCUT rows ({prov['static_rows']}); SWITCHED = edge-diode rows with the slot at 90 deg to the diode axis, real OFF/ON alpha, beta and Delta, filtered by me to alpha_OFF >= 0.9 Np/m and sqrt(alpha_ON/alpha_OFF) <= 0.5: **{prov['switched_rows']} rows** {prov['switched_ids']} (alpha_OFF 1.28-2.96, alpha_ON 0.21, **beta_OFF 1.20-1.30, beta_ON 1.017**). The edge-diode scan has no resonance proxy. Architecture, weights (rho_d about -0.3), rotations and V0-V3 as in the island follow-up. Switched cells lie in straight lines of consecutive cells along x or y (length 10-24, no overlap, no end-to-end contact of same-axis lines), each line one series string of n-1 diodes (0.95 V each at 10 mA); placement a = random compatible lines, b = greedy + remove/replace local search on the ideal-cell projection objective (exactly N_sw cells). Sweep: sides 200 / 220 mm, N_sw 135 / 240, taper 0.5/0.8/1.0, ring exclusion off/on.

## Gates
**MEETS_ALL: {len(meets)} of {len(R)} rows** ({', '.join(f"{int(float(r['side_mm']))} mm N{r['N_sw']} {r['placement']} {r['variant']}" for r in meets) or 'none'}); the misses are SLL (OFF-state sidelobes of the string structure, -9 to -14 dB) and, at N 240, balance / suppression.
{t2}
Best V3 design per (side, N, placement):
{t1}
## Beta-shift effect on beam A and the slot-angle caveat
{t3}
The switched cells' beta drops by 0.19-0.20 from OFF to ON (1.20-1.30 to 1.017); a ray crossing a string of 24 cells accumulates up to 24 x k0 x 0.19 x 11.6 mm = about 370 deg of phase difference between the states, so the static-cell phase error is large (rms 19-65 deg) and G_A,ON changes by up to 1 dB.
The table row (slot at 90 deg to the diode axis) is used for every switched cell, but the hologram sets the slot axis per cell: the mean angle between slot and string is about 40-45 deg and only about 20 % of the switched cells are within 22.5 deg of 90 deg; with a theta-aware alpha (rows theta = 0 / 45 / 90, same q and W) the mean switched amplitude would be only about 0.65-0.70 of the value used here, and about 20 % of the switched cells would fall below half amplitude. **The results are therefore optimistic for B; this is the main limit of the model.**

## Limits
Cell values are periodic-lattice values (every cell has a diode to its +x neighbour, disc r 5.30, gap 1.00); string ends and neighbours of a different kind are not modelled; the edge-diode rows are not resonance-checked; no bias supply, no HFSS check; the rows' sheet resolution is 6 x 2 mesh cells. Strings of 10-24 cells need 8.6-21.9 V at 10 mA (all strings individually, more than one string cannot be in series unless adjacent). Diode count N_sw - n_lines: 124-128 for N 135 (within the 150 on hand), 221-227 for N 240 (above 150). Stopped here.
"""
(D / "BOARD_EDGE.md").write_text(md, encoding="utf-8", newline="\n")
print(len(meets), len(R))
