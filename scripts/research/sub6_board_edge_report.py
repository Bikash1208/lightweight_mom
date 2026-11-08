"""Writes results/sub6_board_edge/BOARD_EDGE.md from points_edge.csv (free segments) and points_edge_chords.csv (EDGE_TO_EDGE + free segments at the same N) (MODEL_ONLY)."""
import csv
import json
from collections import Counter
from pathlib import Path

D = Path(__file__).resolve().parents[2] / "results" / "sub6_board_edge"
R = list(csv.DictReader(open(D / "points_edge.csv")))
RC = list(csv.DictReader(open(D / "points_edge_chords.csv")))
prov = json.load(open(D / "provenance_run_edge.json"))
g = lambda r, k: float(r[k])  # noqa: E731
f = lambda r, k, d=2: f"{g(r, k):.{d}f}"  # noqa: E731
n = len(R) // 4
KEYS = lambda rows: sorted({(r["side_mm"], r["N_sw"]) for r in rows}, key=lambda k: (-float(k[0]), int(float(k[1]))))  # noqa: E731


def best_row(rows, key, pl):
    rs = [r for r in rows if r["variant"] == "V3" and (r["side_mm"], r["N_sw"]) == key and r["placement"] == pl]
    return max(rs, key=lambda r: g(r, "margin"))


t1 = "| side | N_sw | place | best V3 config (taper / excl / w) | OFF G_A / G_B / bal | OFF AR_A / AR_B (-3 dB) | ON G_A / B supp / AR_A | SLL off / on | eta off / on | lines (x) len | string V min-max | diodes | gate (margin) |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|\n"
for key in KEYS(R):
    for pl in ("a", "b"):
        b = best_row(R, key, pl)
        t1 += (f"| {int(float(key[0]))} | {int(float(key[1]))} | {pl} | {float(b['taper']):g} / {'yes' if int(float(b['excl'])) else 'no'} / {float(b['w']):g} | {f(b, 'off_GA', 1)} / {f(b, 'off_GB', 1)} / {f(b, 'bal')} | {f(b, 'off_ARA3')} / {f(b, 'off_ARB3')} | "
               f"{f(b, 'on_GA', 1)} / {f(b, 'supp', 1)} / {f(b, 'on_ARA3')} | {f(b, 'sll_off', 1)} / {f(b, 'sll_on', 1)} | {f(b, 'eta_off')} / {f(b, 'eta_on')} | {b['n_lines']} ({b['n_lines_x']}) {b['line_len_min']}-{b['line_len_max']} | "
               f"{f(b, 'string_V_min', 1)}-{f(b, 'string_V_max', 1)} V | {b['n_diodes_in_strings']} | {b['gate']} ({f(b, 'margin')}) |\n")
t2 = "| variant | designs | MEETS_ALL | mean margin | best margin |\n|---|---|---|---|---|\n"
for v in ("V0", "V1", "V2", "V3"):
    rs = [r for r in R if r["variant"] == v]
    t2 += f"| {v} | {len(rs)} | {sum(r['gate'] == 'MEETS_ALL' for r in rs)} | {sum(g(r, 'margin') for r in rs) / len(rs):.2f} | {max(g(r, 'margin') for r in rs):.2f} |\n"
t3 = "| side | N_sw | place | G_A ON, beta shift off / on | cost (dB) | sw d_beta (ON-OFF) | static dphi rms / max (deg) | slot-to-string angle: mean (deg), fraction within 7.5 deg of 90 / of 0, 45 or 90 | hist 15-deg bins (0-90) |\n|---|---|---|---|---|---|---|---|---|\n"
for key in KEYS(R):
    for pl in ("a", "b"):
        b = best_row(R, key, pl)
        t3 += (f"| {int(float(key[0]))} | {int(float(key[1]))} | {pl} | {f(b, 'on_GA_beta_shift_off', 1)} / {f(b, 'on_GA_beta_shift_on', 1)} | {f(b, 'beta_shift_cost_dB')} | {f(b, 'sw_dbeta_on_minus_off', 3)} | {f(b, 'static_dphi_rms_deg', 1)} / {f(b, 'static_dphi_max_deg', 0)} | "
               f"{f(b, 'theta_rel_mean_deg', 0)}, {f(b, 'theta_frac_within_tol_of_90')} / {f(b, 'theta_frac_within_tol_of_0_45_90')} | {b['theta_hist_15deg_bins']} |\n")
V3c = [r for r in RC if r["variant"] == "V3"]
ch = [r for r in V3c if r["placement"] in ("cx", "cy", "c")]
fr = [r for r in V3c if r["placement"] in ("a", "b")]
t4 = "| side | target N | chords | actual N | lines (x) | chord len min-max | longest string (diodes, V) | total diodes | best gate (margin) | free segments at the same N: best gate (margin) |\n|---|---|---|---|---|---|---|---|---|---|\n"
for side in ("200.0000", "220.0000"):
    for tg in ("135", "240"):
        for pl, name in (("cx", "x only"), ("cy", "y only"), ("c", "x+y")):
            rs = [r for r in ch if r["side_mm"] == side and r["n_target"] == tg and r["placement"] == pl]
            b = max(rs, key=lambda r: g(r, "margin"))
            fs = [r for r in fr if r["side_mm"] == side and int(float(r["N_sw"])) == int(float(b["N_sw_live"])) and r["excl"] == b["excl"]]
            fb = max(fs, key=lambda r: g(r, "margin"))
            t4 += (f"| {int(float(side))} | {tg} | {name} | {int(float(b['N_sw_live']))} | {b['n_lines']} ({b['n_lines_x']}) | {b['line_len_min']}-{b['line_len_max']} | {int(float(b['line_len_max'])) - 1} diodes, {f(b, 'string_V_max', 1)} V | {b['n_diodes_in_strings']} | "
                   f"{b['gate']} ({f(b, 'margin')}) | {fb['gate']} ({f(fb, 'margin')}, place {fb['placement']}) |\n")
meets = [r for r in R if r["gate"] == "MEETS_ALL"]
need = [g(r, "theta_frac_need_angle_table") for r in V3c]
md = f"""# SUB6 board model with EDGE-DIODE switched strings (stage: sub6-board-edge, + EDGE_TO_EDGE amendment)

**MODEL_ONLY / PENDING_HFSS_CELL_CHECK / NOT_HFSS_VALIDATED. NEED_ANGLE_TABLE (see below): every gate on this page is provisional, because the 90-deg table row is applied to cells whose hologram angle to the string is not in {{0, 45, 90}} deg.**
Code `scripts/research/sub6_board_edge.py`, `sub6_board_edge2.py` (chords); they import `sub6_board_island2.py` and `sub6_board_model4.py` unchanged (V3 method as the island follow-up). Data `points_edge.csv` ({len(R)} rows = {n} designs x V0-V3), `points_edge_chords.csv`, `provenance_run_edge.json`, input `input_edge_diode.json/csv` (pushed commit 34c4a37 of compute/sub6-edge-diode-cell). Launch fraction {prov['launch_fraction']}.

## NEED_ANGLE_TABLE
The edge-diode table has the slot at theta = 0, 45 and 90 deg to the diode axis; only the 90-deg rows are used here. The hologram sets the slot axis of every cell, so the angle to its string takes all values: for the switched cells of the V3 designs the angle is within 7.5 deg of 0, 45 or 90 deg for only {100 * (1 - max(need)):.0f}-{100 * (1 - min(need)):.0f} % of the cells (mean over designs {100 * (1 - sum(need) / len(need)):.0f} %); within 7.5 deg of 90 deg only about 20 % (histograms in the table below). The remaining cells need an angle table (alpha, beta, Delta of the edge cell versus the slot angle) from the cell session; **nothing was interpolated or guessed, and the earlier "theta-aware amplitude factor" (0.65-0.70) is withdrawn: it was an interpolation between three rows.** Until the angle table exists the B-beam numbers are not valid for those cells.

## Free segments (10-24 cells, series strings along x or y)
STATIC = NOCUT rows ({prov['static_rows']}); SWITCHED = edge-diode theta-90 rows filtered by me to alpha_OFF >= 0.9 and sqrt(alpha_ON/alpha_OFF) <= 0.5: {prov['switched_rows']} rows {prov['switched_ids']} (beta_OFF 1.20-1.30, beta_ON 1.017). Edge rows carry no resonance proxy (resonance_flip = null: not computed, not "no resonance").
MEETS_ALL: {len(meets)} of {len(R)} rows. {t2}
Best V3 design per (side, N, placement):
{t1}
### Beta-shift effect on beam A and slot-to-string angles
{t3}
The switched cells' beta falls by 0.19-0.20 from OFF to ON (1.20-1.30 to 1.017); a ray over a 24-cell string accumulates about 370 deg of phase difference between the states; static-cell phase error rms 19-65 deg; cost to G_A,ON 0.3-1.0 dB (one design -0.7 dB).

## EDGE_TO_EDGE: every switched line is a full chord of the hexagon (both string ends at the rim)
Chords = complete lattice rows (x) or columns (y) from rim to rim; the row j = 0 and column i = 0 are interrupted by the feed clearance and are not chords; with ring exclusion on, no chord within 2 rows/columns of the feed. Chords are chosen greedily on the ideal-cell projection objective until the cell count is close to the target (135 / 240); x and y chords may not share a cell (a crossing is not a series string); y-only and x+y gave the same set whenever the greedy picks only y chords. Same V3 evaluation.
{t4}
Chord designs MEETS_ALL: {sum(r['gate'] == 'MEETS_ALL' for r in ch)} of {len(ch)}; free segments at the same N: {sum(r['gate'] == 'MEETS_ALL' for r in fr)} of {len(fr)}. **Longest strings: 28 diodes = 26.6 V (x chord, 200 mm), 32 diodes = 30.4 V (y chord, 200 mm), up to 36 diodes = 34.2 V (y chord, 220 mm)**; a full chord is thus always too long for a 24-cell string, the required voltage is reported regardless. Chord total diodes 120-224 (150 on hand: chords at N about 130 fit, N about 225 does not).

## Limits
NEED_ANGLE_TABLE (above); cell values are periodic-lattice values (every cell has a diode and a 1.00 mm gap, disc r 5.30), string ends and neighbours are not modelled; no HFSS check; no bias supply modelled; the chord selection and the row filter are mine. Stopped here.
"""
(D / "BOARD_EDGE.md").write_text(md, encoding="utf-8", newline="\n")
print(len(meets), len(R), sum(r["gate"] == "MEETS_ALL" for r in ch), len(ch))
