"""Writes results/sub6_board_island/BOARD_ISLAND.md from points_island.csv (MODEL_ONLY)."""
import csv
import json
from collections import defaultdict
from pathlib import Path

D = Path(__file__).resolve().parents[2] / "results" / "sub6_board_island"
R = list(csv.DictReader(open(D / "points_island.csv")))
prov = json.load(open(D / "provenance_run_island.json"))
g = lambda r, k: float(r[k])  # noqa: E731
f = lambda r, k, d=2: f"{g(r, k):.{d}f}"  # noqa: E731
n = len(R)
m3 = sum(r["gate_3dB_beam_AR"] == "MEETS_ALL" for r in R)
mp = sum(r["gate"] == "MEETS_ALL" for r in R)
groups = defaultdict(list)
for r in R:
    groups[(int(g(r, "side_mm")), int(g(r, "N_sw")))].append(r)
t1 = "| side | N_sw (f) | best pl / taper / ring-excl | OFF G_A / G_B / bal | OFF AR3dB A / B | ON G_A / B supp | ON AR3dB | SLL off / on | eta off / on | rotated | MEETS_ALL (-3 dB AR / peak AR) of 12 | best gate (margin) |\n|---|---|---|---|---|---|---|---|---|---|---|---|\n"
for k in sorted(groups, key=lambda k: (-k[0], k[1])):
    rs = groups[k]
    r = max(rs, key=lambda r: g(r, "worst_margin_3dB_beam_db"))
    t1 += (f"| {k[0]} | {k[1]} ({f(r, 'f_sw')}) | {r['placement']} / {float(r['taper']):g} / {'yes' if int(g(r, 'excl')) else 'no'} | {f(r, 'off_GA', 1)} / {f(r, 'off_GB', 1)} / {f(r, 'off_balance')} | {f(r, 'off_AR_A_3dB')} / {f(r, 'off_AR_B_3dB')} | "
           f"{f(r, 'on_GA', 1)} / {f(r, 'on_suppression', 1)} | {f(r, 'on_AR_A_3dB')} | {f(r, 'sll_off', 1)} / {f(r, 'sll_on', 1)} | {f(r, 'eta_off')} / {f(r, 'eta_on')} | {r['flip_count']} | "
           f"{sum(x['gate_3dB_beam_AR'] == 'MEETS_ALL' for x in rs)} / {sum(x['gate'] == 'MEETS_ALL' for x in rs)} | {r['gate_3dB_beam_AR']} ({f(r, 'worst_margin_3dB_beam_db')}) |\n")
t2 = "| subset | points | ON supp: unrotated / rotated | ON AR_A 3dB beam | OFF AR_B 3dB beam | SLL_on | MEETS_ALL (-3 dB AR) |\n|---|---|---|---|---|---|---|\n"
subs = [("all", lambda r: True), ("ring exclusion on", lambda r: int(g(r, "excl")) == 1), ("ring exclusion off", lambda r: int(g(r, "excl")) == 0), ("taper 0.5", lambda r: abs(g(r, "taper") - 0.5) < 1e-9),
        ("taper 0.8", lambda r: abs(g(r, "taper") - 0.8) < 1e-9), ("taper 1.0", lambda r: abs(g(r, "taper") - 1.0) < 1e-9), ("placement a", lambda r: r["placement"] == "a"), ("placement b", lambda r: r["placement"] == "b"),
        ("N_sw <= 135 (diodes on hand)", lambda r: int(g(r, "N_sw")) <= 135), ("N_sw >= 180 (ideal)", lambda r: int(g(r, "N_sw")) >= 180)]
mean = lambda rs, k: sum(g(r, k) for r in rs) / len(rs)  # noqa: E731
for name, fn in subs:
    rs = [r for r in R if fn(r)]
    t2 += (f"| {name} | {len(rs)} | {mean(rs, 'preflip_on_supp'):.1f} / {mean(rs, 'on_suppression'):.1f} | {mean(rs, 'on_AR_A_3dB'):.2f} | {mean(rs, 'off_AR_B_3dB'):.2f} | {mean(rs, 'sll_on'):.1f} | "
           f"{sum(r['gate_3dB_beam_AR'] == 'MEETS_ALL' for r in rs)} |\n")
best = max(R, key=lambda r: g(r, "worst_margin_3dB_beam_db"))
md = f"""# SUB6 board model with ISLAND switched cells (stage: sub6-board-island)

**MODEL_ONLY / PENDING_HFSS_CELL_CHECK / NOT_HFSS_VALIDATED.** Code `scripts/research/sub6_board_model4.py` (derived from `sub6_board_model2.py`, same Phase 26 engine, real per-state alpha/beta, joint 90-degree rotations, taper, ring exclusion, placements a/b); data `points_island.csv` ({n} points), `layout_*_t0.8.csv`, `provenance_run_island.json`, inputs `input_island_rows.json/csv` (copied from compute/sub6-island-cell 6d53878). Launch fraction {prov['launch_fraction']} (standard SMA pin, sub6-launcher). {prov['seconds']:.0f} s on {prov['workers']} workers.

## Architecture and cells
STATIC = NOCUT rows ({prov['static_rows']} rows, filters as in round 1). SWITCHED = island rows of the island scan **without the resonance flag** (16 of 48), further filtered by me to alpha_OFF >= 0.9 Np/m and sqrt(alpha_ON/alpha_OFF) <= 0.5 (near-silent in ON): **{prov['switched_rows']} rows**, alpha_OFF 0.96-2.70 Np/m, beta_OFF 1.08-1.11, beta_ON 1.00-1.05 (each cell's own beta in each state). Static cells radiate beam A in both states; island cells radiate beam B in OFF. General weights, rho_d = mean signed ratio J_ON/J_OFF of the table (about -0.33), c = -rho_d/(1-rho_d): static J = (u_A + c u_B)/(1-f), switched J_OFF = u_B/((1-rho_d) f), J_ON = rho_i J_OFF with rho_i = sign(Im Delta_ON/Im Delta_OFF) sqrt(alpha_ON/alpha_OFF) of the chosen row (my choice: magnitude from alpha, sign from the reactance polarity; the complex phase of Delta_ON/Delta_OFF is not used). The island cells have no A_perp (mirror symmetric), so the rotations act on the Delta term only.
Sweep: N_sw 110/120/135 (diodes on hand) and 180/240/300 (ideal) x sides 180/200/220 mm x placement a/b x taper 0.5/0.8/1.0 x ring exclusion off/on = {n} points.

## Gates (balance <= 3 dB, AR <= 3 dB over the -3 dB beam, B suppression >= 12 dB, SLL <= -13 dB both states)
**MEETS_ALL on the stricter -3 dB-beam AR gate: {m3} of {n}** (peak-AR gate: {mp} of {n}). Best point: {int(g(best, 'side_mm'))} mm, N_sw {best['N_sw']}, placement {best['placement']}, taper {float(best['taper']):g}, ring exclusion {'yes' if int(g(best, 'excl')) else 'no'}, margin {f(best, 'worst_margin_3dB_beam_db')} dB.

{t1}
{t2}
Rotation lowers B suppression (cap 17 dB set in the optimiser; pre-rotation values are in the table) and buys AR. Gains are referred to accepted power x model eta x launch 0.656.

## Limits
Switched table is only {prov['switched_rows']} rows (alpha_req is nearly uniform, so one or two rows are used; amplitude mismatch rms 0-8 % at the checked points). rho_d is a table mean, the per-site spread of rho_i leaves a residual B in ON. Island cells: mesh 64 with 0.3-0.6 mm clearances of 2-3 mesh cells (qualitative), no HFSS check, no bias network, floating-island coupling to neighbours not modelled, the OFF-state beta of island cells (1.08-1.11) differs from static cells (1.05-1.2). The AR over the -3 dB beam is from the power formula of `metrics`. Stopped here.
"""
(D / "BOARD_ISLAND.md").write_text(md, encoding="utf-8", newline="\n")
print(n, m3, mp)
