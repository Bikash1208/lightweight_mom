# Diodes only at key areas — MODEL_ONLY

Launch zone: the model reserves the centre feed cell and regularizes ρ at 5.8 mm; it does not solve where the launch becomes a bound surface wave. Thus 5.8 mm is a model sweep floor, not a measured physical launch-zone radius.
Maximum swept gate inner radius with both peak drops ≥10 dB for a perfect gate is 60 mm; no 75 mm configuration passes both peak-drop gates.
The existing model assigns independent A/RHCP and B/LHCP phase ramps by construction; physical CP is unvalidated.

Power obeys dP/dr=−2αP; α=0 passes power outward. The same anchor α profile applies in both halves: 1.11/1.45/1.83 Np/m at r=0/75/150 mm. R=150 mm, 11.6-mm pitch, 30-mm dither, eight seeds, 2° angular grid. A-only and B-only are computed separately.

## Study 1 — two semicircular gates

Inner radius 5.8/10/20/30/40/50/60/75 mm; 1–5 cell rows; blocked-state transmitted power 0/1/2/5/10%; one diode per gate cell. Half-only covers assigned A cells at x≥0 and B cells at x≤0. Extended span adds assigned A cells down to x=−30 mm and B cells up to x=+30 mm. Pass-state gate cells either retain anchor α or have α=0. ABSORBING and REFLECTIVE have identical forward fields at fixed power transmission; the removed power is counted as absorbed or reflected respectively. The redistributed bound directs it into the open geometric half at the gate; for ABSORBING this is only a mathematical upper bound.

| r inner mm | Rows | T % | Block / pass / span / bound | Diodes max | A-only: B peak/power drop; A gain; SLL dB | B-only: A peak/power drop; B gain; SLL dB | BOTH A/B cost dB | Both states |
|---:|---:|---:|---|---:|---|---|---|---|
| 5.8 | 2 | 0 | REFLECTIVE / anchor_like / half / lost | 18 | 22.04/12.06; -2.38; -14.62 | 21.86/11.13; -2.53; -14.22 | 0.00/0.00 | FAIL |
| 5.8 | 2 | 0 | REFLECTIVE / anchor_like / half / redistributed | 18 | 19.78/7.67; 1.57; -14.56 | 19.27/7.92; 1.41; -14.21 | 0.00/0.00 | PASS |
| 5.8 | 2 | 1 | REFLECTIVE / anchor_like / half / redistributed | 18 | 14.16/7.47; 1.71; -14.93 | 14.05/7.70; 1.57; -14.74 | 0.00/0.00 | PASS |
| 5.8 | 2 | 0 | REFLECTIVE / nonleaking / half / redistributed | 18 | 19.33/7.14; 1.46; -15.26 | 18.75/7.41; 1.30; -15.46 | 0.39/0.38 | PASS |
| 10 | 2 | 0 | REFLECTIVE / anchor_like / half / redistributed | 21 | 21.65/7.41; 1.60; -14.84 | 21.12/7.33; 1.52; -15.10 | 0.00/0.00 | PASS |
| 10 | 2 | 1 | REFLECTIVE / anchor_like / half / redistributed | 21 | 15.05/7.21; 1.74; -15.26 | 14.88/7.14; 1.66; -15.43 | 0.00/0.00 | PASS |
| 10 | 2 | 2 | REFLECTIVE / anchor_like / half / redistributed | 21 | 13.32/7.03; 1.78; -15.40 | 13.19/6.96; 1.70; -15.56 | 0.00/0.00 | PASS |

Every one-beam case, including suppressed-hand wide-angle level and SLL, is in `study1_gates_rows.csv` (51200 seed rows); `study1_gates_summary.csv` (6400) and `study1_joint.csv` (3200) aggregate the cases. PASS requires ≥10 dB peak drop and no surviving-beam gain decrease for every seed. Power drop is the integrated suppressed-hand radiated conversion change; gain is referenced to equal launched power. No lost-power-bound configuration meets the full requirement for both one-beam states; the full passes require the ideal redistribution bound.

## Study 2 — selected anchor cells with α=0

Inner/outer rings select target-hand cells by radius. Axis patch selects cells nearest (±75,0) mm; upper/lower boundary patches select cells nearest (±20,±75) mm. Power released by switched cells remains on their radial rays.

| Suppressed | Placement | Fraction of target cells for 10 dB (all seeds) | Fraction of whole board | Fraction meeting 10 dB and no survivor loss |
|---|---|---:|---:|---:|
| A | inner_ring | 0.702 | 0.350 | 0.702 |
| A | outer_ring | 0.802 | 0.400 | 0.802 |
| A | axis_patch | 0.702 | 0.350 | 0.702 |
| A | upper_boundary_patch | 0.702 | 0.350 | 0.702 |
| A | lower_boundary_patch | 0.702 | 0.350 | 0.702 |
| B | inner_ring | 0.702 | 0.352 | 0.702 |
| B | outer_ring | 0.802 | 0.402 | 0.802 |
| B | axis_patch | 0.702 | 0.352 | 0.702 |
| B | upper_boundary_patch | 0.702 | 0.352 | 0.702 |
| B | lower_boundary_patch | 0.702 | 0.352 | 0.702 |

At the 0.702 target-half threshold, independently switching either beam would require about 180 A-cell plus 182 B-cell diodes (about 362 total). The outer-ring threshold requires about 414 total. These selected-cell schemes exceed the earlier 200-diode preference.
`study2_selected_rows.csv` (1600 seed rows), `study2_selected_summary.csv` (200) and `study2_thresholds.csv` provide the full placement/fraction tables.

FLAG: This is a coherent ideal-CP array-factor model, not HFSS. Gate reflection phase, coupling into other modes, feed mismatch, diode impedance, hardware loss, and the physical surface-wave launch radius are unknown. Blocked-power redistribution is an upper bound, not a predicted circuit response. Some broad gates leave no open-half ray to receive redistributed power and are marked infeasible. Reported power conservation is radial bookkeeping, not an electromagnetic scattering solution.

HANDOFF: MODEL_ONLY | both one-beam states and BOTH evaluated | full R150 aperture | no HFSS | next: validate a realizable gate or selected-cell pattern.
