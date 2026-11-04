# Item 2 — architecture study (idealised cells, MODEL_ONLY)

Script `scripts/research/sidequest_arch_study.py` (about 2.7 min on 4 cores) -> `arch_study.json`. 150 mm hexagon (circle of equal
area as comparison), 10.54 GHz, 6.4 mm lattice, one centre probe, launch 0.62, beta = 1.065 k0, eta_d = 0.9 per ray.
Beams: 2 states (25, 0)/(25, 180), pols RHCP/LHCP; 4 states (25, 0/90/180/270), pols RHCP/LHCP/LP-H/LP-V.
Pol convention: local (theta_hat, phi_hat) frame as `copol_component` (LP-H = theta_hat, LP-V = phi_hat, repo E_R/E_L).

**Model limits (read first).** Cells are ideal tangential-current point sources with the polarisation imposed by hand, so XPD/AR
are ideal by construction (> 60 dB / 0 dB) and say nothing about real cells: polarisation purity is an item-3 question. The gate is
only an isolation number (dB per strip), not a cell. No TM/TE split vs ray/beam angle, no edge diffraction, no mutual coupling.
Gain = D x 0.62 x f_az x eta_conv x kappa (referred to accepted power; no HFSS offset applied). Absolute gains are upper-bound-like.

## Result table (hexagon, gate 3 rows at r_g = 13 mm for A)
| Option | diodes | gain per state [dBi] | SLL [dB] | ap. eff. | isolation / contrast for off-beam <= -15 dB |
|---|---|---|---|---|---|
| A, 2 quadrants | 68 | 19.4 / 19.4 | -13.3 | 0.72 | 16 dB per gate |
| A, 4 quadrants | 68 | 13.6 / 13.3 / 13.7 / 13.2 | -12.9 to -14.2 | 0.73 | 17 dB per gate |
| B, 2 quadrants (full) | 710 (50 % of 1420 cells) | 12.2 / 12.2 | -11.1 | 0.26 | 18 dB (amplitude ratio of radiating/non-radiating) |
| B, 4 quadrants (full) | 752 (53 %) | 6.6 | -6.2 to -10.9, **-3.1 (LP-V)** | 0.30 | 40 dB (SLL-floor limited) |
| B, 2 q, keep 50 % / 35 % | 356 / 248 | 9.9 / 8.4 | -9.9 / -8.3 | | |
| B, 4 q, keep 50 % / 35 % | 376 / 262 | 4.4 / 3.5 | -4.9 / -5.1 (LP-V -2.0) | | |

A gate-row sweep (2 / 3 / 4 rows) changes only the diode count (36 / 68 / 100) and the lost aperture area (ap. eff. 0.70-0.75);
the isolation per row count is NOT computed here (needs the gate cell, item 3). Circle vs hexagon: within 0.2 dB, same trends.

## Findings
1. **A beats B on every figure of merit in this model**: 68 diodes vs 710-750, +7 dB gain (2q) and +7 dB (4q), SLL 13 dB vs 3-11 dB,
   and a leakage requirement of about 16-17 dB isolation per gate. B needs a modulated-cell radiating-power contrast of 18 dB (2q) to
   40 dB (4q) because its 1-bit hologram has a high harmonic floor. B also needs ~2x the local leakage constant of A's cells
   (`alpha_mult` 1.8-2.0, a cell-design requirement not yet checked) and thinning below 100 % lowers gain/SLL further.
2. **Centre feed costs 10 log(nq) in power plus the same in area**: the quadrant gets only 1/nq of the launched surface wave and 1/nq of
   the aperture. 2q: 19.4 dBi. 4q: 13.5 dBi (-5.9 dB vs 2q). A single 150 mm hexagon, all area, is 23.0 dBi (HFSS).
   An offset feed (feed inside the active quadrant) is not possible for several quadrants with ONE probe; not evaluated.
3. **Default state with all diodes OFF is possible for A** with two gate-cell types (normally-pass for quadrant 1, normally-block for
   the others), each quadrant's gate a separately biased group (4 groups for 4q). This forces the cell to be designed in both polarities.
4. **Blocked power goes back to the feed in A.** With a reflective block, up to 1 - 1/nq of the launched surface wave returns to the
   probe: coherent S11 perturbation amplitude bound sqrt(0.62 (1 - 1/nq)) = 0.56 (2q) / 0.68 (4q), against the HFSS baseline
   |S11| = 0.18 (-15 dB). Worst case S11 ~ -2 to -3 dB; it depends on the gate ring radius (reflection phase), so r_g must be tuned,
   or the block state made lossy (absorbing). B has no such return path (inactive quadrants just pass the wave to the edge), but
   that unconverted power is lost to grazing/edge radiation (Phase 30 saw 48-68 % of port power there).
5. Remaining unknowns that decide A: block-state stopband depth with the real MADP-000907 model (ON 2-3.3 ohm, OFF 28 fF), rows
   needed for >= 17 dB, insertion loss of the pass state, tolerance +-20 %.

## Recommendation (for approval)
Architecture **A**, **2 quadrants first** (19.4 dBi, 68 diodes with 3-row gates, 36 with 2 rows, 100 with 4 rows); 4 quadrants as a
second step only if 13.5 dBi is acceptable. Required contrast per gate: >= 17 dB isolation (design for 20 dB), IL <= 1 dB in the
pass state, block state either lossy or ring radius tuned so |S11| stays <= -10 dB. Sparse fraction: 3 gate rows (about 5 % of the
1420 cells). B/C only if mixed beams per quadrant (several states per quadrant) are needed later.
