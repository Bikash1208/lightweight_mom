# SUB6 board model with real cell tables (stage: sub6-board-model)

**MODEL_ONLY / PENDING_HFSS_CELL_CHECK / NOT_HFSS_VALIDATED.** Every number is a model number for cells whose tensors come from `results/sub6_slotdisc/rows_RT5880.jsonl` (mesh 64, MoM); no cell has been checked in HFSS yet. Job `jobs/sub6-board-model.md`; branch `compute/sub6-board-model`; data `points.csv`, `layout_<side>mm_N<n>_<a|b>.csv` (one per point: site, x, y, STATIC/SWITCHED, row id, slot axis, alpha and beta OFF/ON, required alpha, ON-OFF phase), `provenance_run.json`; code `scripts/research/sub6_board_model.py` (+ `sub6_board_req.py` copied unmodified from compute/sub6-board-requirements 7694971).

## Setup
5.8 GHz, RT5880 3.175 mm on PEC, p = 11.6 mm, hexagon side 200 mm (758 cells) and 150 mm (442 cells), centre probe, **launch fraction 0.6263** (`src/dispersion/pin_launch.py` at 5.8 GHz; the Phase 26 pin scaled geometrically by the substrate-height ratio 3.175/1.575, as in stage 1b), RHCP beams A (25, 0) and B (25, 180). Phase 26 `march_rays` / `cell_design` unmodified; no frozen code edited.

**Architecture.** OFF: static cells radiate A+B (one linear current, amplitude via the cell alpha, axis psi). ON: static cells unchanged, switched cells radiate -B (near-silent in OFF: alpha_OFF <= 0.1 Np/m). Group-normalised weights: static (uA+uB)/(1-f), switched -uB/f, f = N_sw/N.
**Cells.** Static = NOCUT rows, switched = notch rows (OFF_BAR50 / ON_BAR50 with their own alpha and beta). Of 80 rows 16 are resonance-flagged (proxy) and excluded; the static filter beta_NOCUT <= 1.2 leaves **64 static rows**; the switched filter alpha_OFF <= 0.1, alpha_ON >= 0.3, beta_ON <= 1.2 leaves **54 switched rows** (the filters are mine). Each site takes the row nearest in sqrt(alpha) to the required alpha. Rows differ in gap and W, so neighbouring cells differ in disc size; coupling between different cells is not modelled.
**Real beta.** Per state, along every ray, the cell under the ray gives its own alpha and beta (static beta_NOCUT 1.05-1.2; switched OFF 1.04-1.06, ON 1.06-1.19). The slot axes are designed with the OFF phase; ON keeps the axes and takes the ON wave phase.
**Placement.** (a) low-discrepancy interleave; (b) the ideal-cell projection optimiser of stage 1b (rho = -1000, the limit of this architecture), then evaluated with real cells. **(b) returned the same result as (a) at 6 of 8 points**: the optimiser found no improvement there.
**Sweep.** N_sw 150/200/250/300 at 200 mm (f 0.19-0.39); 87/117/146/175 at 150 mm (same fractions). 9 s for 16 points (the 16 tasks limit the pool to 16 workers; the machine has 64).

## Gates (targets: balance <= 3 dB, AR <= 3 dB (OFF A, B and ON A), B suppression >= 12 dB, SLL <= -13 dB, both states)
**MEETS_ALL: 0 of 16 points.** Balance and suppression are met everywhere. The misses are AR in ON (beam A built from the static lattice with the switched sites removed) and SLL in ON. Best margin -1.15 dB (200 mm, N 200, placement a, MISSES_AR).

| side mm | N_sw (f) | pl | OFF G_A | OFF G_B | balance | OFF AR A/B | ON G_A | B supp | ON AR_A | SLL off/on | gate (margin) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 200 | 150 (0.20) | a | 11.0 | 10.7 | 0.30 | 2.61/2.19 | 15.2 | 23.6 | 5.50 | -14.6/-12.5 | MISSES_AR+SLL (-2.50) |
| 200 | 200 (0.26) | a | 13.3 | 13.0 | 0.32 | 2.38/1.80 | 17.0 | 25.4 | 4.15 | -15.0/-13.1 | MISSES_AR (-1.15) |
| 200 | 250 (0.33) | a | 15.3 | 15.0 | 0.32 | 2.02/1.48 | 18.4 | 25.8 | 2.80 | -16.3/-11.6 | MISSES_SLL (-1.39) |
| 200 | 300 (0.40) | a | 15.7 | 14.0 | 1.76 | 2.11/1.62 | 18.0 | 22.9 | 2.64 | -15.8/-11.7 | MISSES_SLL (-1.29) |
| 200 | 150 (0.20) | b | 11.0 | 10.7 | 0.30 | 2.59/2.00 | 15.4 | 25.0 | 5.18 | -14.6/-12.6 | MISSES_AR+SLL (-2.18) |
| 200 | 200 (0.26) | b | 13.3 | 13.0 | 0.32 | 2.38/1.80 | 17.0 | 25.4 | 4.15 | -15.0/-13.1 | MISSES_AR (-1.15) |
| 200 | 250 (0.33) | b | 15.3 | 14.9 | 0.32 | 2.18/1.58 | 18.2 | 24.3 | 3.12 | -15.9/-11.3 | MISSES_AR+SLL (-1.73) |
| 200 | 300 (0.40) | b | 15.7 | 14.0 | 1.76 | 2.11/1.62 | 18.0 | 22.9 | 2.64 | -15.8/-11.7 | MISSES_SLL (-1.29) |
| 150 | 87 (0.20) | a | 7.8 | 7.5 | 0.26 | 2.87/2.79 | 11.7 | 21.1 | 6.14 | -15.0/-9.4 | MISSES_AR+SLL (-3.58) |
| 150 | 117 (0.26) | a | 10.3 | 9.9 | 0.39 | 2.31/2.52 | 13.9 | 25.8 | 3.44 | -16.0/-10.5 | MISSES_AR+SLL (-2.54) |
| 150 | 146 (0.33) | a | 12.4 | 12.3 | 0.11 | 2.27/2.61 | 15.6 | 24.2 | 2.77 | -16.0/-9.8 | MISSES_SLL (-3.18) |
| 150 | 175 (0.40) | a | 12.3 | 11.8 | 0.43 | 1.69/2.20 | 14.8 | 25.5 | 1.50 | -14.5/-8.7 | MISSES_SLL (-4.29) |
| 150 | 87 (0.20) | b | 7.8 | 7.5 | 0.26 | 2.87/2.79 | 11.7 | 21.1 | 6.14 | -15.0/-9.4 | MISSES_AR+SLL (-3.58) |
| 150 | 117 (0.26) | b | 10.3 | 9.9 | 0.39 | 2.31/2.52 | 13.9 | 25.8 | 3.44 | -16.0/-10.5 | MISSES_AR+SLL (-2.54) |
| 150 | 146 (0.33) | b | 12.4 | 12.3 | 0.11 | 2.27/2.61 | 15.6 | 24.2 | 2.77 | -16.0/-9.8 | MISSES_SLL (-3.18) |
| 150 | 175 (0.40) | b | 12.3 | 11.8 | 0.43 | 1.69/2.20 | 14.8 | 25.5 | 1.50 | -14.5/-8.7 | MISSES_SLL (-4.29) |

Gain is referred to the accepted power x the model conversion eta, so the absolute G_A are model numbers; ON G_A is higher than OFF because the switched cells add conversion in ON. B suppression is the coherent cancellation of discrete real cells (switched amplitude error rms <= 0.4 %, the table is dense) and is the most optimistic number here.

## Beam-A phase error from the switched cells' beta change (ON - OFF)
| side mm | N_sw | pl | phase error ON-OFF at static sites: rms / max / mean deg | switched d_beta (ON-OFF) | beam-A gain cost (dB) | ON AR_A with / without | ON SLL with / without |
|---|---|---|---|---|---|---|---|
| 200 | 150 | a | 8.8 / 37.8 / 6.5 | 0.064 | 0.04 | 5.50 / 5.39 | -12.5 / -12.6 |
| 200 | 200 | a | 11.4 / 41.2 / 8.9 | 0.064 | 0.08 | 4.15 / 4.11 | -13.1 / -13.3 |
| 200 | 250 | a | 14.1 / 49.9 / 11.4 | 0.064 | 0.08 | 2.80 / 2.86 | -11.6 / -11.3 |
| 200 | 300 | a | 16.0 / 56.6 / 12.8 | 0.065 | 0.14 | 2.64 / 2.67 | -11.7 / -11.8 |
| 200 | 150 | b | 9.0 / 33.6 / 6.7 | 0.064 | 0.06 | 5.18 / 5.08 | -12.6 / -11.6 |
| 200 | 200 | b | 11.4 / 41.2 / 8.9 | 0.064 | 0.08 | 4.15 / 4.11 | -13.1 / -13.3 |
| 200 | 250 | b | 13.8 / 50.9 / 10.9 | 0.064 | 0.07 | 3.12 / 3.12 | -11.3 / -11.0 |
| 200 | 300 | b | 16.0 / 56.6 / 12.8 | 0.065 | 0.14 | 2.64 / 2.67 | -11.7 / -11.8 |
| 150 | 87 | a | 7.0 / 23.2 / 5.0 | 0.064 | 0.01 | 6.14 / 5.96 | -9.4 / -9.7 |
| 150 | 117 | a | 8.4 / 26.0 / 6.3 | 0.057 | 0.02 | 3.44 / 3.36 | -10.5 / -10.6 |
| 150 | 146 | a | 13.6 / 39.7 / 10.4 | 0.073 | 0.07 | 2.77 / 2.78 | -9.8 / -9.9 |
| 150 | 175 | a | 12.7 / 36.2 / 9.9 | 0.060 | 0.08 | 1.50 / 1.50 | -8.7 / -9.0 |
| 150 | 87 | b | 7.0 / 23.2 / 5.0 | 0.064 | 0.01 | 6.14 / 5.96 | -9.4 / -9.7 |
| 150 | 117 | b | 8.4 / 26.0 / 6.3 | 0.057 | 0.02 | 3.44 / 3.36 | -10.5 / -10.6 |
| 150 | 146 | b | 13.6 / 39.7 / 10.4 | 0.073 | 0.07 | 2.77 / 2.78 | -9.8 / -9.9 |
| 150 | 175 | b | 12.7 / 36.2 / 9.9 | 0.060 | 0.08 | 1.50 / 1.50 | -8.7 / -9.0 |

The switched cells' beta rises by about 0.06 in ON (mean over the chosen rows). Downstream static cells see a wave phase shifted by up to 57 deg (rms 7-16 deg); the cost is **<= 0.14 dB of beam-A gain and <= 0.3 dB in AR / SLL** (columns "without" = the ON amplitudes with the OFF phase). The beta change is therefore not what makes the ON AR / SLL miss: the miss is present without it.

## Design bookkeeping (placement a)
| side mm | N_sw | alpha0 mean (Np/m) | cap1 | q90 demand | eta OFF / ON | alpha-saturated static / switched sites | switched amplitude error rms % |
|---|---|---|---|---|---|---|---|
| 200 | 150 | 0.33 | 0.35 | 10.4 | 0.11 / 0.28 | 0 / 0 | 2.50 |
| 200 | 200 | 0.61 | 0.68 | 5.4 | 0.18 / 0.39 | 0 / 0 | 0.06 |
| 200 | 250 | 1.07 | 1.18 | 3.1 | 0.30 / 0.51 | 0 / 0 | 2.36 |
| 200 | 300 | 1.07 | 1.18 | 3.1 | 0.30 / 0.47 | 0 / 0 | 0.36 |
| 150 | 87 | 0.33 | 0.36 | 10.3 | 0.08 / 0.21 | 0 / 0 | 0.45 |
| 150 | 117 | 0.61 | 0.70 | 5.2 | 0.14 / 0.30 | 0 / 0 | 0.72 |
| 150 | 146 | 1.07 | 1.20 | 3.1 | 0.23 / 0.41 | 0 / 0 | 0.25 |
| 150 | 175 | 1.07 | 1.19 | 3.1 | 0.23 / 0.37 | 0 / 0 | 0.79 |

alpha0 = mean alpha of the base march (engine, static rows with alpha <= cap1 = min(static, switched max alpha) / q90 of the demand, the alpha cap of stage 1b). The switched cells need 1/f of the amplitude, so the demand q90 is 3.1-10.4 and the base alpha (hence the conversion) is small at small N_sw. No site is alpha-saturated.

## Limits
Cells are mesh-64 MoM, not HFSS-checked (the HFSS eigenmode package is `hfss_exports/sub6_slotdisc_cell/`, 2 of the 80 rows). beta/alpha of a row are those of the infinite periodic array of that row; the ON tensor is the BAR50 model at 5.8 GHz. The A+B current of a static cell is treated as one ideal linear current whose amplitude is the chosen row's alpha-amplitude. The switched OFF polarisation is assumed equal to the ON direction. The placement (b) objective is the ideal-cell one. Stopped here: no cell choice, no board design, no HFSS package.
