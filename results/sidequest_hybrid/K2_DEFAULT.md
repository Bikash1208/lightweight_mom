# K2-default 2-beam board (one diode group, one control line, no inverter) - MODEL_ONLY (2026-11-05)

**Requirement:** default state (all diodes OFF, unpowered) = beam 1 (theta 25, phi 0); circuit ON (all diodes ON) = beam 2 (25, 180) with beam 1 suppressed; linear polarisation; exactly 10.5 GHz; 150 mm hexagon circumradius; RT5880 1.575 mm; 6.4 mm lattice; HFSS-calibrated centre probe; analog static cells. The earlier two-group board `board_K2_LP` (GA/GB, two control states) does not meet this and is **kept in the repo but not used**.
**Nothing here is an HFSS result.** The switched cell is HFSS-validated (C/C 10.670, R/R 11.006 GHz, dY -3.13 mS, alpha 0.844 Np/m, beta/k0 1.0731, back-out at 10.54 GHz, not re-derived at 10.5); the analog static cells, the board and family B are model only.

## Diode model (identical in the model and HFSS)
OFF = open slit, no sheet and no load (HFSS: C/C 10.670 against the MoM open slit 10.662; the lumped 28 fF is ineffective in HFSS, the MoM with a 28 fF load would give 10.525). ON = a resistance only: **HFSS Impedance boundary Rs = R_on x width / gap = 2.44 x 1.0 / 0.4 = 6.1 ohm/sq** (checked on the sheet size of the re-imported STEP); in the MoM a 2.44 ohm lumped load (the measured MA4AGP907 ON impedance 2.44 - j3.84 ohm is not used). Pipeline check on the validated cell: model open 10.662 / ON 11.007 GHz against HFSS 10.670 / 11.006.

## Cell library (`k2_default_lib.json`, MoM sheet tensor, normal incidence, mesh 64, 10.5 GHz, scalar cell = both diodes of a cell on the line)
Normalised m: open slit (OFF) = +1, ON (R 2.44) = -1 (Y 8.69 -> 6.06 mS; swing 2.63 mS in the MoM, 3.13 mS in HFSS).
| static analog cell (no diode) | solid disc | slit half 1.0 | 1.8 | 2.2 | 2.6 | 2.8 (= OFF geometry) | strap 2.2 | strap 2.8 |
|---|---|---|---|---|---|---|---|---|
| m | -1.67 | -1.58 | -1.22 | -0.81 | +0.05 | +1.00 | -1.39 | -1.01 |
Family A (rim cells): ON lowers Y. **Family B (through-slit, 4 islands, diode bridges them): ON raises Y**, but OFF sits at m = -4.2 (Y 1.86 mS) and ON at -1.8 / -1.6 / -0.6 (diode at 0.8 / 1.4 / 2.0 mm): a swing of +2.4 .. +3.6 in m, far below the useful static range. Static levels and pairs use normal-incidence extraction (flagged, as in all rounds).

## Designs (engine `sidequest_coded_k4`, lam 0.28, g 2.0 as in the old design; LP-V and eps_B chosen by the best min-state gain over 2 polarisations x 8 phases)
**The plain per-cell least squares (suppression weight mu = 0) gives only -4.7 dB isolation** (both beams almost equal in both states: the static compromise cells radiate the same field in both states). An explicit wrong-beam penalty (mu = 30, `design_null`: greedy coordinate descent with the exact change of the leaked amplitude of the other beam) is used for every new design.
| design | diodes (switched cells) | gain state 1 / 2 [dBi] | SLL [dB] | ISOLATION: beam 2 in state 1 / beam 1 in state 2 [dBi (dB rel. own peak)] | S11 |
|---|---|---|---|---|---|
| **(i)** family A only, one group, mu 30 | 514 (257 of 1338) | 11.19 / 11.44 | -11.7 / -11.8 | -4.1 (-15.3) / -4.4 (-15.9) | NA (bare probe -15 dB) |
| (i) at mu 0 (no suppression term) | 254 | 13.57 / 13.60 | -4.7 / -4.8 | +8.9 (-4.7) / +8.8 (-4.8) | NA |
| **(ii)** + family-B cells as MoM'd (through-slit), mu 30 | 580 (290) | 12.76 / 12.00 | -9.3 / -12.5 | -4.5 (-17.3) / -4.3 (-16.3) | NA |
| (ii) bound: hypothetical family B (ON raises m by +2, OFF anywhere in the static range) | 864 (432) | 15.09 / 14.78 | -13.8 / -14.4 | -8.4 (-23.5) / -8.9 (-23.7) | NA |
| old two-group `K2_LP`, as built (+-1 statics, 10.54 GHz, no suppression term) | 720 (360) | 15.15 / 15.15 | -13.5 / -13.5 | +1.7 (-13.5) / +1.7 (-13.5) | NA (-15 dB bare) |
| two-group topology at the new settings (analog statics, 10.5 GHz, mu 30) | 824 (412) | 14.78 / 14.78 | -14.2 / -14.0 | -7.7 (-22.5) / -7.7 (-22.5) | NA |
| (i) on the 120 mm hexagon | 336 (168 of 846) | 8.34 / 8.38 | -11.2 / -11.1 | -7.1 (-15.4) / -8.4 (-16.8) | NA |
Peak directions (25, 0) / (25, 180) in every row. XPD (i): 28.4 / 27.0 dB; aperture efficiency (i) 1.5 %. Suppression-weight trade-off of (i) (mu 0 / 10 / 30 / 100 / 300): gain 13.6 / 12.0-12.3 / 11.2-11.4 / 10.8-11.1 / 10.5-10.8 dBi, isolation -4.7 / -12.7 / -15.3 / -17.2 / -17.8 dB (`k2_default_designs.json`). **S11 is not modelled with the hologram: the only number is the bare calibrated probe, -15.0 dB (Phase 26 HFSS).**

## Findings
* One group and no inverter costs 4 dB gain against the two-group board (11.2 against 15.2 dBi) at about 30 % fewer diodes (514 against 720) **if the beam-1 suppression is enforced** (-15 dB isolation); the old board reaches only -13.5 dB isolation without any suppression term and needs two control states.
* Only the cells that go OFF -> ON are switchable (257 of 1338 = 19 %); the rest are static compromise cells. The static cells make both beams appear in both states unless the design pushes them against it, which is what costs gain.
* Family B recovers some of it only if its OFF level can be placed in the useful range: with the through-slit cells as they exist the gain is +1.6 / +0.6 dB (580 diodes); an idealised family-B cell (OFF anywhere, ON +2) would give 15 dBi with -23.5 dB isolation. The HFSS cell check decides whether the real B cell exists (`hfss_exports/sidequest_cellB_through_slit/`): MoM prediction OFF 11.348 -> ON 10.866 GHz (shift -0.48 GHz, ON below OFF as required for Y rising).

## Cost and fallback (`hfss_exports/sidequest_board_K2_default*/hfss_cost.json`; ESTIMATE from scaling only; the Phase 26 150 mm board's tets / RAM / time are NOT_AVAILABLE_FROM_CURRENT_DATA)
150 mm: 8-23 M tets, **82-342 GB (upper end exceeds 256 GB)**, 18-42 h state 1 on 32-64 cores, 6-14 h state 2 with the mesh reused. **120 mm fallback (same concept): 846 cells, 168 switched, 336 diodes, 5-14 M tets, 52-217 GB (fits 256 GB), 10-19 h, 8.34 / 8.38 dBi, same isolation (-15.4 / -16.8 dB).**

## Caveats
* Frequency exactly 10.5 GHz in the engine; the cell constants alpha / beta (HFSS back-out) are from 10.54 GHz, not re-derived. Static cells carry no diode: the same analog levels in both states.
* Discs of the outermost cells are kept 0.2 mm inside the outline (cells with the disc inside the 150 / 120 mm hexagon: 1338 / 846 cells; the old board used 1430 cells by centre, outline 154 mm).
* The wrong-beam penalty is a design choice (mu 30), not an HFSS-validated parameter; gain and isolation trade (table).
* Polarisation LP-V (phi-hat, scalar cell J = m rho-hat E); CP waits for the 4-diode cell's HFSS check (`docs/briefs/SIDEQUEST_4QUADRANT.md`).
* lam 0.28 and g 2.0 are the old design's operating point, not re-optimised (a lower lam would admit more switched cells; not run).

## Lambda scan of design (i) (2026-11-05, 150 mm, LP-V, eps_B 180, mu 30, g 2.0; `k2_default_lambda_scan.json`)
| lam | diodes (switched) | gain [dBi] | SLL [dB] | isolation [dB] |
|---|---|---|---|---|
| 0.28 (current package `board_K2_default`) | 514 (257) | 11.19 / 11.44 | -11.7 / -11.8 | -15.3 / -15.9 |
| 0.20 | 608 (304) | 11.69 / 11.95 | -12.0 / -12.4 | -16.1 / -16.6 |
| 0.14 | 686 (343) | 11.99 / 12.07 | -12.5 / -11.8 | -17.2 / -17.9 |
| 0.10 | 746 (373) | 12.19 / 12.21 | -12.6 / -11.8 | -17.7 / -18.1 |
| 0.05 | 812 (406) | 12.27 / 12.29 | -12.2 / -12.6 | -17.7 / -18.6 |
Gain rises 1.08 dB (min state) at lam 0.05 with isolation <= -17.7 dB: the rule (>= 1 dB at <= -15 dB) is met, so a second package `hfss_exports/sidequest_board_K2_default_lam005/` (812 diodes, all checks PASS, 8.6-24 M tets, 86-359 GB) was built NEXT TO the current one; `board_K2_default` was not touched (it is the one being run in HFSS). Family B stays future work (through-slit cell weak; a better ON-raises-Y topology to be proposed later).
