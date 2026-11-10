# JOB 1 B: small-board HFSS variants (MODEL_ONLY / NOT_HFSS_VALIDATED) - one page
Six STEPs, same RF-only EXPORT_STANDARD as `BOARD_N240_RF_9cell.step` (9 switched cells + the static cells within 17 mm, hexagon, pin, port annulus, 9 diode sheets). Names inside every STEP: METAL_PEC (1), SUBSTRATE, SUBSTRATE_FR4, AIR, PORT, D<string>_<pos> x9. Gates for all six: names survive re-import, BRepCheck VALID, NO_OVERLAP, SHEETS_OK, PORT_OK. V0 has the same geometry as the delivered 9cell STEP (entity-id order of the colour items differs).

| file | the ONE change vs V0 |
|---|---|
| BOARD_V0_baseline.step | none (FR4 0.8 mm, anti-pad dia 1.0 mm, vias at today's spots: island (-1.780, -0.530), disc (-3.630, -3.630) cell-local mm) |
| BOARD_V1_fr4_0p2.step | FR4 0.8 -> 0.2 mm (barrels and via pads follow) |
| BOARD_V2_antipad_0p6.step | anti-pads 1.0 -> 0.6 mm dia (barrel r 0.15 stays) |
| BOARD_V3_vias_lowfield.step | both vias moved to the best low-field spot of JOB 1 A: island (-1.692, +0.121), disc (-4.471, +2.417) cell-local mm |
| BOARD_V4_shielded.step | + ring of 6 ground vias (r 0.15) at radius 0.8 mm around every bias via, FR4 only, from the ground plate down to the FR4 bottom (blind, no pads, fused into the ground lump) |
| BOARD_V5_novias.step | no barrels, no pads, no anti-pads (reference); discs/islands float |

## Manual setup (same as HFSS_SETUP_RF.md, 5.8 GHz)
1. Import each STEP, units mm. Materials: METAL_PEC = pec; SUBSTRATE = RT5880 (eps 2.2, tan 0.0009); SUBSTRATE_FR4 = FR4 (eps 4.4, tan 0.02; use the value of the full-board run); AIR = vacuum (radiation boundary on the 6 faces of AIR, lambda/4 = 12.9 mm from the board is inside the 15 mm box).
2. PORT: wave port or lumped port on the PORT annulus (coax, inner = pin r 0.635, outer r 2.05), 50 ohm, same as the full-board run.
3. ONE Impedance boundary on all 9 D* sheets, variables Rs, Xs (ohm/square, sheet 1.00 x 0.35 mm: ohm/sq = Z x 0.35/1.00). New diode rules, 5.8 GHz:
   - OFF: C 0.11 pF (0.105-0.12), Rp 4 k (3-5 k), series L 0.6 nH: Z = 15.497 - j226.627 ohm -> **Rs 5.424, Xs -79.320 ohm/sq**. Formula: Z = 1/(1/Rp + j 2 pi f C) + j 2 pi f L; sheet impedance = Z x 0.35/1.00. Corrected 2026-11-10 by decision session: former Xs -86.97 omitted the series inductance (+j7.653 ohm/sq). Recompute sensitivity values from the full formula. Prior HFSS runs retain their actual settings; this correction does not relabel them.
   - ON: 3 ohm + 0.6 nH: Z = 3 + j21.87 -> **Rs 1.05, Xs +7.65**.
   For a frequency sweep use the full complex expressions above: OFF Zs(f) = 0.35 [1/(1/Rp + j 2 pi f C) + j 2 pi f 0.6e-9], ON Zs(f) = 0.35 [3 + j 2 pi f 0.6e-9]. Separate real/imaginary parts for Rs/Xs; do not hold the 5.8 GHz numbers constant over a sweep.
4. Setup: 5.8 GHz, max delta S 0.02, then discrete sweep 5.0-6.6 GHz (resonance check, step 0.1), infinite sphere theta 0-180 step 2, phi -180..180 step 2 (Global CS). Run OFF and ON for each of the 6 files.

## What to send back (per file, OFF and ON)
1. S11 (dB) vs f 5.0-6.6 GHz (complex touchstone or csv) and the S11 at 5.8 GHz.
2. Power split at 5.8 GHz, as % of accepted power: radiated into theta 0-90 (up) and into theta 90-180 (down) (Radiated Power from the infinite sphere over the two hemispheres), plus dielectric loss in RT5880 and in FR4 and the ohmic loss in the Impedance sheets (Fields Calculator or the HFSS loss report), so that up + down + losses = 100 %.
3. Power fraction radiated at theta 70-90 (the grazing band), and the fraction at theta 0-60 (the lobe region), same normalisation.
4. Peak RHCP and LHCP directivity (dBi) and their (theta, phi), plus the exported .ffd of OFF and ON (the same name scheme as before) for the correlation tests.
Read-across: V1 vs V0 and V5 vs V0 separate the FR4 drain; V2/V3/V4 vs V0 tell which via measure reduces the downward power and the grazing power.
