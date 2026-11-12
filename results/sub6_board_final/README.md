# Final-board selection (MODEL_ONLY): sweep R = 110-260 mm x N = 60-200, k = 1

**Result: no (R, N) meets the SLL rule, so no winner and `BOARD_FINAL` was NOT built.** Everything below is MODEL_ONLY (jobs1_iter1_study.py engine, rot180, the study's own placement key and cap optimisation, k = 1 only: gathered elements k > 1 are not buildable with the current exporter, one diode per cell).

## Rule applied (S3, replaces step 3)
1. Polarisation reported at each wanted beam peak: RHCP, LHCP, XPD and axial ratio (target AR <= 3 dB, XPD >= 15 dB). The model does give polarisation (E_theta, E_phi of the aperture currents); it is ideal (AR about 1-2 dB) and HFSS on the old board was worse (see baseline), so AR is judged in HFSS.
2. SLL (S3b definition, `sll_def.py`): highest RHCP-or-LHCP directivity over theta 0-90 deg outside the exclusion cones around the actual RHCP peaks, cone radius max(1.5 x HPBW, 10 deg), HPBW = equivalent circular half-power width. OFF: cones around A and B. ON: cone around A only, so any B residual counts as a sidelobe. Levels relative to the strongest RHCP beam of that state. Required <= -15 dB HFSS-equivalent in both states.
3. OFF A/B balance, gain, ON suppression: reported only.
4. Tie-breaks: smallest N, then most converted power.

## Model-to-HFSS offsets (old no-FR4 board, same SLL definition, `sweep/baseline_sll.json`)
Frame: HFSS far fields rotated by 180 deg (`fields_hfss_rot180.npz`) against the model identity prediction (`fields.npz`); A at (25, 0), B at (25, 180).
| State | HFSS SLL | model SLL | offset (HFSS - model) | HFSS D peak | model D peak | directivity offset (model - HFSS) |
|---|---|---|---|---|---|---|
| OFF | -6.3 dB at (28, -36), RHCP | -16.5 dB at (10, -150) | **+10.2 dB** | 17.1 dBi (A), 13.1 dBi (B) | 22.8 / 21.8 dBi | **5.7 dB** (A) |
| ON | -7.5 dB at (22, -40), RHCP | -15.4 dB at (2, 10) | **+7.9 dB** | 18.0 dBi (A), 8.1 dBi (B) | 25.1 dBi (A) | 7.2 dB (A) |
Polarisation on that board (HFSS): OFF A AR 4.3 dB / XPD 12.3 dB, B AR 5.7 / XPD 10.0; ON A AR 2.8 / XPD 15.9. Model: AR 1.0-1.6 dB, XPD 21-25 dB. Applied: SLL_eq = model SLL + 10.2 (OFF) / + 7.9 (ON); gain_eq = model gain - 5.7 dB (reported only). One board only: treat the offsets as rough.

## Outcome
- SLL rule passes: **0 of 35**. Even the raw model SLL never reaches -15 dB in OFF (best -13.6 dB), so with the +10.2 dB offset the best HFSS-equivalent worst-case SLL is -3.4 dB.
- AR <= 3 dB at all wanted peaks in the model: 29 of 35 rows.
- 3 best by the rule (smallest worst-case HFSS-equivalent SLL; all R = 260 mm, placement inner, cap 0.5 Np/m):
| R, N | SLL OFF eq / ON eq (dB) | peak sidelobe OFF / ON (theta, phi) | AR A/B OFF, A ON | OFF balance A-B | gain A/B OFF (model / HFSS-eq) | B suppression ON | converted power OFF/ON | power reaching the edge OFF/ON |
|---|---|---|---|---|---|---|---|---|
| 260, 60 | -3.4 / -7.5 | (22, 206) / (14, 172) | 1.3 / 2.4, 1.5 | +13.0 dB | 18.5 / 5.5 dBi (12.8 / -0.2) | 8.8 dB | 0.253 / 0.181 | 0.747 / 0.819 |
| 260, 80 | -1.9 / -6.1 | (18, 204) / (14, 188) | 2.1 / 1.1, 1.7 | +10.7 dB | 17.9 / 7.2 dBi (12.2 / 1.5) | 8.9 dB | 0.263 / 0.176 | 0.737 / 0.824 |
| 260, 100 | -0.8 / -5.3 | (16, 198) / (14, 188) | 1.9 / 0.7, 1.7 | +9.2 dB | 17.5 / 8.3 dBi (11.8 / 2.6) | 9.3 dB | 0.272 / 0.172 | 0.728 / 0.828 |
"Power reaching the edge" = 1 - converted power (the unconverted remainder; includes ohmic loss).

## Static-cell check (SCREEN_C40B, new stack, `results/sub6_cell_param/compare_c40b/static_cells/`)
Five most-used static cell types (usage counted over all study caps at R = 180 mm); beta/k0 at 5.8 GHz from the lowest-Q surface-wave branch, along (phase parallel to the slot axis) and across, from only two phases per direction (bracketed, 80/88 deg along and 88/100 deg across).
| Cell (table row) | table beta/k0 | HFSS along | diff | HFSS across | diff |
|---|---|---|---|---|---|
| g0.8_q0.8_W1.6 | 1.0683 | 1.0616 | -0.007 | 1.0973 | +0.029 |
| g0.6_q0.8_W1.6 | 1.0743 | 1.0667 | -0.008 | 1.1099 | +0.036 |
| g0.4_q0.7_W3.6 | 1.1018 | 1.0161 | **-0.086** | 1.1365 | +0.035 |
| g0.4_q0.8_W2.4 | 1.1218 | 1.0638 | **-0.058** | 1.1684 | +0.047 |
| g0.6_q0.8_W2.4 | 1.0864 | 1.0668 | -0.020 | 1.1340 | +0.048 |
The table values (all marked `extrapolated` in the table) match "along" for three cells within 0.02 and not for the W 3.6 and W 2.4 / g 0.4 cells. Correction given to the study as an offset on the static-table beta of those two rows: -0.0857 (g0.4_q0.7_W3.6) and -0.0580 (g0.4_q0.8_W2.4), applied as `b_new = b + offset` to the static table before the march (`STATIC_CORR`, `sweep/static_corr.json`); the third borderline cell (-0.0196) was not corrected. The corrected and uncorrected sweeps agree for the best three rows (cap 0.5-0.7 Np/m does not use the corrected rows); both are in `sweep/`.

## Flags
- k = 1 only. The model has no direct-feed radiation, no FR4/bias-layer effect on the static cells (FR4-calibrated table), and a BAR50-02V switched-cell strength (alpha 2.16 Np/m, r_on 0.335), not the BAR90-02ELS.
- The SLL offsets come from one old board (HFSS vs model), +10.2 / +7.9 dB; they dominate the outcome.
- OFF balance is reported and not used for selection, per the rule; gain has no threshold.
- Files: `sweep_ranked_HFSS_equivalent.csv` (all 35 rows with HFSS-equivalent columns and ranks), `best3.json`, `sweep/full_s3b_corr.*` (primary), `sweep/full_s3b.*` (uncorrected), `sweep/baseline_sll.json`, `sweep/static_corr.json`.
