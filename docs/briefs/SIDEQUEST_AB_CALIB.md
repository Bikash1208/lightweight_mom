# SIDEQUEST_AB_CALIB - A/B board: package and 4-state calibration of the board model (2026-11-06). MODEL_ONLY. Step 1 of the A/B track; no design predictions yet.
## 1. Package
`hfss_exports/sidequest_board_k2_lp/` (STEP `board_K2_LP.step`, `cells.csv`, `sheets.csv`, `states.csv`, `HFSS_SETUP.md`; copy in `C:\hmsa_sync\data\sidequest\board_K2_LP`). Stack: RT5880 (eps_r 2.2, tan d 0.0009) 1.575 mm, ground 35 um, patches 35 um, probe r 0.635 mm. Outline: hexagon, circumradius 150 mm (substrate corners (0, +-154), (+-133.4, +-77)). Lattice 6.4 mm, 1430 cells (no cell for rho < 6.5 mm). Cell: disc r 3.0 mm with a cross slit 0.4 x 5.6 mm (rim 0.2). Design frequency 10.54 GHz (HFSS data at 10.5 GHz).
| cell kind | count | meaning |
|---|---|---|
| static open (+1) | 526 | cross slit open, no diode footprint |
| static strap (-1) | 544 | metal in both footprints (permanent ON, R -> 0) |
| switched GA | 180 cells = 360 sheets `GA_<ix>_<iy>_DX|DY` | 2 diodes per cell (0.4 x 1.0 mm sheets at +0.8 mm), ON in state A |
| switched GB | 180 cells = 360 sheets `GB_...` | ON in state B |
Total 360 switched cells, 720 sheets / diodes. Original design (model): both states 15.15 dBi at (25, 0) / (25, 180), SLL -13.5 dB, 77 % of the SW power reaching the edge, launch 0.627 (assumed from the Phase 26 board).
## 2. Calibration set (user HFSS, 10.5 GHz, LP)
| (GA, GB) | gain at (25, 0) | gain at (25, 180) |
|---|---|---|
| ON, OFF | 13.5 | 1.2 |
| ON, ON | 8.9 | 8.9 |
| OFF, ON | 1.4 | 13.3 |
| OFF, OFF | 10.4 | 8.5 |
Also SLL about -5 dB; about 12 % of the power within theta < 30 deg; about 70 % to the horizon / back.
## 3. Model and inputs
SWHB (`sidequest_swhb.py`), cell data = the HFSS-validated K2 library `k2_default_lib.json` (MoM mesh, family A, slit half-length 2.8, raw, NO contrast factor): open 8.69 mS, ON 6.06 mS (+0.06 mS real), strap 6.05 mS (imag parts, Y_xx = Y_yy); OFF = open (the 28.16 fF across the 0.4 mm gap is not added). Env RT5880 1.575 mm, 10.5 GHz. Script `sidequest_ab_calib.py` (`raw`, `grid`), data `results/sidequest_hybrid/ab_calib_grid.json`.
## 4. Result: parameter-free cell, ONE fitted number (launch efficiency)
Raw model, launch efficiency 1: beams 15.9 dBi (ON/OFF, OFF/ON), 11.4 (ON/ON), 11.7 (OFF/OFF); wrong beam 2.5; peak directions (25, 0) / (25, 180) exactly; beta/k0 1.056-1.063; radiated 0.11-0.13, edge power 0.85-0.89 of the SW power.
| state | model G(25,0) / G(25,180) after eta = -2.12 dB | HFSS | error [dB] |
|---|---|---|---|
| ON, OFF | 13.8 / 0.3 | 13.5 / 1.2 | +0.3 / -0.9 |
| ON, ON | 9.3 / 9.3 | 8.9 / 8.9 | +0.4 / +0.4 |
| OFF, ON | 0.3 / 13.8 | 1.4 / 13.3 | -1.1 / +0.5 |
| OFF, OFF | 9.6 / 9.6 | 10.4 / 8.5 | -0.8 / +1.1 |
* **The fitted eta = -2.12 dB (0.61) equals the independent launch value 0.627 (-2.0 dB) of the Phase 26 board**; with eta fixed at 0.627 the errors are +0.4 / -0.8, +0.5 / +0.5, -1.0 / +0.6, -0.7 / +1.2 dB.
* Gain: 6 of 8 numbers within +-1 dB, the other two (a wrong-beam level and the OFF/OFF asymmetry) at 1.1-1.2 dB: **PASS within +-1.2 dB, marginal on the strict +-1 dB** (the model cannot make the 1.9 dB OFF/OFF left-right asymmetry; the model is symmetric there).
* Directions: peaks at (25, 0) and (25, 180) in the states with a beam (grid 2 / 5 deg): within +-2 deg. PASS.
* Isolation: model 13.4 dB (A) / 13.4 (B) against HFSS 12.3 / 11.9: errors +1.2 / +1.6 dB. PASS (+-2).
* Sensitivity: the residuals barely move over contrast factor 0.8-1.1 and common susceptance shift -1..+1 mS (max residual 1.02-1.58 dB; best +0.5 mS: 1.02), so no cell parameter is needed or justified.
* Not matched: SLL (model -7.5 dB in the beam states against about -5 in HFSS: model 2.5 dB optimistic; the model has no rim diffraction / mutual coupling); power fractions: within theta < 30 deg model 8.6 % of the total radiated power (beam + rim) against about 12 %; rim / horizon+back: model 87 % (all of the edge power) against about 70 %, i.e. roughly 80 % of the edge power radiates and the rest is reflected or absorbed.
## 5. Verdict
**Calibration PASSES (gain +-1.2 dB, directions +-2 deg, isolation +-1.6 dB) with the single launch factor 0.627, no cell fit.** Caveats: SLL 2.5 dB optimistic; OFF/OFF asymmetry not reproduced; the rim is modelled as lost (edge power), which is consistent with the HFSS 70 % to the horizon / back. Iterations (a)-(e) are NOT started (instructed stop after step 1). Open input that would tighten eta: the HFSS S11 of this board.

## 6. Addendum (2026-11-06): GA = PEC data
HFSS, GA = PEC: GB OFF -> 13.2 dBi at (26, 2), 2.3 dBi at (25, 180); GB ON -> 8.8 / 8.8 dBi. Model (GA strap, same eta -2.12 dB): GB OFF 13.9 / 0.3 (errors +0.7 / -2.0 dB; direction (25, 0) vs (26, 2): within 2 deg), GB ON 9.4 / 9.3 (+0.6 / +0.5). Beam-2 contrast (GB ON level minus GB OFF level at (25, 180)): HFSS 6.5 dB, model 9.0 dB (**the model is 2.5 dB optimistic on the contrast; the wrong-beam / diffuse floor is under-estimated by about 2 dB**). Main beams within 0.7 dB. The 6-point set (4 states + 2 PEC states) stays within +-1.2 dB on beams, +-2 dB on the floor level.
