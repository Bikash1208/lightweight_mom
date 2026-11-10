> **SUPERSEDED IN PART by `DIAG_FOLLOWUP.md` (identity convention validated on the hex150 board; the 'rot180 + switched swap' reading of section 2 is not supported by the complex fits).**

# DIAGNOSTIC: HFSS full board (BOARD_N240_RF_full_win.step) vs board model. MODEL_ONLY / NOT_HFSS_VALIDATED. No design change.
Inputs: `2_34_ -56_3.ffd` (OFF, Rs 2.34 / Xs -56.3) and `s1.ffd` (ON, 1.05 / +7.65), both found in `C:\Users\Bikash\Desktop\` (NOT in `...\sub6_board_hfss\hfss_results\`, which does not exist; the OFF file name has a space, not a double underscore). Grid theta 0..180 step 2 (91), phi -180..180 step 2 (181), e^{+jwt}, RHCP = (Etheta + j Ephi)/sqrt2 = the project's own `copol_component(.., 'rhcp')`.
Scripts: `diag_step_geometry.py` (item 1), `diag_hfss_vs_model.py` (items 2, 3). Model = the design engine (`final_design('cal')`, mesh-96 cells, `currents` OFF/ON, layered-medium dyad) on the HFSS grid, theta <= 90.

## 1. Geometry: the STEP is the design (no sign, mirror, swap or offset)
All 758 discs, 518 static + 240 switched, measured in the STEP (top face z 0.035: outer circle centre, slot polygon) and matched to the design sites: position error max 3e-14 mm, slot long-axis angle error max 0.054 deg (mean 0.024 deg; estimator resolution), slot centre offset max 0.0002 mm, 240 islands: centre error max 0.0006 mm, axis error max 0.051 deg, 240 diode sheets: centre error max 0.00045 mm (all 240 matched uniquely, also the island side and the diode direction follow the string rotation). Hypotheses on the residuals (all 758 cells): identity max position error 3e-14 mm / axis 0.054 deg; x-mirror 325 mm / 90 deg; y-mirror 394 mm / 90 deg; rot180 394 mm (axes alone cannot tell, they are mod 180); x<->y 377 mm / 90 deg; axis+90 deg 90 deg. => The STEP holds the design exactly; a 12 sample table of 36 cells (12 switched + 24 static spread in azimuth) is in `step_geometry.json`.

## 2. Model vs HFSS (RHCP power pattern correlation over theta <= 60, all phi; beams = RHCP maxima at |phi| < 90 (A side) and >= 90 (B side), theta/phi in deg)
HFSS: OFF one beam at (24, 178) 14.06 dBi (RHCP 13.98 / LHCP 0.06 at (25, 180)); (25, 0): RHCP 3.1 dBi. ON: beam (26, -178) 13.13 dBi; (25, 0): RHCP 8.2 dBi, a lobe (30, 4) 9.4 dBi.
Model (identity): OFF beam A (24, 0) 22.8 dBi, B (24, 180) 21.8 dBi; ON A 25.2 dBi, B -11.6 dBi at (25, 180).
| variant | corr OFF | corr ON | beams (theta/phi) |
|---|---|---|---|
| identity (design as is) | 0.46 | 0.20 | OFF A 24/0 B 24/-180 ; ON A 24/0, B 6/90 (none) |
| x-mirror (x -> -x) | 0.03 | 0.23 | OFF 16/68, 26/-136 ; ON 18/-88 |
| y-mirror (y -> -y) | 0.13 | 0.08 | OFF 26/44, 16/-112 ; ON 18/88 |
| rot 180 | 0.54 | 0.64 | OFF A 24/0 B 24/-180 ; ON only 24/-180 |
| conj(E) (time-convention flip of the field) | 0.05 | 0.16 | scattered |
| conj(cell phase) | 0.25 | 0.23 | scattered |
| x<->y swap | 0.24 | 0.36 | scattered |
| switched cells OFF/ON swapped | -0.01 | 0.56 | OFF A 24/0 ; ON A 24/0 and B 24/-180 |
| both states swapped | -0.02 | 0.57 | same |
| **rot 180 + switched-cell OFF/ON swapped** | **0.73** | **0.64** | **OFF only 24/-180 ; ON 24/0 and 24/-180** |
| mirror x + conj(E) + switched-cell swapped | 0.73 | 0.64 | same as the line above |
| rot 180 + both states swapped | 0.73 | 0.63 | same |
**Variant that matches HFSS (best, equivalence class): beam A mirrored to phi = 180 in both states AND the switched-cell response inverted (the switched B contribution appears in ON, not in OFF).** rot 180, "mirror x + conj(E)" and "identity + conj(E) + conj(cell phase)" give the same pattern (one equivalence class, the data cannot separate them); a pure mirror without the conj is rejected by the polarisation: HFSS is RHCP at (25, 180) with XPD 14 dB, a mirror alone would give LHCP there. Pure x/y-mirror, pure conj and x<->y swap are rejected (corr <= 0.36). Quantitatively the evidence for the hypothesis: HFSS (25, 0) RHCP rises by +5.1 dB from OFF to ON (3.1 -> 8.2) while the model (rotated by 180) has the switched beam there strong in OFF and absent in ON (21.7 -> -11.6 dBi); the beam at phi 180 stays (14.0 -> 13.0, -0.9 dB) as for static cells. **Limit: the best correlation is only 0.73 / 0.64 (complex field correlation 0.1 to 0.2): the sign / mirror error explains the beam positions and the state dependence, NOT the HFSS pattern itself, which is diffuse.** J-direction mirrors / 90-degree rotations of the cell currents change nothing (identical correlations).

## 3. Directivity gap (report only)
HFSS D(RHCP) 14.1 / 13.1 dBi (OFF / ON) vs model 22.8 / 25.2 dBi: -8.7 / -12 dB. Beam widths are similar (theta HPBW 14 / 16 deg vs model 12 deg), but the power is spread: within 10 deg of the main beam HFSS has 12 / 8 % of the power (model 43 / 75 %); SLL -2.1 / -2.2 dB (secondary lobes as strong as the main beam: e.g. OFF lobes (34, -168), (24, 148), (18, -76)); power by theta band, HFSS OFF/ON: 0-10 deg 2 / 1 %, 10-40 deg 39 / 44 % (model 91 %), 40-70 deg 35 / 36 % (model 7 %), 70-90 deg 14 / 10 % (model < 1 %), below the horizon (theta > 90) 10 / 9 % (model 0). Model sensitivities (same engine, one change at a time): scaling all alpha by 0.3 ... 3 changes D by < 0.5 dB (eta only: 0.11 ... 0.69 in OFF); random cell phase errors sigma 30 / 60 / 90 deg give D 22.4 / 20.3 / 16.5 dBi (OFF), 25.0 / 23.3 / 18.8 (ON): D 14 needs about 100+ deg rms phase error, i.e. an almost random aperture phase.
Likely causes, ranked by the evidence: (1) the cell phase / sign convention of the model does not match the physical cells (beam on the wrong side, switched response inverted: demonstrated above) and the tensor-cell phase of the static rows may not be constant (I did not check whether the model carries a row-to-row phase offset): this gives an incoherent aperture, diffuse power, SLL -2 dB; (2) strong edge / grazing radiation (35 % of the power at 40-70 deg, 10-14 % at 70-90 deg, 9-10 % below the horizon): the guided wave is not absorbed or radiated by the cells before the board edge and the 12 x 12 mm FR4 window / ground edge radiate; consistent with alpha of the HFSS cells weaker than the model's table (but scaling alpha in the model leaves D almost unchanged, so alpha alone does not explain D); (3) cell data vs HFSS cell: the model uses periodic supercell tables with 0.6 to 1.0 mm gaps between the discs, the real board has finite-array coupling, the vias of the switched cells (resonance near 5.8 GHz is not excluded) and the lumped Impedance sheet instead of a real diode; (4) feed: only 2 % of the power below theta 10 deg, so direct pin radiation is not the cause. NO fix proposed, NO new STEP.

## Flags
- The ffd location / name differ from the message (see top); sizes and headers are consistent (16471 rows).
- The correlation ceiling 0.73 means no variant reproduces the HFSS pattern; the inverted-switching + 180-degree statement is a qualitative result for the beams, not a calibrated one.
- The model's directivity is computed from the model far field on theta <= 90; HFSS D uses the full sphere (theta up to 180).
