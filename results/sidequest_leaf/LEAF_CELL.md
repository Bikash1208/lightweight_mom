# Leaf cell (Zhao, Pan, Shen, IEEE TAP 74(1) 2026, doi 10.1109/TAP.2025.3621516) - MODEL_ONLY (2026-11-05)
Method: MoM normal-incidence sheet Y (`sheet_tensor`) -> TM0 beta (grounded slab, + PEC lid lambda/4 above) -> X_s = eta alpha/k0 [tanh(alpha H)]; alpha from the real part of the residual. Script `scripts/research/sidequest_leaf_cell.py`, data `v6p9*.json`, `v10p5.json`.
## 1. Validation at 6.9 GHz: INCONCLUSIVE (not validated)
Paper X: 223.3 (OFF) / 100.9 (ON) ohm. Fig. 2 / Table I are not readable from the text, so the leaf is RECONSTRUCTED (lens of r1 6.4 / r2 6.0 mm, tip-to-tip L 6.0, gap 0.5, pads 0.6, a 13.3; diode = their RLC on a 0.5 x 0.6 footprint). Result (OFF / ON, X with lid):
| mesh | 54 | 60 | 66 | 70 | 76 | 84 |
|---|---|---|---|---|---|---|
| OFF [ohm] | 96 | 93 | 90 | 456 | 242 | 143 |
| ON [ohm] | 46 | 47 | 48 | 41 | 42 | 43 |
OFF is NOT mesh-converged (the 0.5 mm gap is not on the 13.3/n grid; OFF sits near the slab pole so X swings 90..456 ohm); ON is stable (41-48, paper 100.9, i.e. 0.45x). Leaf length scan at mesh 54 (L 5.0/5.5/6.0/6.2): OFF 71/80/96/117, ON 35/44/46/47. A 45-deg diode axis gave no load effect (rotated load not engaged; not pursued). Converged run needs pitch 0.1 mm (n = 133, ~35k unknowns): not run. The 5 % criterion is not met and cannot be tested with this geometry guess; the solver is NOT called validated by this paper.
## 2. 10.5 GHz adaptation (a 8.6 = 0.30 lambda0, gap 0.4, pad 0.6, leaf L 3.6 / r 4.1, 3.9; mesh 86 = 0.1 mm pitch, gap and pad edges on the grid)
| substrate | diode model | X OFF / ON (lid) [ohm] | contrast [ohm] | alpha OFF / ON [Np/m] |
|---|---|---|---|---|
| RO3006 1.27 (er 6.15) | SMP1340 (OFF 0.23 pF + 0.45 nH, ON 2 ohm + 0.45 nH) | 361 / 55 | 306 | 1.36 / 0.17 |
| RO3006 1.27 | MA4AGP907 (HFSS-confirmed: OFF open, ON 2.44 ohm) | 109 / 17 | 92 | 0.02 / 0.13 |
| RT5880 1.575 | SMP1340 | 74 / 76 (no contrast) | -2 | 0.01 / 0.04 |
| RT5880 1.575 | MA4AGP907 | 68 / 78 | -10 | 0.01 / 0.04 |
RT5880 (er 2.2) gives no usable reactance swing in this cell. RO3006 works, but the stack changes: the calibrated centre probe (Phase 26, RT5880 1.575) is NOT valid on it (new HFSS probe calibration needed). SMP1340 OFF X 361 ohm is near the slab pole (Y 10 mS), so it is sensitive. 150 mm hexagon: about 740 cells (a 8.6, leaf inside the outline, rho >= 6.5 mm) against 1326 at 6.4; 1 diode per cell, so at most 740 diodes (the K2-default (i) uses 514 for 257 switched cells).
## 3. Board redo with the leaf cell and side-fed launcher: NOT RUN (NA)
Blocked by (a) RT5880 gives no contrast and RO3006 needs a new probe calibration, (b) the engine has no leaf-cell static library or side-fed launcher model, (c) the cell needs HFSS confirmation first.
## 4. HFSS eigenmode package for the leaf cell: NOT BUILT (NA); waits for the substrate decision.

## 1b. Exact Table I run (a 13.3, g 0.5, l 0.6, r1 6.4, r2 6.0; lens leaves along y; stack eps 6.15 1.27 mm; lid lambda/4) - FAIL (`t1_conv.json`, `sidequest_leaf_cell2.py`)
Intercell clearance s (0.8 mm) is not in Table I (assumption). X along x (diode axis), lid, mesh 0.15 / 0.12 / 0.10 mm (g snapped 0.448 / 0.479 / 0.5): OFF 544 / 214 / 440 ohm, ON 52 / 58 / 53 ohm; paper 223.3 / 100.9. Y OFF 21 mS at 0.1 mm = near the slab pole (X unbounded). X along y (across the leaves, no diode in the path) 366-378 ohm in both states. Not converged (OFF), ON half of the paper: not within 10 %.
Likely causes: (1) OFF sits at the slab pole (Y 13-21 mS against B_ind about 13 mS): X = 1/(B_ind - B) swings with 1 % in B or the gap; (2) the normal-incidence Y (k_t = 0) is used for a wave at beta = 1.2-1.4 k0: spatial dispersion of the gap / lens is not modelled, the paper uses the eigenmode at the real phase shift; (3) the intercell clearance and the lens width are not in Table I; (4) lid and reactance definition (X = eta alpha/k0 tanh(alpha H), Pavone convention) not confirmed, the paper may define X on a different plane; (5) 0.1 mm still marginal for a 0.5 mm gap bridge with 0.6 mm pads. The solver is not validated by this paper; step 2 (substrates) was not run as specified (a partial RO3010 run crashed on the alpha step and is discarded).
