# phase28_loop_R4 -- first HFSS run of the loop cell R4 (user, 2026-10-27)

Loop cell R4 (outer 6.6 mm, strip 1.2 mm, one diode 0.6 x 1.2 mm), 9 mm period, RT/duroid 5880 3.175 mm on PEC. Floquet port, normal incidence, de-embedded to the substrate surface, ROT 0, 9 - 11 GHz (0.05 GHz).
Diode as lumped RLC boundaries: R only 2.4 ohm and C only 28.2 fF. Boundaries and port were set by hand (the script failed at AssignMaster; fixed in 9ea6897). Files: `ON.csv`, `OFF.csv` (complex S of Floquet modes 1, 2).

**State labels: swapped in the project.** `OFF.csv` carries the 2.4 ohm boundary (lossy, isotropic) and `ON.csv` the 28.2 fF boundary (open, anisotropic); the loss G(mode 2) = 0.105 mS appears only in OFF.csv and equals the MoM
dissipation of the R state (0.09 mS). Mode 2 is x (along the diode), mode 1 is y. With the physical pairing all sheet-admittance components agree with the MoM within 4.3 % (loss within 17 %):

| 10 GHz | B_x [mS] HFSS / MoM | B_y [mS] HFSS / MoM | G_x [mS] HFSS / MoM |
|---|---|---|---|
| R state (OFF.csv) | 6.544 / 6.276 | 6.545 / 6.270 | 0.105 / 0.089 |
| C state (ON.csv) | 2.326 / 2.230 | 6.419 / 6.237 | 0.002 / 0.002 |

Delta B along the diode 4.22 mS (MoM 4.05), mean b_bar 4.44 mS (MoM 4.25); Delta B 3.50 - 5.43 mS over 9 - 11 GHz. Exact sheet-model leakage with the HFSS tensors 2.86 Np/m at 10 GHz (MoM 2.69).
Extraction and comparison: `scripts/research/phase28_hfss_cell_compare.py`; calibration hook: `src/ris/loop_cell_calibration.py`; results `results/phase28/f10/loop_R4/hfss_vs_mom.json`, `hfss_calibration.json`.
Label: HFSS_CHECKED for this cell at ROT 0, normal incidence, 9 - 11 GHz (one run). Not checked: the directionality at oblique polarisation angles and the rotated cell, the diode with the measured reactance, L_s.
