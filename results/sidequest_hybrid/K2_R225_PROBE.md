# r 2.25 nominal corners, bare-launch gains, probe trim - MODEL_ONLY (2026-11-05); no package built
Rules used: **all gains with the bare launch 0.71** (pad = UNVALIDATED, not used in any number here); gains calibrated to the HFSS state-1 point (offset -1.78 dB); <= 514 diodes, one line. eps_r +-0.15 (Rogers spec, +-0.05 is not obtainable by purchasing: the design must survive +-0.15); etching +-0.03 mm assumed. The dY x 3.13/2.63 factor is a calibration from the baseline HFSS cell, not first principles (A0_mom = without it).
Cell-check decision rule (user): HFSS OFF/OFF within +-0.05 GHz of the predictions (OFF 9.899 / ON 10.356 GHz at r 2.2) -> move to r 2.25; if shifted, re-centre r from the HFSS offset. Prediction sensitivity: d f_OFF / d r = -0.17 GHz per +0.05 mm, d f_OFF / d eps_r = -0.16 GHz per +0.15.

## 1. r 2.25 nominal with corners (`sidequest_tolerance.py`, `k2_tolerance_r225_table.json`; gain state 1 / 2 [dBi])
| r [mm] | eps_r 6.00 | eps_r 6.15 | eps_r 6.30 |
|---|---|---|---|
| 2.22 (-0.03 etch) | A0 3.9, b/k0 1.19, 14.6 / 15.5 | 6.0, 1.22, 15.7 / 16.6 | 9.0, 1.26, 16.2 / 17.1 |
| **2.25 (nominal)** | 4.4, 1.20, 14.9 / 15.9 | **6.7, 1.23, 16.0 / 16.7** | 9.8, 1.27, 16.3 / 17.0 |
| 2.28 (+0.03 etch) | 15.1, 1.29, 16.8 / 15.7 | 19.3, 1.35, 16.4 / 15.4 | 22.2, 1.41, 16.0 / 15.2 |
| 2.30 (+0.05) | 15.4, 1.31, 16.9 / 16.7 | 19.1, 1.36, 17.3 / 16.9 | 21.4, 1.42, 17.0 / 16.4 |
| 2.35 (+0.10) | | | 17.2, 1.585, 17.1 / 15.9 |
Requested corners: **r 2.30, eps_r 6.30: A0 21.4 Np/m (x25 baseline), beta/k0 1.42, gain 17.0 / 16.4 dBi, SW at the edge 4 %, 500 diodes**; **r 2.25, eps_r 6.30: A0 9.75, beta/k0 1.27, 16.3 / 17.0 dBi**.
* **Worst case with etching +-0.03 and eps_r +-0.15 around r 2.25: 14.6 / 15.5 dBi (r 2.22, eps_r 6.00) against the nominal 16.0 / 16.7: -1.4 / -1.2 dB, inside the 2 dB limit.** The best case is 16.8 dBi (state 1). With the larger r +-0.05 (r 2.20, eps_r 6.00) the earlier value 13.8 / 14.6 gives -2.2 dB.
* **Resonance (A0 peak, |Ybar + Y_up + Y_down(kx)| minimum) is approached but not crossed inside the tolerance box**: A0 saturates at r 2.28-2.30 (15-21 Np/m, the MoM A0 is 13.6 at eps_r 6.15 for both r 2.28 and 2.30) and the peak lies at r about 2.30-2.33; at eps_r 6.30 A0 is already past it at r 2.35 (17.2, down from 21.4). The plateau is flat, so gain barely changes there (16.4-17.3 dBi); the price is beta/k0 up to 1.42 (the engine phase constant moves from 1.07 to 1.2-1.4: the hologram period changes; the model has no check of the bound-mode validity at 1.4) and a state imbalance up to 1.1 dB (state 1 above state 2 at r >= 2.28, eps_r 6.00-6.30).
* Edge power (SW reaching the edge): nominal 32 %; 50 % at the low corner; 2-9 % at r >= 2.28. The <= 15 % edge target holds only on the high side.

## 3. Probe trim (thin-pin FD, `sidequest_probe_trim.py`, `k2_probe_trim.json`; MODEL_ONLY; bare launch geometry = cells from 6.5 mm)
Z_in = R + jX of the 50 ohm probe vs pin height t above the surface (a 0.635 mm):
| t [mm] | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|
| RO3006 bare: R / X [ohm] | 9 / -48 | 16 / -18 | 24 / +5 | 35 / +29 | 50 / +57 | 71 / +91 |
| S11 [dB] | -1.6 | -4.9 | -9.1 | -8.8 | -6.1 | -4.2 |
| RT5880 baseline (HFSS: S11 -15 dB at t 6) | S11 -2.4 | -6.1 | -11.5 | **-12.7** | -8.4 | -6.1 |
Calibration check: the model gives -12.7 dB for the baseline (HFSS -15): the X is biased by an unknown few ohm (the coax opening and the feed are not modelled). A thicker pin (a 1.0 mm) lowers X by about 6 ohm at t 6 (R 34, X +23): not enough alone.
**Proposed trim (|X| < 14 ohm at the bare launch):** shortening the pin only moves X to zero at t about 5 mm, where R is 24 ohm (S11 -9 dB): not enough. **Recommended: pin t = 6.5 mm (R about 42 ohm, X about +43 ohm) with a series capacitance C = 1/(w X) = 0.35 pF at the pin base (a small gap / dielectric washer in the centre pin or an annular coupling pad), giving X about 0, S11 about -22 dB; the window |X| < 14 ohm (S11 <= -15 dB) is C = 0.26 .. 0.52 pF.** Alternative: t = 7 mm (R 50, X +57): C 0.27 pF (window 0.22-0.36 pF). Tune in HFSS: sweep t 5.5-7 mm and C 0.2-0.5 pF; the model X offset is unknown (compare the baseline). If the pad is later validated: pad rho_in 2.6 mm with t = 5 mm gives R 43, X +9, S11 -18 dB without an extra part (UNVALIDATED: the pad launch is not trusted).

## 4. Assumptions recorded
eps_r +-0.15 (Rogers RO3006 spec; purchasing cannot do +-0.05) and etching +-0.03 mm are the design tolerances; the diode slit gap +-0.03 mm has no effect (checked on r 2.2). Everything here is MoM + a first-order leakage formula + a 1D engine, not HFSS-checked.
