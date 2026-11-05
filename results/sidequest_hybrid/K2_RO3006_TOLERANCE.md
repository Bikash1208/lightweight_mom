# RO3006 1.27 mm, r 2.2 cell: tolerance, launch pad (finite radius), HFSS cell package - MODEL_ONLY (2026-11-05)
Approved cell: RO3006 1.27 mm, disc r 2.2, rim 0.2, lattice 6.4. Board package NOT built (waits for the user's HFSS state 2 and the cell check). <= 514 diodes, one control line.

## 2. Labelling (A0 calibration)
**The factor 3.13/2.63 applied to the MoM dY is a calibration taken from the baseline HFSS cell (HFSS dY 3.13 mS against MoM 2.63 mS for the RT5880 r 3.0 cell), NOT a first-principles value.** It is applied to every new cell under the assumption that the MoM under-predicts dY by the same 16 % there; without it all A0 values are x0.70 (RO3006 r 2.2: 3.37 instead of 4.77 Np/m, x4.0 instead of x5.6 baseline). The first-order A0 formula itself has no fitted constant (0.836 against HFSS 0.844 on the baseline cell, with the HFSS dY). Gains below use the calibrated A0.

## 1. Tolerance (r 2.2 +-0.05 mm, eps_r 6.15 +-0.15, diode slit gap 0.4 +-0.03 mm; centre + 6 one-at-a-time + 8 vertex corners; `sidequest_tolerance.py`, `k2_tolerance_table.json`)
The disc-disc gap is 6.4 - 2r = 2.0 mm; r +-0.05 moves it by +-0.10 mm, so +-0.03 mm of gap is r +-0.015, inside the r corners. The diode slit gap +-0.03 mm (and Rs = 2.44 x 1.0 / G: 6.6 / 5.7 ohm/sq) has no visible effect (A0 4.72 .. 4.77).
Gain = calibrated to the HFSS state-1 point (offset -1.78 dB), state 1 / 2 [dBi]; bare = launch 0.71, pad = launch 0.96 (the pad is an upper bound, see section 3):
| case | r | eps_r | A0 [Np/m] | beta/k0 | dX [ohm] | gain bare | gain pad | SW at edge |
|---|---|---|---|---|---|---|---|---|
| centre | 2.20 | 6.15 | 4.77 | 1.201 | 108 | 15.2 / 16.1 | 16.5 / 17.4 | 44 % |
| r - | 2.15 | 6.15 | 2.05 | 1.164 | 77 | 13.4 / 12.5 | 14.7 / 13.8 | 64 % |
| r + | 2.25 | 6.15 | 6.67 | 1.232 | 122 | 15.9 / 16.7 | 17.3 / 18.0 | 32 % |
| eps_r - | 2.20 | 6.00 | 3.04 | 1.175 | 89 | 13.8 / 14.6 | 15.1 / 15.9 | 56 % |
| eps_r + | 2.20 | 6.30 | 7.33 | 1.233 | 130 | 16.1 / 16.9 | 17.4 / 18.2 | 29 % |
| **worst corner: r 2.15, eps_r 6.00** | 2.15 | 6.00 | **1.29-1.31** | 1.145 | 63 | **11.7 / 10.8** | **13.0 / 12.1** | 73 % |
| best corner: r 2.25, eps_r 6.30 | 2.25 | 6.30 | 9.6-9.75 | 1.270 | 144 | 16.3 / 17.0 | 17.7 / 18.3 | 20 % |
Ranges over all 15 cases: A0 1.29 .. 9.75 Np/m (x7.6), beta/k0 1.144 .. 1.270, gain (bare) 10.8 .. 17.0 dBi; diodes 500-514 in every case. **Worst case: -4.4 dB (state 1, 15.2 -> 11.7) and -5.3 dB (state 2, 16.1 -> 10.8): more than 2 dB**, driven by r and eps_r together (A0 is a steep function of the loading: A0 falls x3.7 at the low corner). Low-side corners lose gain, high-side corners gain; the gain is flat above A0 about 4 (A0 9.7 -> 16.4 / 17.0 dBi).
**Less resonant / more robust option: move the nominal UP the curve, r 2.25** (A0 6.67, beta/k0 1.23, dX 122 ohm): its low corner (r 2.20, eps_r 6.00, already computed) gives 13.8 / 14.6 dBi bare (15.1 / 15.9 pad) against its own nominal 15.9 / 16.7: a -2.1 dB drop (about the 2 dB limit) and a worst-case level 2.1 dB above the r 2.2 worst case. The r 2.15 option is NOT more robust: it sits on the steep part (A0 2.05, gain 13.4 / 12.5 nominal, 11.7 / 10.8 at its worst). The high-side r 2.30 corners of the r 2.25 option are not computed (A0 about 13 in the uncorrected MoM, beta/k0 1.36). Cheaper alternative: hold r +-0.02 mm and eps_r +-0.05 (RO3006 datasheet tolerance is +-0.15). Both need the HFSS cell check first.

## 3. Launch pad at its real finite radius (axisymmetric FD, `sidequest_launcher_fd.py`, `sidequest_pad_study.py`, `k2_pad_study.json`)
Model: z-pin (6 mm, r 0.635), grounded slab, uniform sheet j*Ybar on rho >= rho_in, bare inside; the FD launch is taken RELATIVE to the uniformly loaded case x the spectral uniform value (the FD energy balance is not closed: pin flux 14 % below the spectral R_in, P_SW agrees to 1 %). R_in, X = thin-pin model with a prescribed sinusoidal current.
| stack | rho_in [mm] | launch (rel.) | R_in [ohm] | X_in [ohm] | S11 untuned / X cancelled [dB] |
|---|---|---|---|---|---|
| RO3006 | 0.64 (loaded to the pin) | 0.96 | 50 | +50 | -7 / -inf |
| RO3006 | **2.6 (first cell ring around the ground opening)** | **0.963** | **43** | **+37** | **-8.6 / -22.8** |
| RO3006 | 4.5 | 0.97 | 34 | +30 | -7.4 / -16 |
| RO3006 | 6.5 (as the present board) | 0.99 | 28 | +28 | -6.7 / -11 |
| RT5880 (baseline stack) | 6.5 (present board) | 0.91 | 40 | +19 | NA |
**The model does NOT pass its own baseline check**: for the present RT5880 board (cells from 6.5 mm) it gives launch 0.91 where the HFSS-based value is 0.627, and X +19 ohm (S11 about -12 dB) against the HFSS bare-probe -15 dB. It also shows almost no dependence of the launch on rho_in (0.96-0.99), which contradicts the bare-slab / HFSS agreement (spectral bare 0.611 against 0.627). So the pad launch (0.96) and any launch >= 0.8 claim are UNVALIDATED; the real launch could stay near the bare value (RO3006 0.71, or 0.63 as measured on RT5880). The reactance is inductive (+37 ohm at 2.6 mm): S11 <= -15 dB needs a pin trim (shorter pin / series C) to bring |X| under about 14 ohm; this is an estimate only. HFSS (bare probe + pad geometry) is the decider.

## 4. HFSS cell package `hfss_exports/sidequest_cell_RO3006_r2p2/` (synced to `C:\hmsa_sync\data\sidequest\cell_RO3006_r2p2\`)
STEP (`METAL_PEC`, `SUBSTRATE`, `AIR`, `DX_SHEET`, `DY_SHEET`), `HFSS_SETUP.md` (1 page), `states.csv`, `predictions.json`, `geometry_checks.json` (all PASS on the re-imported STEP; EXACT_DIMS PASS, 0.001 mm grid), `results_template.csv`. STEP files are gitignored (108 kB).
* Air box: one box from the ground bottom (z -1.305) to the lid (z +30), one face per side wall, Primary/Secondary pairs on all four side faces (they cover the substrate); `METAL_PEC` gets material `pec`; Perfect E lid on the top face.
* States: OFF = sheets unassigned; ON = ONE Impedance boundary on both sheets, **Rs = 2.44 x W/L = 2.44 x 1.0 / 0.4 = 6.1 ohm/sq, W/L = 2.5**, reactance 0.
* Eigenmode at 90 deg along x (period 6.4 mm), y run swapped. **Prediction (MODEL_ONLY, +-0.15 GHz): OFF 9.899 GHz, ON 10.356 GHz (shift +0.457 GHz, ON above OFF).** Sensitivity OFF / ON: r 2.15: 10.038 / 10.434; r 2.25: 9.815 / 10.276; eps_r 6.00: 9.979 / 10.416; eps_r 6.30: 9.819 / 10.295 GHz. The RT5880 version of this pipeline was within 0.01 GHz of HFSS (model 10.662 / 11.007, HFSS 10.670 / 11.006).
* Falsifier: shift below 0.1 GHz or ON below OFF.
