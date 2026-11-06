# SIDEQUEST_VIA_CELL - via-switched GB cell for the CP+LP hybrid (priority 1), RT5880 1.575 mm, one diode per cell (2026-11-06). MODEL_ONLY; HFSS cell package and board package NOT built (stop requested)
Calibration of the board model on the HFSS result of the hybrid package is in `SIDEQUEST_CP_HYBRID.md` (beam 1 -1.8 dB, beam 2 +1.5 dB, balance -3.3 dB, AR x0.5, SLL raw model + 1.0 dB). Sign -1 (beam 2 ON at S_w 6.1) is used throughout. Scripts `sidequest_via_cell.py`, `_design.py`, `_final.py`, `_board.py`, `_step.py`; data `via_design_search.json`, `via_final.json`, `via_board.json`, `via_board_final.json`, `via_step_required.json`, `cellB_rt5880_*.json`.

## 1. Required step (board model, levels symmetric about the tensor mean 7.55 mS = beta/k0 1.062, uncalibrated)
A model balance <= 6.6 dB corresponds to an HFSS balance <= 3.3 dB. Required step: 10 mS already gives a model balance 5.9 dB (180 cells, 52 rotated slots) and 5.9 dB (360 cells); 15 mS gives 2.5 dB (model; about -0.8 dB calibrated). **A step of 12-15 mS is enough; 15 mS leaves margin; 20 mS buys little** (model balance 2.1-2.0 at 36-72 slots). Beam 1 stays 22-24 dBi (model).

## 2. Diode choice: the packaged PINs cannot do it on this stack
Branch (per square cell): ON Z = R + j w L_tot - j/(w C_s); OFF Z = j w L_tot - j/(w C_s) - j/(w (C_j + C_pad)); L_tot = mu0 h (1.98 nH via, dense-lattice mushroom value) + L_package + L_ser; C_pad = 25 fF (pad over the 0.1 mm bias ply, an ESTIMATE). Constraints: step 14-21 mS, mean sheet = 7.55 +- 0.25 mS (beta/k0 1.062 of the tensor cells), both branch series resonances >= 1.5 GHz from 10.5 GHz, minimum sensitivity; knobs C_s, L_ser, patch gap (0.4-1.6).
| diode | R_on / C_j / L_pkg | source | feasible settings |
|---|---|---|---|
| MA4AGP907 chip (user measured) | 2.44 ohm / 28.16 fF / 0 | user | **300** |
| MADP-000907-14020 chip | 5.2 ohm / 25 fF / NA (0 used) | MACOM datasheet Rev V8 via repo `src/ris/diodes.py` | 207 (ON loss 2x) |
| SMP1345-079LF | 1.5 ohm / 0.19 pF / 0.7 nH | Skyworks 200046U via repo `src/ris/diodes.py` (datasheet based) | **0** |
| BAR63-02L | 1.0 ohm / 0.21-0.3 pF / 0.6 nH | Infineon preliminary datasheet (Ls); rf and CT from third-party listings, UNVERIFIED | 0 |
| BAR90-02LRH | 1.3 ohm / 0.19 pF / 0.4 nH | Infineon BAR90 family snippets, UNVERIFIED for this part | 0 |
| BAR50-02V | 4.5 ohm / 0.24-0.4 pF / Ls not found (0.7 assumed) | listings disagree, UNVERIFIED | 0 |
Vendor S-parameter models could not be retrieved here; lumped datasheet models are used. **Why:** with L_tot >= 2.7 nH the OFF branch is above its series resonance when C_j + C_pad is above about 85 fF. Probe with a hypothetical 1.5 ohm / 0.7 nH diode: C_j 28 fF and 40 fF are feasible (margins 4.0 / 3.3 GHz); at 60 fF and above no setting meets even step >= 10 mS together with the mean. A series C_s small enough to keep OFF capacitive also makes ON capacitive and resonant. The packaged parts have 150-400 fF. Options for a packaged diode: a thinner substrate (via L is proportional to h) or accepting a mean-beta mismatch. A chip PIN (MA4AGP907 / MADP-000907 / MA4AGFCP910) works.

## 3. Recommended cell (MA4AGP907 chip PIN, polarity A: ON lowers the sheet susceptance)
RT5880 1.575 (eps_r 2.2, tan d 0.0009), p 6.4, square patch 6.0 (gap 0.4 to the neighbours), central via (diameter free in this lumped model; the +-0.05 mm drill is carried as +-8 % of L_via) through a ground anti-pad to a back island, thin bias layer (about 0.1 mm), MA4AGP907 between island and ground (DC return), **series MIM capacitor C_s = 0.614 pF (the knob), L_ser = 0**, pad capacitance about 25 fF (to be drawn small). Stack: patch / RT5880 / ground / thin bias ply / back layer (3 metal layers, as in Route B).
| quantity | value |
|---|---|
| step (Y_OFF - Y_ON) | **15.0 mS**: branch Y_ON -9.4 mS (+0.22 mS real), Y_OFF +5.6 mS; sheet levels ON +0.10 mS, OFF +15.1 mS (symmetric about the 7.55 mS mean) |
| beta/k0 | mean sheet **1.062** (met by C_s); OFF **1.170**; ON **1.018** (arithmetic mean of ON and OFF 1.094; the board sees the mean-sheet value) |
| convergence (mesh 32 / 64) | beta_mean 1.0598 / 1.0620 (0.2 %), beta_OFF 1.1674 / 1.1704 (0.3 %), beta_ON 1.0179 / 1.0184 (0.05 %); the mesh-96 solve is still running (M2 met on two levels so far) |
| isotropy x vs 45 deg (mesh 64) | beta_mean +0.7 %, beta_OFF +0.9 %, beta_ON +0.15 % (limit 3 %) |
| resonances / stop bands | branch series resonance ON **4.6 GHz (-5.9 GHz)**, OFF **16.2 GHz (+5.7 GHz)**; normal-incidence zeros ON 3.1 / 4.6 / 13.6 GHz, OFF 7.8 / 16.2 GHz (nearest 2.7 GHz); the TM root of the ON, OFF and mean sheets exists and is smooth over 7-14 GHz (OFF beta 1.02 -> 1.35, ON 1.003 -> 1.105): **no stop band within 3.5 GHz** |
| loss | alpha_loss ON 0.13 / OFF 0.12 Np/m (dielectric + diode); diode R_on 2.44 ohm |

Tolerance (change of beta_mean / step), single factor, MoM patch + lumped branch:
| factor | beta_mean | step |
|---|---|---|
| etch +-0.03 mm per edge (gap +-0.06) | -0.005 / +0.009 | 0 (patch only; sheet levels move -0.5 / +0.8 mS) |
| drill +-0.05 mm (via L +-8 %) | -0.007 / +0.007 | +4.8 % / -3.3 % |
| diode C +-20 % | -0.005 / +0.006 | -5.9 % / +6.7 % |
| diode R +-20 % | 0 | 0 |
| eps_r +-0.02 (assumed Rogers RT5880 value, not re-verified) | -0.002 / +0.002 | 0 |
| C_s +-10 % | -0.002 / +0.001 | +1.2 % / -0.7 % |
| pad C +-5 fF | -0.004 / +0.005 | -5.3 % / +5.9 % |
| via-L model (kappa 0.7 instead of 1) | **-0.027** | **+30 %** |
RSS over the real tolerances: beta_mean 0.014 (limit 0.02), step 10 % (alpha +-25 % allows about 12 % in the step); worst-case sum beta_mean 0.031. **The via-inductance model (kappa) is the dominant unknown; the HFSS eigenmode check decides it.** Board-level effect of the same perturbations (same design, perturbed levels): balance moves <= 0.7 dB (Section 4).

## 4. Board prediction with this cell (calibrated model, sign -1, polarity A, one diode per cell; flips = tensor slots rotated by 90 deg)
| diodes | flips | beam 1 OFF / ON [dBi] | beam-1 change | AR OFF / ON [dB] | beam 2 OFF / ON [dBi] | contrast | **balance** | SLL OFF / ON [dB] |
|---|---|---|---|---|---|---|---|---|
| 60 | 80 | 21.8 / 20.8 | -1.0 | 0.1 / 0.3 | 0.3 / 16.7 | 16.4 | 4.1 | -14.8 / -14.5 |
| 60 | 124 | 21.4 / 20.3 | -1.0 | 0.4 / 0.4 | 7.2 / 18.2 | 11.0 | 2.1 | -11.1 / -14.9 |
| **90** (C_s 0.614) | 76 | 22.0 / 20.4 | -1.5 | 0.5 / 0.6 | -6.0 / 17.9 | 23.9 | **2.5** (tolerance 1.9-3.1) | -14.8 / -12.6 |
| 120 | 54 | 22.1 / 20.5 | -1.6 | 0.9 / 0.7 | -6.6 / 18.3 | 24.9 | 2.2 | -15.0 / -12.0 |
| **150** (C_s 0.614) | 52 | 21.9 / 20.6 | -1.4 | 1.2 / 0.9 (worst in beam 1.6) | -6.3 / 18.7 | 25.0 | **1.9** (1.2-2.5) | -16.5 / -14.2 |
| 180 (C_s 0.614) | 43 | 21.9 / 20.7 | -1.2 | 1.3 / 0.9 | -3.9 / 18.9 | 22.9 | 1.8 (1.1-2.4) | -16.3 / -13.9 |
| 270 | 17 | 21.9 / 21.2 | -0.7 | 1.4 / 0.9 | -1.2 / 19.9 | 21.0 | 1.3 | -16.6 / -15.3 |
| 360 | 0 | 21.7 / 21.3 | -0.4 | 1.4 / 0.9 | 0.7 / 20.1 | 19.4 | 1.2 | -16.5 / -15.0 |
(Rows without the C_s 0.614 label use C_s 0.518 pF: step 15.3 mS, mean 7.33 mS.) Targets: balance <= 3 dB, beam 1 >= 20 dBi, AR <= 2 dB, SLL <= -13 dB, beam-1 change <= 1 dB plus the power-sharing minimum (two equal beams must share the power: about -1.8 to -3 dB at balance 0). **Fewest diodes meeting everything: 150 (balance 1.9, SLL -16.5 / -14.2, beam-1 change -1.4). 90 diodes meet balance (2.5) and AR but SLL ON is -12.6 (0.4 dB short) and the change is -1.5.** More diodes reduce the beam-1 change (-0.4 dB at 360 diodes, no rotated slot) and the rotated-slot count (0 at 360).

## 5. Risks
1. The via inductance (kappa) and the pad capacitance are lumped estimates: kappa 0.7 moves beta_mean by -0.027 and the step by +30 %; the pad C is not in the measured diode model. HFSS eigenmode first.
2. The ON sheet level is +0.1 mS (+-1 mS with the drill): beta_ON is close to the light line, and a negative level removes the lowest TM root of the ON cell (the mean-sheet model the board uses is unaffected).
3. Packaged low-cost PIN parts do NOT work on 1.575 mm RT5880 (Section 2); a chip PIN needs flip-chip or wire-bond assembly. The Infineon datasheet values are unverified.
4. The beam-2 OFF level (a cancellation null) is fragile (contrast 17-36 dB under +-20 % cell changes); balance and beam 1 are not.

## 6. Next (not started, per the instruction to stop after the cell design)
HFSS cell package: eigenmode OFF / ON at 2 phases along x, y, 45 deg with the real via, anti-pad, bias ply and series C_s, plus a radiation-top run for alpha; then the board package. Priority 2 (4-quadrant outline) not started.
