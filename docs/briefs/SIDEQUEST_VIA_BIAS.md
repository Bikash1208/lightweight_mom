# SIDEQUEST_VIA_BIAS - the via-switched GB cell WITH the real bias hardware (fix B), fix A assessment, board-level bias (2026-11-06). MODEL_ONLY
Trigger: the 5.6 mm cross slits leave a 0.2 mm rim, so the four quarters of the existing GB cell are DC-connected and its two diodes are DC-shorted (HFSS results stay valid as RF only). Fix B = plain top patch, via through a ground anti-pad, ONE diode on the back, bias through a feed on the back layer (3 Cu layers). Package: `hfss_exports/sidequest_cell_via_RT5880/` (STEP, `HFSS_SETUP.md`, `predictions.json`, `states.csv`, `results_template.csv`, `geometry_checks.json` all PASS, EXACT_DIMS PASS; synced to `C:\hmsa_sync\data\sidequest\cell_via_RT5880`). Scripts `sidequest_via_hw.py`, `_hwsearch.py`, `_hwboard.py`, `_cell_package.py`, `_cell_docs.py`; data `via_hw_search.json`, `via_hw_final.json`, `via_hwboard.json`.

## 1. The hardware that is now in the cell model and in the HFSS geometry
Stack: patch L1 / RT5880 1.575 / ground L2 with a **slot anti-pad** (x -0.55..1.15, y -0.45..1.25 mm, covering the via land, the series-C gap and the diode island) / **bond ply 0.101 mm** (RO4450F-class: eps_r 3.52, tan d 0.004, Rogers datasheet via mirrors; bonding to the PTFE RT5880 must be confirmed, RO3001 film 0.038 mm is the PTFE alternative, Dk not found) / back layer L3 (35 um). L3: via land 0.5 x 0.5, **series-C gap** (0.216 pF knob), **island A** with the pad `A_up`, **diode gap** (MA4AGP907), **cathode pad K with a blind via to L2**, **ballast-resistor gap** (4.7 kohm, 01005-class) to a **bias pad**, a connection to a **bias trunk** (x = 2.5, full pitch, continues into the neighbour cell). Thin-ply capacitance is 0.31 pF/mm2, so every pad over copper costs tens of fF: the island and the land are kept inside the anti-pad (about 10-15 fF each).
## 2. How the hardware changes beta and the step (C_s = 0.614 pF of the ideal cell, unchanged)
| configuration | Y_ON [mS] | Y_OFF [mS] | step | beta_mean (target 1.062) |
|---|---|---|---|---|
| ideal branch (no hardware) | -9.44 + 0.22 real | +5.58 | 15.0 | 1.062 |
| + pads, anti-pad, blind via (C_A 10 fF, C_V 15 fF, 0.15 nH) | -7.37 | +6.55 | 13.9 | 1.081 |
| + ballast 4.7 kohm with 8 fF body (final hardware) | -7.37 | +8.48 (+0.98 real) | 15.8 | **1.094** |
| final hardware, **C_s re-tuned to 0.216 pF** | -10.86 (+0.30) | +7.02 (+0.63 real) | **17.9** | **1.062** |
So the hardware moves beta_mean by +0.032 and the OFF susceptance by +3 mS at fixed C_s; **C_s must fall from 0.614 to 0.216 pF (0.258 pF if the cathode pad floats for series strings)**, which gives a step of 17.9 mS. Resonances of the re-tuned branch: ON 6.8 GHz (-3.7 GHz), OFF 14.8 GHz (+4.3 GHz); beta ON / mean / OFF = 1.015 / 1.062 / 1.19; alpha_loss ON 0.13, **OFF 1.94 Np/m** (the ballast and body capacitance make the OFF cell lossy: Re Y_OFF 0.63 mS; board effect about 0.07 Np/m on the mean).
**Why a resistive feed and not a lambda/4 choke:** a microstrip choke on the 0.1 mm ply has Z0 = 64 ohm; its input admittance is 2-9 mS (as large as the OFF diode's 1.9 mS): at C_s 0.614 pF a 3 / 4 / 4.9 / 6 mm choke gives step 4.8 / 9.2 / 34 / -12 mS and beta_mean 1.03 / 1.05 / 1.21 / 1.01, and any choke longer than 4.5 mm creates an OFF-state |Y| extremum within 1.5 GHz of 10.5 GHz (9.3-11.8 GHz). A stepped-impedance three-section choke (about 13 mm, 4.9 kohm) isolates only within a few percent of 10.5 GHz and puts the OFF parallel resonance inside the 1.5 GHz window. The resistive feed adds a predictable capacitance and loss, no resonance.
**Parasitic budget found by search** (step 12-20 mS, mean at 7.55 mS, Re Y_OFF <= 1.5 mS, both branches >= 1.5 GHz from any |Y| extremum): feasible only for island capacitance C_A <= 15 fF **and** resistor body capacitance C_p <= 8-15 fF, R_b >= 2.2 kohm (4.7 kohm gives Re Y_OFF 0.6 mS), no extra series inductance; with C_A 35 fF and C_p 25 fF (a normal 0402 on a normal pad) nothing is feasible. The layout must therefore keep the island inside the anti-pad and use an 01005-class resistor.
## 3. Tolerance of the final hardware cell (single factor; change of beta_mean / step)
| factor | beta_mean | step |
|---|---|---|
| diode C +-20 % | -0.005 / +0.007 | -5.6 % / +6.3 % |
| drill +-0.05 mm (via L +-8 %) | -0.009 / +0.010 | +5.1 % / -3.1 % |
| C_s +-5 % (a printed gap: etch-limited; a 0201 chip at +-0.05 pF is +-23 %: beta_mean +-0.011, step +-7 %) | -0.003 / +0.002 | +1.8 % / -1.5 % |
| island C_A +-30 % | -0.003 / +0.003 | -3.0 % / +3.1 % |
| resistor body C_p +-30 % | -0.002 / +0.002 | -2.3 % / +2.4 % |
| via-land C_V +-30 % (ply thickness +-15 %) | +-0.007 | +-6 % |
| R_b +-20 %, diode R x8 (1 mA) | <= 0.002 | <= 2 % |
| etch +-0.03 mm per edge (patch gap, from the earlier solves) | -0.005 / +0.009 | 0 |
| via-L model (kappa 0.7) | **-0.035** | **+37 %** |
RSS (diode C, drill, C_s +-20 %, C_A, C_p, C_V, etch): beta_mean **0.020** (limit 0.02: marginal), step **13 %** (alpha +-25 % allows about 12 %: marginal). Worst-case boards (same design, perturbed levels): balance moves <= 0.8 dB. The kappa uncertainty is again the largest: the HFSS eigenmode check decides it.
## 4. Board prediction with the hardware cell (calibrated model, sign -1, one diode per cell, cathode by blind via)
| diodes | flips | beam 1 OFF / ON [dBi] | change | AR OFF / ON [dB] | beam 2 OFF / ON [dBi] | contrast | balance (tolerance range) | SLL OFF / ON [dB] |
|---|---|---|---|---|---|---|---|---|
| 90 | 62 | 22.3 / 20.7 | -1.7 | 1.0 / 0.8 | -7.5 / 18.2 | 25.6 | 2.5 (1.2-3.2) | -15.7 / -12.1 |
| **150** | 44 | 22.2 / 20.9 | -1.3 | 1.4 / 1.1 (worst in beam 1.7) | -1.7 / 19.0 | 20.7 | **1.9** (0.6-2.6) | -16.5 / -13.7 (worst ON -13.0) |
| 180 | 36 | 22.4 / 21.1 | -1.3 | 1.4 / 1.2 | 0.3 / 19.4 | 19.1 | 1.7 (0.4-2.4) | -15.9 / -13.8 |
The hardware cell reproduces the ideal-cell board result (150 diodes: balance 1.9 dB, beam 1 20.9 dBi, AR <= 1.7 dB, SLL <= -13.7 dB); step 17.9 mS. Effect on beam 1 AR: none beyond the model scatter (1.4 / 1.1 dB against 1.2 / 0.9 for the ideal cell). **HFSS cell prediction (lid 30 mm, x, 90 deg): OFF 10.08 GHz, ON 11.38 GHz (shift 1.30 GHz); 110 deg: 11.44 / 13.28; 45 deg (beta p 99 deg): 10.65 / 12.26 GHz; sensitivity to kappa 0.7: +0.1 GHz.**
## 5. Fix A (rim cut + DC-block capacitors + bias vias for the existing 2-diode cell), one paragraph
Cut the 0.2 mm rim at two places (or all four ends) so the quarters are DC-isolated, bridge each cut with a 0201 DC-block capacitor (0.5-1 pF; ESL 0.3 nH gives X = -15..-30 + 20 ohm, comparable to the rim's own reactance), and feed the quarters through bias vias in ground anti-pads to L3 with a choke or resistor per via: per switched cell that is 2 diodes + 2 caps + 2 vias/anti-pads + 2 feed resistors (8 elements, 6 assembled parts) against 1 diode + 1 chip C (or a printed gap) + 1 resistor + 1 blind via (fix B: 3-4 parts); for 180 cells 360 diodes against 150-180. RF impact: the rim is part of the cell's capacitive coupling, so cutting it and adding series capacitors changes Y_OFF and Y_ON of the cell (needs a new MoM + lumped re-evaluation, not done), the ground gets two extra holes under each patch (anti-pad leakage, inductance), and each feed resistor/choke sits across a diode exactly as in section 2 (the same 2-9 mS shunt problem, twice per cell). Above all it keeps the 2.6 mS step: balance stays about 14 dB (SIDEQUEST_CP_STEP3), so fix A repairs the bias but not the beam-2 level; fix B repairs both.
## 6. Board-level bias (chosen diode MA4AGP907; V_F 1.33 V typ / 1.45 V max at 10 mA from the MACOM V4 sheet; low-current V_F and R_on(I) are NOT in the data: 1.2 V and 8 x 2.44 ohm at 1 mA are ESTIMATES, RF-insensitive: diode R x8 moves beta_mean by +0.002, step -1.8 %)
Default = control line unpowered (0 V) = diodes OFF. The 4.7 kohm ballast is also the RF feed, so the string current is set at 1 mA.
| scheme | per cell | string / bus | line voltage | total current | notes |
|---|---|---|---|---|---|
| **P: parallel** (cathode to L2 by the blind via, used in the package) | R_b 4.7 kohm from the bus to the anode | one bus tree, 150 branches | V_F + R_b I = 1.2 + 4.7 = **5.9 V** | 150 x 1 mA = **0.15 A**, 0.9 W | stiff: V_F spread +-0.15 V gives +-0.03 mA (3 %); resistor 4.7 mW (01005 rating about 31 mW) |
| **S: series-parallel** (cathode pad floats, K pad = RF ground; C_s 0.258 pF) | R_b between K(k) and A(k+1) is the cell's own feed | strings of **4 cells** (4 diodes, 4 x 4.7 kohm), 38 strings (152 diodes) from one line | 4 x (1.2 + 4.7) = **23.6 V** (24 V supply, no extra ballast) | 38 x 1 mA = **38 mA**, 0.9 W | last cathode to ground by a via; string of 8: 47 V and 19 mA; the reference 32 x 16 string at 24 V needs a low-resistance feed (the lambda/4 choke), which fails the RF margin (section 2) |
Hot-spot and layout: the bias trunk (0.3 mm) on L3 runs under the same cell row; bus routing for 150 scattered sites uses L3 only, with the ballast right at each island (the resistor is the RF feed, so it must sit at the island, not at the bus).
## 7. Risks
1. Kappa (via inductance) and the pad / resistor capacitance are lumped estimates: kappa 0.7 moves beta_mean by -0.035 and the step by +37 %; the parasitic budget (C_A <= 15 fF, C_p <= 8-15 fF) is tight and unverified.
2. Bond ply: RO4450F on PTFE RT5880 is not confirmed (RO3001 film is thinner, which raises every pad capacitance by 2.5x and breaks the budget).
3. The ballast makes the OFF cell lossy (alpha_loss 1.9 Np/m in a uniform OFF lattice, board effect about 0.07 Np/m).
## 8. Next
User runs the cell package (12 eigenmode runs + 2 radiation-top runs). The re-tune rule for C_s is in `predictions.json`. The board package follows after the cell passes.
