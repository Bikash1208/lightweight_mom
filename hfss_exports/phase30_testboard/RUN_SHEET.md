# RUN_SHEET - Phase 30 test board (short). Full details: `HFSS_TESTBOARD_SETUP.md`. MODEL_ONLY predictions, written before any HFSS run.

Two runs, identical setup, eps_r 11.2, tan d 0.0022, 10.05 GHz, one 1 W lumped port (50 ohm, coax pin):
1. `testboard_active.step` (572 slotted cells)  2. `testboard_control.step` (same board, unslotted cells)

Setup (both runs): import STEP (mm, heal OFF) > SUB = RO3010 eps_r 11.2 / tan d 0.0022; GND, PATCH, PIN = Perfect E > region padding 15 mm on all sides with **Radiation** boundary (say if PML) > Lumped Port on PORT, 1 W > solve **10.05 GHz only**, 15 passes, delta S 0.02 > far field sphere theta 0..180, phi 0..355, 3 deg.

Edge power: cylinder r 43.5 mm, h 20 mm, base z = 0, lateral face = sheet `EDGE`; Calculator: Re(E x H*)/2, integrate over EDGE (outward).

**Export, per board** (fill the table, send the numbers; no Touchstone, no field files):
| quantity | where |
|---|---|
| S11 (dB) at 10.05 GHz | Results > S Parameter |
| accepted power (W) | Fields Calculator or solution data |
| EDGE flux / accepted power (%) | Calculator, as above |
| radiation efficiency | Radiation > Antenna parameters |
| gain, RHCP and LHCP, at theta = 0 (dBi) and peak RHCP (dBi, theta, phi) | far-field report, Realized Gain |
| axial ratio at theta = 0 (dB) | far-field report |
| passes / final delta S | Convergence |

**How to read it** (ACTIVE edge fraction): ~40 % or less = direct MoM supercell is right (cell leaks as modelled); ~80 % = HFSS S1 PML is right (leaks ~4-5x less); CONTROL edge ~100 % and low broadside gain = no hidden direct radiation. Gain upper bound for ACTIVE: 17.1 dBi.
