# phase28_cell_G -- HFSS check of the cell G, all four states, Floquet port

**State: PREDICTION_NOT_HFSS_VALIDATED (MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET, VIA_OMITTED_APPROXIMATION).** Nothing here has an HFSS number yet. Diode: MADP-000907-14020, ON R = 5.2 ohm, OFF C = 0.025 pF, **package inductance L_s = 0 nH assumed (NOT_AVAILABLE_FROM_CURRENT_DATA)**; flip-chip pads, via and bias network are not modelled.

Cell G (Phase 28 item 3 candidate): period 6.4 mm, RT/duroid 5880 (eps_r 2.2, tan d 0.0009) 1.575 mm on a PEC ground; zero-thickness metal at the top surface: a 2.4 mm hub, four arms, four diagonal pieces, and **eight diode necks** (0.4 mm long) bridging hub-arm and hub-diagonal gaps. State k = the two necks of one pair ON (R), the other six OFF (C): axis 0 / 45 / 90 / 135 deg for k = 0 .. 3. The polygons are generated from the same geometry code as the MoM (`phase28_g_geometry.py`; tested against the MoM mask).

## Files
| File | Content |
|---|---|
| `build_cell.py` | AEDT script: builds the cell (state `STATE` at the top), 8 Lumped RLC sheets, periodic walls and a 2-mode Floquet port with scan angle variables `scan_theta` / `scan_phi`, setup 10 GHz, discrete sweep 9.5 - 10.5 GHz (0.25 GHz) |
| `set_state.py` | AEDT script: puts another state on the open model (edit `STATE`) |
| `mom_reference.csv` | MoM prediction: S_TM_TM, S_TE_TE, S_TM_TE, S_TE_TM and r_RR, r_RL, r_LR, r_LL (magnitude dB, phase deg) per state, f = 9.5 .. 10.5 GHz; theta = 0 and theta = 30 deg / phi = 225 deg (10.0 GHz); mesh 64 and 96 |
| `bare_grounded_slab_reference.csv` | analytic reflection of the bare grounded slab at the pattern plane (convention check) |

## Predictions (MoM, mesh 96, normal incidence, 10.0 GHz, pattern plane = top of the substrate)
| state | axis | r_RR [dB] | r_RR phase [deg] | max(opposite-handed) [dB] | r_LL [dB] |
|---|---|---|---|---|---|
| 0 | 0 deg | -3.30 | -27.5 | -6.27 | -3.30 |
| 1 | 45 deg | -3.53 | 69.2 | -4.78 | -3.53 |
| 2 | 90 deg | -3.30 | 152.5 | -6.27 | -3.30 |
| 3 | 135 deg | -3.53 | -110.8 | -4.78 | -3.53 |

Phase steps between neighbouring states (r_RR): 96.8, 83.2, 96.8 deg (a Pancharatnam-Berry cell turned by 45 deg steps wants 90). The tensor-sheet model behind the hologram design uses b_bar, delta and the axis of the same solves (`../../results/phase28/f10/g_library/G_LIBRARY.md`).
Mesh convergence 64 -> 96 over 9.5 - 10.5 GHz: all five frequencies converged (5 % / 10 % criterion).

Oblique prediction (theta 30 deg, phi 225 deg, TE / TM, mesh 64, 10.0 GHz) is in `mom_reference.csv` (rows with theta = 30): this is the incidence of the beam the hologram radiates to; the sheet tensors of the design are extracted at normal incidence only.

## HFSS setup
1. Run `build_cell.py` (Tools > Run Script). It builds state 0. Read the messages: three assignments are NOT yet verified in AEDT (periodic pairs with scan angle, Floquet port, and in `set_state.py` `EditLumpedRLC`); each falls back to a message telling you what to assign by hand.
2. By hand if needed: master/slave on the four side walls with the scan angle; ONE Floquet port on the top face of `Air`, modes 1 and 2 (TE: E along y, TM: E along x at phi = 0), **de-embed 20 mm** to the pattern plane (top of the substrate); `scan_theta` = 0 deg first, then 30 deg with `scan_phi` = 225 deg.
3. Diode boundaries are Lumped RLC on the red sheets `D0 .. D7`, current line along the neck, **one element type each** (ON: R only; OFF: C only; the AEDT-verified form of the 3x3 smoke build).
4. Solve states 0 .. 3 (`set_state.py` between runs), theta = 0 and 30 deg: 8 solves. Export the complex Floquet S11 (modes 1, 2).

## Compare (no modelling)
* **Convention first:** run the bare grounded slab (delete the metal and the RLC sheets, or use `bare_grounded_slab_reference.csv`) at theta = 0: apply any 180 deg TM sign shift HFSS shows to every TM entry. Circular terms as in `phase20_r3_cp_cell`: incident R = (TM + jTE)/sqrt2, reflected R = (TM - jTE)/sqrt2; a metal plate has r_LR = r_RL = -1.
* **Pass (ASSUMED thresholds, taken from the R3 package):** r_RR phase of every state within 20 deg of `mom_reference.csv`, neighbouring-state steps 90 +- 20 deg, |r_RR| within 1 dB of the prediction. **What would be news:** a step far from 90 deg (the pinwheel cell is not a clean rotation in HFSS), |r_RR| more than 2 dB below the prediction (the table already contains a 3 - 4 dB diode loss: this cell is a surface-wave modulator, not a reflectarray cell, and its opposite-handed terms of -5 to -6 dB are predicted, not a fault), or an opposite-handed term more than 3 dB away from the prediction.
* The package inductance decides this cell: the MoM sensitivity (assumed L_s = 0.1 - 0.3 nH) is in `../../results/phase28/f10/g_refine/`. If HFSS includes any package inductance, compare with that file, not with the L_s = 0 prediction.

## Honest limits
Zero-thickness metal; staircased polygons in the MoM against exact ones in HFSS (self-similar rational slopes, converged to the criterion above); no via, no bias network, no pad geometry; the tensor-sheet extraction assumes the local response; MADP L_s unknown.
