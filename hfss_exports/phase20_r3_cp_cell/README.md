# phase20_r3_cp_cell -- HFSS check of the CP-2 octagon cell with MACOM MADP-000907-14020 diodes (Phase 20 R3)

**State: PREDICTION_NOT_HFSS_VALIDATED (MODEL_ONLY), diode values PROVISIONAL_DATASHEET, VIA_OMITTED_APPROXIMATION.** Nothing here has an HFSS number yet; this package
is meant to produce one. It is a **candidate cell, not a design**: the geometry is the optimiser result `results/phase20_r3/OPTIMISATION.md` (CP-2, finalist A), which meets the
phase / amplitude / tolerance targets only within a model uncertainty of order 10 deg and **misses the opposite-handed target** (-13.6 dB at normal incidence against -15 dB).
This same package is also the check for the **LP interim** (LP as 1-bit codes of the CP-2 cell): the Floquet port gives TE and TM, so the linear co/cross numbers come out of the
same eight solves; there is no separate LP geometry (2-bit LP is deferred, see `docs/STATUS.md`).

Diode (MADP-000907-14020, MACOM datasheet Rev V8, DC-0004063, p. 2): ON R_s = 5.2 ohm (typ., 10 mA, 1 GHz), OFF C_T = 0.025 pF (typ., -10 V, 1 MHz; max 0.030 pF = +20 %).
**L_s is NOT_AVAILABLE_FROM_CURRENT_DATA**: the nominal model has L_s = 0 (an assumption); the tolerance file adds L_s = +0.1 nH (ASSUMED_SENSITIVITY_NOT_DEVICE_DATA). The part is a
**flip-chip** device: pads, solder and the real assembly geometry are NOT modelled (the diode is one lumped boundary), and flip-chip assembly is not yet confirmed with fabrication.
Because L_s = 0 each diode is **one single-type boundary** (ON: R only; OFF: C only): the parallel/series question of HFSS Lumped RLC boundaries (R2 package, item 0) does not
enter this package.

## The cell
Period 15.0 x 15.0 mm; RT5880 eps_r 2.2, tan d 0.0009, h 1.575 mm; PEC ground under the substrate. Metal: a regular **octagon** (circumradius 4.601 mm, flat edges
facing 0/45/90/135 deg ...) plus four radial stubs at **0, 45, 90, 135 deg** (width 1.375 mm), each = 0.1 mm overlap into the patch | 0.25 mm pad | **0.75 mm diode span** |
0.25 mm pad | 1.748 mm line. The centre via is omitted (VIA_OMITTED_APPROXIMATION); no bias network. `cell.step` parts (see `step_parts.json`): `Substrate`, `Metal`
(35 um slab: the octagon and the four stubs, five solids), `Ground`, and **`Load_1` ... `Load_4`**: single sheets bridging the diode span of each stub (a rectangle `span x stubw`
rotated with the stub; current line ALONG the stub axis):

| sheet | stub | centre (x, y) [mm] | span x width [mm] | current-line angle [deg from +x] |
|---|---|---|---|---|
| Load_1 | stub at 0 deg | (+4.876, +0.000) | 0.75 x 1.37 | 0 |
| Load_2 | stub at 45 deg | (+3.448, +3.448) | 0.75 x 1.37 | 45 |
| Load_3 | stub at 90 deg | (+0.000, +4.876) | 0.75 x 1.37 | 90 |
| Load_4 | stub at 135 deg | (-3.448, +3.448) | 0.75 x 1.37 | 135 |

## HFSS setup
1. Import `cell.step`. **Perfect E** on `Metal` and `Ground` (the MoM metal has zero thickness; state whether you replaced the 35 um slab by sheets at z = 0).
2. Master/slave on the four side walls, scan angle **theta = 0 and 30 deg**, phi = 0 (plane of incidence xz); ONE **Floquet port** on top of an air box (>= lambda/4 at 8 GHz), modes 1 and 2
   (TE: E along y, TM: E along x at phi = 0). **De-embed** to the pattern plane z = 0.
3. Diode boundaries: **Lumped RLC** on each load sheet, current line along the stub (table above), ONE element type each:
   * diode ON : R = 5.2 ohm (L, C unchecked); diode OFF : C = 0.025 pF (R, L unchecked).
   * **State k** (k = 1..4) = `Load_k` ON and the other three OFF (one-hot).  The four states are the same cell turned by 0 / 45 / 90 / 135 deg in the Pancharatnam-Berry sense.
4. Sweep 8-12 GHz (0.5 GHz steps or interpolating), adaptive at 10 GHz, delta S <= 0.01. Export Floquet-mode complex S11 (modes 1, 2: the four entries TE/TM co and cross), 4 states x 2 angles
   = 8 solves. Optional (tolerance, 10 GHz only, 4 states): C_T = 0.030 pF (datasheet max, +20 %), C_T = 0.020 pF (-20 %), and ON with an extra L_s = 0.1 nH (second element, series).

## What to compare (`mom_reference.csv`, `mom_tolerance_10GHz.csv`)
* Columns: `quantity` is `S_TE_TE`, `S_TM_TM`, `S_TE_TM` (TE received, TM incident), `S_TM_TE`, and the circular terms `r_RR`, `r_RL`, `r_LR`, `r_LL`
  (reflected row, incident column; **propagation-relative handedness**, `docs/theory/CONVENTIONS.md` "Circular basis for the RIS cells": incident R = (TM + jTE)/sqrt2, reflected R = (TM - jTE)/sqrt2
  on the power-wave amplitudes). A metal plate has r_LR = r_RL = -1 and r_RR = r_LL = 0; the PB cell wants r_RR = exp(+j 2 alpha) with alpha = state turn.
* `mom_*` = MoM on the 90 x 90 mesh (pitch 0.167 mm); `mesh60_*` the 60 x 60 mesh; `mesh_delta_*` their difference. **The mesh uncertainty is of order 10 deg**: the octagon diagonals are staircased
  and the convergence order was not usable (see `OPTIMISATION.md`); at 10 GHz the 60 -> 90 phase change of r_RR is up to 6.1 deg here. Do not read a smaller HFSS-MoM difference as agreement to better than that.
* **Convention first:** run the bare grounded slab (no metal, no loads; `bare_grounded_slab_reference.csv`, analytic) at theta = 0, 30 deg and apply any 180 deg TM sign shift HFSS shows to every TM entry
  before comparing; then build r_RR etc. from the corrected TE/TM entries with the formula above. Report what you saw.
* Predictions at 10 GHz, normal incidence, states 1-4: **r_RR phase** -107 / -13 / 87 / 175 deg (a rotating ladder, nominally 90 deg steps), **|r_RR|** -0.90 / -1.05 / -1.10 / -0.90 dB, **opposite-handed max(|r_LR|, |r_RL|)** -16.5 / -15.0 / -13.5 / -19.1 dB.
  The opposite-handed term carries no PB phase (it is the same for every state), so in an array it goes to the specular lobe, not the steered beam (`results/phase20_r3/LOCALITY.md`).
* What would be news: an r_RR ladder step differing from ~90 deg by more than ~20 deg (the model uncertainty is ~10 deg), |r_RR| below -2 dB (diode loss: R_s = 5.2 ohm is the dominant loss), or an
  opposite-handed level much worse than -13 dB. The expected bandwidth is limited by the diode loss: |r_RR| >= -1.5 dB over about 10.0-10.75 GHz (MoM grid), the ladder steps staying within 90 +- 20 deg across 8.5-11.5 GHz.
* LP interim (`results/phase20_r3/lp_interim.json`): with x-polarised incidence the co-LP and cross-LP terms of the four states follow from the same S11 (TM = x at phi = 0).
* Honest limits: zero-thickness metal; staircased octagon in the MoM vs exact in HFSS; omitted via and bias network; lumped diode without flip-chip pad geometry; L_s assumed 0.
