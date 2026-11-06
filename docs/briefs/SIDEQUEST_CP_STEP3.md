# SIDEQUEST_CP_STEP3 - stronger switched GB cell for the CP+LP hybrid, both sign hypotheses (2026-11-06). MODEL_ONLY, no package
Model, tensor background, GB cell amplitude and joint optimisation as in `SIDEQUEST_CP_HYBRID.md`; GB cell = same crossed-slit disc family, free geometry = disc radius r (slit to r - 0.2), MoM pairs OFF (open) / ON (2.44 ohm) from `k2_family_scan_raw*.json` (RT5880 1.575, 10.5 GHz, mesh 64): r 3.0: Y 8.69 / 6.06 mS (step 2.6), r 3.10: 12.12 / 7.52 (4.6), r 3.12: 12.73 / 7.82 (4.9), r 3.15: 16.86 / 8.55 (8.3); r 3.17 flips sign (resonance). Scripts `sidequest_cp_step3.py`, `sidequest_cp_hybrid.py` (`joint_opt2`), data `results/sidequest_hybrid/cp_step3_ar3.json`, `cp_step3_free.json`. Sign hypotheses: -1 (assumed in the package) and +1: each gets its own design (sites + rotated slots), the cell itself is the same.
**Correction found on the way:** the optimiser of the delivered package (`joint_opt`) took `circular_components()[0]` as the beam-1 amplitude; in this project's basis the RHCP co-pol is `E_L` (index 1, `copol_label('rhcp') == 'E_L'`), so its beam-1 floor and the loss budget acted on the cross-pol. The package's reported patterns were computed with the correct `copol_component` and stay valid (22.2 dBi RHCP, AR 2.1 dB), but the design was not optimised against AR. `joint_opt2` (correct index, explicit cross-pol floor XPD >= 16 dB in both GB states, cap on rotated slots) re-optimises the SAME cell (r 3.0, 360 diodes): beam 1 22.5 / 22.3 dBi, **AR 1.0 / 1.0 dB**, contrast 11.7 dB, balance 14.8 dB, SLL -12.2 / -12.0 (with the offset), 38 rotated slots (rel -1). A package v2 with that design is available after your HFSS result.
## Results (180 cells = 360 diodes unless stated; SLL with the +2.5 dB offset; AR peak OFF / ON)
| cell | sign | rotated slots | beam 1 OFF / ON [dBi] | AR OFF / ON | beam 2 OFF / ON [dBi] | contrast | balance (b1 - b2, ON) | SLL OFF / ON |
|---|---|---|---|---|---|---|---|---|
| r 3.0 (step 2.6 mS, as now) | -1 | 38 | 22.5 / 22.3 | 1.0 / 1.0 | -4.1 / 7.5 | 11.7 | 14.8 | -12.2 / -12.0 |
| r 3.0 | +1 | 45 | 22.5 / 22.3 | 0.95 / 0.9 | -3.4 / 7.6 | 11.0 | 14.6 | -11.2 / -11.0 |
| r 3.12 (4.9 mS, disc gap 0.16 mm) | -1 | 80 | 23.2 / 22.0 | 1.0 / 0.25 | -3.0 / 11.7 | 14.6 | 10.3 | -13.0 / -12.2 |
| r 3.12 | +1 | 80 | 23.0 / 22.0 | 1.3 / 0.7 | -1.7 / 11.7 | 13.4 | 10.3 | -13.0 / -13.1 |
| r 3.12, freer (128 slots, loss budget 3 dB) | -1 | 128 | 22.6 / 21.6 | 0.9 / 1.2 | 0.9 / 13.3 | 12.5 | 8.2 | -12.3 / -11.6 |
| r 3.12, freer (160 slots) | +1 | 160 | 22.1 / 21.2 | 1.6 / 1.9 | 3.6 / 14.2 | 10.6 | 7.0 | -8.9 / -11.0 |
| r 3.15 (8.3 mS, disc gap 0.10 mm) | -1 | 79 | 24.2 / 22.1 | 2.55 / 0.2 | -7.2 / 12.2 | 19.4 | 9.9 | -13.8 / -16.1 |
| r 3.15, 131 slots | -1 | 131 | 24.1 / 21.7 | 2.0 / 0.8 | -3.6 / 14.7 | 18.3 | 7.0 | -13.0 / -15.1 |
| r 3.15, 176 slots | -1 | 176 | 23.5 / 21.0 | 0.6 / 1.7 | -0.5 / 16.0 | 16.6 | 5.0 | -8.9 / -12.9 |
| r 3.15, 228 slots | +1 | 228 | 22.5 / 20.3 | 1.0 / 2.6 | 1.2 / 17.1 | 15.9 | **3.2** | -7.5 / -11.3 |
| r 3.15, 270 cells = 540 diodes, 153 slots | -1 | 153 | 23.8 / 20.7 | 2.4 / 0.9 | -4.0 / 15.5 | 19.4 | 5.2 | -9.6 / -14.9 |
| r 3.15, 270 cells = 540 diodes, 205 slots | +1 | 205 | 23.2 / 19.9 | 2.7 / 0.8 | -1.0 / 16.6 | 17.7 | 3.3 | -8.3 / -12.5 |
| r 3.12, 360 cells = 720 diodes | -1 / +1 | 15 / 53 | 22.8 / 21.4 ; 22.6 / 21.1 | 2.4 / 1.0 ; 2.8 / 1.0 | -7.1 / 10.5 ; -4.7 / 11.8 | 17.7 ; 16.6 | 10.9 ; 9.2 | -14.2 / -17.2 ; -11.0 / -16.5 |
## Reading
* **A stronger cell is the lever**: step 2.6 -> 4.9 mS (r 3.12) lowers the balance from 14.8 to 10.3 dB at unchanged diodes (360) and keeps beam 1 at 22-23 dBi, AR <= 1.5 dB, SLL -13 / -12 (both signs behave alike). Step 8.3 mS (r 3.15) gives beam 2 12-17 dBi and contrast 16-19 dB.
* **Balance <= 3 dB is reached only with the aggressive settings** (r 3.15, 205-228 of the 1250 tensor slots rotated by 90 deg, sign +1: balance 3.2-3.3 dB, beam 1 20.3-19.9 dBi in the ON state, beam 2 16.6-17.1 dBi, AR <= 2.7 dB) **and then the OFF-state SLL is only -7.5 / -8.3 dB (target -13)**. Sign -1 reaches 5.0 (360 diodes) / 5.2 dB (540). With SLL <= -13 dB in both states the best balance is about 7 dB (r 3.15, 131 slots, sign -1) or 10 dB (r 3.12, 80 slots).
* **Beam 1 stays >= 17 dBi in every case** (19.9-24.2); beam 1 can be traded down (ON state loses 1.5-3 dB when the balance is pushed).
* **Diodes**: the stronger cell has the SAME topology (2 diodes per cell, 180 cells = 360 diodes); 270 cells (540 diodes) and 360 cells (720) buy 1-2 dB of balance at SLL cost; the diode count is not the lever, the step (mS) is. A step above about 5 mS needs the neighbour gap to shrink below 0.16 mm.
* **AR of beam 1**: <= 1.5 dB at r 3.0-3.12 with <= 80 slots; 2.0-2.8 dB at r 3.15 or >= 130 slots (worst in the -3 dB beam up to 3.9 dB at 270-cell variants): AR grows with the GB OFF-state scatter (r 3.15 OFF is +9 mS against the local mean).
## Risks (model-only)
1. **r 3.15 is a resonant, 0.10 mm-gap cell** (r 3.17 flips the sign of Y): +-0.02 mm in r changes the step by several mS; a gap below 0.15 mm may not be buildable. r 3.12 (gap 0.16, step 4.9) is the practical choice; its balance limit is about 7-10 dB.
2. Local beta of the GB cells rises (r 3.12: Y 12.7 mS against 7.5 mS mean: beta/k0 about 1.12 against 1.06); with 12 % of the cells this shifts the mean beta by about +0.008 (phase +15 deg over the radius, not modelled).
3. 130-230 rotated slots (10-18 % of the tensor cells) mean the Phase 26 hologram is no longer the one HFSS validated; the model ignores the mutual coupling of the rotated cells.
4. The relative sign decides which design is built; both are listed (the cell is identical, only the site list and rotations differ).
## What each HFSS outcome leads to
* Sign -1 confirmed (beam 2 on only with S_w 6.1): build r 3.12, 360 diodes, 80 slots rotated (balance 10.3 dB, SLL -13.0 / -12.2) or r 3.15 with 131 slots (balance 7.0, SLL -13.0 / -15.1) if the 0.10 mm gap is acceptable.
* Sign +1 (beam 2 on with S_w 1e6): build the +1 designs of the same table (r 3.12, 80 slots: balance 10.3, SLL -13.0 / -13.1; aggressive r 3.15, 228 slots: balance 3.2 with SLL -7.5 / -11.3).
* In both cases the balance <= 3 dB with SLL <= -13 dB is NOT reached by this cell family; it needs a step above 8 mS at a buildable gap (another stack / a two-layer cell) - not designed here.
