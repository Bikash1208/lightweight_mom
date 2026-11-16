# RUN_README - BOARD_FINAL_SPREAD_ELS.aedt
**Status: built, UNSOLVED, clean.** `BOARD_FINAL_ELS.aedt` in the same folder is the older coherence-layout project and is **NOT FINAL**; do not solve it.

## Run order (one full-board solve at a time)
1. `ON_ELS_ABS` = TWO_BEAM_ABS (diodes forward, 10 mA): A at theta 25 / phi 0, B at theta 25 / phi 180 (rot180 frame in the model; check which beam sits at which phi in the far field).
2. `OFF_ELS10_ABS` = ONE_BEAM_ABS (-10 V): beam A only.
3. Optional references without the absorber: `ON_ELS`, `OFF_ELS10`.
Solver: direct (Setup1). Expect about 60 GB at pass 1, growing; stop and switch to iterative if memory passes 120 GB.

## Gates, in order
1. Convergence: DeltaS <= 0.02 with at least 6 passes.
2. Accepted power >= 50 % at 5.8 GHz.
3. Outputs: RHCP and LHCP at theta 24 / 26 deg, phi 0 / 180 deg; AR at the beam peaks (target <= 3 dB); XPD >= 15 dB; SLL <= -15 dB in both states (definition S3b: cones max(1.5 HPBW, 10 deg); both beams excluded in the two-beam state, A only in the one-beam state); A-B balance within 3 dB in the two-beam state; B suppression in the one-beam state; power below the board (`Lower` sphere); export .ffd.

## Model expectation (MODEL_ONLY, offsets to HFSS +10.2 dB OFF / +7.9 dB ON on SLL, directivity over-predicted by 5.7-7.2 dB)
Balance +1.6 dB; one-beam SLL about -13.9 dB; two-beam SLL about -16.6 dB. See `readback.md` and `checks.md`.

## Files
STEP `BOARD_FINAL.step`, `manifest.json`, `diode_sheets.json` (note: lists old sheet values; the project uses OFF_ELS10 / ON_ELS), `layout_FINAL2.json`, `readback.md`, `checks.md`.
