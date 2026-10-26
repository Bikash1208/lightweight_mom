# hex150_tensor_cp_single_v2_pad

V2 of `hex150_tensor_cp_single`: **one change** -- a graded launch pad around the pin. Everything else (hexagon, 1366 hologram cells with
r >= 30 mm: slot sizes and angles, board, pin 6.0 mm above the surface = `mono_h` 7.575 mm / `r_pin` 0.635 mm, setup, far-field sphere) is
the current run's, so the HFSS comparison is clean. Predictions below were written before any HFSS run (MODEL_ONLY, HFSS-anchored to the current run).

## What changed
* Cells with centre radius r < 30 mm (64 hologram cells, the feed-cell keep-out of r < 6.5 mm included) are replaced by **68 isotropic cross-slotted
  circles** (D 6.0 mm, two orthogonal slots of length 4.8 mm, width W(r)); the feed cell (r < 3.2 mm) is empty (pin hole). No leakage design inside r < 30 mm.
* Loading grade (mean sheet admittance at kt = 1.06 k0, `results/phase26/v2_pad_design.json`): 7.55 mS for r <= 15 mm
  (W = 2.0 mm), linear in r down to the hologram's local mean 6.90 mS at r = 30 mm (W = 1.31 mm,
  rounded to 0.05 mm per cell).
* **Why cross slots.** The family's isotropic loading is limited: an unslotted D 6.0 circle gives only 4.4 mS -- LESS than the hologram's slotted
  cells (5.4-7.8 mS) -- so a plain-circle pad would unload the pin's surroundings. The cross-slotted circle (bridge >= 0.43 mm like the
  hologram's) reaches 7.55 mS isotropically. The contrast to the hologram mean is therefore small (7.55 vs 6.9 mS): the expected benefit is
  mainly **filling the bare keep-out with loaded cells right up to the pin**, not a strongly higher loading. (The earlier proposal of a 12-16 mS pad
  of unslotted circles is withdrawn: that loading is not available in this family at mesh-stable bridges.)
* The pad's own beta only adds a constant phase to the surface wave beyond r = 30 mm (circular pad): the hologram's slot angles stay valid.

## Predictions (10.5 GHz, same quantities as the current run)
The model cannot resolve a finite pad (its sheet is infinite), so the launch fraction is bracketed: current value 0.627 (bare, calibrated
feed model) ... 0.88 (the infinite-sheet value at the pad loading). Conversion of the arriving power falls 0.54 -> 0.51 (no radiation inside 30 mm); the aperture
directivity changes by -0.08 dB (hole in the aperture), main-beam share of the aperture power 0.95 -> 0.85, RHCP SLL -22.8 -> -19.1 dB.

| Quantity | Current run (HFSS) | V2, launch unchanged (lower) | V2 central | V2, infinite-sheet launch (upper) |
|---|---|---|---|---|
| Launch fraction (model) | 0.62 | 0.63 | 0.75 | 0.88 |
| **Peak RHCP gain [dBi]** | 23.03 | **22.7** | **23.5** | **24.1** |
| Main beam (<= 8 deg) share of radiated power | 0.271 | 0.23 | 0.27 | 0.32 |
| **Grazing + back share** | 0.533 | 0.56 | 0.53 | 0.50 |
| **LHCP share (all angles)** | 0.310 | 0.33 | 0.31 | 0.29 |
| XPD at peak [dB] | 23.0 | ~ 23 (model: 16.8 -> 16.8) | | |
| AR at peak [dB] | 1.2 | <= 2 (model unchanged) | | |
| S11 at 10.54 GHz (model, calibrated feed) | -15.1 | -13.5 est. | -12.5 est. | **-11.5 est.** (R rises 45 -> 59 ohm) |

* **S11:** the calibrated thin-wire model gives -13.5 dB for the current pin (HFSS -15.1) and degrades to -11.5 dB with the fully loaded surface
  (resistance up, reactance unchanged): expect about -10 ... -13 dB; the pin was deliberately not retuned (one change). Resonant retuning of the pin
  to ~4.8 mm above the surface is the known remedy (est. -17 ... -28 dB at the same launch) if the match is the problem.
* Grazing + back and LHCP shares scale with the non-beam power at the HFSS split of the current run (0.73 and 0.43 of the non-beam power): an
  ASSUMPTION; if the pin's direct radiation is more LHCP-rich than the edge power, LHCP falls faster than shown.
* **Possible loss, stated:** the lower column is a *loss* of 0.4 dB (hole in the aperture, 3 % less conversion, no launch gain): V2 is not guaranteed
  to help. The benefit rests on the pad raising the launch fraction; the E-field line of the current run will bound that (pin-direct vs edge split).
* What would falsify the premise: gain within +-0.3 dB of 23.0 dBi with unchanged grazing + back (the loading near the pin does not change the
  direct radiation). What would confirm it: grazing + back <= 0.50 and gain >= 23.5 dBi.

## Files
`hfss_build_board.py` (AEDT script, exact geometry; each cross cell = circle minus two rotated rectangles), `layout.json/.csv/.dxf/.png` (`cross_slot`
column), `predictions.json` (the full design JSON incl. the loading curve), `board.step` (not committed). Comparison:
`python scripts/research/phase26_compare_single_hex.py --gain <csv> --params <txt> --s11 <csv> --ar <csv> --pin <label>`; the V1 fallback
(keep-out filled, no cross cells) is not exported.
