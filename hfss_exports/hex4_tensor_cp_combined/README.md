# hex4_tensor_cp_combined -- 4-hexagon tensor-CP board B (combined beam), Phase 26

**State: PREDICTION_NOT_HFSS_VALIDATED (MODEL_ONLY, SOLVER_SUPERCELL cells).** Generated 2026-10-21 by the advising
session from the committed Phase 26 design `results/phase26/fourhex/fourhex_B_cells.json` (no re-design; same cells the
session's scope used). Board B is the former "Butler-ready" board: Butler scanning on this layout is **shelved**
(quantisation lobes, see `results/phase26/README.md` P); board B is used as the **combined-beam board** — all four
ports driven with equal amplitude and phase through a 1:4 in-phase divider.

Same board as board A (`hex4_tensor_cp_static`): four 150 mm hexagons on RT/duroid 5880 1.575 mm, 652.8 x 537.6 mm,
one shared 6.4 mm lattice, 5684 slotted circles (D 6.0 mm), pins 6.0 mm above the surface (mono_h 7.575 mm), r_pin
0.635 mm, lumped 50-ohm ports. Only the slot angles differ: every hexagon aims at (theta 30, phi 225) RHCP with a
**global phase reference at the board centre**, so equal-phase drive gives one continuous aperture.

## Files
* `hfss_build_board.py` — AEDT script (native geometry, same builder as board A). `DRIVE_PORT = 5` (default here) sets
  all four ports to 0.25 W / 0 deg (combined beam); 1..4 = one port alone. The EditSources call is unverified (wrapped):
  check Results > Edit Sources after the solve.
* `board.step` — the same geometry as STEP (AP214, mm; parts: Substrate_RT5880_1p575 with pin holes, Ground with
  r_out = 2.05 mm coax openings, Patches_boardB (5684 zero-thickness slotted-circle faces at z = 1.575 mm),
  Monopole1..4). 27.9 MB; large, keep out of git. Verified by reading back: 5684 patch faces with the expected areas per
  slot size (4.8x1.6: 3460, 4.6x1.6: 1370, 4.8x2.0: 751, 4.4x1.2: 103), extent 652.8 x 537.6 mm, z 0..7.575 mm.
  For HFSS prefer the AEDT script (native, exact); the STEP is for other tools / fabrication drawings. If you import the
  STEP into HFSS: assign Perfect E to Ground and Patches, RT5880 to the substrate, PEC to the pins, and draw the four
  lumped-port annuli (r_pin..r_out) at the feeds yourself.
* `layout.csv` — x, y, D, slot L, W, angle, owning hexagon (port).

## Predictions (model, s = 1, anchored to the single-hexagon HFSS run; 10.54 GHz)
| drive | beam (theta, phi) | RHCP gain | XPD | AR peak / worst in beam | first lobes |
|---|---|---|---|---|---|
| all four, 0.25 W, 0 deg | (29.9, 225.0) | 28.6 dBi (28.4 at s = 1.7) | 16.9 dB (model; pessimistic vs HFSS) | 2.5 / 3.0 dB | about -11 dB at 5-6 deg (tile period); lobes beyond 20 deg <= -26 dB; strongest LHCP lobe -20.6 dB at (30, 45) |
Each port alone: about 22.8-22.9 dBi toward (29-30, 223-225). Gain over one port: +5.5 dB (ideal +6.0).
Quiet-zone variant (-13.2 dB first lobes at -1.2 dB gain) exists only as a study, not exported.

## What to export from HFSS
S-matrix 8-12 GHz; with DRIVE_PORT = 5: RHCP/LHCP realized gain 1 deg grid at 10.54 GHz, antenna parameters, AR on the
phi = 225 cut. Optionally port 3 alone (DRIVE_PORT = 3) to compare with board A's port 3.
