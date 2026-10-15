# hex4_tensor_cp_static -- 4-hexagon tensor-CP board A (static beams), Phase 26

**State: PREDICTION_NOT_HFSS_VALIDATED (MODEL_ONLY, SOLVER_SUPERCELL cells).** Four 150 mm hexagons on RT/duroid 5880
1.575 mm, layout of `four_port_hex150_gap0` (652.8 x 537.6 mm board), **one shared 6.4 mm cell lattice** (5684 slotted
circles, D 6.0 mm), no boundary row and no gap strip. Every hexagon is the validated single-hexagon design
(`hex150_tensor_cp_single`, HFSS 23.03 dBi) re-aimed: port 1 -> (theta 30, phi 45), port 2 -> (30, 135), port 3 ->
(30, 225), port 4 -> (30, 315), all RHCP. Feed: pin 6.0 mm above the surface (mono_h = 7.575 mm from the ground),
r_pin 0.635 mm, as the single-hexagon run.

## Files
`hfss_build_board.py` (AEDT script, native geometry), `layout.csv/.dxf/.png/.json` (port column = owning hexagon),
`predictions.json`. Raw `.s4p` stays local.

## HFSS setup
* Build with Tools > Run Script. Build time is long (5684 slotted circles, one polyline + subtract each): expect
  hours, several times the single-hexagon build. `ADD_SWEEP = False` solves 10.54 GHz only.
* Boundary and solution as the single run (radiation box 15 mm sides / 30 mm top; Setup1 at 10.54 GHz; far-field
  sphere theta 0-180, phi -180..180, 1 deg).
* **One port driven, the others 50-ohm terminated:** set `DRIVE_PORT = k` (1..4) near the top of the script. It tries
  `EditSources` (an unverified API call, wrapped); if HFSS reports a message instead, use Results > Edit Sources: port k
  1 W / 0 deg, the other three 0 W. The ports are always terminated in 50 ohm (lumped ports); only the far-field
  excitation changes. S-parameters need no source setting.
* Mesh: see the size note in `predictions.json` (`hfss_size`): scale the single-hexagon run by ~4 in cells/tetrahedra
  (4.06x cells, 4.3x board area) and 5-8x in solve time; drive one port at a time for the far fields.

## What to export
1. S11..S44 and the full S matrix, 8-12 GHz (or the sweep range) -- compare the five adjacent-pair couplings
   (1-2, 1-3, 2-3, 1-4, 3-4) with the estimate in `predictions.json`.
2. For each driven port k: RHCP and LHCP realized gain, 1 deg grid (theta 0-90, phi 0-359) at 10.54 GHz and
   the antenna-parameter table; axial ratio on the beam cut.
3. **Check item (explicit): the neighbour LHCP lobe.** Per port the predictions list the level (relative to the
   main RHCP peak) and direction of the largest LHCP lobe of the unconverted wave that crosses into the
   neighbouring hexagons. Read the LHCP gain at those directions and at the mirror directions; state level, direction
   and which neighbour edge faces it.
4. Optional, decisive for the edge model: Mag_E along the shared edge of one driven hexagon, 1 mm above the surface
   (the 4-hex equivalent of the single-hexagon E-line).

## Predictions (port driven alone, anchored to the single-hexagon HFSS result; s = scale of the unconverted edge power)
| port | beam (theta, phi) | model peak | gain dBi (s=1 / s=1.7) | XPD dB | AR peak / worst-in-beam dB | neighbour LHCP lobe (rel. main peak, direction) |
|---|---|---|---|---|---|---|
| 1 | (30, 45) | [29.7, 45.0] | 23.11 / 23.12 | 16.4 / 16.3 | 2.65 / 3.93 (s=1) | -19.9 / -17.6 dB at [30.0, 75.0] (14.9 deg off) |
| 2 | (30, 135) | [29.7, 135.2] | 22.97 / 22.96 | 16.4 / 16.0 | 2.65 / 3.41 (s=1) | -21.7 / -19.4 dB at [21.0, 206.0] (29.9 deg off) |
| 3 | (30, 225) | [29.7, 225.0] | 23.11 / 23.12 | 16.4 / 16.3 | 2.65 / 3.93 (s=1) | -19.9 / -17.6 dB at [30.0, 255.0] (14.9 deg off) |
| 4 | (30, 315) | [29.7, 315.2] | 22.97 / 22.96 | 16.4 / 16.0 | 2.65 / 3.41 (s=1) | -21.7 / -19.4 dB at [21.0, 26.0] (29.9 deg off) |

Edge power of each hexagon (fraction of the launched surface-wave power, design march): see `predictions.json`
`edges`. Coupling estimate (model, uncalibrated and bias-corrected from the 1-bit board): see `coupling`.

## Honest limits
The model's XPD/AR are pessimistic against the single-hexagon HFSS (model 16.9 dB / 2.5 dB vs HFSS 23 dB / 1.3 dB); the
pickup uses the cells' own designed alpha for the oblique/counter wave (ASSUMPTION); the unconverted split between the
free edge and the shared edges follows the design march; the neighbour interaction is not an HFSS result until this
board is run.
