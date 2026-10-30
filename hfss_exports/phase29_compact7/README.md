# phase29_compact7: seven hexagonal tiles side by side, one monopole each (Phase 29, compact option)

Seven R 40 mm hexagons (vertex at +/-y, flat edges shared, no gap) around a centre tile, one centre monopole + lumped 50 ohm port per tile (7 ports), rotating-slot
tensor-CP cells (circle D 6 mm, slot L x W and angle per cell) on one shared 6.4 mm lattice, RT/duroid 5880 1.575 mm, 10.54 GHz, 676 cells. Port 1 = centre tile (broadside beam);
ports 2..7 = outer tiles at azimuth 0, 60, ..., 300 deg (centre of the tile seen from the board centre), each beam at theta 30 deg pointing outward. Continuous scanning = complex weights on
2-3 ports at a time (`weights_table.csv`), theta 0-45 deg in 3 deg steps, phi in 15 deg steps. **MODEL_ONLY / NOT_HFSS_VALIDATED**: all numbers below are predictions; inter-tile
coupling is NOT modelled; gain is anchored only to the single-hexagon HFSS run (23.03 dBi at 150 mm), so tile size and tilt are extrapolations; S11 and port coupling are
NOT_AVAILABLE_FROM_CURRENT_DATA.

## Files
| File | Content |
|---|---|
| `predictions.json` | everything below, machine readable (per port alone, three weighted examples, scan statistics, check items) |
| `pattern_cuts.csv` | RHCP/LHCP gain (dBi) vs signed theta (-90..90 deg) in the plane of each port's beam (`portK_alone`) and of the three weighted examples (`weighted_*`) |
| `weights_table.csv` | per direction: for ports 1..7 the incident power (W, sum = 1 W) and phase (deg), method, active ports, predicted peak gain / direction / pointing error |
| `build_board.py` | AEDT (IronPython) script: substrate, ground, 7 pins + lumped ports, 676 slotted circles, radiation box, setup 10.54 GHz, far-field sphere |
| `set_weights.py` | AEDT script: sets the port sources (EditSources) from one table row (edit THETA_DEG / PHI_DEG at the top); prints the values if the call fails |
| `layout.csv`, `layout.json`, `layout.png` | every cell (centre, slot, tile/port), tile hexagons, feeds |

## Run order
1. AEDT: Tools > Run Script > `build_board.py` (a few minutes). Check the messages; it saves `phase29_compact7.aedt`. `ADD_SWEEP = True` at the top adds a 10.0-11.1 GHz sweep (slower; not needed for the weights).
2. Validate, then Analyze (Setup1, 10.54 GHz).
3. **Read the S-matrix first**: S11..S77 and the coupling between neighbouring ports (centre-outer, outer-outer). If a neighbouring |S_ij| is above about -15 dB the uncoupled weights are unreliable (note it in the comparison).
4. Port alone: Results > Edit Sources: port k = 1 W, all others 0 W, then the far field (Radiation > Far Fields, Realized Gain, RHCP/LHCP, Axial Ratio). Compare with `pattern_cuts.csv` (`portK_alone`).
5. Weighted cases: run `set_weights.py` after setting THETA_DEG / PHI_DEG (or type the row of `weights_table.csv` by hand), re-run the far field; compare with the weighted examples below.

## How to apply the port weights in HFSS (Edit Sources)
Results (project tree) > Edit Sources... > type **Power** (not Voltage), then for each Port1..Port7: **Magnitude = `P<k>_power_W` from the table in W**, **Phase = `P<k>_phase_deg` in deg**; ports with 0 W stay terminated in 50 ohm. The table's weights are field weights
w_k with |w_k|^2 = incident power (sum = 1 W) and phase = angle(w_k), exp(+j w t) (the project convention, also HFSS's); the phase reference of all seven ports is the identical port plane, so a common offset cancels.
If the beams come out mirrored in phi, the phase sign convention differs: negate all phases. The model assumes uncoupled ports: with strong coupling the realised excitation differs (read the HFSS "port incident / accepted" values). `set_weights.py` automates this
(EditSources call NOT verified against an AEDT version; wrapped, it prints the values to type if it fails).

## Predictions (MODEL_ONLY, 10.54 GHz)
Port alone (1 W accepted; RHCP realised gain; SLL = highest RHCP outside 1.5 x HPBW about the peak, relative):

| Port | Feed (mm) | Design beam (theta, phi) | Peak gain (dBi) | Peak (theta, phi) | HPBW (deg) | SLL (dB) | AR peak (dB) |
|---|---|---|---|---|---|---|---|
| P1 | (0.0, 0.0) | (0, 0) | 7.2 | (0, 205) | 23 | -19.8 | 0.7 |
| P2 | (69.3, 0.0) | (30, 0) | 7.9 | (26, 0) | 23 | -15.4 | 2.8 |
| P3 | (34.6, 60.0) | (30, 60) | 7.5 | (26, 60) | 24 | -14.7 | 3.1 |
| P4 | (-34.6, 60.0) | (30, 120) | 7.5 | (26, 119) | 24 | -15.0 | 0.8 |
| P5 | (-69.3, 0.0) | (30, 180) | 7.9 | (26, 180) | 23 | -15.4 | 2.8 |
| P6 | (-34.6, -60.0) | (30, 240) | 7.5 | (26, 240) | 24 | -14.7 | 3.1 |
| P7 | (34.6, -60.0) | (30, 300) | 7.5 | (26, 299) | 24 | -15.0 | 0.8 |

Weighted examples (sum of port powers = 1 W):

| Case | Active ports | Peak gain (dBi) | Peak (theta, phi) | Pointing error (deg) | HPBW (deg) | AR peak (dB) |
|---|---|---|---|---|---|---|
| broadside | [1, 2, 5] | 7.3 | (0, 244) | 0.5 | 23 | 0.7 |
| theta15_phi0 | [1, 2, 4] | 6.8 | (16, 1) | 1.0 | 24 | 1.3 |
| theta30_phi30 | [2, 3, 4] | 5.9 | (26, 28) | 4.3 | 20 | 2.1 |

Scan over the table (tuned weights, <= 3 active ports): theta <= 30 deg: worst pointing error 5.0 deg (median 1.5), gain dip 2.1 dB, minimum peak gain 5.9 dBi
(criterion: pointing <= 8 deg, dip <= 3 dB: MET). 30-45 deg band, as is: worst pointing error 21.3 deg, only 15 % of directions within 8 deg,
gain dip 2.6 dB, minimum peak gain 5.4 dBi. Gains are LOW by design (small tiles, weak slot leakage: about 6-8 dBi). The tilted tiles peak near theta 26 deg (not 30) in the model, so the theta 30 deg rows of the table are reached by the combination of ports, not by one tile.

What would falsify the model: a port alone whose peak is more than ~10 deg off its design direction or more than 5 dB from the table; a weighted example whose peak is more than ~15 deg off the target; neighbouring-port coupling above -10 dB.
The weights are a starting point: with coupling in HFSS they should be re-tuned from the measured port patterns (the 7 single-port patterns are all that is needed, the model is linear).

## STEP import (compact7_board.step, MODEL_ONLY)

`compact7_board.step` (gitignored; copy with `ports_and_excitations.csv` in `C:\hmsa_sync\data\phase29_compact7\`), generated by `scripts/research/phase29_compact7_step.py` from `build_board.py` (676 cells, 7 pins). Bodies: Substrate (solid, 1.575 mm, pin holes), Ground (sheet z = 0, coax opening r 2.05 mm per feed), Metal_pattern (676 slotted circles, D 6.0 mm, sheets at z = 1.575), Pin_1..7 (r 0.635 mm, z 0..7.575 mm), Port_1..7 (annulus sheets r 0.635..2.05 mm at z = 0). Setup: the STEP section of `../phase29_twofeed/README.md` (ports 1..7 instead of 2), 10.54 GHz single frequency. Excitation per port: `ports_and_excitations.csv`.
