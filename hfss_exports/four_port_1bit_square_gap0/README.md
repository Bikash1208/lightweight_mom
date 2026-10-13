# 4-port 1-bit coded HMSA, square cells only (gap 0 mm)

10.54 GHz. Rogers RT/duroid 5880, 1.575 mm, PEC ground. 6.4 mm cells. The board is 281.6 × 281.6 mm, centred at (0, 0).
It is a 2×2 of 140.8 mm sub-apertures (22 × 22 whole cells each), 0 mm apart, with one monopole at each sub-aperture centre.

**What changed from `four_port_1bit_gap0/`:**
- **Cells.** Only plain squares now: state 0 is empty and state 1 is a 5.85 × 5.85 mm copper square. There are no annular slots.
- **Sub-aperture size.** 140.8 mm instead of 142.2 mm. 142.2 mm is not a whole number of 6.4 mm cells, so the edges fell off-grid; 22 cells is the nearest exact size.
- **Feeds.** The 4 cells touching each feed point are kept empty so the monopole pin never touches a patch.
- **Unchanged:** beams, polarisation, gap, and the synthesis method (coded 1-bit, model β).

| Port | Monopole (x, y) mm | Target beam | Predicted beam (θ, φ) | Directivity dBi | Cross-pol (LHCP) dB |
|---|---|---|---|---|---|
| 1 | (70.4, 70.4) | (30°, 45°) | 30.0°, 44.4° | 21.5 | -18.3 |
| 2 | (-70.4, 70.4) | (30°, 135°) | 29.0°, 133.2° | 21.3 | -17.1 |
| 3 | (-70.4, -70.4) | (30°, 225°) | 29.0°, 226.8° | 21.4 | -17.7 |
| 4 | (70.4, -70.4) | (30°, 315°) | 29.0°, 315.6° | 21.5 | -16.6 |

All beams are RHCP (IEEE, +z). In `port_patterns.png` the label "LHCP power" is the solver's internal E_L basis, which is IEEE RHCP.

## Cell pair
These are MoM values at the patch plane, mesh-converged with 0.2/0.1 mm Richardson extrapolation:

| State | Cell | Phase |
|---|---|---|
| 0 | Empty (bare grounded substrate) | +138.1° (exact) |
| 1 | 5.85 mm square | ≈ −42° |

The two states are ≈ 180° apart. The square sits close to resonance: the phase changes by about 2.5° per 0.01 mm of side length, and the remaining mesh uncertainty is about ±5°. Run one HFSS unit cell before the full board:
- Master/slave boundaries, 6.4 mm period, Floquet port, PEC ground, 10.54 GHz.
- Compare the 5.85 mm square with the empty cell.
- If the separation is far from 180°, nudge L in 0.05 mm steps. 5.80 mm gives about 166°; 5.90 mm moves past 180°.

## Files
- **`hfss_build_board.py`** (recommended): a native AEDT script. In HFSS, run Tools > Run Script. It builds the substrate, ground, all patches (exact mm values, united into one body "Patches" with Perfect E), the monopoles with 50 Ω coax lumped ports, the air box with a radiation boundary, a 10.54 GHz setup, a 10–11 GHz sweep and the far-field sphere. There is no import step, so there are no 1.99999 values.
  - Design variables: `mono_h` (default 7.1 mm, total pin height from the ground; set it to your reference monopole), `r_pin` 0.635 mm, `r_out` 2.05 mm, `pad_xy`/`pad_bot` 15 mm, `pad_top` 30 mm.
  - At the top of the script, set `BUILD_MONOPOLES = False` for geometry only.
- **`board.step`**: the same geometry as STEP, with bodies Substrate, Ground and Patches. Every vertex is written at its exact coordinate (no transforms); the check found no values like 1.9999 or 2.0000001, and every point lies on a 0.001 mm grid. Assign Perfect E to Ground and Patches, set the substrate material to RT/duroid 5880, and add your own feeds at the coordinates above.
- **`layout.csv`**: every patch (port, row, column, x/y min/max in mm).
- **`design_report.json`**: the full numbers, including the cell pair, β and the isolation-gap flags.
- **`board_layout.png`** and **`port_patterns.png`**: the layout and the predicted pattern for each port.

## HFSS run
1. Analyse all ports in one solve. S21, S31 and S41 give the isolation.
2. Look at the far field one port at a time (Edit Sources).
3. Compare with `four_port_1bit_square_gap28`. The prediction is that isolation differs by only about 1–2 dB, because a bare gap does not stop the TM0 wave.
