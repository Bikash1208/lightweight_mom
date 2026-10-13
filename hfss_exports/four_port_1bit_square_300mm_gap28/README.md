# 4-port 1-bit coded HMSA, square cells only, 300.8 mm sub-apertures (gap 28 mm)

10.54 GHz. Rogers RT/duroid 5880, 1.575 mm, PEC ground. 6.4 mm cells. The board is 629.6 × 629.6 mm, centred at (0, 0).
It is a 2×2 of 300.8 mm sub-apertures (47 × 47 whole cells each, about 10.6 λ0), 28 mm apart, with one monopole at each sub-aperture centre.

This replaces `four_port_1bit_300mm_gap28/` (same aperture and sub-aperture sizes, gap and beams) with these changes:
- **Cells.** Only plain squares: state 0 is empty and state 1 is a 5.85 × 5.85 mm copper square. There are no annular slots.
- **Feed cell.** Each feed sits at the centre of a cell, and that cell is kept empty.
- **Exact dimensions.** Every coordinate is exact. Use `hfss_build_board.py` (native AEDT) or the cleaned `board.step`.

| Port | Monopole (x, y) mm | Target beam | Predicted beam (θ, φ) | Directivity dBi | Cross-pol (LHCP) dB |
|---|---|---|---|---|---|
| 1 | (164.4, 164.4) | (30°, 45°) | 30.0°, 44.4° | 27.7 | -17.2 |
| 2 | (-164.4, 164.4) | (30°, 135°) | 30.0°, 133.2° | 27.3 | -18.3 |
| 3 | (-164.4, -164.4) | (30°, 225°) | 30.0°, 226.8° | 27.3 | -18.9 |
| 4 | (164.4, -164.4) | (30°, 315°) | 30.0°, 315.6° | 27.7 | -16.8 |

All beams are RHCP (IEEE, +z). In `port_patterns.png` the label "LHCP power" is the solver's internal E_L basis, which is IEEE RHCP.

## Cell pair
These are MoM values at the patch plane, mesh-converged with 0.2/0.1 mm Richardson extrapolation:

| State | Cell | Phase |
|---|---|---|
| 0 | Empty (bare grounded substrate) | +138.1° (exact) |
| 1 | 5.85 mm square | ≈ −42° |

The two states are ≈ 180° apart. The square sits close to resonance: the phase changes by about 2.5° per 0.01 mm of side length, and the remaining mesh uncertainty is about ±5°. Run one HFSS unit cell before the full board:
- Master/slave boundaries, 6.4 mm period, Floquet port, PEC ground, 10.54 GHz.
- If the separation is far from 180°, nudge L in 0.05 mm steps. 5.80 mm gives about 166°.

## Files
- **`hfss_build_board.py`** (recommended): a native AEDT script. In HFSS, run Tools > Run Script. It builds the substrate, ground, all 4416 patches (exact mm values, united into one body "Patches" with Perfect E), the monopoles with 50 Ω coax lumped ports, the air box with a radiation boundary, a 10.54 GHz setup, a 10–11 GHz sweep and the far-field sphere.
  - Design variables: `mono_h` (default 7.1 mm, total pin height from the ground; set it to your monopole), `r_pin` 0.635 mm, `r_out` 2.05 mm, `pad_xy`/`pad_bot` 15 mm, `pad_top` 30 mm.
  - At the top of the script, set `BUILD_MONOPOLES = False` for geometry only.
  - Creating 4416 rectangles takes a few minutes.
- **`board.step`**: the same geometry as STEP, with bodies Substrate (z 0–1.575 mm), Ground (z = 0) and Patches (z = 1.575 mm). Every vertex is written at its exact coordinate (no transforms). The check found no values like 1.9999 or 2.0000001; all points lie on a 0.001 mm grid.
- **`layout.csv`**: every patch (port, row, column, x/y min/max in mm).
- **`design_report.json`**, **`board_layout.png`** and **`port_patterns.png`**.

## HFSS notes
The board is about 22 λ0 across, so a full-wave solve is large. A practical sequence:
1. Run the unit cell first.
2. Solve with the FE-BI or IE-region option, or a mixed-order mesh, with one port excited at a time (Edit Sources).
3. S21, S31 and S41 give the isolation.
4. Compare with `four_port_1bit_square_300mm_gap0`.
