# HFSS gate test board — setup (one page). MODEL_ONLY predictions in `predictions.json`; nothing here is HFSS-validated.

**Purpose.** Measure what the 3-row dog-bone gate really does on RT5880: pass (OFF) and block (ON). The model predicts the ON state blocks only ~3 dB (see `predictions.json`).

**Import** `gate_test_board.step` (mm, centre probe at the origin, +z up). Named parts: `Substrate` (RT5880, eps_r 2.2, tan d 0.0009, z -1.575..0), `Ground` (35 um PEC, z -1.610..-1.575, hole r 2.05),
`Probe_Pin` (PEC, r 0.635, z -1.610..+6.0), `Port1` (sheet, annulus r 0.635..2.05 at z -1.610), `Metal_REF` (83 square patches 5.8 mm = hologram-mean sheet, b = 7.5 mS), `Metal_TRANS`
(20 patches 5.4 mm, b = 5.4 mS), `Metal_GATE` (35 dog-bone cells, 70 solids), sheets `D_ON_01..35` and `D_OFF_01..35` (identical positions, 0.4 x 0.8 mm, z +0.0175). All metal/ground/pin = PEC (35 um solids).
Board 83.2 x 70.4 mm (3 mm over the 80 mm guideline: the 5-row gate stack needs 41.6 mm from the probe). Gate on the **-x half only**; +x half is plain reference surface (front = +x, back = -x). `layout.png`, `layout.csv`, `sites.csv` (axis x/y per site).

**Feed / boundary.** Lumped port 50 ohm on `Port1`, integration line radial (pin -> r 2.05), same feed as `hex150_tensor_cp_single`. Radiation boundary on an air box padded 15 mm (sides, bottom) / 30 mm (top).

**Diode sites (two groups).** Assign a *Lumped RLC* boundary to every sheet of the active group; current axis = `sites.csv` column `current_axis` (integration line across the 0.4 mm gap, along that axis):
- **State OFF (pass):** `D_OFF_*` = C **28.16 fF** (no R, no L); delete or suppress `D_ON_*`.
- **State ON (block):** `D_ON_*` = R **2.44 ohm** (no L, no C; the measured ON residual X = -3.9 ohm is not represented, as in the verified Phase 28 package); delete or suppress `D_OFF_*`.
(MA4AGP907 measured values at 10.54 GHz, `src/ris/diodes_measured.py`, 10 mA / 0 V.) Run both states with identical mesh.

**Solve** 10.54 GHz, adaptive passes to delta S 0.01; interpolating sweep 10.0-11.1 GHz, 0.05 GHz steps (the model says any band effect, if it exists, sits at 10.2-10.95 GHz).

**Export (both states).**
1. S11 vs frequency (`.s1p` or CSV).
2. Far field at 10.54 GHz: realised gain (total) in the xz plane, phi = 0 (front, +x) and phi = 180 (back, -x), theta 0..180 deg, 1 deg step (CSV). Front/back ratio = gain(theta, 0)/gain(theta, 180) at theta = 60, 80, 90 deg and its maximum over theta (the pin's direct radiation is symmetric and dilutes it: also report the ratio at theta = 90).
3. E-field on the surface: Mag_E and complex Ez on the line y = 0, z = +1 mm, x = -38 ... +38 mm, 1 mm step (CSV); and Mag_E on the plane z = +1 mm over the whole board, 1 mm grid.

**Read-out.** Gate transmission = |E(x = -32 mm)| / |E(x = +32 mm)| (OFF should be 0 +- 1 dB; ON is the isolation, model: -3 to -5 dB leaky block; -15 dB or better would mean the original cell claim holds). Compare with `predictions.json` falsifiers.
Caveats: the reference/transition cells are plain square patches (equivalent sheets, not the slotted circles of the hologram); the board edges radiate (no edge treatment).
