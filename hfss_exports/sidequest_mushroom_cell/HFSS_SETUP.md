# HFSS: switchable via-to-ground (mushroom-type) gate cell — eigenmode (one page). MODEL_ONLY predictions in `predictions.json`; nothing is HFSS-validated.

**Files.** `cell_ON.step` (sheet `D_ON`), `cell_OFF.step` (sheet `D_OFF`): identical geometry, mm, +z up, patch plane z = 0. Parts: `Substrate` (RT5880, eps_r 2.2, tan d 0.0009, z -1.575..0, via hole),
`Ground` (35 um, z -1.610..-1.575), `Via` (r 0.2, z -1.575..0), `Pad` (0.8 x 0.8, 35 um), `Patch` (5.6 x 5.6 with a 1.6 x 1.6 hole, 35 um; 0.8 mm gap to the next cell), `Air` (6.4 x 6.4, z 0..30, metal cut out), one sheet.
The MA4AGP907 sits **at the via top**: it bridges the 0.4 mm gap between `Pad` (on the via) and `Patch`, current along +x (sheet x 0.4..0.8, y -0.4..0.4 at z +0.0175).

**Assign.** PEC: `Patch`, `Pad`, `Via`, `Ground` and the **top face of `Air` (the PEC lid at z = +30 mm)**. Substrate and Air as materials. Lumped RLC on the sheet (integration line along +x across the 0.4 mm gap):
**ON** (`cell_ON`): R = **2.44 ohm** only. **OFF** (`cell_OFF`): C = **28.16 fF** only (MA4AGP907 measured, 10 mA / 0 V at 10.54 GHz; the ON residual X = -3.9 ohm is not represented, as in the verified Phase 28 package).
Periodic: master/slave on the x faces with phase = **0, 30, 60, 90, 120, 150, 180 deg**, on the y faces 0 deg (linked-phase sweep or 7 solves).

**Solve.** Eigenmode, minimum frequency 2 GHz, 8 modes, delta f 0.1 %, mesh seeded 1 mm in the substrate / 0.2 mm at the gap and via (about 30-60 k tets, minutes per phase). ON is lossy (R): HFSS gives complex f, compare Re(f) (and Q).

**Identify the modes.** The surface-wave (TM-like) mode lies just *below* the light line f = c*phase/(2 pi a) = 3.90 / 7.81 / 11.71 / 15.61 GHz at 30 / 60 / 90 / 120 deg, with Ez concentrated at the via/pad; modes at or above it are PEC-lid (parallel-plate) modes, Ey-dominated modes are TE-like. Phase 0: the band starts at 0 Hz.

**Prediction (model, +-5 %)** for the TM surface band:

| phase [deg] | light line [GHz] | ON: TM surface band [GHz] | OFF: TM surface band [GHz] |
|---|---|---|---|
| 0 | 0.00 | - | - |
| 30 | 3.90 | 3.79 | 3.84 |
| 60 | 7.81 | 7.22 | 7.60 |
| 90 | 11.71 | 9.01 | 10.76 |
| 120 | 15.61 | 9.65 | 13.07 |
| 150 | 19.52 | 9.96 | 15.13 |
| 180 | 23.42 | 10.15 | - |

**ON stopband** (no TM band between): lower edge = f(180 deg) = **10.15 GHz**, upper edge = via-patch resonance = **10.97 GHz** (target 10.54 +- 0.3 = 10.24-10.84: covered). Read the upper edge from the ON band flattening toward 10.98 at 150-180 deg.
**OFF:** beta/k0 at 10.54 GHz = **1.079** (hologram average 1.065); in HFSS beta = phase/a at the OFF band crossing f = 10.54 GHz (model: about 88 deg, f(90 deg) = 10.76 GHz).
**Falsifiers:** ON f(180 deg) below 9.7 or above 10.6 GHz, or a TM band crossing 10.54 GHz in ON for any phase; OFF beta/k0 outside 1.04-1.12. Tuning knob: patch side (5.2 / 5.6 / 6.0 mm moves the ON resonance 12.4 / 11.0 / 9.5 GHz).
