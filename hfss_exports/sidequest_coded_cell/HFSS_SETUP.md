# HFSS: 1-bit coded-hologram cell, ON / OFF eigenmode (one page). NOT SENT. MODEL_ONLY predictions in `predictions.json`; nothing is HFSS-validated.

**Files.** `cell_ON.step` (sheet `D_ON`), `cell_OFF.step` (sheet `D_OFF`): identical geometry, mm, +z up, patch plane z = 0, every coordinate on a 0.001 mm grid (`EXACT_DIMS` = PASS, see below).
Bodies (5 solids + 1 sheet each): `Substrate` (RT5880, eps_r 2.2, tan d 0.0009, 6.4 x 6.4, z -1.575..0), `Ground` (35 um, z -1.610..-1.575), `Patch_L`, `Patch_R` (2.8 x 4.8, 35 um, x = -3.0..-0.2 and 0.2..3.0, y = -2.4..2.4),
`Air` (6.4 x 6.4, z 0..30, metal cut out), sheet `D_ON` / `D_OFF` (0.4 x 1.2, x -0.2..0.2, y -0.6..0.6, z +0.0175) bridging the 0.4 mm gap. The split patch is 6.0 mm long: 0.4 mm gap to the next cell.

**Assign.** PEC: `Patch_L`, `Patch_R`, `Ground`, and the **top face of `Air` (PEC lid, z = +30 mm)**. Substrate / Air as materials. Lumped RLC on the sheet, integration line along **+x** across the gap:
**ON** (`cell_ON`): R = **2.44 ohm** only. **OFF** (`cell_OFF`): C = **28.16 fF** only (MA4AGP907 measured, 10 mA / 0 V at 10.54 GHz; the ON residual X = -3.9 ohm is not represented, as in the verified Phase 28 package).
Periodic: master/slave on the x faces with phase = **90 and 120 deg**, on the y faces 0 deg (two solves per state, four in all).

**Solve.** Eigenmode, minimum frequency 5 GHz, 8 modes, delta f 0.1 %, mesh seeded 1 mm in the substrate / 0.2 mm at the gap (about 30-60 k tets). ON is lossy (R): HFSS reports complex f, compare Re(f).

**Identify the mode.** The surface-wave (TM0-like) band lies just *below* the light line f = c*phase/(2 pi a) = 11.71 GHz (90 deg) and 15.61 GHz (120 deg), Ez concentrated under the patch and across the gap; modes at or above the line are PEC-lid (parallel-plate) modes, Ey-dominated modes are TE-like.
In HFSS beta = phase/a, so beta/k0 = phase/(a k0(f)).

**Prediction** (uniform-sheet TM0 model from the MoM Y_xx(f) of this cell, mesh 32; +-5 % in frequency):

| phase [deg] | light line [GHz] | OFF: f [GHz] (beta/k0) | ON: f [GHz] (beta/k0) |
|---|---|---|---|
| 90 | 11.71 | 11.32 (1.035) | 10.62 (1.103) |
| 120 | 15.61 | > 14.5 (outside the Y(f) table) | 12.68 (1.232) |

Reading: the diode state moves the 90-deg band by **0.70 GHz** (ON lower). At the design frequency 10.54 GHz the ON/OFF contrast in beta is the difference of the two columns at the same f; the hologram mean sheet (beta/k0 = 1.054) sits between them.
**Falsifiers:** the 90-deg ON-OFF split below 0.35 GHz or above 1.4 GHz; ON f(90 deg) outside 10.1-11.1 GHz; any ON band crossing from below to above 10.54 GHz between 90 and 120 deg without the predicted ordering (ON < OFF at both phases).
**Why this cell:** MoM sheet pair OFF 0.002 + 2.985j mS, ON 1.199 + 10.593j mS (mesh 64), Y_yy 3.885j mS static; alpha_beam 4.0 Np/m.

**EXACT_DIMS** (re-opened both STEP files): 5 solids + 1 sheet, bounding boxes and volumes equal the intended values to < 1e-6, 0 inexact literals, 8 snapped tokens (max 2e-16 mm), no merged/collapsed vertices.
