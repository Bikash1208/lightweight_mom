# HFSS: crossed two-diode rim cell, 4 states, eigenmode 90 / 120 deg along x and y (one page). NOT SENT. MODEL_ONLY predictions in `predictions.json`.

**Files** (identical geometry, mm, +z up, patch plane z = 0, all coordinates on a 0.001 mm grid; `EXACT_DIMS` = PASS for all four): `cell_xOFF_yOFF.step`, `cell_xON_yOFF.step`, `cell_xOFF_yON.step`, `cell_xON_yON.step`.
Bodies (4 solids + 2 sheets each): `Substrate` (RT5880, eps_r 2.2, tan d 0.0009, 6.4 x 6.4, z -1.575..0), `Ground` (35 um, z -1.610..-1.575), `Patch` (circle r = 3.000, 35 um, minus two 0.4 mm slits
along y (x -0.2..0.2) and along x (y -0.2..0.2), each 5.6 mm long: **rim 0.200 mm** at every slit end, so the patch is one connected body), `Air` (6.4 x 6.4, z 0..30, patch cut out).
Sheets: `Dx_<state>` 0.4 x 1.2 (x -0.2..0.2, y 1.0..2.2, z +0.0175), integration line **along x**; `Dy_<state>` 1.2 x 0.4 (x 1.0..2.2, y -0.2..0.2), integration line **along y**. `<state>` = OFF or ON, set by the file name.

**Assign.** PEC: `Patch`, `Ground`, **top face of `Air` (PEC lid, z = +30 mm)**. Lumped RLC: ON sheet R = **2.44 ohm** only; OFF sheet C = **28.16 fF** only (MA4AGP907 measured at 10.54 GHz; ON residual X = -3.9 ohm not represented).
**Periodic, x propagation:** master/slave on the x faces, phase **90 and 120 deg**, y faces 0 deg. **y propagation:** the same with the faces swapped. (8 solves per state in all; the x/y symmetry means cell_xON_yOFF along y = cell_xOFF_yON along x, a free consistency check.)

**Solve.** Eigenmode, min frequency 5 GHz, 8 modes, delta f 0.1 %, mesh 1 mm in the substrate, 0.1 mm at the slits/rims and sheets. ON is lossy: compare Re(f).
**Mode.** TM0 surface band just below the light line f = c*phase/(2 pi a) = 11.71 GHz (90 deg), 15.61 GHz (120 deg); modes at or above it are PEC-lid modes, Ey-dominated modes are TE-like. beta/k0 = phase/(a k0(f)).

**Prediction** (uniform-sheet TM0 from the MoM Y(f) of each state, mesh 64, **+-0.15 GHz**: the mesh-96 dY is 6 % smaller than mesh 64, so splits are, if anything, over-predicted). Frequencies [GHz] of the TM0 band, x propagation:
| state (x diode / y diode) | 90 deg | 120 deg |
|---|---|---|
| OFF / OFF | 10.914 | > 14.5 (outside the Y(f) table) |
| ON / OFF | 11.033 | 13.397 |
| OFF / ON | 10.906 | > 14.5 |
| ON / ON | 11.041 | 13.383 |
y propagation: the same table with the diode roles swapped (the y diode is the on-axis one; OFF/ON along y = 11.033 GHz, ON/OFF along y = 10.906 GHz).
Reading: **the diode on the propagation axis raises the 90-deg band by +0.12 GHz** (ON removes slit capacitance, Y lower); the crossed diode moves it by < 0.01 GHz (independence); ON/ON ~ ON/OFF (additive). Flat band near the operating point (beta/k0 at 10.54 GHz ~1.07).
**Falsifiers:** on-axis ON-OFF split at 90 deg below +0.04 or above +0.25 GHz (wrong sign = wrong diode orientation); crossed-diode shift above 0.05 GHz; ON/ON more than 0.05 GHz from ON/OFF. The rim is the sensitive feature: a rim of >= 0.25 mm would cut the split by ~35-45 % (CODED_RIM.md).
