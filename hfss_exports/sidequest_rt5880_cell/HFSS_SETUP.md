# HFSS: RT5880 coded cell (rim disc, diodes at 0.8 mm), 4 states, eigenmode 90 deg along x then y (one page). Predictions MODEL_ONLY: `predictions.json`.

**Files** `cell_xOFF_yOFF.step`, `cell_xON_yOFF.step`, `cell_xOFF_yON.step`, `cell_xON_yON.step` (x diode first). mm, +z up, patch plane z = 0, every coordinate on a 0.001 mm grid (`EXACT_DIMS` PASS: 4 solids + 2 sheets each, bounding boxes/volumes < 1e-5, 0 inexact literals).
Bodies: `Substrate` (RT5880 eps_r 2.2, tan d 0.0009, 6.4 x 6.4, z -1.575..0), `Ground` (35 um, below), `Patch` (circle r 3.000, 35 um, minus two 0.4 mm slits 5.6 mm long along y and x; rim 0.2 mm; one connected body),
`Air` (6.4 x 6.4, z 0..30, patch cut out), sheets `Dx_<state>` 0.4 x 1.2 at (0, +0.8) [integration line along **x**] and `Dy_<state>` 1.2 x 0.4 at (+0.8, 0) [line along **y**], z +0.0175.

**Assign (nothing is implicit).** PEC: `Patch`, `Ground`, and the **top face of `Air` (the PEC lid, z = +30 mm): assign it explicitly**. Unassigned outer faces of the vacuum region default to PEC, so the **side faces must be assigned**
(Lattice Pair) or the cell becomes a PEC box. Lumped RLC on each sheet: ON **R = 2.44 ohm**, OFF **C = 28.16 fF** (MA4AGP907 measured, 10.54 GHz). The state is set by the file (ON/OFF in the sheet name).
**Phase variables.** Project > Project Variables: `$ph_x`, `$ph_y` (deg). Boundaries > Lattice Pair 1: Primary/Secondary on the x = -3.2 / +3.2 mm side faces of `Substrate` and `Air` (Phase Delay = `$ph_x`); Lattice Pair 2: y = -3.2 / +3.2 faces, `$ph_y`.
Run **x**: `$ph_x = 90deg`, `$ph_y = 0deg`; run **y**: `$ph_x = 0deg`, `$ph_y = 90deg` (a 2-point parametric sweep). 4 files x 2 runs = **8 runs**.
**Eigenmode.** Minimum frequency 5 GHz, 8 modes, delta f 0.1 %, mesh 1 mm in the substrate, 0.1 mm on the slits and sheets. The bound mode is the TM0 band just **below the light line 11.711 GHz** (Ez under the patch, field decaying upward; modes above it are PEC-lid modes).
ON is lossy: take Re f (and Q = Re f / (2 Im f) if you want the loss).

**Prediction per run** (Bloch tensor of the MoM at kt = beta = 90 deg / 6.4 mm = 245.4 rad/m, transverse resonance; about +-0.15 GHz, mesh spread 3 %). Y in mS (susceptance); Y_axis = the entry along the propagation axis:
| state (x/y diode) | axis | f [GHz] | beta/k0 | Y_axis | Y_other | Y_xy |
|---|---|---|---|---|---|---|
| OFF_OFF | x | 10.527 | 1.113 | 11.33 | 6.64 | +0.09 |
| OFF_OFF | y | 10.527 | 1.113 | 11.33 | 6.64 | +0.09 |
| ON_OFF | x | 11.001 | 1.065 | 7.09 | 8.40 | -0.35 |
| ON_OFF | y | 10.528 | 1.112 | 11.32 | 4.41 | -0.23 |
| OFF_ON | x | 10.528 | 1.112 | 11.32 | 4.41 | -0.23 |
| OFF_ON | y | 11.001 | 1.065 | 7.09 | 8.40 | -0.33 |
| ON_ON | x | 11.027 | 1.062 | 6.83 | 4.84 | +0.66 |
| ON_ON | y | 11.027 | 1.062 | 6.83 | 4.84 | +0.66 |
Reading: the on-axis diode ON raises the band by **+0.47 GHz** (Y 11.33 -> 7.09 mS); the crossed diode moves it < 0.002 GHz; ON/ON is within 0.03 GHz of one diode ON (additive; Y_xy < 0.7 mS). **Falsifiers:** on-axis split below +0.25 or above +0.8 GHz; crossed-diode shift above 0.05 GHz; ON/ON more than 0.1 GHz from ON/OFF.

**Back out alpha_beam.** Fill `results.csv` (`state,axis,f_GHz[,Q]`, state = OFF_OFF / ON_OFF / OFF_ON / ON_ON, axis x|y) and run `python scripts/research/hfss_alpha_backout.py results.csv`. Formulas (beta = 245.4 rad/m for every run):
Y_s = j B_s, B_s = - Im[Y_air(beta, f) + Y_slab(beta, f)] (transverse resonance); beta/k0 = beta c / (2 pi f); dY = B_s(ON on that axis) - B_s(OFF/OFF); Y_m = mean of all B_s; beta_m from Im[jY_m + Y_air + Y_slab] = 0 at 10.54 GHz;
**alpha_beam = leakage_alpha_first_order(beta_m, jY_m, c1 = |dY|/pi, theta 25 deg, zeta 0)** (alpha ~ c1^2); alpha_abs ~ (d beta_m/dB) G_ON/2 with G_ON from Q. Round trip on the predictions: dY = -4.24 mS, beta_m 1.081 k0, alpha_beam **1.43 Np/m** (model range 1.5-1.9; mesh-96 1.74).
