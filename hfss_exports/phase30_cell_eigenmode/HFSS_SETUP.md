# HFSS_SETUP - Phase 30 cell eigenmode check (manual, from STEP; MODEL_ONLY predictions below were written before any HFSS run)

**Files** (mm, origin at the cell centre on the ground plane; z up): `E1_strong_M0p64.step`, `E2_longthin_M0p62.step`, `E3_mid_M0p27.step`, `E4_weak_M0p12.step`, `S1_supercell_leaky.step`.
Bodies in each: **SUB** (solid, 3.2 x 3.2 x 0.635, z 0..0.635; S1: 22.4 x 3.2), **GND** (sheet, z = 0), **PATCH** (sheet at z = 0.635, 72-vertex n = 4 super-ellipse 3.05 mm, gap 0.15 mm to the cell edge, rectangular slot cut out; S1: 7 sheets), **AIR** (solid, same footprint, z 0.635..60.635: **height 60.0 mm = 94.5 x the substrate**; lid is 6x = 3.81 mm minimum, but 60 mm keeps the lid field ~2e-5 so the predictions hold - do not shorten).
Slots (L x W, along x): E1 1.89 x 0.58, E2 2.29 x 0.10, E3 1.89 x 0.10, E4 1.53 x 0.10. S1 = E1 x 7 cells (x = -9.6 ... +9.6), slot angle 0, 25.7, 51.4 ... 154.3 deg; outline fixed, only the slot turns.

**Import**: Modeler > Import (STEP), units mm, heal OFF. Check the four body names. New project per file, Solution type **Eigenmode**.

**Materials**: SUB = new `RO3010`, eps_r **10.2** (run 1) and **11.2** (run 2, E1 only; edit the material), tan delta 0.0022, no conductivity. AIR = vacuum. PATCH, GND = **Perfect E** (lossless; copper not needed - conductor loss is ~3 % of the radiation in the model).

**Boundaries** (unit cells E1-E4; cell side = 3.2 mm, faces at x = +-1.6, y = +-1.6; each plane has TWO faces: SUB side + AIR side - select both):
1. Perfect E on PATCH and GND.
2. Perfect E on the top face of AIR (z = 60.635) = PEC lid.
3. Primary/Secondary pair **X**: Primary = faces at x = -1.6, Secondary = faces at x = +1.6; U vector on both along +y (start y = -1.6, end y = +1.6); Phase delay = **ph_x** (design variable, deg).
4. Pair **Y**: Primary y = -1.6, Secondary y = +1.6; U along +x; Phase delay = **ph_y**.
   Direction x (along the slot): ph_x = phase, ph_y = 0. Direction y: ph_x = 0, ph_y = phase. (Sign of the phase does not matter.)
   **S1** (supercell 22.4 x 3.2): same, faces at x = +-11.2 and y = +-1.6; the top of AIR is a **Radiation** boundary (no PEC lid); ph_x = 6.74 deg, ph_y = 0 for eps_r 10.2 (computed from beta at 10.54 GHz). **For eps_r 11.2 (design value, operating point ~10.1 GHz) use ph_x = 0.4 deg, ph_y = 0**: broadside needs beta = 2 pi/22.4 mm = 280.5 rad/m, which the model gives at 10.10 GHz (beta/k0 1.326), i.e. 360.4 deg per supercell. Expected: f_r ~10.07 GHz, Q 45 (first-order) to 91 (Floquet); Q = f_r/(2 Im f); alpha/k0 = c/(2 Q v_g), v_g from the 5-point sweep (model 0.160 c).

**Eigenmode setup**: Minimum frequency **8 GHz**; Number of modes **6** (S1: **10**); Max passes **15**; Min passes 2; Min converged passes 2; **Max delta frequency per pass 0.1 %**; Converge on real frequency only (S1: leave off - complex); basis order first. Mesh: default adaptive, no extra seeding. Bound mode = the lowest mode whose E_z peaks at the patch edges, decays up into the air and is nodal at the ground (field plot of |E|; other modes are expected to be lid modes above ~13 GHz, not verified). If two bound modes appear (x/y polarisations), report both as f1, f2.

**Runs** (phase = beta(10.54 GHz) x 3.2 mm from the model; predicted bound-mode f, GHz, mesh 48 / mesh 64):
| # | cell | eps_r | direction | ph (deg) | model f |
|---|---|---|---|---|---|
| 1,2 | E1 | 10.2 | x, y | 52.7 | x 10.90/11.14, y 10.26/10.43 |
| 3,4 | E2 | 10.2 | x, y | 51.9 | x 10.91/11.16, y 10.25/10.53 |
| 5,6 | E3 | 10.2 | x, y | 49.0 | x 10.72/10.94, y 10.38/10.64 |
| 7,8 | E4 | 10.2 | x, y | 47.7 | x 10.62/10.84, y 10.46/10.69 |
| 9,10 | E1 | 11.2 | x, y | 63.6 | NA (about 0.8-1 GHz below runs 1,2) |
| 11,12 | E1 | 10.2 | x | 40 and 65 | NA (gives v_g) |
| 13 | S1 | 10.2 | x | 6.74 | Q 49 (first order) ... 91 (Floquet); Im f 0.109 ... 0.058 GHz |
Pass: f_x, f_y within 0.3 GHz of the mesh-64 values and the x/y split within 30 %; S1 Q between ~40 and ~100 (Q > 200 = alpha far below the model).

**Copy back to the reviewing agent** (one row per run; f in GHz; for S1 give Re f and Im f, Q = Re f / (2 Im f); passes = adaptive passes used; also say if any run did not converge):
| cell | eps_r | direction | phase (deg) | f1 | f2 | Q | passes |
|---|---|---|---|---|---|---|---|
| E1 | 10.2 | x | 52.7 | | | | |
(for the other runs give the same columns; Q only for S1, otherwise NA)
