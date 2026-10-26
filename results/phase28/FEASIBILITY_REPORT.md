# Phase 28 item 2 - feasibility of a coded monopole-fed hologram (stop point)

**Labels:** `MODEL_ONLY`, `IDEALISED_CELLS`, `NOT_HFSS_VALIDATED`, `SUPERCELL_TABLE_6P4MM_10P54GHZ_ASSUMED_AT_ALL_PITCHES`. **Generated from** `results/phase28/feasibility_results.json` only.

Forward model: the Phase 26 tensor-CP model unchanged (supercell table, ray march, aperture far field on the grounded slab). Frequency 10.54 GHz, the only frequency of the supercell table; 10.0 GHz is `NOT_AVAILABLE_FROM_CURRENT_DATA (the supercell table exists at 10.54 GHz only)`. Cells are idealised (a state is a tensor or scalar value, no geometry); the table, computed for 6.4 mm cells, is assumed at every pitch, so the pitch sweep changes only how the hologram is sampled. Gain is referred to the accepted power with the Phase 26 bare launch fraction 0.592 (the calibrated launch fraction 0.62 of the four-hexagon work would add 0.20 dB to every gain here, and nothing to any loss).

## Anchor: no quantisation reproduces the Phase 26 single-hexagon prediction

Source: results/phase26/tensor_cp_hexagon.json (selected design: hand +1, eta_design 0.9, taper 0.5). Same grid, same code path; the finer refinement of the peak is reported separately.

| quantity | this pipeline | Phase 26 | abs diff | tolerance | within |
|---|---|---|---|---|---|
| d_ap_rhcp_dbi | 28.477911 | 28.477911 | 7.1e-15 | 1e-06 | yes |
| gain_launch_bare_dbi | 23.511963 | 23.511963 | 7.1e-15 | 1e-06 | yes |
| ar_at_peak_db | 2.504165 | 2.504165 | 0.0e+00 | 1e-06 | yes |
| xpd_at_peak_db | 16.883500 | 16.883500 | 0.0e+00 | 1e-06 | yes |
| rhcp_sll_db | -22.773995 | -22.773995 | 7.1e-15 | 1e-06 | yes |
| eta_conv | 0.538373 | 0.538373 | 0.0e+00 | 1e-12 | yes |

Refined peak directivity 28.565 dBi (0.087 dB above the coarse-grid peak), pointing error 0.44 degrees.

## What the model got wrong in the one HFSS run

Source: results/phase26/hex150_tensor_cp_single_hfss_comparison.json (one HFSS run, one board, 10.5 GHz). The prediction is the continuous Phase 26 design; the HFSS run is one board at one frequency.

| quantity | HFSS | model |
|---|---|---|
| peak RHCP gain [dBi] | 23.03 | 23.51 |
| RHCP sidelobe level [dB] | -18.5 | -22.8 |
| XPD at peak [dB] | 23.0 | 16.9 |
| AR at peak [dB] | 1.23 | 2.50 |
| LHCP share of the radiated power | 0.310 | 0.008 |

The model sidelobe levels below are therefore optimistic by about the difference above, and the HFSS LHCP share includes the pin's direct radiation, which the model does not radiate. The bit-depth comparison uses the model's relative numbers; absolute sidelobe levels and cross-polar levels are not predictions.

## Switching contrast, first

**Definition.** The contrast is what the cell must be able to change so that the coded hologram can leak as strongly as the continuous one. Four measures, all at the design frequency:

- `alpha_cap`: the largest local leakage constant [Np/m] the coded surface can realise. The continuous Phase 26 design asks for a profile alpha(x) up to the supercell table's maximum; a coded surface that cannot reach it converts less of the launched surface wave.
- Tensor (CP, rotating axis) cell: the fixed anisotropy `delta` [mS] of `Y = Ybar I + j delta R(a) diag(+,-) R(a)^T`, with alpha equal to a constant times delta squared, from the Phase 26 closed form (`leakage_alpha_first_order_tensor`). `dX_aniso` [ohm] is the difference of the TM surface reactances `X = Im(1/(Y + Y_slab,TM(beta)))` seen along the two principal axes; the switch rotates the axis and does not change it. The largest state-to-state change of the sheet tensor (axes perpendicular) is twice delta.
- Scalar cell (also the amplitude control of the `lpx` bound), states `Ybar +- j dB/2` (`dB` the full swing between the extreme levels): the first harmonic of the quantised surface is c1, one quarter of dB times g (g the first-harmonic gain of the level set), and alpha is a constant times c1 squared (`leakage_alpha_first_order`, TM and TE parts, the larger `dB` reported). `dX_switch` [ohm] is the TM surface-reactance difference of the two extreme states; `d_beta` the difference of the TM0 phase constants of the two uniform states.
- Modulation depth: `delta/|Y_down|` (tensor) or half of dB over `|Y_down|` (scalar), `Y_down = Ybar + Y_slab,TM(beta)`, the Phase 26 `required_map` definition. That map asks for far larger values because it targets a much higher conversion of the launched wave; the continuous Phase 26 design itself converts only the eta_conv of the scan table below.

**Method.** The supercell table's leakage constants are capped at `alpha_cap` and the single-hexagon CP design is re-marched; `alpha_cap` is the value at which the gain falls 1.0 dB below the uncapped design (table maximum 4.259 Np/m, minimum 0.422 Np/m). It is converted to physical units with the Phase 26 closed forms at the mean `beta/k0` of the design's cells; the scalar case is cross-checked with the exact truncated-Floquet complex beta. This row is the contrast-only requirement: it ignores the quantisation loss, which is added in the combined table below.

| beam theta / phi [deg] | march points outside the table | alpha_cap [Np/m] | loss re-check [dB] | beta/k0 mean | tensor delta [mS] | tensor dX_aniso [ohm] | tensor d_beta/k0 | tensor depth | LP binary dB [mS] | LP dX_switch [ohm] | LP d_beta/k0 | LP four-level dB [mS] | LP four-level dX_switch [ohm] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 / 0 | 0.00 inside the table (closed form indicative only) | 1.787 | 0.999 | 1.0593 | 2.054 | 76.8 | 0.0425 | 0.270 | 4.563 | 86.9 | 0.0472 | 5.498 | 109.6 |
| 30 / 225 | 0.22 TABLE_EXTRAPOLATED | 2.308 | 0.996 | 1.0621 | 2.296 | 92.3 | 0.0494 | 0.310 | 5.448 | 114.5 | 0.0585 | 6.564 | 148.4 |
| 50 / 225 | 0.57 TABLE_EXTRAPOLATED | 2.412 | 0.995 | 1.0623 | 2.337 | 94.7 | 0.0504 | 0.316 | 7.224 | 172.9 | 0.0772 | 8.704 | 242.5 |

Exact-Floquet cross-check of the scalar closed form (binary sheet at the closed-form full swing found for each beam; TM, collinear):

| beam theta / phi [deg] | alpha closed form [Np/m] | alpha Floquet exact [Np/m] | ratio exact/closed | note |
|---|---|---|---|---|
| 0 / 0 | 1.787 | n/a | n/a | NOT_AVAILABLE: the Floquet balance has no finite root at broadside (open stop band) |
| 30 / 225 | 2.308 | 2.343 | 1.015 | - |
| 50 / 225 | 2.412 | 2.867 | 1.189 | - |

A ratio away from one means the closed form under- or overstates the leakage for that swing (finite modulation depth); the contrast numbers of that row carry that uncertainty.

**Combined requirement: 2-bit coded surface.** The quantisation loss adds to the contrast loss. Each capped design is also quantised to four slot axes; the table gives the smallest `alpha_cap` whose total loss against the uncapped continuous design stays within the tolerance, converted to the tensor and scalar contrast (four-level swing). The 2-bit loss of the uncapped design is the first column. The cap is interpolated linearly in loss between the grid points of the cap sweep, so it can be slightly conservative (a larger cap than strictly needed).

| beam theta / phi [deg] | 2-bit loss, uncapped [dB] | total loss tolerance [dB] | alpha_cap [Np/m] | tensor delta [mS] | tensor dX_aniso [ohm] | scalar four-level dX_switch [ohm] | caveat |
|---|---|---|---|---|---|---|---|
| 0 / 0 | 0.93 | 1.0 | 2.581 | 2.469 | 95.6 | 141.2 | inside the table; closed form indicative only |
| 0 / 0 | 0.93 | 1.5 | 2.063 | 2.207 | 83.5 | 120.6 | inside the table; closed form indicative only |
| 0 / 0 | 0.93 | 2.0 | 1.745 | 2.030 | 75.8 | 107.9 | inside the table; closed form indicative only |
| 30 / 225 | 0.95 | 1.0 | 4.010 | 3.027 | 132.0 | 238.4 | TABLE_EXTRAPOLATED |
| 30 / 225 | 0.95 | 1.5 | 2.737 | 2.500 | 102.6 | 169.3 | TABLE_EXTRAPOLATED |
| 30 / 225 | 0.95 | 2.0 | 2.263 | 2.273 | 91.2 | 146.2 | TABLE_EXTRAPOLATED |
| 50 / 225 | 0.91 | 1.0 | 3.966 | 2.997 | 130.7 | 471.1 | TABLE_EXTRAPOLATED; Floquet exact/closed ratio 1.189 |
| 50 / 225 | 0.91 | 1.5 | 2.853 | 2.542 | 105.1 | 291.9 | TABLE_EXTRAPOLATED; Floquet exact/closed ratio 1.189 |
| 50 / 225 | 0.91 | 2.0 | 2.336 | 2.300 | 92.8 | 234.7 | TABLE_EXTRAPOLATED; Floquet exact/closed ratio 1.189 |

For scale only (the cell-level contrast depends on how the diode is embedded, which is item 3): the terminal reactance swing of the diode models at the design frequency.

| diode (datasheet, PROVISIONAL) | R_on [ohm] | X_on [ohm] | X_off [ohm] | X_off - X_on [ohm] |
|---|---|---|---|---|
| MADP-000907-14020 | 5.20 | 0.0 | -604.0 | -604.0 |
| SMP1345-079LF | 1.50 | 46.4 | -33.1 | -79.5 |
| SMP1345-079LF_1mA | 3.50 | 46.4 | -33.1 | -79.5 |

MADP-000907-14020 package inductance is NOT_AVAILABLE from the datasheet (L_s = 0 here); datasheet values are PROVISIONAL.

**Table range caveat.** The supercell table covers a limited range of modulation wavenumbers; rows with `TABLE_EXTRAPOLATED` above (and in the scan table below) have a large share of march points outside it, and their contrast and scan numbers are extrapolations.

## Design point: 150 mm hexagon, beam theta 30 / phi 225 degrees

This beam has 0.22 of its march points outside the supercell table: `TABLE_EXTRAPOLATED`; the design-point numbers share that caveat.

Quantisation loss is relative to the continuous design of the same family and pitch. Two definitions are given because they differ for a phase-only quantisation:

- *fixed source*: the same cell currents, so the peak field scales as the pattern's directivity times its radiated power. A phase-only quantisation (`cp`) leaves every current magnitude unchanged, and the quantisation error is radiated into the non-radiating part of the spectrum; that power does not contribute to the beam (where it goes is not modelled), and this is the loss used as the headline for `cp`.
- *renormalised*: the directivity of the quantised pattern alone, i.e. its own radiated power renormalised to the launched power. It hides the lost power. It is the energy-consistent one for an amplitude quantisation (`lp`, `lpx`), where the currents change magnitude and would change the leakage, and is the headline for those.

`cp` = rotating-slot tensor, RHCP (the Phase 26 family); `lp` = scalar surface, current along the radial direction; `lpx` = linear polarisation with IDEAL polarisation control (current along x in every cell), an upper bound no cell has been shown to realise. AR and XPD beyond the bound are shown as bounds, never as numbers. A negative renormalised loss is a quantised pattern whose directivity is higher than the continuous one on this sampling; a negative fixed-source loss for `lp` and `lpx` is not a gain: the plus and minus one levels have larger magnitude than the continuous modulation, so the radiated power rises (radiated power change column).

| pitch [mm] | variant | gain, headline [dBi] | loss, headline [dB] | loss, fixed source [dB] | loss, renormalised [dB] | radiated power change [dB] | SLL [dB] | pointing error [deg] | XPD [dB] | AR [dB] |
|---|---|---|---|---|---|---|---|---|---|---|
| 4.0 | cp_cont | 23.57 | 0.00 | 0.00 | 0.00 | 0.00 | -22.6 | 0.44 | 17.00 | 2.47 |
| 4.0 | cp_1bit | 19.67 | 3.89 | 3.89 | 3.00 | -0.90 | -19.0 | 0.44 | 15.71 | 2.87 |
| 4.0 | cp_2bit | 22.66 | 0.91 | 0.91 | 0.06 | -0.85 | -22.5 | 0.44 | 17.12 | 2.44 |
| 4.0 | cp_3bit | 23.34 | 0.22 | 0.22 | 0.02 | -0.20 | -22.5 | 0.44 | 16.89 | 2.50 |
| 4.0 | lp_cont | 18.53 | 0.00 | 0.00 | 0.00 | 0.00 | -2.5 | 5.06 | - | - |
| 4.0 | lp_1bit | 18.48 | 0.04 | -2.11 | 0.04 | 2.16 | -2.5 | 5.06 | - | - |
| 4.0 | lp_2bit | 18.55 | -0.03 | -0.48 | -0.03 | 0.46 | -2.5 | 5.06 | - | - |
| 4.0 | lpx_cont | 23.68 | 0.00 | 0.00 | 0.00 | 0.00 | -21.6 | 0.44 | - | - |
| 4.0 | lpx_1bit | 23.61 | 0.07 | -2.10 | 0.07 | 2.17 | -21.0 | 0.44 | - | - |
| 4.0 | lpx_2bit | 23.67 | 0.01 | -0.45 | 0.01 | 0.45 | -21.6 | 0.44 | - | - |
| 5.0 | cp_cont | 23.52 | 0.00 | 0.00 | 0.00 | 0.00 | -22.5 | 0.44 | 16.93 | 2.49 |
| 5.0 | cp_1bit | 19.62 | 3.90 | 3.90 | 3.07 | -0.83 | -18.7 | 0.44 | 15.76 | 2.86 |
| 5.0 | cp_2bit | 22.60 | 0.92 | 0.92 | 0.16 | -0.76 | -23.0 | 0.44 | 17.01 | 2.47 |
| 5.0 | cp_3bit | 23.31 | 0.22 | 0.22 | 0.05 | -0.17 | -22.4 | 0.44 | 16.77 | 2.54 |
| 5.0 | lp_cont | 18.48 | 0.00 | 0.00 | 0.00 | 0.00 | -2.6 | 5.17 | - | - |
| 5.0 | lp_1bit | 18.22 | 0.26 | -1.99 | 0.26 | 2.25 | -2.1 | 5.16 | - | - |
| 5.0 | lp_2bit | 18.48 | 0.00 | -0.54 | 0.00 | 0.54 | -2.7 | 5.17 | - | - |
| 5.0 | lpx_cont | 23.63 | 0.00 | 0.00 | 0.00 | 0.00 | -21.7 | 0.44 | - | - |
| 5.0 | lpx_1bit | 23.47 | 0.16 | -2.09 | 0.16 | 2.25 | -20.0 | 0.44 | - | - |
| 5.0 | lpx_2bit | 23.61 | 0.02 | -0.50 | 0.02 | 0.53 | -21.9 | 0.44 | - | - |
| 6.4 | cp_cont | 23.60 | 0.00 | 0.00 | 0.00 | 0.00 | -22.9 | 0.44 | 17.17 | 2.42 |
| 6.4 | cp_1bit | 19.78 | 3.82 | 3.82 | 3.12 | -0.70 | -18.2 | 0.44 | 15.57 | 2.92 |
| 6.4 | cp_2bit | 22.65 | 0.95 | 0.95 | 0.20 | -0.75 | -22.8 | 0.61 | 17.22 | 2.41 |
| 6.4 | cp_3bit | 23.38 | 0.22 | 0.22 | 0.06 | -0.17 | -22.3 | 0.44 | 17.27 | 2.39 |
| 6.4 | lp_cont | 18.60 | 0.00 | 0.00 | 0.00 | 0.00 | -2.6 | 5.06 | - | - |
| 6.4 | lp_1bit | 18.45 | 0.16 | -2.10 | 0.16 | 2.25 | -2.6 | 5.06 | - | - |
| 6.4 | lp_2bit | 18.48 | 0.13 | -0.42 | 0.13 | 0.55 | -2.5 | 5.06 | - | - |
| 6.4 | lpx_cont | 23.70 | 0.00 | 0.00 | 0.00 | 0.00 | -22.1 | 0.44 | - | - |
| 6.4 | lpx_1bit | 23.51 | 0.19 | -2.04 | 0.19 | 2.23 | -20.7 | 0.44 | - | - |
| 6.4 | lpx_2bit | 23.63 | 0.07 | -0.47 | 0.07 | 0.54 | -22.0 | 0.44 | - | - |
| 8.0 | cp_cont | 23.58 | 0.00 | 0.00 | 0.00 | 0.00 | -22.4 | 0.44 | 16.92 | 2.49 |
| 8.0 | cp_1bit | 19.71 | 3.87 | 3.87 | 3.17 | -0.70 | -20.7 | 0.44 | 15.01 | 3.12 |
| 8.0 | cp_2bit | 22.68 | 0.90 | 0.90 | 0.29 | -0.61 | -21.4 | 0.44 | 16.95 | 2.49 |
| 8.0 | cp_3bit | 23.35 | 0.23 | 0.23 | 0.10 | -0.13 | -22.5 | 0.44 | 17.05 | 2.46 |
| 8.0 | lp_cont | 18.53 | 0.00 | 0.00 | 0.00 | 0.00 | -1.4 | 5.21 | - | - |
| 8.0 | lp_1bit | 18.21 | 0.32 | -2.14 | 0.32 | 2.47 | -2.6 | 5.17 | - | - |
| 8.0 | lp_2bit | 18.48 | 0.05 | -0.60 | 0.05 | 0.66 | -2.5 | 5.06 | - | - |
| 8.0 | lpx_cont | 23.68 | 0.00 | 0.00 | 0.00 | 0.00 | -21.9 | 0.44 | - | - |
| 8.0 | lpx_1bit | 23.37 | 0.31 | -2.13 | 0.31 | 2.44 | -21.2 | 0.44 | - | - |
| 8.0 | lpx_2bit | 23.62 | 0.06 | -0.58 | 0.06 | 0.64 | -21.6 | 0.44 | - | - |

Peak search: 0.25 degree direction-cosine refinement around the coarse-grid peak; the sidelobe level uses the strongest coarse-grid sidelobe refined the same way. The headline gain is the continuous design's gain (aperture directivity plus ten times the base-ten logarithm of launch fraction x eta_conv) minus the headline loss; eta_conv of this design is given in the scan table.

## Finding: a scalar surface on a monopole cannot form a linear-polarised pencil beam

With the induced current along the surface wave's own (radial) direction the continuous scalar design has a pointing error of 5.06 degrees and a strongest sidelobe at -2.6 dB: the aperture field is radially polarised, so the contributions cancel in the beam direction (a ring or split beam, not a pencil). Gain 18.60 dBi against 23.70 dBi for the same hologram with ideal polarisation control (`lpx`) and 23.60 dBi for the CP tensor design. The `lp` rows above are therefore NOT a candidate; the bit-depth conclusion for LP rests on `lpx`, which needs a cell that sets the current direction (not designed).

## Steering: pitch 6.4 mm, 150 mm hexagon, CP tensor

Headline gain by beam polar angle (azimuth 225 degrees; broadside has no azimuth). Scan loss is relative to broadside of the same variant; a negative value means the gain is higher than at broadside (the conversion efficiency rises with the steering angle in this model). `outside` is the share of the design's march points whose modulation wavenumber lies outside the supercell table; the table supports the beam up to the polar angle 20 degrees (no point outside), and the rows with a larger share are `TABLE_EXTRAPOLATED`.

| theta0 [deg] | eta_conv | outside the table | cp continuous [dBi] | cp 3-bit [dBi] | cp 2-bit [dBi] | cp 1-bit [dBi] | 2-bit scan loss [dB] | 2-bit SLL [dB] | 2-bit pointing [deg] | 2-bit gain over the three azimuths [dBi] |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 0.453 | 0.00 inside the table | 23.60 | 23.37 | 22.67 | 19.71 | 0.00 | -22.8 | 0.00 | 22.67 to 22.67 |
| 10 | 0.471 | 0.00 inside the table | 23.66 | 23.45 | 22.75 | 19.70 | -0.07 | -22.4 | 0.18 | 22.74 to 22.75 |
| 20 | 0.519 | 0.00 inside the table | 23.86 | 23.64 | 22.95 | 19.93 | -0.28 | -22.7 | 0.25 | 22.94 to 22.95 |
| 30 | 0.538 | 0.22 TABLE_EXTRAPOLATED | 23.60 | 23.38 | 22.65 | 19.78 | 0.02 | -22.8 | 0.61 | 22.65 to 22.70 |
| 40 | 0.548 | 0.48 TABLE_EXTRAPOLATED | 22.99 | 22.76 | 22.08 | 19.13 | 0.59 | -21.6 | 0.79 | 22.08 to 22.09 |
| 50 | 0.554 | 0.57 TABLE_EXTRAPOLATED | 21.99 | 21.76 | 21.07 | 18.07 | 1.60 | -19.1 | 1.01 | 21.07 to 21.07 |

## Pitch and aperture

Headline quantisation loss at the design point for the CP family, with the renormalised 2-bit value for comparison:

| pitch [mm] | cells | cp continuous [dBi] | cp 3-bit loss [dB] | cp 2-bit loss [dB] | cp 2-bit loss, renormalised [dB] | cp 1-bit loss [dB] | lpx 1-bit loss, renormalised [dB] | lpx 2-bit loss, renormalised [dB] |
|---|---|---|---|---|---|---|---|---|
| 4.0 | 3646 | 23.57 | 0.22 | 0.91 | 0.06 | 3.89 | 0.07 | 0.01 |
| 5.0 | 2306 | 23.52 | 0.22 | 0.92 | 0.16 | 3.90 | 0.16 | 0.02 |
| 6.4 | 1430 | 23.60 | 0.22 | 0.95 | 0.20 | 3.82 | 0.19 | 0.07 |
| 8.0 | 924 | 23.58 | 0.23 | 0.90 | 0.29 | 3.87 | 0.31 | 0.06 |

The larger reference aperture (hexagon of circumradius 225 mm), pitch 6.4 mm:

| circumradius [mm] | theta0 [deg] | cells | eta_conv | outside the table | cp continuous [dBi] | cp 2-bit [dBi] | cp 1-bit [dBi] | cp 2-bit SLL [dB] |
|---|---|---|---|---|---|---|---|---|
| 225 | 0 | 3210 | 0.574 | 0.00 inside the table | 28.17 | 27.25 | 24.26 | -25.8 |
| 225 | 30 | 3210 | 0.655 | 0.22 TABLE_EXTRAPOLATED | 27.93 | 27.00 | 24.10 | -24.1 |
| 225 | 50 | 3210 | 0.669 | 0.57 TABLE_EXTRAPOLATED | 26.24 | 25.32 | 22.30 | -23.2 |

## Control hardware (150 mm hexagon, beam theta 30 / phi 225)

Counting rules (one method for every row): diodes per cell = bits (one diode per bit, an ASSUMPTION; the topology is item 3); bit = 1 is the forward-biased diode at 10 mA, the share of forward bits is counted from the actual coded states; DC power = (forward diodes) x I_F x V_rail, V_rail = 3.3 V (the supply power, series resistor included; the diode alone dissipates 0.89 V x I_F). The same formula is evaluated at 1 mA. Independent bias lines = diodes (arbitrary patterns need an independent state per diode). Sharing a line among several diodes constrains the pattern: the columns give the line count only, they are NOT simulated, and no row/column or per-ring scheme is evaluated (a per-ring scheme would also not carry the azimuth-dependent slot axes of this family). Shift registers: eight- or sixteen-output devices in chains of eight clocked at 10 MHz; FPGA pins = chains plus the shared clock, latch and enable lines (generic counts, not checked against a specific FPGA board); output current per pin is the diode forward current.

| pitch [mm] | variant | cells | diodes | bias lines | lines at 2 diodes per line | lines at 4 diodes per line | eight-output chips | sixteen-output chips | FPGA pins (eight-output) | update [us] | forward share | DC current [A] | DC power expected, 10 mA [W] | DC power expected, 1 mA [W] | DC power all on, 10 mA [W] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 4.0 | cp_1bit | 3646 | 3646 | 3646 | 1823 | 912 | 456 | 228 | 60 | 6.4 | 0.493 | 18.0 | 59.4 | 5.94 | 120.3 |
| 4.0 | cp_2bit | 3646 | 7292 | 7292 | 3646 | 1823 | 912 | 456 | 117 | 6.4 | 0.495 | 36.1 | 119.2 | 11.92 | 240.6 |
| 4.0 | cp_3bit | 3646 | 10938 | 10938 | 5469 | 2735 | 1368 | 684 | 174 | 6.4 | 0.497 | 54.4 | 179.4 | 17.94 | 361.0 |
| 5.0 | cp_1bit | 2306 | 2306 | 2306 | 1153 | 577 | 289 | 145 | 40 | 6.4 | 0.492 | 11.3 | 37.5 | 3.75 | 76.1 |
| 5.0 | cp_2bit | 2306 | 4612 | 4612 | 2306 | 1153 | 577 | 289 | 76 | 6.4 | 0.495 | 22.8 | 75.3 | 7.53 | 152.2 |
| 5.0 | cp_3bit | 2306 | 6918 | 6918 | 3459 | 1730 | 865 | 433 | 112 | 6.4 | 0.493 | 34.1 | 112.5 | 11.25 | 228.3 |
| 6.4 | cp_1bit | 1430 | 1430 | 1430 | 715 | 358 | 179 | 90 | 26 | 6.4 | 0.498 | 7.1 | 23.5 | 2.35 | 47.2 |
| 6.4 | cp_2bit | 1430 | 2860 | 2860 | 1430 | 715 | 358 | 179 | 48 | 6.4 | 0.493 | 14.1 | 46.6 | 4.66 | 94.4 |
| 6.4 | cp_3bit | 1430 | 4290 | 4290 | 2145 | 1073 | 537 | 269 | 71 | 6.4 | 0.505 | 21.7 | 71.4 | 7.14 | 141.6 |
| 8.0 | cp_1bit | 924 | 924 | 924 | 462 | 231 | 116 | 58 | 18 | 6.4 | 0.497 | 4.6 | 15.1 | 1.51 | 30.5 |
| 8.0 | cp_2bit | 924 | 1848 | 1848 | 924 | 462 | 231 | 116 | 32 | 6.4 | 0.490 | 9.1 | 29.9 | 2.99 | 61.0 |
| 8.0 | cp_3bit | 924 | 2772 | 2772 | 1386 | 693 | 347 | 174 | 47 | 6.4 | 0.492 | 13.7 | 45.0 | 4.50 | 91.5 |

Driver current per output: 10 mA (a sixteen-output device with every output on draws 160 mA, an eight-output device 80 mA): the output stage must source or sink this per pin; whether common logic-level shift registers are specified for it is NOT checked here (datasheets to be checked in item 3), a current-driver class may be needed. The LP rows are identical in count to the CP rows of the same bit depth (`lp` bias tables are in the JSON).

For `lpx` the renormalised loss is a lower bound on the penalty: it hides the power moved out of the beam, and an equal-current-power estimate is not computed here.

The lower-bias column is only a power figure; the diode series resistance rises at low bias (item 3 must trade it against insertion loss).

## Recommendation and trade-offs

- **Polarisation scope: CP, the rotating-slot tensor family.** It is the family the existing forward model covers, with a single HFSS check of the continuous hexagon (one board, one run); the scalar LP surface does not form a pencil beam (finding above), and the `lpx` bound needs a cell that controls the current direction, which no design exists for.
- **Bit depth: 2-bit as the baseline, 3-bit as the fallback.** At pitch 6.4 mm the 2-bit loss is 0.95 dB (fixed source; the renormalised figure 0.20 dB understates it) against 3.82 dB for 1-bit and 0.22 dB for 3-bit. The 2-bit loss uses almost all of the 1.0 dB budget, so the coded surface is within 1.0 dB of the continuous one only if the cell delivers essentially the full table contrast (combined table above); 3-bit costs one more diode per cell (4290 diodes against 2860).
- **Pitch: 6.4 mm** (the pitch of the supercell table and of the HFSS-checked hexagon). The 2-bit headline loss is flat in pitch (from 0.90 to 0.95 dB over the four pitches), so in this model pitch is not a loss lever and the choice is set by the diode count and by what item 3 can build; 8.0 mm is the alternative worth testing with a real cell (1848 diodes against 2860).
- **Diode count at the recommendation: 2860 diodes (1430 cells x 2), 2860 independent bias lines**, 179 sixteen-output shift registers, 26 FPGA pins at the stated chain length; expected DC power 46.6 W at 10 mA per forward diode (4.66 W at 1 mA). The bias power, not the diode count, is the sharpest cost of this architecture.
- **Required contrast (beam 30 degrees, contrast loss alone within 1.0 dB):** tensor anisotropy 2.296 mS, i.e. a principal-axis reactance difference of 92.3 ohm, present in every state, with the axis rotating between 4 orientations; the equivalent scalar swing would be 148.4 ohm. These are sheet-level numbers; whether a diode-loaded subwavelength cell reaches them is the question of item 3. With the 2-bit loss included the requirement is larger: a total loss of 1.5 dB needs a tensor anisotropy of 2.500 mS (102.6 ohm), and a total of 1.0 dB needs 3.027 mS (132.0 ohm), close to the table maximum of the leakage constant.
- **Scan range:** the supercell table supports the beam up to the polar angle 20 degrees; the design point has 0.22 of its march points outside the table and the larger angles more (scan table), so the gain at those angles is `TABLE_EXTRAPOLATED`.
- **Reliability of the sidelobe and cross-polar numbers:** the one HFSS run has a sidelobe level of -18.5 dB against the model's -22.8 dB; treat the model's absolute sidelobe levels as optimistic.

## What is not available or not established

- 10.0 GHz: `NOT_AVAILABLE_FROM_CURRENT_DATA (the supercell table exists at 10.54 GHz only)`; a table at that frequency is a server job, not run.
- Everything is `MODEL_ONLY`: idealised cells, the table assumed at all pitches, no coupling between neighbouring states, no diode loss, the unconverted surface wave and the feed's direct radiation assumed outside the beam (Phase 26 assumptions).
- Bias-line sharing schemes (row/column, per-ring) and the fit to a specific FPGA board are not evaluated; only generic counts are given.
- The literature table (item 1) does not support any of the contrast, bit-depth or bias statements here; they are model results.

## Stop point

No cell has been designed. Approve or change the pitch, bit depth, polarisation scope and diode count above before item 3 starts.
