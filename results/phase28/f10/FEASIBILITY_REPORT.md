# Phase 28 item 2 - feasibility of a coded monopole-fed hologram (stop point)

**Labels:** `MODEL_ONLY`, `IDEALISED_CELLS`, `NOT_HFSS_VALIDATED`, `SUPERCELL_TABLE_6P4MM_10P0GHZ_ASSUMED_AT_ALL_PITCHES`. **Generated from** `results/phase28/feasibility_results.json` only.

Forward model: the Phase 26 tensor-CP model unchanged (supercell table, ray march, aperture far field on the grounded slab). Frequency 10.00 GHz. 10 GHz is the PI's operating frequency (2026-10-26); the 10.54 GHz Phase 26 table and the item-2 sweep at that frequency are history. The launch fraction is the 10.54 GHz Phase 26 value (NOT recomputed at this frequency): absolute gains carry that assumption, losses and contrast do not. Cells are idealised (a state is a tensor or scalar value, no geometry); the table, computed for 6.4 mm cells, is assumed at every pitch, so the pitch sweep changes only how the hologram is sampled. Gain is referred to the accepted power with the Phase 26 bare launch fraction 0.592 (the calibrated launch fraction 0.62 of the four-hexagon work would add 0.20 dB to every gain here, and nothing to any loss).

## Anchor: no quantisation reproduces the Phase 26 single-hexagon prediction

Source: results/phase26/tensor_cp_hexagon.json (10.54 GHz). The committed Phase 26 prediction is at 10.54 GHz; this code path was anchored to it in the 10.54 GHz sweep (results/phase28/feasibility_results.json). There is nothing to compare at this frequency.

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

**Method.** The supercell table's leakage constants are capped at `alpha_cap` and the single-hexagon CP design is re-marched; `alpha_cap` is the value at which the gain falls 1.0 dB below the uncapped design (table maximum 4.123 Np/m, minimum 0.007 Np/m). It is converted to physical units with the Phase 26 closed forms at the mean `beta/k0` of the design's cells; the scalar case is cross-checked with the exact truncated-Floquet complex beta. This row is the contrast-only requirement: it ignores the quantisation loss, which is added in the combined table below.

| beam theta / phi [deg] | march points outside the table | alpha_cap [Np/m] | loss re-check [dB] | beta/k0 mean | tensor delta [mS] | tensor dX_aniso [ohm] | tensor d_beta/k0 | tensor depth | LP binary dB [mS] | LP dX_switch [ohm] | LP d_beta/k0 | LP four-level dB [mS] | LP four-level dX_switch [ohm] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 / 0 | 0.00 inside the table (closed form indicative only) | 1.899 | 0.997 | 1.0511 | 2.335 | 75.6 | 0.0416 | 0.285 | 5.188 | 85.8 | 0.0463 | 6.250 | 108.8 |
| 30 / 225 | 0.00 inside the table | 2.089 | 0.979 | 1.0527 | 2.444 | 82.6 | 0.0448 | 0.303 | 6.020 | 107.3 | 0.0553 | 7.252 | 139.5 |
| 50 / 225 | 0.01 inside the table | 2.178 | 0.982 | 1.0541 | 2.427 | 84.4 | 0.0456 | 0.305 | 8.030 | 169.9 | 0.0755 | 9.674 | 241.9 |

Exact-Floquet cross-check of the scalar closed form (binary sheet at the closed-form full swing found for each beam; TM, collinear):

| beam theta / phi [deg] | alpha closed form [Np/m] | alpha Floquet exact [Np/m] | ratio exact/closed | note |
|---|---|---|---|---|
| 0 / 0 | 1.899 | n/a | n/a | NOT_AVAILABLE: the Floquet balance has no finite root at broadside (open stop band) |
| 30 / 225 | 2.089 | 2.200 | 1.053 | - |
| 50 / 225 | 2.178 | 2.765 | 1.269 | - |

A ratio away from one means the closed form under- or overstates the leakage for that swing (finite modulation depth); the contrast numbers of that row carry that uncertainty.

**Combined requirement: 2-bit coded surface.** The quantisation loss adds to the contrast loss. Each capped design is also quantised to four slot axes; the table gives the smallest `alpha_cap` whose total loss against the uncapped continuous design stays within the tolerance, converted to the tensor and scalar contrast (four-level swing). The 2-bit loss of the uncapped design is the first column. The cap is interpolated linearly in loss between the grid points of the cap sweep, so it can be slightly conservative (a larger cap than strictly needed).

| beam theta / phi [deg] | 2-bit loss, uncapped [dB] | total loss tolerance [dB] | alpha_cap [Np/m] | tensor delta [mS] | tensor dX_aniso [ohm] | scalar four-level dX_switch [ohm] | caveat |
|---|---|---|---|---|---|---|---|
| 0 / 0 | 0.91 | 1.0 | 2.598 | 2.732 | 91.4 | 135.7 | inside the table; closed form indicative only |
| 0 / 0 | 0.91 | 1.5 | 2.177 | 2.500 | 82.0 | 119.4 | inside the table; closed form indicative only |
| 0 / 0 | 0.91 | 2.0 | 1.845 | 2.302 | 74.3 | 106.7 | inside the table; closed form indicative only |
| 30 / 225 | 0.90 | 1.0 | 3.661 | 3.236 | 118.4 | 228.1 | inside the table |
| 30 / 225 | 0.90 | 1.5 | 2.443 | 2.643 | 90.9 | 157.6 | inside the table |
| 30 / 225 | 0.90 | 2.0 | 2.017 | 2.402 | 80.9 | 135.8 | inside the table |
| 50 / 225 | 0.85 | 1.0 | 3.616 | 3.127 | 116.7 | 507.5 | inside the table; Floquet exact/closed ratio 1.269 |
| 50 / 225 | 0.85 | 1.5 | 2.540 | 2.621 | 92.8 | 289.3 | inside the table; Floquet exact/closed ratio 1.269 |
| 50 / 225 | 0.85 | 2.0 | 2.067 | 2.364 | 81.8 | 228.7 | inside the table; Floquet exact/closed ratio 1.269 |

For scale only (the cell-level contrast depends on how the diode is embedded, which is item 3): the terminal reactance swing of the diode models at the design frequency.

| diode (datasheet, PROVISIONAL) | R_on [ohm] | X_on [ohm] | X_off [ohm] | X_off - X_on [ohm] |
|---|---|---|---|---|
| MADP-000907-14020 | 5.20 | 0.0 | -636.6 | -636.6 |
| SMP1345-079LF | 1.50 | 44.0 | -39.8 | -83.8 |
| SMP1345-079LF_1mA | 3.50 | 44.0 | -39.8 | -83.8 |

MADP-000907-14020 package inductance is NOT_AVAILABLE from the datasheet (L_s = 0 here); datasheet values are PROVISIONAL.

**Table range caveat.** The supercell table covers a limited range of modulation wavenumbers; rows with `TABLE_EXTRAPOLATED` above (and in the scan table below) have a large share of march points outside it, and their contrast and scan numbers are extrapolations.

## Design point: 150 mm hexagon, beam theta 30 / phi 225 degrees

This beam has 0.00 of its march points outside the supercell table: `inside the table`; the design-point numbers share that caveat.

Quantisation loss is relative to the continuous design of the same family and pitch. Two definitions are given because they differ for a phase-only quantisation:

- *fixed source*: the same cell currents, so the peak field scales as the pattern's directivity times its radiated power. A phase-only quantisation (`cp`) leaves every current magnitude unchanged, and the quantisation error is radiated into the non-radiating part of the spectrum; that power does not contribute to the beam (where it goes is not modelled), and this is the loss used as the headline for `cp`.
- *renormalised*: the directivity of the quantised pattern alone, i.e. its own radiated power renormalised to the launched power. It hides the lost power. It is the energy-consistent one for an amplitude quantisation (`lp`, `lpx`), where the currents change magnitude and would change the leakage, and is the headline for those.

`cp` = rotating-slot tensor, RHCP (the Phase 26 family); `lp` = scalar surface, current along the radial direction; `lpx` = linear polarisation with IDEAL polarisation control (current along x in every cell), an upper bound no cell has been shown to realise. AR and XPD beyond the bound are shown as bounds, never as numbers. A negative renormalised loss is a quantised pattern whose directivity is higher than the continuous one on this sampling; a negative fixed-source loss for `lp` and `lpx` is not a gain: the plus and minus one levels have larger magnitude than the continuous modulation, so the radiated power rises (radiated power change column).

| pitch [mm] | variant | gain, headline [dBi] | loss, headline [dB] | loss, fixed source [dB] | loss, renormalised [dB] | radiated power change [dB] | SLL [dB] | pointing error [deg] | XPD [dB] | AR [dB] |
|---|---|---|---|---|---|---|---|---|---|---|
| 4.0 | cp_cont | 22.83 | 0.00 | 0.00 | 0.00 | 0.00 | -21.7 | 0.44 | 17.01 | 2.47 |
| 4.0 | cp_1bit | 18.93 | 3.90 | 3.90 | 2.99 | -0.91 | -19.0 | 0.44 | 15.29 | 3.02 |
| 4.0 | cp_2bit | 21.92 | 0.91 | 0.91 | 0.06 | -0.85 | -21.3 | 0.44 | 17.06 | 2.45 |
| 4.0 | cp_3bit | 22.60 | 0.23 | 0.23 | 0.02 | -0.21 | -21.7 | 0.44 | 16.96 | 2.48 |
| 4.0 | lp_cont | 17.71 | 0.00 | 0.00 | 0.00 | 0.00 | -1.1 | 5.39 | - | - |
| 4.0 | lp_1bit | 17.71 | 0.00 | -2.12 | 0.00 | 2.12 | -1.2 | 5.39 | - | - |
| 4.0 | lp_2bit | 17.73 | -0.02 | -0.55 | -0.02 | 0.53 | -1.2 | 5.39 | - | - |
| 4.0 | lpx_cont | 22.91 | 0.00 | 0.00 | 0.00 | 0.00 | -21.9 | 0.44 | - | - |
| 4.0 | lpx_1bit | 22.84 | 0.07 | -2.09 | 0.07 | 2.16 | -21.8 | 0.44 | - | - |
| 4.0 | lpx_2bit | 22.89 | 0.02 | -0.54 | 0.02 | 0.56 | -21.5 | 0.44 | - | - |
| 5.0 | cp_cont | 22.79 | 0.00 | 0.00 | 0.00 | 0.00 | -21.7 | 0.44 | 17.13 | 2.43 |
| 5.0 | cp_1bit | 18.90 | 3.89 | 3.89 | 3.02 | -0.87 | -18.2 | 0.44 | 15.35 | 3.00 |
| 5.0 | cp_2bit | 21.87 | 0.92 | 0.92 | 0.10 | -0.81 | -21.9 | 0.44 | 17.46 | 2.34 |
| 5.0 | cp_3bit | 22.57 | 0.22 | 0.22 | 0.02 | -0.19 | -21.7 | 0.44 | 16.93 | 2.49 |
| 5.0 | lp_cont | 17.69 | 0.00 | 0.00 | 0.00 | 0.00 | -1.1 | 5.53 | - | - |
| 5.0 | lp_1bit | 17.60 | 0.09 | -2.10 | 0.09 | 2.19 | -1.1 | 5.39 | - | - |
| 5.0 | lp_2bit | 17.68 | 0.02 | -0.45 | 0.02 | 0.46 | -0.6 | 5.38 | - | - |
| 5.0 | lpx_cont | 22.89 | 0.00 | 0.00 | 0.00 | 0.00 | -21.7 | 0.44 | - | - |
| 5.0 | lpx_1bit | 22.79 | 0.10 | -2.10 | 0.10 | 2.20 | -21.6 | 0.44 | - | - |
| 5.0 | lpx_2bit | 22.88 | 0.01 | -0.44 | 0.01 | 0.45 | -21.7 | 0.44 | - | - |
| 6.4 | cp_cont | 22.85 | 0.00 | 0.00 | 0.00 | 0.00 | -21.8 | 0.44 | 16.97 | 2.48 |
| 6.4 | cp_1bit | 18.96 | 3.90 | 3.90 | 3.13 | -0.77 | -18.3 | 0.44 | 14.93 | 3.15 |
| 6.4 | cp_2bit | 21.95 | 0.90 | 0.90 | 0.18 | -0.72 | -23.1 | 0.62 | 16.98 | 2.48 |
| 6.4 | cp_3bit | 22.62 | 0.23 | 0.23 | 0.03 | -0.20 | -22.3 | 0.44 | 17.00 | 2.47 |
| 6.4 | lp_cont | 17.77 | 0.00 | 0.00 | 0.00 | 0.00 | -1.2 | 5.29 | - | - |
| 6.4 | lp_1bit | 17.63 | 0.14 | -2.10 | 0.14 | 2.24 | -1.2 | 5.39 | - | - |
| 6.4 | lp_2bit | 17.72 | 0.04 | -0.47 | 0.04 | 0.52 | -1.1 | 5.39 | - | - |
| 6.4 | lpx_cont | 22.95 | 0.00 | 0.00 | 0.00 | 0.00 | -22.2 | 0.44 | - | - |
| 6.4 | lpx_1bit | 22.76 | 0.18 | -2.14 | 0.18 | 2.32 | -21.2 | 0.61 | - | - |
| 6.4 | lpx_2bit | 22.91 | 0.03 | -0.50 | 0.03 | 0.54 | -21.9 | 0.44 | - | - |
| 8.0 | cp_cont | 22.85 | 0.00 | 0.00 | 0.00 | 0.00 | -21.5 | 0.44 | 16.73 | 2.55 |
| 8.0 | cp_1bit | 18.92 | 3.93 | 3.93 | 3.18 | -0.75 | -18.7 | 0.62 | 14.97 | 3.14 |
| 8.0 | cp_2bit | 21.94 | 0.91 | 0.91 | 0.34 | -0.57 | -20.4 | 0.44 | 15.99 | 2.78 |
| 8.0 | cp_3bit | 22.62 | 0.23 | 0.23 | 0.07 | -0.16 | -21.8 | 0.44 | 17.49 | 2.33 |
| 8.0 | lp_cont | 17.74 | 0.00 | 0.00 | 0.00 | 0.00 | -1.2 | 5.29 | - | - |
| 8.0 | lp_1bit | 17.54 | 0.20 | -2.20 | 0.20 | 2.40 | -1.4 | 5.45 | - | - |
| 8.0 | lp_2bit | 17.66 | 0.08 | -0.66 | 0.08 | 0.74 | -1.3 | 5.45 | - | - |
| 8.0 | lpx_cont | 22.94 | 0.00 | 0.00 | 0.00 | 0.00 | -22.0 | 0.44 | - | - |
| 8.0 | lpx_1bit | 22.62 | 0.33 | -2.06 | 0.33 | 2.39 | -22.0 | 0.44 | - | - |
| 8.0 | lpx_2bit | 22.84 | 0.10 | -0.64 | 0.10 | 0.74 | -22.4 | 0.44 | - | - |

Peak search: 0.25 degree direction-cosine refinement around the coarse-grid peak; the sidelobe level uses the strongest coarse-grid sidelobe refined the same way. The headline gain is the continuous design's gain (aperture directivity plus ten times the base-ten logarithm of launch fraction x eta_conv) minus the headline loss; eta_conv of this design is given in the scan table.

## Finding: a scalar surface on a monopole cannot form a linear-polarised pencil beam

With the induced current along the surface wave's own (radial) direction the continuous scalar design has a pointing error of 5.29 degrees and a strongest sidelobe at -1.2 dB: the aperture field is radially polarised, so the contributions cancel in the beam direction (a ring or split beam, not a pencil). Gain 17.77 dBi against 22.95 dBi for the same hologram with ideal polarisation control (`lpx`) and 22.85 dBi for the CP tensor design. The `lp` rows above are therefore NOT a candidate; the bit-depth conclusion for LP rests on `lpx`, which needs a cell that sets the current direction (not designed).

## Steering: pitch 6.4 mm, 150 mm hexagon, CP tensor

Headline gain by beam polar angle (azimuth 225 degrees; broadside has no azimuth). Scan loss is relative to broadside of the same variant; a negative value means the gain is higher than at broadside (the conversion efficiency rises with the steering angle in this model). `outside` is the share of the design's march points whose modulation wavenumber lies outside the supercell table; the table supports the beam up to the polar angle 50 degrees (no point far from every column of the table), and the rows with a larger share are `TABLE_EXTRAPOLATED`. Counting only the march points past the largest tabulated wavenumber, the supported angle is 40 degrees; above it a small share of the march points (the `outside` column) lies past the end of the table.

For this run the share also counts the march points whose modulation wavenumber is farther than fifteen percent from every column of the table (a gap between two columns, not only the ends of the table): the larger of the two shares is printed.

| theta0 [deg] | eta_conv | outside the table | cp continuous [dBi] | cp 3-bit [dBi] | cp 2-bit [dBi] | cp 1-bit [dBi] | 2-bit scan loss [dB] | 2-bit SLL [dB] | 2-bit pointing [deg] | 2-bit gain over the three azimuths [dBi] |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 0.473 | 0.00 inside the table | 23.34 | 23.11 | 22.43 | 19.44 | 0.00 | -23.9 | 0.00 | 22.43 to 22.43 |
| 10 | 0.443 | 0.00 inside the table | 22.94 | 22.72 | 22.12 | 18.92 | 0.31 | -22.5 | 0.18 | 22.02 to 22.12 |
| 20 | 0.472 | 0.00 inside the table | 22.98 | 22.76 | 22.06 | 19.07 | 0.37 | -22.2 | 0.38 | 22.06 to 22.06 |
| 30 | 0.504 | 0.00 inside the table | 22.85 | 22.62 | 21.95 | 18.96 | 0.48 | -23.1 | 0.62 | 21.95 to 21.96 |
| 40 | 0.514 | 0.00 inside the table | 22.25 | 22.03 | 21.35 | 18.27 | 1.08 | -21.2 | 0.79 | 21.35 to 21.35 |
| 50 | 0.514 | 0.01 inside the table | 21.22 | 21.00 | 20.37 | 17.27 | 2.06 | -19.3 | 1.55 | 20.32 to 20.37 |

## Pitch and aperture

Headline quantisation loss at the design point for the CP family, with the renormalised 2-bit value for comparison:

| pitch [mm] | cells | cp continuous [dBi] | cp 3-bit loss [dB] | cp 2-bit loss [dB] | cp 2-bit loss, renormalised [dB] | cp 1-bit loss [dB] | lpx 1-bit loss, renormalised [dB] | lpx 2-bit loss, renormalised [dB] |
|---|---|---|---|---|---|---|---|---|
| 4.0 | 3646 | 22.83 | 0.23 | 0.91 | 0.06 | 3.90 | 0.07 | 0.02 |
| 5.0 | 2306 | 22.79 | 0.22 | 0.92 | 0.10 | 3.89 | 0.10 | 0.01 |
| 6.4 | 1430 | 22.85 | 0.23 | 0.90 | 0.18 | 3.90 | 0.18 | 0.03 |
| 8.0 | 924 | 22.85 | 0.23 | 0.91 | 0.34 | 3.93 | 0.33 | 0.10 |

The larger reference aperture (hexagon of circumradius 225 mm), pitch 6.4 mm:

| circumradius [mm] | theta0 [deg] | cells | eta_conv | outside the table | cp continuous [dBi] | cp 2-bit [dBi] | cp 1-bit [dBi] | cp 2-bit SLL [dB] |
|---|---|---|---|---|---|---|---|---|
| 225 | 0 | 3210 | 0.595 | 0.00 inside the table | 27.86 | 26.94 | 23.94 | -27.6 |
| 225 | 30 | 3210 | 0.620 | 0.00 inside the table | 27.23 | 26.33 | 23.34 | -25.2 |
| 225 | 50 | 3210 | 0.629 | 0.01 inside the table | 25.53 | 24.66 | 21.61 | -22.3 |

## Control hardware (150 mm hexagon, beam theta 30 / phi 225)

Counting rules (one method for every row): diodes per cell = bits (one diode per bit, an ASSUMPTION; the topology is item 3); bit = 1 is the forward-biased diode at 10 mA, the share of forward bits is counted from the actual coded states; DC power = (forward diodes) x I_F x V_rail, V_rail = 3.3 V (the supply power, series resistor included; the diode alone dissipates 0.89 V x I_F). The same formula is evaluated at 1 mA. Independent bias lines = diodes (arbitrary patterns need an independent state per diode). Sharing a line among several diodes constrains the pattern: the columns give the line count only, they are NOT simulated, and no row/column or per-ring scheme is evaluated (a per-ring scheme would also not carry the azimuth-dependent slot axes of this family). Shift registers: eight- or sixteen-output devices in chains of eight clocked at 10 MHz; FPGA pins = chains plus the shared clock, latch and enable lines (generic counts, not checked against a specific FPGA board); output current per pin is the diode forward current.

| pitch [mm] | variant | cells | diodes | bias lines | lines at 2 diodes per line | lines at 4 diodes per line | eight-output chips | sixteen-output chips | FPGA pins (eight-output) | update [us] | forward share | DC current [A] | DC power expected, 10 mA [W] | DC power expected, 1 mA [W] | DC power all on, 10 mA [W] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 4.0 | cp_1bit | 3646 | 3646 | 3646 | 1823 | 912 | 456 | 228 | 60 | 6.4 | 0.491 | 17.9 | 59.1 | 5.91 | 120.3 |
| 4.0 | cp_2bit | 3646 | 7292 | 7292 | 3646 | 1823 | 912 | 456 | 117 | 6.4 | 0.499 | 36.4 | 120.0 | 12.00 | 240.6 |
| 4.0 | cp_3bit | 3646 | 10938 | 10938 | 5469 | 2735 | 1368 | 684 | 174 | 6.4 | 0.499 | 54.6 | 180.1 | 18.01 | 361.0 |
| 5.0 | cp_1bit | 2306 | 2306 | 2306 | 1153 | 577 | 289 | 145 | 40 | 6.4 | 0.486 | 11.2 | 37.0 | 3.70 | 76.1 |
| 5.0 | cp_2bit | 2306 | 4612 | 4612 | 2306 | 1153 | 577 | 289 | 76 | 6.4 | 0.495 | 22.8 | 75.4 | 7.54 | 152.2 |
| 5.0 | cp_3bit | 2306 | 6918 | 6918 | 3459 | 1730 | 865 | 433 | 112 | 6.4 | 0.496 | 34.3 | 113.2 | 11.32 | 228.3 |
| 6.4 | cp_1bit | 1430 | 1430 | 1430 | 715 | 358 | 179 | 90 | 26 | 6.4 | 0.492 | 7.0 | 23.2 | 2.32 | 47.2 |
| 6.4 | cp_2bit | 1430 | 2860 | 2860 | 1430 | 715 | 358 | 179 | 48 | 6.4 | 0.498 | 14.2 | 47.0 | 4.70 | 94.4 |
| 6.4 | cp_3bit | 1430 | 4290 | 4290 | 2145 | 1073 | 537 | 269 | 71 | 6.4 | 0.498 | 21.4 | 70.6 | 7.06 | 141.6 |
| 8.0 | cp_1bit | 924 | 924 | 924 | 462 | 231 | 116 | 58 | 18 | 6.4 | 0.500 | 4.6 | 15.2 | 1.52 | 30.5 |
| 8.0 | cp_2bit | 924 | 1848 | 1848 | 924 | 462 | 231 | 116 | 32 | 6.4 | 0.486 | 9.0 | 29.6 | 2.96 | 61.0 |
| 8.0 | cp_3bit | 924 | 2772 | 2772 | 1386 | 693 | 347 | 174 | 47 | 6.4 | 0.499 | 13.8 | 45.6 | 4.56 | 91.5 |

Driver current per output: 10 mA (a sixteen-output device with every output on draws 160 mA, an eight-output device 80 mA): the output stage must source or sink this per pin; whether common logic-level shift registers are specified for it is NOT checked here (datasheets to be checked in item 3), a current-driver class may be needed. The LP rows are identical in count to the CP rows of the same bit depth (`lp` bias tables are in the JSON).

For `lpx` the renormalised loss is a lower bound on the penalty: it hides the power moved out of the beam, and an equal-current-power estimate is not computed here.

The lower-bias column is only a power figure; the diode series resistance rises at low bias (item 3 must trade it against insertion loss).

## Recommendation and trade-offs

- **Polarisation scope: CP, the rotating-slot tensor family.** It is the family the existing forward model covers, with a single HFSS check of the continuous hexagon (one board, one run); the scalar LP surface does not form a pencil beam (finding above), and the `lpx` bound needs a cell that controls the current direction, which no design exists for.
- **Bit depth: 2-bit as the baseline, 3-bit as the fallback.** At pitch 6.4 mm the 2-bit loss is 0.90 dB (fixed source; the renormalised figure 0.18 dB understates it) against 3.90 dB for 1-bit and 0.23 dB for 3-bit. The 2-bit loss uses almost all of the 1.0 dB budget, so the coded surface is within 1.0 dB of the continuous one only if the cell delivers essentially the full table contrast (combined table above); 3-bit costs one more diode per cell (4290 diodes against 2860).
- **Pitch: 6.4 mm** (the pitch of the supercell table and of the HFSS-checked hexagon). The 2-bit headline loss is flat in pitch (from 0.90 to 0.92 dB over the four pitches), so in this model pitch is not a loss lever and the choice is set by the diode count and by what item 3 can build; 8.0 mm is the alternative worth testing with a real cell (1848 diodes against 2860).
- **Diode count at the recommendation: 2860 diodes (1430 cells x 2), 2860 independent bias lines**, 179 sixteen-output shift registers, 26 FPGA pins at the stated chain length; expected DC power 47.0 W at 10 mA per forward diode (4.70 W at 1 mA). The bias power, not the diode count, is the sharpest cost of this architecture.
- **Required contrast (beam 30 degrees, contrast loss alone within 1.0 dB):** tensor anisotropy 2.444 mS, i.e. a principal-axis reactance difference of 82.6 ohm, present in every state, with the axis rotating between 4 orientations; the equivalent scalar swing would be 139.5 ohm. These are sheet-level numbers; whether a diode-loaded subwavelength cell reaches them is the question of item 3. With the 2-bit loss included the requirement is larger: a total loss of 1.5 dB needs a tensor anisotropy of 2.643 mS (90.9 ohm), and a total of 1.0 dB needs 3.236 mS (118.4 ohm), close to the table maximum of the leakage constant.
- **Scan range:** the supercell table supports the beam up to the polar angle 50 degrees; the design point has 0.00 of its march points outside the table and the larger angles more (scan table), so the gain at those angles is `TABLE_EXTRAPOLATED`.
- **Reliability of the sidelobe and cross-polar numbers:** the one HFSS run has a sidelobe level of -18.5 dB against the model's -22.8 dB; treat the model's absolute sidelobe levels as optimistic.

## What is not available or not established

- 10 GHz is the PI's operating frequency (2026-10-26); the 10.54 GHz Phase 26 table and the item-2 sweep at that frequency are history. The launch fraction is the 10.54 GHz Phase 26 value (NOT recomputed at this frequency): absolute gains carry that assumption, losses and contrast do not.
- Everything is `MODEL_ONLY`: idealised cells, the table assumed at all pitches, no coupling between neighbouring states, no diode loss, the unconverted surface wave and the feed's direct radiation assumed outside the beam (Phase 26 assumptions).
- Bias-line sharing schemes (row/column, per-ring) and the fit to a specific FPGA board are not evaluated; only generic counts are given.
- The literature table (item 1) does not support any of the contrast, bit-depth or bias statements here; they are model results.

## Stop point

No cell has been designed. Approve or change the pitch, bit depth, polarisation scope and diode count above before item 3 starts.
