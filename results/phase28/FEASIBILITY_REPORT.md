# Phase 28 item 2 - feasibility of a coded monopole-fed hologram (stop point)

**Labels:** `MODEL_ONLY`, `IDEALISED_CELLS`, `NOT_HFSS_VALIDATED`, `SUPERCELL_TABLE_6P4MM_10P54GHZ_ASSUMED_AT_ALL_PITCHES`. **Generated from** `results/phase28/feasibility_results.json` only.

Forward model: the Phase 26 tensor-CP model unchanged (supercell table, ray march, aperture far field on the grounded slab). Frequency 10.54 GHz, the only frequency of the supercell table; 10.0 GHz is `NOT_AVAILABLE_FROM_CURRENT_DATA (the supercell table exists at 10.54 GHz only)`. Cells are idealised (a state is a tensor or scalar value, no geometry); the table, computed for 6.4 mm cells, is assumed at every pitch, so the pitch sweep changes only how the hologram is sampled. Gain is referred to the accepted power with the Phase 26 bare launch fraction 0.592; quantisation changes the aperture directivity only (the conversion efficiency is held at the continuous design's value, an assumption).

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

## Switching contrast, first

**Definition.** The contrast is what the cell must be able to change so that the coded hologram can leak as strongly as the continuous one. Four measures, all at the design frequency:

- `alpha_cap`: the largest local leakage constant [Np/m] the coded surface can realise. The continuous Phase 26 design asks for a profile alpha(x) up to the supercell table's maximum; a coded surface that cannot reach it converts less of the launched surface wave.
- Tensor (CP, rotating axis) cell: the fixed anisotropy `delta` [mS] of `Y = Ybar I + j delta R(a) diag(+,-) R(a)^T`, with alpha equal to a constant times delta squared, from the Phase 26 closed form (`leakage_alpha_first_order_tensor`). `dX_aniso` [ohm] is the difference of the TM surface reactances `X = Im(1/(Y + Y_slab,TM(beta)))` seen along the two principal axes; the switch rotates the axis and does not change it. The largest state-to-state change of the sheet tensor (axes perpendicular) is twice delta.
- Scalar cell (also the amplitude control of the `lpx` bound), states `Ybar +- j dB/2` (`dB` the full swing between the extreme levels): the first harmonic of the quantised surface is c1, one quarter of dB times g (g the first-harmonic gain of the level set), and alpha is a constant times c1 squared (`leakage_alpha_first_order`, TM and TE parts, the larger `dB` reported). `dX_switch` [ohm] is the TM surface-reactance difference of the two extreme states; `d_beta` the difference of the TM0 phase constants of the two uniform states.
- Modulation depth: `delta/|Y_down|` (tensor) or half of dB over `|Y_down|` (scalar), `Y_down = Ybar + Y_slab,TM(beta)`, the Phase 26 `required_map` definition.

**Method.** The supercell table's leakage constants are capped at `alpha_cap` and the single-hexagon CP design is re-marched; `alpha_cap` is the value at which the gain falls 1.0 dB below the uncapped design (table maximum 4.259 Np/m, minimum 0.422 Np/m). It is converted to physical units with the Phase 26 closed forms at the mean `beta/k0` of the design's cells; the scalar case is cross-checked with the exact truncated-Floquet complex beta.

| beam theta / phi [deg] | alpha_cap [Np/m] | loss re-check [dB] | beta/k0 mean | tensor delta [mS] | tensor dX_aniso [ohm] | tensor d_beta/k0 | tensor depth | LP binary dB [mS] | LP dX_switch [ohm] | LP d_beta/k0 | LP four-level dB [mS] | LP four-level dX_switch [ohm] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 / 0 | 1.787 | 0.999 | 1.0593 | 2.054 | 76.8 | 0.0425 | 0.270 | 4.563 | 86.9 | 0.0472 | 5.498 | 109.6 |
| 30 / 225 | 2.308 | 0.996 | 1.0621 | 2.296 | 92.3 | 0.0494 | 0.310 | 5.448 | 114.5 | 0.0585 | 6.564 | 148.4 |
| 50 / 225 | 2.412 | 0.995 | 1.0623 | 2.337 | 94.7 | 0.0504 | 0.316 | 7.224 | 172.9 | 0.0772 | 8.704 | 242.5 |

Exact-Floquet cross-check of the scalar closed form (binary sheet at the closed-form full swing found for each beam; TM, collinear):

| beam theta / phi [deg] | alpha closed form [Np/m] | alpha Floquet exact [Np/m] | ratio exact/closed | note |
|---|---|---|---|---|
| 0 / 0 | 1.787 | n/a | n/a | NOT_AVAILABLE: the Floquet balance has no finite root at broadside (open stop band) |
| 30 / 225 | 2.308 | 2.343 | 1.015 | - |
| 50 / 225 | 2.412 | 2.867 | 1.189 | - |

For scale only (the cell-level contrast depends on how the diode is embedded, which is item 3): the terminal reactance swing of the diode models at the design frequency.

| diode (datasheet, PROVISIONAL) | R_on [ohm] | X_on [ohm] | X_off [ohm] | X_off - X_on [ohm] |
|---|---|---|---|---|
| MADP-000907-14020 | 5.20 | 0.0 | -604.0 | -604.0 |
| SMP1345-079LF | 1.50 | 46.4 | -33.1 | -79.5 |
| SMP1345-079LF_1mA | 3.50 | 46.4 | -33.1 | -79.5 |

MADP-000907-14020 package inductance is NOT_AVAILABLE from the datasheet (L_s = 0 here); datasheet values are PROVISIONAL.

**Table range caveat.** The supercell table covers a limited range of modulation wavenumbers; the share of the design's march points outside it is given in the scan table below and the contrast rows for steering beyond the table are `TABLE_EXTRAPOLATED`.

## Design point: 150 mm hexagon, beam theta 30 / phi 225 degrees

Quantisation loss is relative to the continuous design of the same family and pitch. `cp` = rotating-slot tensor, RHCP (the Phase 26 family); `lp` = scalar surface, current along the radial direction; `lpx` = linear polarisation with IDEAL polarisation control (current along x in every cell), an upper bound no cell has been shown to realise. AR and XPD beyond the bound are shown as bounds, never as numbers.

| pitch [mm] | variant | gain [dBi] | quantisation loss [dB] | SLL [dB] | pointing error [deg] | XPD [dB] | AR [dB] |
|---|---|---|---|---|---|---|---|
| 4.0 | cp_cont | 23.57 | 0.00 | -22.6 | 0.44 | 17.00 | 2.47 |
| 4.0 | cp_1bit | 20.57 | 3.00 | -19.0 | 0.44 | 15.71 | 2.87 |
| 4.0 | cp_2bit | 23.51 | 0.06 | -22.5 | 0.44 | 17.12 | 2.44 |
| 4.0 | lp_cont | 18.53 | 0.00 | -2.5 | 5.06 | - | - |
| 4.0 | lp_1bit | 18.48 | 0.04 | -2.5 | 5.06 | - | - |
| 4.0 | lp_2bit | 18.55 | -0.03 | -2.5 | 5.06 | - | - |
| 4.0 | lpx_cont | 23.68 | 0.00 | -21.6 | 0.44 | - | - |
| 4.0 | lpx_1bit | 23.61 | 0.07 | -21.0 | 0.44 | - | - |
| 4.0 | lpx_2bit | 23.67 | 0.01 | -21.6 | 0.44 | - | - |
| 5.0 | cp_cont | 23.52 | 0.00 | -22.5 | 0.44 | 16.93 | 2.49 |
| 5.0 | cp_1bit | 20.45 | 3.07 | -18.7 | 0.44 | 15.76 | 2.86 |
| 5.0 | cp_2bit | 23.36 | 0.16 | -23.0 | 0.44 | 17.01 | 2.47 |
| 5.0 | lp_cont | 18.48 | 0.00 | -2.6 | 5.17 | - | - |
| 5.0 | lp_1bit | 18.22 | 0.26 | -2.1 | 5.16 | - | - |
| 5.0 | lp_2bit | 18.48 | 0.00 | -2.7 | 5.17 | - | - |
| 5.0 | lpx_cont | 23.63 | 0.00 | -21.7 | 0.44 | - | - |
| 5.0 | lpx_1bit | 23.47 | 0.16 | -20.0 | 0.44 | - | - |
| 5.0 | lpx_2bit | 23.61 | 0.02 | -21.9 | 0.44 | - | - |
| 6.4 | cp_cont | 23.60 | 0.00 | -22.9 | 0.44 | 17.17 | 2.42 |
| 6.4 | cp_1bit | 20.48 | 3.12 | -18.2 | 0.44 | 15.57 | 2.92 |
| 6.4 | cp_2bit | 23.40 | 0.20 | -22.8 | 0.61 | 17.22 | 2.41 |
| 6.4 | lp_cont | 18.60 | 0.00 | -2.6 | 5.06 | - | - |
| 6.4 | lp_1bit | 18.45 | 0.16 | -2.6 | 5.06 | - | - |
| 6.4 | lp_2bit | 18.48 | 0.13 | -2.5 | 5.06 | - | - |
| 6.4 | lpx_cont | 23.70 | 0.00 | -22.1 | 0.44 | - | - |
| 6.4 | lpx_1bit | 23.51 | 0.19 | -20.7 | 0.44 | - | - |
| 6.4 | lpx_2bit | 23.63 | 0.07 | -22.0 | 0.44 | - | - |
| 8.0 | cp_cont | 23.58 | 0.00 | -22.4 | 0.44 | 16.92 | 2.49 |
| 8.0 | cp_1bit | 20.41 | 3.17 | -20.7 | 0.44 | 15.01 | 3.12 |
| 8.0 | cp_2bit | 23.29 | 0.29 | -21.4 | 0.44 | 16.95 | 2.49 |
| 8.0 | lp_cont | 18.53 | 0.00 | -1.4 | 5.21 | - | - |
| 8.0 | lp_1bit | 18.21 | 0.32 | -2.6 | 5.17 | - | - |
| 8.0 | lp_2bit | 18.48 | 0.05 | -2.5 | 5.06 | - | - |
| 8.0 | lpx_cont | 23.68 | 0.00 | -21.9 | 0.44 | - | - |
| 8.0 | lpx_1bit | 23.37 | 0.31 | -21.2 | 0.44 | - | - |
| 8.0 | lpx_2bit | 23.62 | 0.06 | -21.6 | 0.44 | - | - |

Peak search: 0.25 degree direction-cosine refinement around the coarse-grid peak; the sidelobe level uses the strongest coarse-grid sidelobe refined the same way. Gain is the aperture directivity plus ten times the base-ten logarithm of (launch fraction x eta_conv); eta_conv of this design is given in the scan table.

## Finding: a scalar surface on a monopole cannot form a linear-polarised pencil beam

With the induced current along the surface wave's own (radial) direction the continuous scalar design has a pointing error of 5.06 degrees and a strongest sidelobe at -2.6 dB: the aperture field is radially polarised, so the contributions cancel in the beam direction (a ring or split beam, not a pencil). Gain 18.60 dBi against 23.70 dBi for the same hologram with ideal polarisation control (`lpx`) and 23.60 dBi for the CP tensor design. The `lp` rows above are therefore NOT a candidate; the bit-depth conclusion for LP rests on `lpx`, which needs a cell that sets the current direction (not designed).

## Steering: pitch 6.4 mm, 150 mm hexagon, CP tensor

Gain by beam polar angle (azimuth 225 degrees; broadside has no azimuth). Scan loss is relative to broadside of the same variant; a negative value means the gain is higher than at broadside (the conversion efficiency rises with the steering angle in this model). `table outside` is the share of the design's march points whose modulation wavenumber lies outside the supercell table: results with a large share are `TABLE_EXTRAPOLATED`.

| theta0 [deg] | eta_conv | table outside | cp cont [dBi] | cp 2-bit [dBi] | cp 1-bit [dBi] | 2-bit scan loss [dB] | 2-bit SLL [dB] | 2-bit pointing [deg] | 2-bit gain over phi0 = 0, 45, 225 [dBi] |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 0.453 | 0.00 | 23.60 | 23.34 | 20.53 | 0.00 | -22.8 | 0.00 | 23.34 to 23.34 |
| 10 | 0.471 | 0.00 | 23.66 | 23.50 | 20.59 | -0.16 | -22.4 | 0.18 | 23.50 to 23.52 |
| 20 | 0.519 | 0.00 | 23.86 | 23.72 | 20.79 | -0.38 | -22.7 | 0.25 | 23.72 to 23.74 |
| 30 | 0.538 | 0.22 | 23.60 | 23.40 | 20.48 | -0.06 | -22.8 | 0.61 | 23.40 to 23.40 |
| 40 | 0.548 | 0.48 | 22.99 | 22.75 | 19.84 | 0.59 | -21.6 | 0.79 | 22.62 to 22.75 |
| 50 | 0.554 | 0.57 | 21.99 | 21.62 | 18.69 | 1.72 | -19.1 | 1.01 | 21.44 to 21.62 |

## Pitch and aperture

| pitch [mm] | cells | cp continuous [dBi] | cp 2-bit loss [dB] | cp 1-bit loss [dB] | lpx 1-bit loss [dB] | lpx 2-bit loss [dB] |
|---|---|---|---|---|---|---|
| 4.0 | 3646 | 23.57 | 0.06 | 3.00 | 0.07 | 0.01 |
| 5.0 | 2306 | 23.52 | 0.16 | 3.07 | 0.16 | 0.02 |
| 6.4 | 1430 | 23.60 | 0.20 | 3.12 | 0.19 | 0.07 |
| 8.0 | 924 | 23.58 | 0.29 | 3.17 | 0.31 | 0.06 |

The larger reference aperture (hexagon of the larger circumradius), pitch 6.4 mm:

| circumradius [mm] | theta0 [deg] | cells | eta_conv | table outside | cp cont [dBi] | cp 2-bit [dBi] | cp 1-bit [dBi] | cp 2-bit SLL [dB] |
|---|---|---|---|---|---|---|---|---|
| 225 | 0 | 3210 | 0.574 | 0.00 | 28.17 | 27.97 | 24.96 | -25.8 |
| 225 | 30 | 3210 | 0.655 | 0.22 | 27.93 | 27.73 | 24.76 | -24.1 |
| 225 | 50 | 3210 | 0.669 | 0.57 | 26.24 | 25.88 | 22.87 | -23.2 |

## Control hardware (150 mm hexagon, beam theta 30 / phi 225)

Counting rules (one method for every row): diodes per cell = bits (one diode per bit, an ASSUMPTION; the topology is item 3); bit = 1 is the forward-biased diode at 10 mA, the share of forward bits is counted from the actual coded states; DC power = (forward diodes) x I_F x V_rail, V_rail = 3.3 V (the supply power, series resistor included; the diode alone dissipates 0.89 V x I_F). The same formula is evaluated at 1 mA. Independent bias lines = diodes (arbitrary patterns need an independent state per diode); sharing a line among several diodes constrains the pattern (the line-count columns below give the count only; four cells per line is the hologram sampled at twice the pitch, see the pitch table). Shift registers: eight- or sixteen-output devices in chains of eight clocked at 10 MHz; FPGA pins = chains plus the shared clock, latch and enable lines; output current per pin is the diode forward current.

| pitch [mm] | variant | cells | diodes | bias lines | lines, 2 per line | lines, 4 per line | eight-output chips | sixteen-output chips | FPGA pins (eight-output) | update [us] | forward share | DC current [A] | DC power expected, 10 mA [W] | DC power expected, 1 mA [W] | DC power all on, 10 mA [W] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 4.0 | cp_1bit | 3646 | 3646 | 3646 | 1823 | 912 | 456 | 228 | 60 | 6.4 | 0.493 | 18.0 | 59.4 | 5.94 | 120.3 |
| 4.0 | cp_2bit | 3646 | 7292 | 7292 | 3646 | 1823 | 912 | 456 | 117 | 6.4 | 0.495 | 36.1 | 119.2 | 11.92 | 240.6 |
| 5.0 | cp_1bit | 2306 | 2306 | 2306 | 1153 | 577 | 289 | 145 | 40 | 6.4 | 0.492 | 11.3 | 37.5 | 3.75 | 76.1 |
| 5.0 | cp_2bit | 2306 | 4612 | 4612 | 2306 | 1153 | 577 | 289 | 76 | 6.4 | 0.495 | 22.8 | 75.3 | 7.53 | 152.2 |
| 6.4 | cp_1bit | 1430 | 1430 | 1430 | 715 | 358 | 179 | 90 | 26 | 6.4 | 0.498 | 7.1 | 23.5 | 2.35 | 47.2 |
| 6.4 | cp_2bit | 1430 | 2860 | 2860 | 1430 | 715 | 358 | 179 | 48 | 6.4 | 0.493 | 14.1 | 46.6 | 4.66 | 94.4 |
| 8.0 | cp_1bit | 924 | 924 | 924 | 462 | 231 | 116 | 58 | 18 | 6.4 | 0.497 | 4.6 | 15.1 | 1.51 | 30.5 |
| 8.0 | cp_2bit | 924 | 1848 | 1848 | 924 | 462 | 231 | 116 | 32 | 6.4 | 0.490 | 9.1 | 29.9 | 2.99 | 61.0 |

Driver current per output: 10 mA (a sixteen-output device with every output on draws 160 mA, an eight-output device 80 mA): the output stage must source or sink this per pin; whether common logic-level shift registers are specified for it is NOT checked here (datasheets to be checked in item 3), a current-driver class may be needed. The LP rows are identical in count to the CP rows of the same bit depth (`lp` and `lpx` bias tables are in the JSON).

The lower-bias column is only a power figure; the diode series resistance rises at low bias (item 3 must trade it against insertion loss).

## Recommendation and trade-offs

- **Polarisation scope: CP, the rotating-slot tensor family.** It is the family the existing forward model covers, with a single HFSS check of the continuous hexagon (one board, one run); the scalar LP surface does not form a pencil beam (finding above), and the `lpx` bound needs a cell that controls the current direction, which no design exists for.
- **Bit depth: 2-bit (four slot axes, equally spaced).** At pitch 6.4 mm the 2-bit loss is 0.20 dB against 3.12 dB for 1-bit; 1-bit costs one diode per cell fewer but loses that much gain and its sidelobe level degrades to -18.2 dB (2-bit: -22.8 dB).
- **Pitch: 6.4 mm** (the pitch of the supercell table and of the HFSS-checked hexagon). In this idealised model the pitch hardly matters (2-bit loss from 0.06 dB at 4.0 mm to 0.29 dB at 8.0 mm), so the choice is set by diode count and by what item 3 can build; 8.0 mm is the alternative worth testing with a real cell (1848 diodes against 2860).
- **Diode count at the recommendation: 2860 diodes (1430 cells x 2), 2860 independent bias lines**, 179 sixteen-output shift registers, 26 FPGA pins at the stated chain length; expected DC power 46.6 W at 10 mA per forward diode (4.66 W at 1 mA). The bias power, not the diode count, is the sharpest cost of this architecture.
- **Required contrast for the recommendation (beam 30 degrees):** tensor anisotropy 2.296 mS, i.e. a principal-axis reactance difference of 92.3 ohm, present in every state, with the axis rotating between four orientations; the equivalent scalar swing would be 148.4 ohm. These are sheet-level numbers; whether a diode-loaded subwavelength cell reaches them is the question of item 3.
- **Scan range:** the model is supported by the supercell table up to about the 30 degree beam (the share of table-extrapolated march points is in the scan table); beyond it the share of table-extrapolated march points rises (scan table) and the scan losses are `TABLE_EXTRAPOLATED`.

## What is not available or not established

- 10.0 GHz: `NOT_AVAILABLE_FROM_CURRENT_DATA (the supercell table exists at 10.54 GHz only)`; a table at that frequency is a server job, not run.
- Everything is `MODEL_ONLY`: idealised cells, the table assumed at all pitches, no coupling between neighbouring states, no diode loss, the unconverted surface wave and the feed's direct radiation assumed outside the beam (Phase 26 assumptions).
- The literature table (item 1) does not support any of the contrast, bit-depth or bias statements here; they are model results.

## Stop point

No cell has been designed. Approve or change the pitch, bit depth, polarisation scope and diode count above before item 3 starts.
