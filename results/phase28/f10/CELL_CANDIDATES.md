# Phase 28 item 3 - diode-switched CP cell candidates at 10.0 GHz (stop point)

**Labels:** `MODEL_ONLY`, `NOT_HFSS_VALIDATED`, `PROVISIONAL_DATASHEET`, `VIA_OMITTED_APPROXIMATION`. Generated from the characterisation JSON files only.

Method: the normal-incidence tensor sheet of each state of the cell (`coded_cell_eval`: the Phase 26 B extraction `Y_s` from the 2 x 2 reflection matrix, the lumped diode loads of the Phase 20 R1 / R2 solver extension acting inside the MoM), its principal susceptances (mean loading and anisotropy), the TM0 phase constant of the mean sheet and the reactance difference of the two principal axes (the contrast measure of the item-2 report). This is a screening of the sheet model, not a dispersion (alpha, beta) of a row of cells; that is the state library of item 4.

## What was tried

| family | pitch [mm] | samples | valid | best smallest-state dX_aniso own [ohm] | best smallest-state alpha at the beam [Np/m] |
|---|---|---|---|---|---|
| octagon_stub | 6.4 | 600 | 88 | 2.56 | 0.000 |
| octagon_stub | 8.0 | 600 | 213 | 2.57 | 0.000 |
| wedge_screen | 6.4 | 72 | 72 | 22.37 | 0.077 |
| wedge_screen | 8.0 | 72 | 72 | 128.54 | 3.564 |
| wedge_balanced | 6.4 | 400 | 400 | 134.49 | 1.208 |
| wedge_balanced | 8.0 | 400 | 400 | 1360.67 | 17.868 |
| wedge_refined | 6.4 | 1200 | 1200 | 153.27 | 0.854 |
| wedge_notch | 6.4 | 400 | 400 | 555.92 | 1.683 |
| wedge_notch | 8.0 | 400 | 400 | 3154.15 | 3.936 |
| manhattan_bowtie | 6.4 | 30 | 30 | 10.84 | 0.018 |
| manhattan_notch | 6.4 | 30 | 30 | 14.51 | 0.061 |
| manhattan_bowtie | 8.0 | 110 | 110 | 41.37 | 0.387 |
| manhattan_notch | 8.0 | 110 | 110 | 272.99 | 0.983 |
| rational_bowtie | 6.4 | 94 | 94 | 47.57 | 0.496 |
| rational_notch | 6.4 | 94 | 94 | 134.04 | 1.277 |
| rational_notch_random | 6.4 | 1200 | 601 | 125.45 | 1.096 |
| rational_bowtie | 8.0 | 12 | 12 | 59.07 | 0.808 |
| rational_notch | 8.0 | 168 | 168 | 11225.81 | 2.363 |
| rational_notch_random | 8.0 | 1200 | 898 | 2306.52 | 1.969 |

- octagon_stub: the repo's CP-2 octagon cell with four diode stubs (ris.cells.Cp2Cell) put into the 6.4 and 8.0 mm cells: the stubs cannot be long enough inside the cell and the reactance difference stays at a few ohm. Not viable.
- slot star (scripts/research/phase28_star_negative.py, results/phase28/f10/search/star_cell_negative.json.gz): four through-slots with a diode wall across each arm. With ideal switches the cell with the walls closed has less than one percent of the anisotropy of the same cell with the perpendicular slot really shortened: a diode wall across a slot does not close the slot. Not viable.
- wedge (SectorHubCell, a hub and eight wedge sectors with radial edges, eight diode necks): the strongest anisotropy and the highest leakage, but NOT mesh-convergent: the wedge edges have irrational slopes and the rooftop mesh rasterises them differently at every mesh, so the leakage at the beam moves by tens of percent between consecutive meshes (candidates A, B, C in the raw files). The screening numbers of the wedge families are not to be used.
- manhattan (ManhattanHubCell, every edge axis-aligned): converges, but the pieces leave too much of the cell empty; the leakage at the beam stays far below the requirement at both pitches.
- rational (RationalSectorCell, wedges whose edges have rational slopes 0, 1/2, 1, 2 so that the mesh represents them by a self-similar staircase): converges and keeps the metal fill; the candidates D to J come from this family. The search score is the smaller leakage of the first two states at the coarse screening mesh; `n_balanced` counts samples whose two states differ by less than fifteen percent in loading and anisotropy.

## Summary of the candidates

| candidate | pitch [mm] | mode | diodes per cell | forward diodes per state | converged | smallest dX_aniso own [ohm] | smallest dX_aniso design [ohm] | meets the target in either definition | alpha at the beam over the states [Np/m] | delta over the states [mS] | worst diode-tolerance change in alpha [%] | DC power at the forward current [W] | DC power at 1 mA [W] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| B_wedge_bowtie_6p4mm | 6.4 | bowtie | 8 | 2 | NO | 43.6 | 91.7 | NO | 0.37 to 0.75 | 2.36 to 2.88 | 69.1 | 94 | 9.4 |
| C_wedge_notch_6p4mm | 6.4 | notch | 8 | 6 | NO | 73.3 | 69.7 | NO | 1.34 to 2.16 | 1.87 to 2.40 | 17.7 | 283 | 28.3 |
| D_rational_notch_8mm | 8.0 | notch | 8 | 6 | NO | 80.2 | 74.1 | NO | 1.48 to 1.93 | 1.97 to 2.25 | 4.8 | 183 | 18.3 |
| E_rational_bowtie_8mm | 8.0 | bowtie | 8 | 2 | NO | 41.6 | 85.4 | NO | 0.35 to 0.40 | 2.22 to 2.32 | 4.8 | 61 | 6.1 |
| F_rational_notch_6p4mm | 6.4 | notch | 8 | 6 | yes | 76.3 | 67.6 | NO | 1.23 to 1.59 | 1.81 to 2.04 | 4.7 | 283 | 28.3 |
| G_rational_bowtie_6p4mm | 6.4 | bowtie | 8 | 2 | yes | 43.9 | 88.6 | NO | 0.40 to 0.51 | 2.29 to 2.47 | 4.6 | 94 | 9.4 |
| I_rational_notch_tuned_8mm | 8.0 | notch | 8 | 6 | NO | 321.8 | 173.0 | yes | 0.12 to 3.72 | 3.73 to 8.49 | 171.4 | 183 | 18.3 |
| J_rational_notch_wide_8mm | 8.0 | notch | 8 | 6 | yes | 231.5 | 120.0 | yes | 1.75 to 1.75 | 2.91 to 3.22 | 4.7 | 183 | 18.3 |

`converged` is the mesh criterion below; a candidate that is not converged is listed for the record and its numbers are not to be used. `mode`: `bowtie` = the two necks of the state's axis conduct, `notch` = every neck conducts except the two of the perpendicular axis. The 103 ohm target and 132 ohm stretch goal are the values approved from the 10.54 GHz item-2 sweep (combined 1.5 dB and 1 dB totals); at 10.0 GHz the same item-2 totals need lower anisotropies, and the leakage the surface must reach is `alpha_cap`; the 10.0 GHz values are in `FEASIBILITY_REPORT.md` and `FREQUENCY_COMPARISON.md`. `meets the target` refers to the anisotropy only; compare the leakage column with `alpha_cap` separately.

## Candidate B_wedge_bowtie_6p4mm

Hub-and-sector cell (`coded_cell_candidates`, the family is in the table): a hub and eight sectors, each joined to the hub by a neck that carries one diode. Mode `bowtie`: state k = the two necks on axis k conduct (the hub and the two opposite sectors form a bow-tie along that axis), the other six are open. Mode `notch`: every neck conducts except the two of the perpendicular axis (the metal fill stays high, two opposite sectors float). The four states are the same cell turned in 45 degree steps, so the sheet tensor rotates with the state.

| family | period [mm] | diode kit | hub_radius [mm] | radial_gap [mm] | neck_width [mm] | sector_gap [mm] | radius [mm] | radius_diag [mm] | edge_margin [mm] |
|---|---|---|---|---|---|---|---|---|---|
| sector_hub | 6.4 | madp | 0.80 | 0.80 | 0.40 | 0.20 | 3.06 | 3.37 | 0.10 |

### Mesh convergence (normal-incidence tensor sheet, 10.0 GHz)

| mesh | state | mean loading b [mS] | anisotropy delta [mS] | beta/k0 | dX_aniso own loading [ohm] | dX_aniso design loading [ohm] | alpha at the beam [Np/m] | largest conductance [mS] | principal axis [deg] |
|---|---|---|---|---|---|---|---|---|---|
| 32 | 0 | 4.833 | 2.938 | 1.0336 | 62.2 | 121.5 | 0.809 | 0.822 | 0.0 |
| 32 | 1 | 4.923 | 3.041 | 1.0341 | 65.9 | 127.4 | 0.915 | 0.850 | 45.0 |
| 32 | 2 | 4.833 | 2.938 | 1.0336 | 62.2 | 121.5 | 0.809 | 0.822 | 90.0 |
| 32 | 3 | 4.923 | 3.041 | 1.0341 | 65.9 | 127.4 | 0.915 | 0.850 | 135.0 |
| 48 | 0 | 5.407 | 3.398 | 1.0370 | 82.8 | 149.5 | 1.519 | 1.606 | 180.0 |
| 48 | 1 | 4.968 | 2.958 | 1.0344 | 64.3 | 122.6 | 0.888 | 0.710 | 45.0 |
| 48 | 2 | 5.408 | 3.397 | 1.0370 | 82.8 | 149.4 | 1.518 | 1.606 | 90.0 |
| 48 | 3 | 4.967 | 2.958 | 1.0344 | 64.3 | 122.6 | 0.888 | 0.710 | 135.0 |
| 64 | 0 | 4.232 | 2.360 | 1.0305 | 43.6 | 91.7 | 0.369 | 0.598 | 180.0 |
| 64 | 1 | 4.761 | 2.879 | 1.0332 | 60.0 | 118.2 | 0.745 | 0.938 | 45.0 |
| 64 | 2 | 4.232 | 2.360 | 1.0305 | 43.6 | 91.7 | 0.369 | 0.598 | 90.0 |
| 64 | 3 | 4.761 | 2.879 | 1.0332 | 60.0 | 118.2 | 0.745 | 0.938 | 135.0 |

Convergence is judged on the relative change between mesh 48 and mesh 64, largest over the four states, against the tolerance of the table (an ASSUMED criterion: 0.05 on the anisotropy, the loading and the reactance difference; the leakage at the beam grows as the square of the anisotropy, so its tolerance is twice that):

| metric | largest relative change | tolerance | converged |
|---|---|---|---|
| b_bar_ms | 0.2779 | 0.05 | NO |
| delta_ms | 0.4401 | 0.05 | NO |
| dx_aniso_ohm | 0.8988 | 0.05 | NO |
| dx_aniso_at_design_loading_ohm | 0.6295 | 0.05 | NO |
| alpha_beam_np_per_m | 3.1188 | 0.10 | NO |

All metrics converged: NO - the numbers of this candidate are not to be used until they are.

Floquet cutoff: the same four states at mesh 48 with the cutoff 2 times the mesh number (the converged setting of the Phase 20 R4 audit for lumped loads) against the cutoff equal to the mesh number: largest relative change b_bar_ms 0.0021, delta_ms 0.0052, dx_aniso_ohm 0.0081, dx_aniso_at_design_loading_ohm 0.0071, alpha_beam_np_per_m 0.0166; cutoff converged: yes.

### Targets and state balance (mesh 64)

| dX_aniso definition | smallest over the four states [ohm] | meets 103 ohm | meets 132 ohm (stretch) |
|---|---|---|---|
| own mean loading (the cell's beta) | 43.6 | NO | NO |
| design mean loading (beta/k0 of the item-2 contrast derivation) | 91.7 | NO | NO |

The same anisotropy gives a different reactance difference on a different mean loading (the surface wave is more or less bound), so both definitions are shown. The leakage the cell delivers to a beam at the design polar angle (`alpha at the beam`) is the quantity the hologram needs: compare it with the contrast requirement of the item-2 feasibility sweep at this frequency.

| quantity over the four states | smallest | largest | relative spread |
|---|---|---|---|
| b_bar_ms | 4.232 | 4.761 | 0.111 |
| delta_ms | 2.360 | 2.879 | 0.180 |
| dx_aniso_ohm | 43.626 | 59.955 | 0.272 |
| dx_aniso_at_design_loading_ohm | 91.746 | 118.164 | 0.224 |
| alpha_beam_np_per_m | 0.369 | 0.745 | 0.505 |
| beta_over_k0 | 1.031 | 1.033 | 0.003 |

### Diode tolerance (plus and minus twenty percent of the diode model parameters, states 0 and 1, mesh 64)

ASSUMED_SENSITIVITY_NOT_DEVICE_DATA: the scale factors multiply the series resistance, the package inductance and the junction capacitance of the kit; the MADP package inductance is NOT_AVAILABLE_FROM_CURRENT_DATA and is zero in the nominal model, so its scale has no effect there.

| scale | state | delta [mS] | change in delta [%] | dX_aniso own [ohm] | change [%] | alpha at the beam [Np/m] | change [%] |
|---|---|---|---|---|---|---|---|
| {"C_T": 0.8} | 0 | 1.657 | -29.8 | 26.1 | -40.3 | 0.114 | -69.0 |
| {"C_T": 0.8} | 1 | 2.839 | -1.4 | 57.2 | -4.6 | 0.654 | -12.2 |
| {"C_T": 1.2} | 0 | 2.601 | 10.2 | 52.0 | 19.1 | 0.558 | 51.4 |
| {"C_T": 1.2} | 1 | 2.955 | 2.7 | 64.2 | 7.1 | 0.888 | 19.2 |
| {"L_s": 0.8} | 0 | 2.360 | 0.0 | 43.6 | 0.0 | 0.369 | 0.0 |
| {"L_s": 0.8} | 1 | 2.879 | 0.0 | 60.0 | 0.0 | 0.745 | 0.0 |
| {"L_s": 1.2} | 0 | 2.360 | 0.0 | 43.6 | 0.0 | 0.369 | 0.0 |
| {"L_s": 1.2} | 1 | 2.879 | 0.0 | 60.0 | 0.0 | 0.745 | 0.0 |
| {"R_s": 0.8} | 0 | 2.366 | 0.3 | 43.8 | 0.4 | 0.372 | 0.9 |
| {"R_s": 0.8} | 1 | 2.899 | 0.7 | 60.7 | 1.2 | 0.765 | 2.7 |
| {"R_s": 1.2} | 0 | 2.350 | -0.4 | 43.4 | -0.6 | 0.364 | -1.3 |
| {"R_s": 1.2} | 1 | 2.852 | -0.9 | 59.0 | -1.5 | 0.720 | -3.3 |
| {"C_T": 0.8, "L_s": 1.2, "R_s": 1.2} | 0 | 1.656 | -29.8 | 26.0 | -40.3 | 0.114 | -69.1 |
| {"C_T": 0.8, "L_s": 1.2, "R_s": 1.2} | 1 | 2.818 | -2.1 | 56.5 | -5.8 | 0.637 | -14.5 |
| {"C_T": 1.2, "L_s": 0.8, "R_s": 0.8} | 0 | 2.611 | 10.6 | 52.3 | 19.8 | 0.565 | 53.3 |
| {"C_T": 1.2, "L_s": 0.8, "R_s": 0.8} | 1 | 2.978 | 3.5 | 65.1 | 8.6 | 0.914 | 22.7 |

Worst relative change over all tolerance cases: delta 29.8 %, dX_aniso (own loading) 40.3 %, dX_aniso (design loading) 33.5 %, alpha at the beam 69.1 %.

Package inductance of the MADP model (NOT_AVAILABLE_FROM_CURRENT_DATA, zero in the nominal model): ASSUMED values as a sensitivity, not device data.

| package inductance [nH] | state | delta [mS] | change [%] | dX_aniso own [ohm] | change [%] | alpha at the beam [Np/m] | change [%] |
|---|---|---|---|---|---|---|---|
| 0.3 | 0 | 4.184 | 77.3 | 127.0 | 191.2 | 3.388 | 818.8 |
| 0.3 | 1 | 6.049 | 110.1 | 517.5 | 763.1 | 13.859 | 1760.3 |
| 0.7 | 0 | 5.063 | 114.6 | 222.8 | 410.6 | 7.678 | 1982.0 |
| 0.7 | 1 | 9.277 | 222.3 | 48.6 | -18.9 | 0.068 | -90.9 |

### Contrast bandwidth (four states, mesh 64)

| frequency [GHz] | smallest dX_aniso own [ohm] | smallest dX_aniso design [ohm] | smallest alpha at the beam [Np/m] |
|---|---|---|---|
| 9.00 | 25.7 | 76.3 | 0.091 |
| 9.25 | 36.5 | 96.0 | 0.212 |
| 9.50 | 54.2 | 124.7 | 0.540 |
| 9.75 | 25.8 | 63.4 | 0.108 |
| 10.00 | 43.6 | 91.7 | 0.369 |
| 10.25 | 58.5 | 108.9 | 0.740 |
| 10.50 | 76.9 | 126.3 | 1.387 |
| 10.75 | 102.6 | 147.0 | 2.562 |
| 11.00 | 142.0 | 173.8 | 4.646 |

- smallest design-loading dX_aniso at or above the target: threshold 103.0 ohm is not met at the design frequency.
- smallest own-loading dX_aniso within 70 percent of its 10.0 GHz value (threshold 30.5 ohm): from 10.00 to 11.00 GHz, width 1.00 GHz (grid step 0.25 GHz).

### Control hardware (150 mm hexagon, same counting rules as the item-2 report)

| diodes per cell | forward diodes per state | cells | diodes | forward diodes | forward current [mA] | DC power at that current [W] | DC power at 1 mA [W] |
|---|---|---|---|---|---|---|---|
| 8 | 2 | 1430 | 11440 | 2860 | 10 | 94.4 | 9.44 |

DC power = forward diodes x forward current x supply rail (3.3 V). Not modelled: the bias network (a connection to the hub and to every sector pad through chokes) and the diode package footprint; both are `NOT_AVAILABLE_FROM_CURRENT_DATA` here.

## Candidate C_wedge_notch_6p4mm

Hub-and-sector cell (`coded_cell_candidates`, the family is in the table): a hub and eight sectors, each joined to the hub by a neck that carries one diode. Mode `bowtie`: state k = the two necks on axis k conduct (the hub and the two opposite sectors form a bow-tie along that axis), the other six are open. Mode `notch`: every neck conducts except the two of the perpendicular axis (the metal fill stays high, two opposite sectors float). The four states are the same cell turned in 45 degree steps, so the sheet tensor rotates with the state.

| family | period [mm] | diode kit | hub_radius [mm] | radial_gap [mm] | neck_width [mm] | sector_gap [mm] | radius [mm] | radius_diag [mm] | edge_margin [mm] |
|---|---|---|---|---|---|---|---|---|---|
| sector_hub | 6.4 | madp | 0.40 | 0.80 | 0.40 | 0.40 | 2.98 | 3.30 | 0.10 |

### Mesh convergence (normal-incidence tensor sheet, 10.0 GHz)

| mesh | state | mean loading b [mS] | anisotropy delta [mS] | beta/k0 | dX_aniso own loading [ohm] | dX_aniso design loading [ohm] | alpha at the beam [Np/m] | largest conductance [mS] | principal axis [deg] |
|---|---|---|---|---|---|---|---|---|---|
| 32 | 0 | 8.543 | 1.875 | 1.0632 | 74.3 | 70.1 | 1.353 | 1.088 | 180.0 |
| 32 | 1 | 7.676 | 2.158 | 1.0545 | 74.0 | 82.4 | 1.696 | 1.643 | 45.0 |
| 32 | 2 | 8.543 | 1.875 | 1.0632 | 74.3 | 70.1 | 1.353 | 1.088 | 90.0 |
| 32 | 3 | 7.676 | 2.158 | 1.0545 | 74.0 | 82.4 | 1.696 | 1.643 | 135.0 |
| 48 | 0 | 7.998 | 1.747 | 1.0576 | 61.9 | 64.8 | 1.162 | 1.134 | 180.0 |
| 48 | 1 | 7.398 | 1.803 | 1.0519 | 57.3 | 67.1 | 1.113 | 1.009 | 45.0 |
| 48 | 2 | 8.006 | 1.743 | 1.0577 | 61.9 | 64.6 | 1.157 | 1.136 | 90.0 |
| 48 | 3 | 7.397 | 1.803 | 1.0519 | 57.3 | 67.1 | 1.112 | 1.008 | 135.0 |
| 64 | 0 | 8.504 | 1.865 | 1.0628 | 73.3 | 69.7 | 1.341 | 1.229 | 0.0 |
| 64 | 1 | 7.836 | 2.402 | 1.0560 | 86.7 | 93.8 | 2.156 | 2.642 | 45.0 |
| 64 | 2 | 8.504 | 1.865 | 1.0628 | 73.3 | 69.7 | 1.341 | 1.229 | 90.0 |
| 64 | 3 | 7.836 | 2.402 | 1.0560 | 86.7 | 93.8 | 2.156 | 2.642 | 135.0 |

Convergence is judged on the relative change between mesh 48 and mesh 64, largest over the four states, against the tolerance of the table (an ASSUMED criterion: 0.05 on the anisotropy, the loading and the reactance difference; the leakage at the beam grows as the square of the anisotropy, so its tolerance is twice that):

| metric | largest relative change | tolerance | converged |
|---|---|---|---|
| b_bar_ms | 0.0595 | 0.05 | NO |
| delta_ms | 0.2492 | 0.05 | NO |
| dx_aniso_ohm | 0.3383 | 0.05 | NO |
| dx_aniso_at_design_loading_ohm | 0.2843 | 0.05 | NO |
| alpha_beam_np_per_m | 0.4839 | 0.10 | NO |

All metrics converged: NO - the numbers of this candidate are not to be used until they are.

Floquet cutoff: the same four states at mesh 48 with the cutoff 2 times the mesh number (the converged setting of the Phase 20 R4 audit for lumped loads) against the cutoff equal to the mesh number: largest relative change b_bar_ms 0.0025, delta_ms 0.0100, dx_aniso_ohm 0.0117, dx_aniso_at_design_loading_ohm 0.0112, alpha_beam_np_per_m 0.0209; cutoff converged: yes.

### Targets and state balance (mesh 64)

| dX_aniso definition | smallest over the four states [ohm] | meets 103 ohm | meets 132 ohm (stretch) |
|---|---|---|---|
| own mean loading (the cell's beta) | 73.3 | NO | NO |
| design mean loading (beta/k0 of the item-2 contrast derivation) | 69.7 | NO | NO |

The same anisotropy gives a different reactance difference on a different mean loading (the surface wave is more or less bound), so both definitions are shown. The leakage the cell delivers to a beam at the design polar angle (`alpha at the beam`) is the quantity the hologram needs: compare it with the contrast requirement of the item-2 feasibility sweep at this frequency.

| quantity over the four states | smallest | largest | relative spread |
|---|---|---|---|
| b_bar_ms | 7.836 | 8.504 | 0.079 |
| delta_ms | 1.865 | 2.402 | 0.223 |
| dx_aniso_ohm | 73.276 | 86.658 | 0.154 |
| dx_aniso_at_design_loading_ohm | 69.694 | 93.754 | 0.257 |
| alpha_beam_np_per_m | 1.341 | 2.156 | 0.378 |
| beta_over_k0 | 1.056 | 1.063 | 0.006 |

### Diode tolerance (plus and minus twenty percent of the diode model parameters, states 0 and 1, mesh 64)

ASSUMED_SENSITIVITY_NOT_DEVICE_DATA: the scale factors multiply the series resistance, the package inductance and the junction capacitance of the kit; the MADP package inductance is NOT_AVAILABLE_FROM_CURRENT_DATA and is zero in the nominal model, so its scale has no effect there.

| scale | state | delta [mS] | change in delta [%] | dX_aniso own [ohm] | change [%] | alpha at the beam [Np/m] | change [%] |
|---|---|---|---|---|---|---|---|
| {"C_T": 0.8} | 0 | 1.892 | 1.4 | 74.1 | 1.1 | 1.380 | 2.9 |
| {"C_T": 0.8} | 1 | 2.581 | 7.4 | 98.0 | 13.1 | 2.538 | 17.7 |
| {"C_T": 1.2} | 0 | 1.840 | -1.4 | 72.5 | -1.0 | 1.303 | -2.8 |
| {"C_T": 1.2} | 1 | 2.278 | -5.1 | 79.9 | -7.8 | 1.912 | -11.3 |
| {"L_s": 0.8} | 0 | 1.865 | 0.0 | 73.3 | 0.0 | 1.341 | 0.0 |
| {"L_s": 0.8} | 1 | 2.402 | 0.0 | 86.7 | 0.0 | 2.156 | 0.0 |
| {"L_s": 1.2} | 0 | 1.865 | 0.0 | 73.3 | 0.0 | 1.341 | 0.0 |
| {"L_s": 1.2} | 1 | 2.402 | 0.0 | 86.7 | 0.0 | 2.156 | 0.0 |
| {"R_s": 0.8} | 0 | 1.869 | 0.2 | 73.7 | 0.5 | 1.346 | 0.4 |
| {"R_s": 0.8} | 1 | 2.483 | 3.4 | 92.6 | 6.9 | 2.340 | 8.5 |
| {"R_s": 1.2} | 0 | 1.860 | -0.3 | 72.7 | -0.8 | 1.334 | -0.5 |
| {"R_s": 1.2} | 1 | 2.352 | -2.1 | 83.1 | -4.2 | 2.040 | -5.4 |
| {"C_T": 0.8, "L_s": 1.2, "R_s": 1.2} | 0 | 1.886 | 1.1 | 73.5 | 0.3 | 1.373 | 2.4 |
| {"C_T": 0.8, "L_s": 1.2, "R_s": 1.2} | 1 | 2.494 | 3.8 | 91.3 | 5.4 | 2.334 | 8.3 |
| {"C_T": 1.2, "L_s": 0.8, "R_s": 0.8} | 0 | 1.844 | -1.1 | 72.9 | -0.5 | 1.309 | -2.4 |
| {"C_T": 1.2, "L_s": 0.8, "R_s": 0.8} | 1 | 2.308 | -3.9 | 82.0 | -5.4 | 1.978 | -8.2 |

Worst relative change over all tolerance cases: delta 7.4 %, dX_aniso (own loading) 13.1 %, dX_aniso (design loading) 9.3 %, alpha at the beam 17.7 %.

Package inductance of the MADP model (NOT_AVAILABLE_FROM_CURRENT_DATA, zero in the nominal model): ASSUMED values as a sensitivity, not device data.

| package inductance [nH] | state | delta [mS] | change [%] | dX_aniso own [ohm] | change [%] | alpha at the beam [Np/m] | change [%] |
|---|---|---|---|---|---|---|---|
| 0.3 | 0 | 2.685 | 44.0 | 327.9 | 347.5 | 0.661 | -50.7 |
| 0.3 | 1 | 4.079 | 69.8 | 472.7 | 445.5 | 3.496 | 62.2 |
| 0.7 | 0 | 16.770 | 799.1 | 1508.4 | 1958.6 | 0.783 | -41.6 |
| 0.7 | 1 | 5.182 | 115.8 | 3058.8 | 3429.7 | 4.450 | 106.5 |

### Contrast bandwidth (four states, mesh 64)

| frequency [GHz] | smallest dX_aniso own [ohm] | smallest dX_aniso design [ohm] | smallest alpha at the beam [Np/m] |
|---|---|---|---|
| 9.00 | 24.9 | 46.3 | 0.176 |
| 9.25 | 31.7 | 50.8 | 0.333 |
| 9.50 | 41.1 | 56.1 | 0.602 |
| 9.75 | 54.3 | 62.3 | 0.972 |
| 10.00 | 73.3 | 69.7 | 1.341 |
| 10.25 | 101.6 | 78.8 | 1.600 |
| 10.50 | 145.7 | 90.2 | 1.706 |
| 10.75 | 207.8 | 104.9 | 1.680 |
| 11.00 | 292.1 | 124.7 | 1.570 |

- smallest design-loading dX_aniso at or above the target: threshold 103.0 ohm is not met at the design frequency.
- smallest own-loading dX_aniso within 70 percent of its 10.0 GHz value (threshold 51.3 ohm): from 9.75 to 11.00 GHz, width 1.25 GHz (grid step 0.25 GHz).

### Control hardware (150 mm hexagon, same counting rules as the item-2 report)

| diodes per cell | forward diodes per state | cells | diodes | forward diodes | forward current [mA] | DC power at that current [W] | DC power at 1 mA [W] |
|---|---|---|---|---|---|---|---|
| 8 | 6 | 1430 | 11440 | 8580 | 10 | 283.1 | 28.31 |

DC power = forward diodes x forward current x supply rail (3.3 V). Not modelled: the bias network (a connection to the hub and to every sector pad through chokes) and the diode package footprint; both are `NOT_AVAILABLE_FROM_CURRENT_DATA` here.

## Candidate D_rational_notch_8mm

Hub-and-sector cell (`coded_cell_candidates`, the family is in the table): a hub and eight sectors, each joined to the hub by a neck that carries one diode. Mode `bowtie`: state k = the two necks on axis k conduct (the hub and the two opposite sectors form a bow-tie along that axis), the other six are open. Mode `notch`: every neck conducts except the two of the perpendicular axis (the metal fill stays high, two opposite sectors float). The four states are the same cell turned in 45 degree steps, so the sheet tensor rotates with the state.

| family | period [mm] | diode kit | hub_half [mm] | neck_gap [mm] | neck_width [mm] | gap_t [mm] | arm_outer [mm] | diag_outer [mm] | clip [mm] |
|---|---|---|---|---|---|---|---|---|---|
| rational_sector_hub | 8.0 | madp | 1.60 | 0.40 | 0.80 | 0.40 | 3.60 | 3.20 | 3.80 |

### Mesh convergence (normal-incidence tensor sheet, 10.0 GHz)

| mesh | state | mean loading b [mS] | anisotropy delta [mS] | beta/k0 | dX_aniso own loading [ohm] | dX_aniso design loading [ohm] | alpha at the beam [Np/m] | largest conductance [mS] | principal axis [deg] |
|---|---|---|---|---|---|---|---|---|---|
| 40 | 0 | 8.659 | 2.261 | 1.0645 | 94.5 | 87.1 | 1.954 | 0.898 | 2.0 |
| 40 | 1 | 8.656 | 1.968 | 1.0644 | 80.2 | 74.1 | 1.482 | 1.037 | 42.7 |
| 40 | 2 | 8.659 | 2.261 | 1.0645 | 94.5 | 87.1 | 1.954 | 0.898 | 92.0 |
| 40 | 3 | 8.656 | 1.968 | 1.0644 | 80.2 | 74.1 | 1.482 | 1.037 | 132.7 |
| 60 | 0 | 8.705 | 2.192 | 1.0650 | 91.9 | 84.0 | 1.832 | 1.074 | 1.7 |
| 60 | 1 | 8.534 | 2.055 | 1.0631 | 82.4 | 77.9 | 1.625 | 1.206 | 43.1 |
| 60 | 2 | 8.705 | 2.192 | 1.0650 | 91.9 | 84.0 | 1.832 | 1.074 | 91.7 |
| 60 | 3 | 8.534 | 2.055 | 1.0631 | 82.4 | 77.9 | 1.625 | 1.206 | 133.1 |
| 80 | 0 | 8.671 | 2.249 | 1.0646 | 94.2 | 86.6 | 1.932 | 0.903 | 1.9 |
| 80 | 1 | 8.655 | 1.969 | 1.0644 | 80.2 | 74.1 | 1.483 | 1.047 | 42.8 |
| 80 | 2 | 8.671 | 2.249 | 1.0646 | 94.2 | 86.6 | 1.932 | 0.903 | 91.9 |
| 80 | 3 | 8.655 | 1.969 | 1.0644 | 80.2 | 74.1 | 1.483 | 1.047 | 132.8 |

Convergence is judged on the relative change between mesh 60 and mesh 80, largest over the four states, against the tolerance of the table (an ASSUMED criterion: 0.05 on the anisotropy, the loading and the reactance difference; the leakage at the beam grows as the square of the anisotropy, so its tolerance is twice that):

| metric | largest relative change | tolerance | converged |
|---|---|---|---|
| b_bar_ms | 0.0141 | 0.05 | yes |
| delta_ms | 0.0436 | 0.05 | yes |
| dx_aniso_ohm | 0.0269 | 0.05 | yes |
| dx_aniso_at_design_loading_ohm | 0.0504 | 0.05 | NO |
| alpha_beam_np_per_m | 0.0962 | 0.10 | yes |

All metrics converged: NO - the numbers of this candidate are not to be used until they are.

Floquet cutoff: the same four states at mesh 60 with the cutoff 2 times the mesh number (the converged setting of the Phase 20 R4 audit for lumped loads) against the cutoff equal to the mesh number: largest relative change b_bar_ms 0.0013, delta_ms 0.0016, dx_aniso_ohm 0.0032, dx_aniso_at_design_loading_ohm 0.0019, alpha_beam_np_per_m 0.0028; cutoff converged: yes.

### Targets and state balance (mesh 80)

| dX_aniso definition | smallest over the four states [ohm] | meets 103 ohm | meets 132 ohm (stretch) |
|---|---|---|---|
| own mean loading (the cell's beta) | 80.2 | NO | NO |
| design mean loading (beta/k0 of the item-2 contrast derivation) | 74.1 | NO | NO |

The same anisotropy gives a different reactance difference on a different mean loading (the surface wave is more or less bound), so both definitions are shown. The leakage the cell delivers to a beam at the design polar angle (`alpha at the beam`) is the quantity the hologram needs: compare it with the contrast requirement of the item-2 feasibility sweep at this frequency.

| quantity over the four states | smallest | largest | relative spread |
|---|---|---|---|
| b_bar_ms | 8.655 | 8.671 | 0.002 |
| delta_ms | 1.969 | 2.249 | 0.124 |
| dx_aniso_ohm | 80.217 | 94.157 | 0.148 |
| dx_aniso_at_design_loading_ohm | 74.113 | 86.564 | 0.144 |
| alpha_beam_np_per_m | 1.483 | 1.932 | 0.232 |
| beta_over_k0 | 1.064 | 1.065 | 0.000 |

### Diode tolerance (plus and minus twenty percent of the diode model parameters, states 0 and 1, mesh 80)

ASSUMED_SENSITIVITY_NOT_DEVICE_DATA: the scale factors multiply the series resistance, the package inductance and the junction capacitance of the kit; the MADP package inductance is NOT_AVAILABLE_FROM_CURRENT_DATA and is zero in the nominal model, so its scale has no effect there.

| scale | state | delta [mS] | change in delta [%] | dX_aniso own [ohm] | change [%] | alpha at the beam [Np/m] | change [%] |
|---|---|---|---|---|---|---|---|
| {"C_T": 0.8} | 0 | 2.280 | 1.4 | 95.2 | 1.1 | 1.990 | 3.0 |
| {"C_T": 0.8} | 1 | 2.008 | 2.0 | 81.5 | 1.5 | 1.546 | 4.2 |
| {"C_T": 1.2} | 0 | 2.218 | -1.4 | 93.2 | -1.0 | 1.876 | -2.9 |
| {"C_T": 1.2} | 1 | 1.931 | -1.9 | 79.0 | -1.5 | 1.423 | -4.0 |
| {"L_s": 0.8} | 0 | 2.249 | 0.0 | 94.2 | 0.0 | 1.932 | 0.0 |
| {"L_s": 0.8} | 1 | 1.969 | 0.0 | 80.2 | 0.0 | 1.483 | 0.0 |
| {"L_s": 1.2} | 0 | 2.249 | 0.0 | 94.2 | 0.0 | 1.932 | 0.0 |
| {"L_s": 1.2} | 1 | 1.969 | 0.0 | 80.2 | 0.0 | 1.483 | 0.0 |
| {"R_s": 0.8} | 0 | 2.261 | 0.6 | 95.0 | 0.8 | 1.953 | 1.1 |
| {"R_s": 0.8} | 1 | 1.962 | -0.3 | 79.9 | -0.4 | 1.473 | -0.6 |
| {"R_s": 1.2} | 0 | 2.234 | -0.7 | 93.1 | -1.1 | 1.908 | -1.2 |
| {"R_s": 1.2} | 1 | 1.974 | 0.3 | 80.4 | 0.2 | 1.491 | 0.5 |
| {"C_T": 0.8, "L_s": 1.2, "R_s": 1.2} | 0 | 2.265 | 0.7 | 94.1 | -0.0 | 1.966 | 1.8 |
| {"C_T": 0.8, "L_s": 1.2, "R_s": 1.2} | 1 | 2.013 | 2.2 | 81.6 | 1.7 | 1.554 | 4.8 |
| {"C_T": 1.2, "L_s": 0.8, "R_s": 0.8} | 0 | 2.231 | -0.8 | 94.0 | -0.2 | 1.896 | -1.8 |
| {"C_T": 1.2, "L_s": 0.8, "R_s": 0.8} | 1 | 1.925 | -2.2 | 78.7 | -1.9 | 1.414 | -4.6 |

Worst relative change over all tolerance cases: delta 2.2 %, dX_aniso (own loading) 1.9 %, dX_aniso (design loading) 2.6 %, alpha at the beam 4.8 %.

Package inductance of the MADP model (NOT_AVAILABLE_FROM_CURRENT_DATA, zero in the nominal model): ASSUMED values as a sensitivity, not device data.

| package inductance [nH] | state | delta [mS] | change [%] | dX_aniso own [ohm] | change [%] | alpha at the beam [Np/m] | change [%] |
|---|---|---|---|---|---|---|---|
| 0.3 | 0 | 3.509 | 56.1 | 368.6 | 291.5 | 2.207 | 14.2 |
| 0.3 | 1 | 2.769 | 40.6 | 257.3 | 220.7 | 1.176 | -20.7 |
| 0.7 | 0 | 10.104 | 349.3 | 223.6 | 137.5 | 0.699 | -63.8 |
| 0.7 | 1 | 5.141 | 161.1 | 637.1 | 694.3 | 0.083 | -94.4 |

### Contrast bandwidth (four states, mesh 80)

| frequency [GHz] | smallest dX_aniso own [ohm] | smallest dX_aniso design [ohm] | smallest alpha at the beam [Np/m] |
|---|---|---|---|
| 9.00 | 29.1 | 51.7 | 0.258 |
| 9.25 | 36.7 | 56.2 | 0.463 |
| 9.50 | 47.0 | 61.3 | 0.775 |
| 9.75 | 60.9 | 67.2 | 1.150 |
| 10.00 | 80.2 | 74.1 | 1.483 |
| 10.25 | 107.9 | 82.3 | 1.698 |
| 10.50 | 148.8 | 92.1 | 1.776 |
| 10.75 | 212.6 | 104.3 | 1.741 |
| 11.00 | 321.8 | 119.8 | 1.636 |

- smallest design-loading dX_aniso at or above the target: threshold 103.0 ohm is not met at the design frequency.
- smallest own-loading dX_aniso within 70 percent of its 10.0 GHz value (threshold 56.2 ohm): from 9.75 to 11.00 GHz, width 1.25 GHz (grid step 0.25 GHz).

### Control hardware (150 mm hexagon, same counting rules as the item-2 report)

| diodes per cell | forward diodes per state | cells | diodes | forward diodes | forward current [mA] | DC power at that current [W] | DC power at 1 mA [W] |
|---|---|---|---|---|---|---|---|
| 8 | 6 | 924 | 7392 | 5544 | 10 | 183.0 | 18.30 |

DC power = forward diodes x forward current x supply rail (3.3 V). Not modelled: the bias network (a connection to the hub and to every sector pad through chokes) and the diode package footprint; both are `NOT_AVAILABLE_FROM_CURRENT_DATA` here.

## Candidate E_rational_bowtie_8mm

Hub-and-sector cell (`coded_cell_candidates`, the family is in the table): a hub and eight sectors, each joined to the hub by a neck that carries one diode. Mode `bowtie`: state k = the two necks on axis k conduct (the hub and the two opposite sectors form a bow-tie along that axis), the other six are open. Mode `notch`: every neck conducts except the two of the perpendicular axis (the metal fill stays high, two opposite sectors float). The four states are the same cell turned in 45 degree steps, so the sheet tensor rotates with the state.

| family | period [mm] | diode kit | hub_half [mm] | neck_gap [mm] | neck_width [mm] | gap_t [mm] | arm_outer [mm] | diag_outer [mm] | clip [mm] |
|---|---|---|---|---|---|---|---|---|---|
| rational_sector_hub | 8.0 | madp | 1.60 | 0.40 | 0.80 | 0.40 | 3.60 | 3.20 | 3.80 |

### Mesh convergence (normal-incidence tensor sheet, 10.0 GHz)

| mesh | state | mean loading b [mS] | anisotropy delta [mS] | beta/k0 | dX_aniso own loading [ohm] | dX_aniso design loading [ohm] | alpha at the beam [Np/m] | largest conductance [mS] | principal axis [deg] |
|---|---|---|---|---|---|---|---|---|---|
| 40 | 0 | 4.421 | 2.336 | 1.0314 | 44.6 | 90.6 | 0.403 | 0.566 | 0.2 |
| 40 | 1 | 4.316 | 2.225 | 1.0309 | 41.5 | 85.5 | 0.344 | 0.395 | 44.8 |
| 40 | 2 | 4.421 | 2.336 | 1.0314 | 44.6 | 90.6 | 0.403 | 0.566 | 90.2 |
| 40 | 3 | 4.316 | 2.225 | 1.0309 | 41.5 | 85.5 | 0.344 | 0.395 | 134.8 |
| 60 | 0 | 4.287 | 2.258 | 1.0308 | 41.9 | 87.0 | 0.348 | 0.544 | 0.0 |
| 60 | 1 | 4.352 | 2.317 | 1.0311 | 43.6 | 89.7 | 0.381 | 0.603 | 44.8 |
| 60 | 2 | 4.287 | 2.258 | 1.0308 | 41.9 | 87.0 | 0.348 | 0.544 | 90.0 |
| 60 | 3 | 4.352 | 2.317 | 1.0311 | 43.6 | 89.7 | 0.381 | 0.603 | 134.8 |
| 80 | 0 | 4.437 | 2.321 | 1.0315 | 44.4 | 89.9 | 0.402 | 0.570 | 0.2 |
| 80 | 1 | 4.345 | 2.223 | 1.0311 | 41.6 | 85.4 | 0.349 | 0.393 | 44.8 |
| 80 | 2 | 4.437 | 2.321 | 1.0315 | 44.4 | 89.9 | 0.402 | 0.570 | 90.2 |
| 80 | 3 | 4.345 | 2.223 | 1.0311 | 41.6 | 85.4 | 0.349 | 0.393 | 134.8 |

Convergence is judged on the relative change between mesh 60 and mesh 80, largest over the four states, against the tolerance of the table (an ASSUMED criterion: 0.05 on the anisotropy, the loading and the reactance difference; the leakage at the beam grows as the square of the anisotropy, so its tolerance is twice that):

| metric | largest relative change | tolerance | converged |
|---|---|---|---|
| b_bar_ms | 0.0339 | 0.05 | yes |
| delta_ms | 0.0423 | 0.05 | yes |
| dx_aniso_ohm | 0.0548 | 0.05 | NO |
| dx_aniso_at_design_loading_ohm | 0.0509 | 0.05 | NO |
| alpha_beam_np_per_m | 0.1325 | 0.10 | NO |

All metrics converged: NO - the numbers of this candidate are not to be used until they are.

Floquet cutoff: the same four states at mesh 60 with the cutoff 2 times the mesh number (the converged setting of the Phase 20 R4 audit for lumped loads) against the cutoff equal to the mesh number: largest relative change b_bar_ms 0.0021, delta_ms 0.0020, dx_aniso_ohm 0.0037, dx_aniso_at_design_loading_ohm 0.0024, alpha_beam_np_per_m 0.0091; cutoff converged: yes.

### Targets and state balance (mesh 80)

| dX_aniso definition | smallest over the four states [ohm] | meets 103 ohm | meets 132 ohm (stretch) |
|---|---|---|---|
| own mean loading (the cell's beta) | 41.6 | NO | NO |
| design mean loading (beta/k0 of the item-2 contrast derivation) | 85.4 | NO | NO |

The same anisotropy gives a different reactance difference on a different mean loading (the surface wave is more or less bound), so both definitions are shown. The leakage the cell delivers to a beam at the design polar angle (`alpha at the beam`) is the quantity the hologram needs: compare it with the contrast requirement of the item-2 feasibility sweep at this frequency.

| quantity over the four states | smallest | largest | relative spread |
|---|---|---|---|
| b_bar_ms | 4.345 | 4.437 | 0.021 |
| delta_ms | 2.223 | 2.321 | 0.043 |
| dx_aniso_ohm | 41.648 | 44.381 | 0.062 |
| dx_aniso_at_design_loading_ohm | 85.380 | 89.950 | 0.051 |
| alpha_beam_np_per_m | 0.349 | 0.402 | 0.131 |
| beta_over_k0 | 1.031 | 1.032 | 0.000 |

### Diode tolerance (plus and minus twenty percent of the diode model parameters, states 0 and 1, mesh 80)

ASSUMED_SENSITIVITY_NOT_DEVICE_DATA: the scale factors multiply the series resistance, the package inductance and the junction capacitance of the kit; the MADP package inductance is NOT_AVAILABLE_FROM_CURRENT_DATA and is zero in the nominal model, so its scale has no effect there.

| scale | state | delta [mS] | change in delta [%] | dX_aniso own [ohm] | change [%] | alpha at the beam [Np/m] | change [%] |
|---|---|---|---|---|---|---|---|
| {"C_T": 0.8} | 0 | 2.359 | 1.6 | 44.3 | -0.3 | 0.387 | -3.5 |
| {"C_T": 0.8} | 1 | 2.268 | 2.0 | 41.8 | 0.3 | 0.341 | -2.3 |
| {"C_T": 1.2} | 0 | 2.286 | -1.5 | 44.5 | 0.2 | 0.416 | 3.5 |
| {"C_T": 1.2} | 1 | 2.179 | -2.0 | 41.5 | -0.3 | 0.357 | 2.3 |
| {"L_s": 0.8} | 0 | 2.321 | 0.0 | 44.4 | 0.0 | 0.402 | 0.0 |
| {"L_s": 0.8} | 1 | 2.223 | 0.0 | 41.6 | 0.0 | 0.349 | 0.0 |
| {"L_s": 1.2} | 0 | 2.321 | 0.0 | 44.4 | 0.0 | 0.402 | 0.0 |
| {"L_s": 1.2} | 1 | 2.223 | 0.0 | 41.6 | 0.0 | 0.349 | 0.0 |
| {"R_s": 0.8} | 0 | 2.328 | 0.3 | 44.6 | 0.4 | 0.405 | 0.9 |
| {"R_s": 0.8} | 1 | 2.228 | 0.2 | 41.8 | 0.4 | 0.352 | 0.8 |
| {"R_s": 1.2} | 0 | 2.312 | -0.4 | 44.1 | -0.6 | 0.396 | -1.3 |
| {"R_s": 1.2} | 1 | 2.216 | -0.3 | 41.5 | -0.4 | 0.346 | -1.0 |
| {"C_T": 0.8, "L_s": 1.2, "R_s": 1.2} | 0 | 2.349 | 1.2 | 44.0 | -0.9 | 0.382 | -4.8 |
| {"C_T": 0.8, "L_s": 1.2, "R_s": 1.2} | 1 | 2.261 | 1.7 | 41.6 | -0.1 | 0.338 | -3.2 |
| {"C_T": 1.2, "L_s": 0.8, "R_s": 0.8} | 0 | 2.292 | -1.3 | 44.7 | 0.6 | 0.419 | 4.4 |
| {"C_T": 1.2, "L_s": 0.8, "R_s": 0.8} | 1 | 2.185 | -1.7 | 41.7 | 0.0 | 0.360 | 3.1 |

Worst relative change over all tolerance cases: delta 2.0 %, dX_aniso (own loading) 0.9 %, dX_aniso (design loading) 2.4 %, alpha at the beam 4.8 %.

Package inductance of the MADP model (NOT_AVAILABLE_FROM_CURRENT_DATA, zero in the nominal model): ASSUMED values as a sensitivity, not device data.

| package inductance [nH] | state | delta [mS] | change [%] | dX_aniso own [ohm] | change [%] | alpha at the beam [Np/m] | change [%] |
|---|---|---|---|---|---|---|---|
| 0.3 | 0 | 3.776 | 62.7 | 105.6 | 138.0 | 2.522 | 528.2 |
| 0.3 | 1 | 3.170 | 42.6 | 74.6 | 79.1 | 1.253 | 258.9 |
| 0.7 | 0 | 11.776 | 407.3 | 205.5 | 363.0 | 13.376 | 3231.0 |
| 0.7 | 1 | 6.695 | 201.2 | 2171.6 | 5114.1 | 16.828 | 4721.1 |

### Contrast bandwidth (four states, mesh 80)

| frequency [GHz] | smallest dX_aniso own [ohm] | smallest dX_aniso design [ohm] | smallest alpha at the beam [Np/m] |
|---|---|---|---|
| 9.00 | 20.1 | 60.5 | 0.055 |
| 9.25 | 23.9 | 65.5 | 0.085 |
| 9.50 | 28.5 | 71.2 | 0.134 |
| 9.75 | 34.3 | 77.8 | 0.215 |
| 10.00 | 41.6 | 85.4 | 0.349 |
| 10.25 | 51.1 | 94.3 | 0.579 |
| 10.50 | 63.7 | 105.1 | 0.981 |
| 10.75 | 80.7 | 118.0 | 1.673 |
| 11.00 | 101.5 | 131.1 | 2.710 |

- smallest design-loading dX_aniso at or above the target: threshold 103.0 ohm is not met at the design frequency.
- smallest own-loading dX_aniso within 70 percent of its 10.0 GHz value (threshold 29.2 ohm): from 9.75 to 11.00 GHz, width 1.25 GHz (grid step 0.25 GHz).

### Control hardware (150 mm hexagon, same counting rules as the item-2 report)

| diodes per cell | forward diodes per state | cells | diodes | forward diodes | forward current [mA] | DC power at that current [W] | DC power at 1 mA [W] |
|---|---|---|---|---|---|---|---|
| 8 | 2 | 924 | 7392 | 1848 | 10 | 61.0 | 6.10 |

DC power = forward diodes x forward current x supply rail (3.3 V). Not modelled: the bias network (a connection to the hub and to every sector pad through chokes) and the diode package footprint; both are `NOT_AVAILABLE_FROM_CURRENT_DATA` here.

## Candidate F_rational_notch_6p4mm

Hub-and-sector cell (`coded_cell_candidates`, the family is in the table): a hub and eight sectors, each joined to the hub by a neck that carries one diode. Mode `bowtie`: state k = the two necks on axis k conduct (the hub and the two opposite sectors form a bow-tie along that axis), the other six are open. Mode `notch`: every neck conducts except the two of the perpendicular axis (the metal fill stays high, two opposite sectors float). The four states are the same cell turned in 45 degree steps, so the sheet tensor rotates with the state.

| family | period [mm] | diode kit | hub_half [mm] | neck_gap [mm] | neck_width [mm] | gap_t [mm] | arm_outer [mm] | diag_outer [mm] | clip [mm] |
|---|---|---|---|---|---|---|---|---|---|
| rational_sector_hub | 6.4 | madp | 1.20 | 0.40 | 0.40 | 0.20 | 3.00 | 2.80 | 3.00 |

### Mesh convergence (normal-incidence tensor sheet, 10.0 GHz)

| mesh | state | mean loading b [mS] | anisotropy delta [mS] | beta/k0 | dX_aniso own loading [ohm] | dX_aniso design loading [ohm] | alpha at the beam [Np/m] | largest conductance [mS] | principal axis [deg] |
|---|---|---|---|---|---|---|---|---|---|
| 32 | 0 | 8.547 | 1.906 | 1.0633 | 75.7 | 71.4 | 1.399 | 0.757 | 3.1 |
| 32 | 1 | 8.558 | 1.822 | 1.0634 | 72.1 | 67.9 | 1.277 | 0.808 | 41.9 |
| 32 | 2 | 8.547 | 1.906 | 1.0633 | 75.7 | 71.4 | 1.399 | 0.757 | 93.1 |
| 32 | 3 | 8.558 | 1.822 | 1.0634 | 72.1 | 67.9 | 1.277 | 0.808 | 131.9 |
| 64 | 0 | 8.745 | 2.044 | 1.0654 | 85.3 | 77.4 | 1.588 | 0.832 | 3.0 |
| 64 | 1 | 8.888 | 1.815 | 1.0670 | 76.3 | 67.6 | 1.235 | 0.890 | 41.8 |
| 64 | 2 | 8.745 | 2.044 | 1.0654 | 85.3 | 77.4 | 1.588 | 0.832 | 93.0 |
| 64 | 3 | 8.888 | 1.815 | 1.0670 | 76.3 | 67.6 | 1.235 | 0.890 | 131.8 |
| 96 | 0 | 8.732 | 2.028 | 1.0653 | 84.2 | 76.7 | 1.564 | 0.836 | 2.9 |
| 96 | 1 | 8.853 | 1.819 | 1.0666 | 76.0 | 67.8 | 1.245 | 0.895 | 41.9 |
| 96 | 2 | 8.732 | 2.028 | 1.0653 | 84.2 | 76.7 | 1.564 | 0.836 | 92.9 |
| 96 | 3 | 8.853 | 1.819 | 1.0666 | 76.0 | 67.8 | 1.245 | 0.895 | 131.9 |

Convergence is judged on the relative change between mesh 64 and mesh 96, largest over the four states, against the tolerance of the table (an ASSUMED criterion: 0.05 on the anisotropy, the loading and the reactance difference; the leakage at the beam grows as the square of the anisotropy, so its tolerance is twice that):

| metric | largest relative change | tolerance | converged |
|---|---|---|---|
| b_bar_ms | 0.0040 | 0.05 | yes |
| delta_ms | 0.0083 | 0.05 | yes |
| dx_aniso_ohm | 0.0123 | 0.05 | yes |
| dx_aniso_at_design_loading_ohm | 0.0096 | 0.05 | yes |
| alpha_beam_np_per_m | 0.0155 | 0.10 | yes |

All metrics converged: yes.

Floquet cutoff: the same four states at mesh 64 with the cutoff 2 times the mesh number (the converged setting of the Phase 20 R4 audit for lumped loads) against the cutoff equal to the mesh number: largest relative change b_bar_ms 0.0017, delta_ms 0.0019, dx_aniso_ohm 0.0044, dx_aniso_at_design_loading_ohm 0.0022, alpha_beam_np_per_m 0.0029; cutoff converged: yes.

### Targets and state balance (mesh 64)

| dX_aniso definition | smallest over the four states [ohm] | meets 103 ohm | meets 132 ohm (stretch) |
|---|---|---|---|
| own mean loading (the cell's beta) | 76.3 | NO | NO |
| design mean loading (beta/k0 of the item-2 contrast derivation) | 67.6 | NO | NO |

The same anisotropy gives a different reactance difference on a different mean loading (the surface wave is more or less bound), so both definitions are shown. The leakage the cell delivers to a beam at the design polar angle (`alpha at the beam`) is the quantity the hologram needs: compare it with the contrast requirement of the item-2 feasibility sweep at this frequency.

| quantity over the four states | smallest | largest | relative spread |
|---|---|---|---|
| b_bar_ms | 8.745 | 8.888 | 0.016 |
| delta_ms | 1.815 | 2.044 | 0.112 |
| dx_aniso_ohm | 76.327 | 85.255 | 0.105 |
| dx_aniso_at_design_loading_ohm | 67.581 | 77.395 | 0.127 |
| alpha_beam_np_per_m | 1.235 | 1.588 | 0.222 |
| beta_over_k0 | 1.065 | 1.067 | 0.001 |

### Diode tolerance (plus and minus twenty percent of the diode model parameters, states 0 and 1, mesh 64)

ASSUMED_SENSITIVITY_NOT_DEVICE_DATA: the scale factors multiply the series resistance, the package inductance and the junction capacitance of the kit; the MADP package inductance is NOT_AVAILABLE_FROM_CURRENT_DATA and is zero in the nominal model, so its scale has no effect there.

| scale | state | delta [mS] | change in delta [%] | dX_aniso own [ohm] | change [%] | alpha at the beam [Np/m] | change [%] |
|---|---|---|---|---|---|---|---|
| {"C_T": 0.8} | 0 | 2.073 | 1.4 | 86.2 | 1.1 | 1.637 | 3.1 |
| {"C_T": 0.8} | 1 | 1.850 | 1.9 | 77.5 | 1.5 | 1.287 | 4.2 |
| {"C_T": 1.2} | 0 | 2.016 | -1.4 | 84.3 | -1.1 | 1.541 | -3.0 |
| {"C_T": 1.2} | 1 | 1.782 | -1.8 | 75.2 | -1.5 | 1.185 | -4.0 |
| {"L_s": 0.8} | 0 | 2.044 | 0.0 | 85.3 | 0.0 | 1.588 | 0.0 |
| {"L_s": 0.8} | 1 | 1.815 | 0.0 | 76.3 | 0.0 | 1.235 | 0.0 |
| {"L_s": 1.2} | 0 | 2.044 | 0.0 | 85.3 | 0.0 | 1.588 | 0.0 |
| {"L_s": 1.2} | 1 | 1.815 | 0.0 | 76.3 | 0.0 | 1.235 | 0.0 |
| {"R_s": 0.8} | 0 | 2.052 | 0.4 | 85.8 | 0.7 | 1.599 | 0.7 |
| {"R_s": 0.8} | 1 | 1.812 | -0.2 | 76.3 | -0.0 | 1.229 | -0.4 |
| {"R_s": 1.2} | 0 | 2.035 | -0.5 | 84.5 | -0.9 | 1.575 | -0.8 |
| {"R_s": 1.2} | 1 | 1.817 | 0.1 | 76.2 | -0.2 | 1.240 | 0.4 |
| {"C_T": 0.8, "L_s": 1.2, "R_s": 1.2} | 0 | 2.064 | 0.9 | 85.4 | 0.2 | 1.624 | 2.3 |
| {"C_T": 0.8, "L_s": 1.2, "R_s": 1.2} | 1 | 1.852 | 2.0 | 77.3 | 1.3 | 1.293 | 4.7 |
| {"C_T": 1.2, "L_s": 0.8, "R_s": 0.8} | 0 | 2.024 | -1.0 | 84.9 | -0.4 | 1.552 | -2.3 |
| {"C_T": 1.2, "L_s": 0.8, "R_s": 0.8} | 1 | 1.779 | -2.0 | 75.2 | -1.5 | 1.180 | -4.4 |

Worst relative change over all tolerance cases: delta 2.0 %, dX_aniso (own loading) 1.5 %, dX_aniso (design loading) 2.3 %, alpha at the beam 4.7 %.

Package inductance of the MADP model (NOT_AVAILABLE_FROM_CURRENT_DATA, zero in the nominal model): ASSUMED values as a sensitivity, not device data.

| package inductance [nH] | state | delta [mS] | change [%] | dX_aniso own [ohm] | change [%] | alpha at the beam [Np/m] | change [%] |
|---|---|---|---|---|---|---|---|
| 0.3 | 0 | 3.094 | 51.3 | 278.4 | 226.5 | 1.766 | 11.2 |
| 0.3 | 1 | 2.278 | 25.5 | 194.4 | 154.7 | 0.762 | -38.3 |
| 0.7 | 0 | 6.139 | 200.3 | 504.1 | 491.2 | 0.492 | -69.0 |
| 0.7 | 1 | 4.102 | 126.0 | 1414.8 | 1753.6 | 0.083 | -93.3 |

### Contrast bandwidth (four states, mesh 64)

| frequency [GHz] | smallest dX_aniso own [ohm] | smallest dX_aniso design [ohm] | smallest alpha at the beam [Np/m] |
|---|---|---|---|
| 9.00 | 29.7 | 49.9 | 0.288 |
| 9.25 | 37.0 | 53.6 | 0.488 |
| 9.50 | 46.6 | 57.7 | 0.753 |
| 9.75 | 59.3 | 62.3 | 1.025 |
| 10.00 | 76.3 | 67.6 | 1.235 |
| 10.25 | 99.5 | 73.6 | 1.347 |
| 10.50 | 131.8 | 80.6 | 1.365 |
| 10.75 | 177.8 | 88.8 | 1.313 |
| 11.00 | 246.5 | 98.7 | 1.219 |

- smallest design-loading dX_aniso at or above the target: threshold 103.0 ohm is not met at the design frequency.
- smallest own-loading dX_aniso within 70 percent of its 10.0 GHz value (threshold 53.4 ohm): from 9.75 to 11.00 GHz, width 1.25 GHz (grid step 0.25 GHz).

### Control hardware (150 mm hexagon, same counting rules as the item-2 report)

| diodes per cell | forward diodes per state | cells | diodes | forward diodes | forward current [mA] | DC power at that current [W] | DC power at 1 mA [W] |
|---|---|---|---|---|---|---|---|
| 8 | 6 | 1430 | 11440 | 8580 | 10 | 283.1 | 28.31 |

DC power = forward diodes x forward current x supply rail (3.3 V). Not modelled: the bias network (a connection to the hub and to every sector pad through chokes) and the diode package footprint; both are `NOT_AVAILABLE_FROM_CURRENT_DATA` here.

## Candidate G_rational_bowtie_6p4mm

Hub-and-sector cell (`coded_cell_candidates`, the family is in the table): a hub and eight sectors, each joined to the hub by a neck that carries one diode. Mode `bowtie`: state k = the two necks on axis k conduct (the hub and the two opposite sectors form a bow-tie along that axis), the other six are open. Mode `notch`: every neck conducts except the two of the perpendicular axis (the metal fill stays high, two opposite sectors float). The four states are the same cell turned in 45 degree steps, so the sheet tensor rotates with the state.

| family | period [mm] | diode kit | hub_half [mm] | neck_gap [mm] | neck_width [mm] | gap_t [mm] | arm_outer [mm] | diag_outer [mm] | clip [mm] |
|---|---|---|---|---|---|---|---|---|---|
| rational_sector_hub | 6.4 | madp | 1.20 | 0.40 | 0.40 | 0.20 | 3.00 | 2.80 | 3.00 |

### Mesh convergence (normal-incidence tensor sheet, 10.0 GHz)

| mesh | state | mean loading b [mS] | anisotropy delta [mS] | beta/k0 | dX_aniso own loading [ohm] | dX_aniso design loading [ohm] | alpha at the beam [Np/m] | largest conductance [mS] | principal axis [deg] |
|---|---|---|---|---|---|---|---|---|---|
| 32 | 0 | 4.426 | 2.327 | 1.0315 | 44.4 | 90.2 | 0.401 | 0.516 | 0.3 |
| 32 | 1 | 4.390 | 2.290 | 1.0313 | 43.4 | 88.5 | 0.380 | 0.404 | 44.9 |
| 32 | 2 | 4.426 | 2.327 | 1.0315 | 44.4 | 90.2 | 0.401 | 0.516 | 90.3 |
| 32 | 3 | 4.390 | 2.290 | 1.0313 | 43.4 | 88.5 | 0.380 | 0.404 | 134.9 |
| 64 | 0 | 4.631 | 2.469 | 1.0325 | 49.1 | 97.0 | 0.508 | 0.597 | 0.3 |
| 64 | 1 | 4.457 | 2.293 | 1.0316 | 43.9 | 88.6 | 0.396 | 0.406 | 44.8 |
| 64 | 2 | 4.631 | 2.469 | 1.0325 | 49.1 | 97.0 | 0.508 | 0.597 | 90.3 |
| 64 | 3 | 4.457 | 2.293 | 1.0316 | 43.9 | 88.6 | 0.396 | 0.406 | 134.8 |
| 96 | 0 | 4.616 | 2.446 | 1.0324 | 48.5 | 95.9 | 0.494 | 0.593 | 0.3 |
| 96 | 1 | 4.463 | 2.291 | 1.0317 | 43.9 | 88.5 | 0.397 | 0.413 | 44.8 |
| 96 | 2 | 4.616 | 2.446 | 1.0324 | 48.5 | 95.9 | 0.494 | 0.593 | 90.3 |
| 96 | 3 | 4.463 | 2.291 | 1.0317 | 43.9 | 88.5 | 0.397 | 0.413 | 134.8 |

Convergence is judged on the relative change between mesh 64 and mesh 96, largest over the four states, against the tolerance of the table (an ASSUMED criterion: 0.05 on the anisotropy, the loading and the reactance difference; the leakage at the beam grows as the square of the anisotropy, so its tolerance is twice that):

| metric | largest relative change | tolerance | converged |
|---|---|---|---|
| b_bar_ms | 0.0034 | 0.05 | yes |
| delta_ms | 0.0095 | 0.05 | yes |
| dx_aniso_ohm | 0.0134 | 0.05 | yes |
| dx_aniso_at_design_loading_ohm | 0.0118 | 0.05 | yes |
| alpha_beam_np_per_m | 0.0284 | 0.10 | yes |

All metrics converged: yes.

Floquet cutoff: the same four states at mesh 64 with the cutoff 2 times the mesh number (the converged setting of the Phase 20 R4 audit for lumped loads) against the cutoff equal to the mesh number: largest relative change b_bar_ms 0.0020, delta_ms 0.0024, dx_aniso_ohm 0.0043, dx_aniso_at_design_loading_ohm 0.0030, alpha_beam_np_per_m 0.0103; cutoff converged: yes.

### Targets and state balance (mesh 64)

| dX_aniso definition | smallest over the four states [ohm] | meets 103 ohm | meets 132 ohm (stretch) |
|---|---|---|---|
| own mean loading (the cell's beta) | 43.9 | NO | NO |
| design mean loading (beta/k0 of the item-2 contrast derivation) | 88.6 | NO | NO |

The same anisotropy gives a different reactance difference on a different mean loading (the surface wave is more or less bound), so both definitions are shown. The leakage the cell delivers to a beam at the design polar angle (`alpha at the beam`) is the quantity the hologram needs: compare it with the contrast requirement of the item-2 feasibility sweep at this frequency.

| quantity over the four states | smallest | largest | relative spread |
|---|---|---|---|
| b_bar_ms | 4.457 | 4.631 | 0.038 |
| delta_ms | 2.293 | 2.469 | 0.071 |
| dx_aniso_ohm | 43.928 | 49.146 | 0.106 |
| dx_aniso_at_design_loading_ohm | 88.625 | 96.987 | 0.086 |
| alpha_beam_np_per_m | 0.396 | 0.508 | 0.220 |
| beta_over_k0 | 1.032 | 1.033 | 0.001 |

### Diode tolerance (plus and minus twenty percent of the diode model parameters, states 0 and 1, mesh 64)

ASSUMED_SENSITIVITY_NOT_DEVICE_DATA: the scale factors multiply the series resistance, the package inductance and the junction capacitance of the kit; the MADP package inductance is NOT_AVAILABLE_FROM_CURRENT_DATA and is zero in the nominal model, so its scale has no effect there.

| scale | state | delta [mS] | change in delta [%] | dX_aniso own [ohm] | change [%] | alpha at the beam [Np/m] | change [%] |
|---|---|---|---|---|---|---|---|
| {"C_T": 0.8} | 0 | 2.514 | 1.8 | 49.2 | 0.0 | 0.493 | -3.0 |
| {"C_T": 0.8} | 1 | 2.345 | 2.2 | 44.2 | 0.6 | 0.389 | -1.9 |
| {"C_T": 1.2} | 0 | 2.425 | -1.8 | 49.1 | -0.1 | 0.523 | 3.0 |
| {"C_T": 1.2} | 1 | 2.244 | -2.1 | 43.7 | -0.6 | 0.403 | 1.8 |
| {"L_s": 0.8} | 0 | 2.469 | 0.0 | 49.1 | 0.0 | 0.508 | 0.0 |
| {"L_s": 0.8} | 1 | 2.293 | 0.0 | 43.9 | 0.0 | 0.396 | 0.0 |
| {"L_s": 1.2} | 0 | 2.469 | 0.0 | 49.1 | 0.0 | 0.508 | 0.0 |
| {"L_s": 1.2} | 1 | 2.293 | 0.0 | 43.9 | 0.0 | 0.396 | 0.0 |
| {"R_s": 0.8} | 0 | 2.478 | 0.4 | 49.4 | 0.5 | 0.514 | 1.2 |
| {"R_s": 0.8} | 1 | 2.299 | 0.2 | 44.1 | 0.4 | 0.399 | 0.8 |
| {"R_s": 1.2} | 0 | 2.457 | -0.5 | 48.8 | -0.7 | 0.500 | -1.6 |
| {"R_s": 1.2} | 1 | 2.286 | -0.3 | 43.7 | -0.4 | 0.392 | -1.0 |
| {"C_T": 0.8, "L_s": 1.2, "R_s": 1.2} | 0 | 2.503 | 1.4 | 48.8 | -0.7 | 0.485 | -4.6 |
| {"C_T": 0.8, "L_s": 1.2, "R_s": 1.2} | 1 | 2.338 | 1.9 | 44.0 | 0.1 | 0.385 | -2.9 |
| {"C_T": 1.2, "L_s": 0.8, "R_s": 0.8} | 0 | 2.434 | -1.4 | 49.4 | 0.5 | 0.530 | 4.2 |
| {"C_T": 1.2, "L_s": 0.8, "R_s": 0.8} | 1 | 2.250 | -1.9 | 43.8 | -0.2 | 0.407 | 2.7 |

Worst relative change over all tolerance cases: delta 2.2 %, dX_aniso (own loading) 0.7 %, dX_aniso (design loading) 2.7 %, alpha at the beam 4.6 %.

Package inductance of the MADP model (NOT_AVAILABLE_FROM_CURRENT_DATA, zero in the nominal model): ASSUMED values as a sensitivity, not device data.

| package inductance [nH] | state | delta [mS] | change [%] | dX_aniso own [ohm] | change [%] | alpha at the beam [Np/m] | change [%] |
|---|---|---|---|---|---|---|---|
| 0.3 | 0 | 4.055 | 64.2 | 125.5 | 155.3 | 3.473 | 583.6 |
| 0.3 | 1 | 3.275 | 42.8 | 79.9 | 81.8 | 1.457 | 267.7 |
| 0.7 | 0 | 12.307 | 398.5 | 191.0 | 288.6 | 12.361 | 2332.8 |
| 0.7 | 1 | 6.949 | 203.0 | 25506.3 | 57963.9 | 17.444 | 4303.1 |

### Contrast bandwidth (four states, mesh 64)

| frequency [GHz] | smallest dX_aniso own [ohm] | smallest dX_aniso design [ohm] | smallest alpha at the beam [Np/m] |
|---|---|---|---|
| 9.00 | 22.2 | 65.4 | 0.069 |
| 9.25 | 26.1 | 70.3 | 0.105 |
| 9.50 | 30.9 | 75.7 | 0.162 |
| 9.75 | 36.7 | 81.7 | 0.252 |
| 10.00 | 43.9 | 88.6 | 0.396 |
| 10.25 | 53.0 | 96.5 | 0.631 |
| 10.50 | 64.6 | 105.8 | 1.017 |
| 10.75 | 79.9 | 116.7 | 1.653 |
| 11.00 | 100.7 | 129.8 | 2.682 |

- smallest design-loading dX_aniso at or above the target: threshold 103.0 ohm is not met at the design frequency.
- smallest own-loading dX_aniso within 70 percent of its 10.0 GHz value (threshold 30.7 ohm): from 9.50 to 11.00 GHz, width 1.50 GHz (grid step 0.25 GHz).

### Control hardware (150 mm hexagon, same counting rules as the item-2 report)

| diodes per cell | forward diodes per state | cells | diodes | forward diodes | forward current [mA] | DC power at that current [W] | DC power at 1 mA [W] |
|---|---|---|---|---|---|---|---|
| 8 | 2 | 1430 | 11440 | 2860 | 10 | 94.4 | 9.44 |

DC power = forward diodes x forward current x supply rail (3.3 V). Not modelled: the bias network (a connection to the hub and to every sector pad through chokes) and the diode package footprint; both are `NOT_AVAILABLE_FROM_CURRENT_DATA` here.

## Candidate I_rational_notch_tuned_8mm

Hub-and-sector cell (`coded_cell_candidates`, the family is in the table): a hub and eight sectors, each joined to the hub by a neck that carries one diode. Mode `bowtie`: state k = the two necks on axis k conduct (the hub and the two opposite sectors form a bow-tie along that axis), the other six are open. Mode `notch`: every neck conducts except the two of the perpendicular axis (the metal fill stays high, two opposite sectors float). The four states are the same cell turned in 45 degree steps, so the sheet tensor rotates with the state.

| family | period [mm] | diode kit | hub_half [mm] | neck_gap [mm] | neck_width [mm] | neck_width_diag [mm] | gap_t [mm] | arm_outer [mm] | diag_outer [mm] | clip [mm] |
|---|---|---|---|---|---|---|---|---|---|---|
| rational_sector_hub | 8.0 | madp | 1.40 | 0.60 | 0.40 | 0.40 | 0.20 | 3.40 | 3.80 | 3.80 |

### Mesh convergence (normal-incidence tensor sheet, 10.0 GHz)

| mesh | state | mean loading b [mS] | anisotropy delta [mS] | beta/k0 | dX_aniso own loading [ohm] | dX_aniso design loading [ohm] | alpha at the beam [Np/m] | largest conductance [mS] | principal axis [deg] |
|---|---|---|---|---|---|---|---|---|---|
| 40 | 0 | 11.292 | 3.029 | 1.0980 | 241.5 | 126.7 | 1.969 | 83.187 | 21.4 |
| 40 | 1 | 10.137 | 2.702 | 1.0823 | 159.7 | 108.7 | 2.159 | 73.542 | 48.3 |
| 40 | 2 | 11.292 | 3.029 | 1.0980 | 241.5 | 126.7 | 1.969 | 83.187 | 111.4 |
| 40 | 3 | 10.137 | 2.702 | 1.0823 | 159.7 | 108.7 | 2.159 | 73.542 | 138.3 |
| 80 | 0 | 24.029 | 3.729 | 1.2539 | 4330.3 | 173.0 | 0.123 | 76.905 | 22.2 |
| 80 | 1 | 16.170 | 8.485 | 1.1686 | 321.8 | 1131.3 | 3.723 | 68.023 | 43.4 |
| 80 | 2 | 24.029 | 3.729 | 1.2539 | 4330.3 | 173.0 | 0.123 | 76.905 | 112.2 |
| 80 | 3 | 16.170 | 8.485 | 1.1686 | 321.8 | 1131.3 | 3.723 | 68.023 | 133.4 |

Convergence is judged on the relative change between mesh 40 and mesh 80, largest over the four states, against the tolerance of the table (an ASSUMED criterion: 0.05 on the anisotropy, the loading and the reactance difference; the leakage at the beam grows as the square of the anisotropy, so its tolerance is twice that):

| metric | largest relative change | tolerance | converged |
|---|---|---|---|
| b_bar_ms | 0.5301 | 0.05 | NO |
| delta_ms | 0.6816 | 0.05 | NO |
| dx_aniso_ohm | 0.9442 | 0.05 | NO |
| dx_aniso_at_design_loading_ohm | 0.9039 | 0.05 | NO |
| alpha_beam_np_per_m | 14.9903 | 0.10 | NO |

All metrics converged: NO - the numbers of this candidate are not to be used until they are.

Floquet cutoff: the same four states at mesh 40 with the cutoff 2 times the mesh number (the converged setting of the Phase 20 R4 audit for lumped loads) against the cutoff equal to the mesh number: largest relative change b_bar_ms 0.1378, delta_ms 0.2691, dx_aniso_ohm 0.5408, dx_aniso_at_design_loading_ohm 0.3398, alpha_beam_np_per_m 0.4016; cutoff converged: NO.

### Targets and state balance (mesh 80)

| dX_aniso definition | smallest over the four states [ohm] | meets 103 ohm | meets 132 ohm (stretch) |
|---|---|---|---|
| own mean loading (the cell's beta) | 321.8 | yes | yes |
| design mean loading (beta/k0 of the item-2 contrast derivation) | 173.0 | yes | yes |

The same anisotropy gives a different reactance difference on a different mean loading (the surface wave is more or less bound), so both definitions are shown. The leakage the cell delivers to a beam at the design polar angle (`alpha at the beam`) is the quantity the hologram needs: compare it with the contrast requirement of the item-2 feasibility sweep at this frequency.

| quantity over the four states | smallest | largest | relative spread |
|---|---|---|---|
| b_bar_ms | 16.170 | 24.029 | 0.327 |
| delta_ms | 3.729 | 8.485 | 0.561 |
| dx_aniso_ohm | 321.779 | 4330.298 | 0.926 |
| dx_aniso_at_design_loading_ohm | 173.008 | 1131.293 | 0.847 |
| alpha_beam_np_per_m | 0.123 | 3.723 | 0.967 |
| beta_over_k0 | 1.169 | 1.254 | 0.068 |

### Diode tolerance (plus and minus twenty percent of the diode model parameters, states 0 and 1, mesh 80)

ASSUMED_SENSITIVITY_NOT_DEVICE_DATA: the scale factors multiply the series resistance, the package inductance and the junction capacitance of the kit; the MADP package inductance is NOT_AVAILABLE_FROM_CURRENT_DATA and is zero in the nominal model, so its scale has no effect there.

| scale | state | delta [mS] | change in delta [%] | dX_aniso own [ohm] | change [%] | alpha at the beam [Np/m] | change [%] |
|---|---|---|---|---|---|---|---|
| {"C_T": 0.8} | 0 | 3.799 | 1.9 | 3406.3 | -21.3 | 0.127 | 3.2 |
| {"C_T": 0.8} | 1 | 8.550 | 0.8 | 318.2 | -1.1 | 3.845 | 3.3 |
| {"C_T": 1.2} | 0 | 3.660 | -1.8 | 5796.2 | 33.9 | 0.119 | -3.1 |
| {"C_T": 1.2} | 1 | 8.421 | -0.8 | 325.4 | 1.1 | 3.604 | -3.2 |
| {"L_s": 0.8} | 0 | 3.729 | 0.0 | 4330.3 | 0.0 | 0.123 | 0.0 |
| {"L_s": 0.8} | 1 | 8.485 | 0.0 | 321.8 | 0.0 | 3.723 | 0.0 |
| {"L_s": 1.2} | 0 | 3.729 | 0.0 | 4330.3 | 0.0 | 0.123 | 0.0 |
| {"L_s": 1.2} | 1 | 8.485 | 0.0 | 321.8 | 0.0 | 3.723 | 0.0 |
| {"R_s": 0.8} | 0 | 4.730 | 26.9 | 739.9 | -82.9 | 0.043 | -65.2 |
| {"R_s": 0.8} | 1 | 13.109 | 54.5 | 166.2 | -48.3 | 2.894 | -22.3 |
| {"R_s": 1.2} | 0 | 3.282 | -12.0 | 1170.4 | -73.0 | 0.323 | 162.7 |
| {"R_s": 1.2} | 1 | 5.840 | -31.2 | 1343.2 | 317.4 | 3.783 | 1.6 |
| {"C_T": 0.8, "L_s": 1.2, "R_s": 1.2} | 0 | 3.342 | -10.4 | 1285.0 | -70.3 | 0.334 | 171.4 |
| {"C_T": 0.8, "L_s": 1.2, "R_s": 1.2} | 1 | 5.905 | -30.4 | 1276.1 | 296.6 | 3.941 | 5.9 |
| {"C_T": 1.2, "L_s": 0.8, "R_s": 0.8} | 0 | 4.649 | 24.7 | 773.6 | -82.1 | 0.042 | -66.1 |
| {"C_T": 1.2, "L_s": 0.8, "R_s": 0.8} | 1 | 13.045 | 53.7 | 167.1 | -48.1 | 2.826 | -24.1 |

Worst relative change over all tolerance cases: delta 54.5 %, dX_aniso (own loading) 317.4 %, dX_aniso (design loading) 79.8 %, alpha at the beam 171.4 %.

Package inductance of the MADP model (NOT_AVAILABLE_FROM_CURRENT_DATA, zero in the nominal model): ASSUMED values as a sensitivity, not device data.

| package inductance [nH] | state | delta [mS] | change [%] | dX_aniso own [ohm] | change [%] | alpha at the beam [Np/m] | change [%] |
|---|---|---|---|---|---|---|---|
| 0.3 | 0 | 6.255 | 67.8 | 15.6 | -99.6 | 0.006 | -95.2 |
| 0.3 | 1 | 10.372 | 22.2 | 343.5 | 6.7 | 1.891 | -49.2 |
| 0.7 | 0 | 9.703 | 160.2 | 21.0 | -99.5 | 0.009 | -93.1 |
| 0.7 | 1 | 7.109 | -16.2 | 11.4 | -96.5 | 0.002 | -99.9 |

### Contrast bandwidth (four states, mesh 80)

| frequency [GHz] | smallest dX_aniso own [ohm] | smallest dX_aniso design [ohm] | smallest alpha at the beam [Np/m] |
|---|---|---|---|
| 9.00 | 249.3 | 50.3 | 0.014 |
| 9.25 | 180.8 | 63.8 | 0.009 |
| 9.50 | 136.9 | 91.6 | 0.007 |
| 9.75 | 120.1 | 142.1 | 0.009 |
| 10.00 | 321.8 | 173.0 | 0.123 |
| 10.25 | 2.7 | 64.5 | 0.000 |
| 10.50 | 3.6 | 123.8 | 0.000 |
| 10.75 | 4.5 | 139.3 | 0.000 |
| 11.00 | 5.9 | 174.0 | 0.001 |

- smallest design-loading dX_aniso at or above the target (threshold 103.0 ohm): from 9.75 to 10.00 GHz, width 0.25 GHz (grid step 0.25 GHz).
- smallest own-loading dX_aniso within 70 percent of its 10.0 GHz value (threshold 225.2 ohm): from 10.00 to 10.00 GHz, width 0.00 GHz (grid step 0.25 GHz).

### Control hardware (150 mm hexagon, same counting rules as the item-2 report)

| diodes per cell | forward diodes per state | cells | diodes | forward diodes | forward current [mA] | DC power at that current [W] | DC power at 1 mA [W] |
|---|---|---|---|---|---|---|---|
| 8 | 6 | 924 | 7392 | 5544 | 10 | 183.0 | 18.30 |

DC power = forward diodes x forward current x supply rail (3.3 V). Not modelled: the bias network (a connection to the hub and to every sector pad through chokes) and the diode package footprint; both are `NOT_AVAILABLE_FROM_CURRENT_DATA` here.

## Candidate J_rational_notch_wide_8mm

Hub-and-sector cell (`coded_cell_candidates`, the family is in the table): a hub and eight sectors, each joined to the hub by a neck that carries one diode. Mode `bowtie`: state k = the two necks on axis k conduct (the hub and the two opposite sectors form a bow-tie along that axis), the other six are open. Mode `notch`: every neck conducts except the two of the perpendicular axis (the metal fill stays high, two opposite sectors float). The four states are the same cell turned in 45 degree steps, so the sheet tensor rotates with the state.

| family | period [mm] | diode kit | hub_half [mm] | neck_gap [mm] | neck_width [mm] | neck_width_diag [mm] | gap_t [mm] | arm_outer [mm] | diag_outer [mm] | clip [mm] |
|---|---|---|---|---|---|---|---|---|---|---|
| rational_sector_hub | 8.0 | madp | 1.40 | 0.60 | 0.80 | 0.40 | 0.20 | 3.60 | 3.20 | 3.80 |

### Mesh convergence (normal-incidence tensor sheet, 10.0 GHz)

| mesh | state | mean loading b [mS] | anisotropy delta [mS] | beta/k0 | dX_aniso own loading [ohm] | dX_aniso design loading [ohm] | alpha at the beam [Np/m] | largest conductance [mS] | principal axis [deg] |
|---|---|---|---|---|---|---|---|---|---|
| 40 | 0 | 12.154 | 3.241 | 1.1105 | 328.8 | 139.4 | 1.739 | 2.129 | 6.3 |
| 40 | 1 | 11.476 | 2.932 | 1.1007 | 237.2 | 121.1 | 1.747 | 2.326 | 38.4 |
| 40 | 2 | 12.154 | 3.241 | 1.1105 | 328.8 | 139.4 | 1.739 | 2.129 | 96.3 |
| 40 | 3 | 11.476 | 2.932 | 1.1007 | 237.2 | 121.1 | 1.747 | 2.326 | 128.4 |
| 80 | 0 | 12.088 | 3.218 | 1.1095 | 319.5 | 138.0 | 1.749 | 2.124 | 6.1 |
| 80 | 1 | 11.415 | 2.912 | 1.0998 | 231.5 | 120.0 | 1.754 | 2.334 | 38.5 |
| 80 | 2 | 12.088 | 3.218 | 1.1095 | 319.5 | 138.0 | 1.749 | 2.124 | 96.1 |
| 80 | 3 | 11.415 | 2.912 | 1.0998 | 231.5 | 120.0 | 1.754 | 2.334 | 128.5 |

Convergence is judged on the relative change between mesh 40 and mesh 80, largest over the four states, against the tolerance of the table (an ASSUMED criterion: 0.05 on the anisotropy, the loading and the reactance difference; the leakage at the beam grows as the square of the anisotropy, so its tolerance is twice that):

| metric | largest relative change | tolerance | converged |
|---|---|---|---|
| b_bar_ms | 0.0055 | 0.05 | yes |
| delta_ms | 0.0072 | 0.05 | yes |
| dx_aniso_ohm | 0.0291 | 0.05 | yes |
| dx_aniso_at_design_loading_ohm | 0.0104 | 0.05 | yes |
| alpha_beam_np_per_m | 0.0057 | 0.10 | yes |

All metrics converged: yes.

Floquet cutoff: the same four states at mesh 40 with the cutoff 2 times the mesh number (the converged setting of the Phase 20 R4 audit for lumped loads) against the cutoff equal to the mesh number: largest relative change b_bar_ms 0.0021, delta_ms 0.0032, dx_aniso_ohm 0.0090, dx_aniso_at_design_loading_ohm 0.0046, alpha_beam_np_per_m 0.0117; cutoff converged: yes.

### Targets and state balance (mesh 80)

| dX_aniso definition | smallest over the four states [ohm] | meets 103 ohm | meets 132 ohm (stretch) |
|---|---|---|---|
| own mean loading (the cell's beta) | 231.5 | yes | yes |
| design mean loading (beta/k0 of the item-2 contrast derivation) | 120.0 | yes | NO |

The same anisotropy gives a different reactance difference on a different mean loading (the surface wave is more or less bound), so both definitions are shown. The leakage the cell delivers to a beam at the design polar angle (`alpha at the beam`) is the quantity the hologram needs: compare it with the contrast requirement of the item-2 feasibility sweep at this frequency.

| quantity over the four states | smallest | largest | relative spread |
|---|---|---|---|
| b_bar_ms | 11.415 | 12.088 | 0.056 |
| delta_ms | 2.912 | 3.218 | 0.095 |
| dx_aniso_ohm | 231.493 | 319.492 | 0.275 |
| dx_aniso_at_design_loading_ohm | 119.988 | 137.966 | 0.130 |
| alpha_beam_np_per_m | 1.749 | 1.754 | 0.003 |
| beta_over_k0 | 1.100 | 1.110 | 0.009 |

### Diode tolerance (plus and minus twenty percent of the diode model parameters, states 0 and 1, mesh 80)

ASSUMED_SENSITIVITY_NOT_DEVICE_DATA: the scale factors multiply the series resistance, the package inductance and the junction capacitance of the kit; the MADP package inductance is NOT_AVAILABLE_FROM_CURRENT_DATA and is zero in the nominal model, so its scale has no effect there.

| scale | state | delta [mS] | change in delta [%] | dX_aniso own [ohm] | change [%] | alpha at the beam [Np/m] | change [%] |
|---|---|---|---|---|---|---|---|
| {"C_T": 0.8} | 0 | 3.246 | 0.9 | 322.9 | 1.1 | 1.796 | 2.7 |
| {"C_T": 0.8} | 1 | 2.946 | 1.2 | 234.5 | 1.3 | 1.814 | 3.4 |
| {"C_T": 1.2} | 0 | 3.191 | -0.8 | 316.2 | -1.0 | 1.704 | -2.6 |
| {"C_T": 1.2} | 1 | 2.878 | -1.2 | 228.6 | -1.3 | 1.696 | -3.3 |
| {"L_s": 0.8} | 0 | 3.218 | 0.0 | 319.5 | 0.0 | 1.749 | 0.0 |
| {"L_s": 0.8} | 1 | 2.912 | 0.0 | 231.5 | 0.0 | 1.754 | 0.0 |
| {"L_s": 1.2} | 0 | 3.218 | 0.0 | 319.5 | 0.0 | 1.749 | 0.0 |
| {"L_s": 1.2} | 1 | 2.912 | 0.0 | 231.5 | 0.0 | 1.754 | 0.0 |
| {"R_s": 0.8} | 0 | 3.247 | 0.9 | 329.1 | 3.0 | 1.753 | 0.2 |
| {"R_s": 0.8} | 1 | 2.911 | -0.0 | 232.7 | 0.5 | 1.740 | -0.8 |
| {"R_s": 1.2} | 0 | 3.183 | -1.1 | 307.8 | -3.6 | 1.748 | -0.1 |
| {"R_s": 1.2} | 1 | 2.909 | -0.1 | 229.0 | -1.1 | 1.775 | 1.2 |
| {"C_T": 0.8, "L_s": 1.2, "R_s": 1.2} | 0 | 3.210 | -0.2 | 311.1 | -2.6 | 1.795 | 2.6 |
| {"C_T": 0.8, "L_s": 1.2, "R_s": 1.2} | 1 | 2.943 | 1.1 | 231.9 | 0.2 | 1.836 | 4.7 |
| {"C_T": 1.2, "L_s": 0.8, "R_s": 0.8} | 0 | 3.220 | 0.1 | 325.7 | 1.9 | 1.708 | -2.3 |
| {"C_T": 1.2, "L_s": 0.8, "R_s": 0.8} | 1 | 2.878 | -1.2 | 229.7 | -0.8 | 1.683 | -4.0 |

Worst relative change over all tolerance cases: delta 1.2 %, dX_aniso (own loading) 3.6 %, dX_aniso (design loading) 1.6 %, alpha at the beam 4.7 %.

Package inductance of the MADP model (NOT_AVAILABLE_FROM_CURRENT_DATA, zero in the nominal model): ASSUMED values as a sensitivity, not device data.

| package inductance [nH] | state | delta [mS] | change [%] | dX_aniso own [ohm] | change [%] | alpha at the beam [Np/m] | change [%] |
|---|---|---|---|---|---|---|---|
| 0.3 | 0 | 6.342 | 97.1 | 471.4 | 47.5 | 0.523 | -70.1 |
| 0.3 | 1 | 5.716 | 96.3 | 603.2 | 160.6 | 0.486 | -72.3 |
| 0.7 | 0 | 11.580 | 259.8 | 11.2 | -96.5 | 0.002 | -99.9 |
| 0.7 | 1 | 12.207 | 319.3 | 14.9 | -93.6 | 0.003 | -99.8 |

### Contrast bandwidth (four states, mesh 80)

| frequency [GHz] | smallest dX_aniso own [ohm] | smallest dX_aniso design [ohm] | smallest alpha at the beam [Np/m] |
|---|---|---|---|
| 9.00 | 48.5 | 66.6 | 0.816 |
| 9.25 | 67.0 | 75.4 | 1.251 |
| 9.50 | 96.0 | 86.4 | 1.592 |
| 9.75 | 143.8 | 100.7 | 1.755 |
| 10.00 | 231.5 | 120.0 | 1.749 |
| 10.25 | 429.2 | 147.7 | 1.469 |
| 10.50 | 1249.6 | 191.3 | 1.175 |
| 10.75 | 906.2 | 271.2 | 0.905 |
| 11.00 | 454.4 | 469.3 | 0.678 |

- smallest design-loading dX_aniso at or above the target (threshold 103.0 ohm): from 10.00 to 11.00 GHz, width 1.00 GHz (grid step 0.25 GHz).
- smallest own-loading dX_aniso within 70 percent of its 10.0 GHz value (threshold 162.0 ohm): from 10.00 to 11.00 GHz, width 1.00 GHz (grid step 0.25 GHz).

### Control hardware (150 mm hexagon, same counting rules as the item-2 report)

| diodes per cell | forward diodes per state | cells | diodes | forward diodes | forward current [mA] | DC power at that current [W] | DC power at 1 mA [W] |
|---|---|---|---|---|---|---|---|
| 8 | 6 | 924 | 7392 | 5544 | 10 | 183.0 | 18.30 |

DC power = forward diodes x forward current x supply rail (3.3 V). Not modelled: the bias network (a connection to the hub and to every sector pad through chokes) and the diode package footprint; both are `NOT_AVAILABLE_FROM_CURRENT_DATA` here.

## Sensitivity to the package inductance of the diode (candidate J_rational_sector_notch_wide_8mm)

ASSUMED_SENSITIVITY_NOT_DEVICE_DATA. The MADP model has no package inductance in the repo (NOT_AVAILABLE_FROM_CURRENT_DATA, zero in the nominal model); the ON state of a diode is a series resistance plus a reactance that grows with the inductance (about six ohm per tenth of a nanohenry at 10.0 GHz). The documented SMP1345-079LF model (package inductance 0.7 nH) is shown for reference. Mesh 40, states 0 and 1.

| diode model | package inductance [nH] | state | delta [mS] | dX_aniso own [ohm] | dX_aniso design [ohm] | alpha at the beam [Np/m] |
|---|---|---|---|---|---|---|
| madp | 0.00 | 0 | 3.24 | 328.8 | 139.4 | 1.739 |
| madp | 0.00 | 1 | 2.93 | 237.2 | 121.1 | 1.747 |
| madp | 0.05 | 0 | 3.52 | 497.1 | 157.4 | 1.515 |
| madp | 0.05 | 1 | 3.24 | 351.5 | 139.3 | 1.578 |
| madp | 0.10 | 0 | 3.86 | 938.3 | 183.0 | 1.288 |
| madp | 0.10 | 1 | 3.59 | 594.4 | 162.5 | 1.358 |
| madp | 0.15 | 0 | 4.28 | 4879.7 | 221.3 | 1.066 |
| madp | 0.15 | 1 | 3.99 | 1452.0 | 194.1 | 1.119 |
| madp | 0.20 | 0 | 4.82 | 1649.4 | 284.7 | 0.856 |
| madp | 0.20 | 1 | 4.47 | 5280.4 | 241.3 | 0.882 |
| madp | 0.30 | 0 | 6.40 | 459.8 | 799.5 | 0.502 |
| madp | 0.30 | 1 | 5.82 | 572.8 | 504.0 | 0.481 |
| madp | 0.50 | 0 | 14.04 | 150.4 | 200.5 | 0.671 |
| madp | 0.50 | 1 | 17.60 | 118.2 | 139.3 | 2.286 |
| madp | 0.70 | 0 | 11.69 | 11.1 | 293.3 | 0.001 |
| madp | 0.70 | 1 | 12.89 | 15.3 | 236.4 | 0.003 |
| smp | documented 0.7 nH | 0 | 31.25 | 13.3 | 68.0 | 0.001 |
| smp | documented 0.7 nH | 1 | 48.77 | 24.6 | 42.0 | 0.002 |

