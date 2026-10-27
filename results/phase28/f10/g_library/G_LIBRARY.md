# G state library, 9.5 - 10.5 GHz (Phase 28, HFSS-first path)

Labels: MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET, VIA_OMITTED_APPROXIMATION. **Diode package inductance L_s = 0 nH** (`NOT_AVAILABLE_FROM_CURRENT_DATA`; zero is an assumption, the sensitivity is in `../g_refine/`).

Cell G: 6.4 mm, rational-sector hub, 'bowtie' mode, eight MADP-000907-14020 diode sites (ON = 5.2 ohm, OFF = 0.025 pF), RT/duroid 5880 1.575 mm on PEC. Four states = slot axis 0, 45, 90, 135 deg (two diodes forward per state).
The sheet tensor is the normal-incidence one of `coded_cell_eval` (the Phase 26 B screening method); `alpha_beam` is the leakage of a row of cells to a beam at 30 deg (closed form), not a dispersion.

## Which numbers were reused and which were recomputed
The item-3 characterisation (`../cells/cell_G_6p4_bowtie.json`) stored only derived metrics (no Y_s / reflection matrices), at mesh 96 only for 10.0 GHz. This library recomputes every state at **mesh 64 and mesh 96 at all five frequencies** (40 solves, Floquet cutoff = mesh), plus the oblique solve; the 10.0 GHz values agree with the item-3 file (same code).

## Mesh convergence (mesh 64 -> mesh 96, maximum over the four states; criterion 5 %, 10 % for alpha, an ASSUMED criterion as in item 3)
| f [GHz] | b_bar | delta | dX own | dX design | alpha_beam | converged |
|---|---|---|---|---|---|---|
| 9.5 | 0.25 % | 0.82 % | 1.05 % | 0.97 % | 2.23 % | yes |
| 9.75 | 0.29 % | 0.88 % | 1.18 % | 1.07 % | 2.51 % | yes |
| 10 | 0.34 % | 0.95 % | 1.34 % | 1.18 % | 2.84 % | yes |
| 10.25 | 0.39 % | 1.03 % | 1.54 % | 1.31 % | 3.21 % | yes |
| 10.5 | 0.45 % | 1.11 % | 1.78 % | 1.47 % | 3.62 % | yes |

## Library (mesh 96; sheet tensor Y_s in the file `g_library.json`, rows `kind: normal`)
| f [GHz] | state | axis [deg] | b_bar [mS] | delta [mS] | g_max [mS] | beta/k0 | dX own [ohm] | dX design [ohm] | alpha_beam [Np/m] |
|---|---|---|---|---|---|---|---|---|---|
| 9.5 | 0 | 0.26 | 4.212 | 2.161 | 0.467 | 1.0266 | 34.4 | 82.6 | 0.205 |
| 9.5 | 1 | 44.84 | 4.058 | 2.006 | 0.322 | 1.0260 | 31.0 | 75.7 | 0.163 |
| 9.5 | 2 | 90.26 | 4.212 | 2.161 | 0.467 | 1.0266 | 34.4 | 82.6 | 0.205 |
| 9.5 | 3 | 134.84 | 4.058 | 2.006 | 0.322 | 1.0260 | 31.0 | 75.7 | 0.163 |
| 9.75 | 0 | 0.27 | 4.408 | 2.297 | 0.526 | 1.0294 | 40.7 | 88.8 | 0.317 |
| 9.75 | 1 | 44.83 | 4.254 | 2.142 | 0.364 | 1.0287 | 36.8 | 81.7 | 0.253 |
| 9.75 | 2 | 90.27 | 4.408 | 2.297 | 0.526 | 1.0294 | 40.7 | 88.8 | 0.317 |
| 9.75 | 3 | 134.83 | 4.254 | 2.142 | 0.364 | 1.0287 | 36.8 | 81.7 | 0.253 |
| 10 | 0 | 0.29 | 4.616 | 2.446 | 0.593 | 1.0324 | 48.5 | 95.9 | 0.494 |
| 10 | 1 | 44.82 | 4.463 | 2.291 | 0.413 | 1.0317 | 43.9 | 88.5 | 0.397 |
| 10 | 2 | 90.29 | 4.616 | 2.446 | 0.593 | 1.0324 | 48.5 | 95.9 | 0.494 |
| 10 | 3 | 134.82 | 4.463 | 2.291 | 0.413 | 1.0317 | 43.9 | 88.5 | 0.397 |
| 10.25 | 0 | 0.30 | 4.838 | 2.608 | 0.671 | 1.0359 | 58.3 | 103.9 | 0.777 |
| 10.25 | 1 | 44.82 | 4.687 | 2.456 | 0.470 | 1.0350 | 52.9 | 96.4 | 0.631 |
| 10.25 | 2 | 90.30 | 4.838 | 2.608 | 0.671 | 1.0359 | 58.3 | 103.9 | 0.777 |
| 10.25 | 3 | 134.82 | 4.687 | 2.456 | 0.470 | 1.0350 | 52.9 | 96.4 | 0.631 |
| 10.5 | 0 | 0.32 | 5.076 | 2.785 | 0.762 | 1.0399 | 70.7 | 113.1 | 1.230 |
| 10.5 | 1 | 44.81 | 4.930 | 2.638 | 0.538 | 1.0389 | 64.5 | 105.4 | 1.013 |
| 10.5 | 2 | 90.32 | 5.076 | 2.785 | 0.762 | 1.0399 | 70.7 | 113.1 | 1.230 |
| 10.5 | 3 | 134.81 | 4.930 | 2.638 | 0.538 | 1.0389 | 64.5 | 105.4 | 1.013 |

`dX design` is the anisotropy at the item-2 design loading (beta/k0 = 1.06); the 10.0 GHz item-2 requirements are 90.9 ohm (1.5 dB total) and 118.4 ohm (1 dB total), and the leakage the surface has to reach is 2.09 - 3.66 Np/m (`../FEASIBILITY_REPORT.md`).
