# phase28_coded_small -- PLUMBING TEST: 70 mm coded hexagon with G cells, centre monopole

**PLUMBING TEST, NOT A BEAM TEST.** The questions this board answers: does the model build, are the 816 diode states applied by `set_coding.py`, does it solve (S11, far field). It is not expected to show the designed beam: G leaks only about 0.45 Np/m, so only about 3 % of the surface wave converts on this board.

State: PREDICTION_NOT_HFSS_VALIDATED (MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET, VIA_OMITTED_APPROXIMATION). Diode package inductance 0 nH assumed (NOT_AVAILABLE_FROM_CURRENT_DATA). 10.0 GHz. Predictions below were written before any HFSS run.

## Read this first: what this board can and cannot show
G leaks only alpha = 0.45 Np/m at 10.0 GHz (mean of its four states; the surface has to reach 2.09 - 3.66 Np/m). On a board of this size the aperture therefore converts **2.6 %** of the surface wave into the beam (the full 10.0 GHz supercell table would give 16.5 %), so the predicted RHCP gain at the beam is **-0.6 dBi** for 30/225 (about 7.4 dBi if G could deliver the table's leakage). The monopole's own direct radiation (about 0.4 of the accepted power on the 150 mm board, not modelled here) is expected to be comparable to or larger than the aperture beam. **Expect a weak or hidden beam.** The pass criteria below are therefore PLUMBING criteria first (build, boundaries, solve, coding switch); the beam criteria are kept as an optional extra, not the purpose.

**Size decision (made unattended, to be confirmed):** '70 mm hexagon' is taken as 70 mm *across flats* (circumradius 40.41 mm, vertex at +-y, as the Phase 26 hexagon), which gives **102 cells and 816 diode boundaries** (limit ~1,200). Taking 70 mm as the circumradius would give 2448 boundaries; the largest circumradius within the limit is 50.4 mm (1200 boundaries).

## Files
| File | Content |
|---|---|
| `build_board.py` | AEDT script, builds the geometry ONCE: substrate/ground (board + 5 mm margin), monopole (Phase 26 launch: pin 6.0 mm above the substrate = `mono_h` 7.575 mm, `r_pin` 0.635 mm, coax annulus `r_out` 2.05 mm, lumped port 50 ohm), 9 metal sheets and 8 diode sheets per cell, one Lumped RLC boundary per diode sheet (initial coding: 30/225), radiation box, setup 10 GHz + interpolating sweep 9.5 - 10.5 GHz, far-field sphere theta 0 - 180, phi -180 .. 180, 1 deg |
| `set_coding.py` | AEDT script: sets all diode states from a coding CSV (edit `CODING_CSV`); the model is not rebuilt |
| `layout.csv` | cell_id, x, y [mm] (lattice 6.4 mm, node at the feed, feed keep-out 6.5 mm) |
| `codings/coding_beam_000.csv`, `coding_beam_030_225.csv`, `coding_beam_045_225.csv` | cell_id, x, y, state 0 .. 3, slot axis [deg] (state k = two necks ON, axis 45 k deg) |
| `mom_reference/predictions_*.json`, `cut_*.csv`, `summary.json` | the prediction per coding; RHCP/LHCP gain cut in the beam plane (phi0 and phi0 + 180), 1 deg |

## Run order
1. **First** the single-cell check `../phase28_cell_G/` (it tests the diode boundaries, the geometry and the cell's own states with a Floquet port in minutes).
2. Run `build_board.py` in AEDT (building 1734 sheets takes a few minutes). Check `Port1`, `Rad1`, that `D0_0_RLC .. D101_7_RLC` exist (816; 204 with R = 5.2 ohm, the rest C = 0.025 pF), validate.
3. Analyze (setup 10 GHz + sweep). Export (per coding): realised gain RHCP and LHCP and axial ratio over the sphere (1 deg CSV), the S11 sweep, the antenna-parameters table (accepted/radiated power, efficiency); optional: Mag_E and complex Ez along rho = 10 .. 35 mm at phi = 225 deg, 1 mm above the board.
4. Edit `CODING_CSV` in `set_coding.py` to the next coding (full path), run it on the open model (watch the message: how many assignments went through `EditLumpedRLC` and how many via delete + assign), Analyze again, export. Repeat for the three codings. `EditLumpedRLC` / `DeleteBoundaries` are not yet run in AEDT.

## Predictions (MODEL_ONLY, 10.0 GHz, RHCP; fixed launch fraction 0.592 of the 10.54 GHz bare slab; G-limited conversion)
| Coding (theta / phi) | Peak (theta, phi) | Pointing error | Gain [dBi] | AR at peak [dB] | XPD at peak [dB] | HPBW (2 cuts) [deg] | SLL > 2 HPBW away [dB] | RHCP share of radiated | 2-bit quantisation loss [dB] |
|---|---|---|---|---|---|---|---|---|---|
| 0 deg (broadside) | 0.0, 0.0 | 0.0 deg | -0.07 | 1.66 | 20.42 | 23.0 / 23.0 | -21.6 | 0.95 | 0.85 |
| 30 deg / 225 deg | 25.6, 228.0 | 4.6 deg | -0.56 | 2.71 | 16.22 | 24.5 / 22.0 | -18.0 | 0.91 | 0.78 |
| 45 deg / 225 deg | 37.8, 224.3 | 7.2 deg | -1.58 | 4.29 | 12.31 | 25.5 / 22.0 | -15.5 | 0.84 | 1.03 |

If G could deliver the leakage of the 10.0 GHz supercell table (not the case): gain 0 deg (broadside) 7.4 dBi, 30 deg / 225 deg 7.4 dBi, 45 deg / 225 deg 6.6 dBi. Aperture directivity D_ap: 0 deg (broadside) 18.0 dBi, 30 deg / 225 deg 17.6 dBi, 45 deg / 225 deg 16.5 dBi.
The SLL of `analyse` (15 deg exclusion) is not used: at a 23 deg half-power width it is inside the main-lobe skirt. State counts per coding are in the predictions files.

**Known model-versus-HFSS gaps (Phase 26, `results/phase26/hex150_tensor_cp_single_hfss_comparison.md`, 150 mm board, 10.54 GHz, one HFSS run):** the model's RHCP sidelobe level was about 4 dB optimistic (-22.8 dB predicted, -18.5 dB in HFSS); the cross-polar (LHCP) power was strongly underestimated (0.31 of the radiated power in HFSS against 0.008 in the model, because the model gives no pattern to the pin's direct radiation); the peak gain agreed within 0.5 dBi (23.0 against 23.5), the pointing exactly, XPD at the peak was better in HFSS (23.0 against 16.9 dB). Those gaps are for a board whose aperture beam dominated; here the pin's share is relatively larger.

## Pass criteria (ASSUMED, stated before the run)
* **PLUMBING PASS:** the build finishes without the RLC-failure message; the boundary list shows the expected R / C counts; `set_coding.py` completes for all three codings; every coding solves (converged, S11 exported); the three far fields are not identical (the states reach the solver).
Optional beam criteria (only meaningful if the beam is not hidden):
* **PASS (per coding):** the RHCP gain maximum in the half space is within half a beamwidth of the target (about 12 deg for 0 deg (broadside), 11 deg for 30 deg / 225 deg, 11 deg for 45 deg / 225 deg) **and** RHCP exceeds LHCP at that maximum.
* **INCONCLUSIVE (expected possibility):** no RHCP maximum near the target because the monopole's direct radiation dominates. Then compare the RHCP patterns of the three codings with each other: the common (pin) part cancels in the difference, and the prediction is that the three differ in the sense of `mom_reference/cut_*.csv` (beam moves from broadside toward 225 deg as theta0 grows). Report the sidelobe/beam positions you see.
* **FAIL:** LHCP dominant at the maximum for all codings (a handedness or state-order error: report it, it is the most useful result), or the maximum moves the wrong way (away from phi 225 deg) when theta0 grows.

## Honest limits
The prediction is the idealised aperture model of Phase 26 (local response, ideal state tensors, no feed scattering, no edge diffraction at 70 mm); the leakage of the four G states differs by up to 25 % and is taken equal; 0.592 is the 10.54 GHz bare-slab launch fraction, not recomputed for 10.0 GHz or for this board; the table minimum and maximum are those of the 10.0 GHz table; flip-chip pads, bias lines, via and the connector are not modelled; the diode package inductance is unknown and G's contrast disappears for a large one (see the cell package).
