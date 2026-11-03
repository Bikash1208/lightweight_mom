# phase30_testboard - single-feed test board (MODEL_ONLY, NOT_HFSS_VALIDATED, not yet run in AEDT)

Purpose: separate (a) cell leakage from (b) the monopole's direct radiation and (c) edge diffraction. The Phase 29 HFSS boards put 48-68 %
of the power at theta 60-90 deg with ~98 % E_theta; the model cannot tell which of (a)/(b)/(c) causes it (no launch-fraction or
diffraction model in the repo).  Three runs, same outline, same feed, same substrate:

| Run | Script | Cells | What it isolates |
|---|---|---|---|
| ACTIVE | `build_board_active.py` | 572 squircle cells, slot length/width/angle from the ray-march design (taper: M_t 0.40 at 4-12 mm to 0.62 beyond 20 mm), broadside RHCP | everything |
| CONTROL | `build_board_control.py` | same 572 positions, every cell the unslotted squircle (isotropic, no leakage) | (b)+(c) only: monopole direct radiation + edge diffraction of the un-leaked wave |
| BARE (optional, edit CONTROL: delete the CELLS list) | -- | grounded RO3010, no patches | (b) + edge diffraction of the bare-slab TM0 wave (beta ~ k0, weakly bound) |

Board: 90 mm circle of 3.2 mm square-lattice cells (feed at a lattice-cell corner, 4-fold symmetric), RO3010 0.635 mm on a 100 x 100 mm
substrate/ground, coax-probe monopole (pin r 0.635 mm, 6.0 mm above the substrate, opening r 2.05 mm), lumped 50 ohm port, 10.54 GHz,
radiation box padded 15 mm laterally / 30 mm above (the Phase 29 builders' padding; 15 mm = 0.53 lambda).  Cell: full-fill super-ellipse n = 4,
major 3.05 mm (gap 0.15 mm), slot along the cell axis; the OUTLINE is fixed and only the slot turns (the tensor rotation model holds to
1-3 % for this construction, not for turning the whole cell: `results/phase30/refine/rotation_check.json`).
Files: `layout.csv` (x, y, slot L, W, angle, alpha, beta, M_t per cell), `predictions.json`.

## Predictions (written before any HFSS run)
| Quantity | ACTIVE | CONTROL |
|---|---|---|
| Surface-wave power reaching the edge | **40 %** (first-order alpha, table max 13.1 Np/m) ... **58 %** (alpha x 0.53 = the Floquet/first-order ratio of the strongest cell) | ~100 % (alpha_diel + alpha_cond ~0.2 + 0.4 Np/m over 45 mm: ~6 % absorbed, estimated from the E1 cell losses) |
| Radiating-aperture directivity, RHCP, broadside | 19.7 dBi (aperture current only; flat-aperture limit 19.9 dBi: treat as optimistic, ideal model) | NA |
| Gain upper bound (launch fraction 1) | 17.4 dBi (edge 40 %) ... 15.9 dBi (edge 58 %) | NA |
| AR at peak / RHCP SLL outside 24 deg / LHCP-RHCP | 1.1 dB / -22 dB / -24 dB (aperture current only) | NA |
| Radiated power in theta > 60 deg from the aperture current | 0.05 % | NA |
| Monopole direct radiation, edge diffraction, launch fraction | **NA** (not modelled) | **NA** |

## Decision rule (read from the HFSS far field: integrate the radiated power over bands; Gain = realized RHCP)
1. P60 = fraction of radiated power at theta > 60 deg, for ACTIVE and CONTROL.
2. If P60(CONTROL) >= 0.8 x P60(ACTIVE): the near-surface power is NOT a leakage-cell problem (monopole direct radiation and/or edge
   diffraction); a better cell cannot fix it - change the launcher (shorter/loaded pin, annular launcher) or add an edge treatment.
3. If P60(ACTIVE) - P60(CONTROL) is large and the active board's RHCP gain is within 3 dB of the 15.9-17.4 dBi band: the model's
   leakage physics is right; the remaining loss is the edge fraction -> larger board or edge taper.
4. If ACTIVE gain < 13 dBi and P60(CONTROL) is small: the cell leakage/phase model is wrong (alpha too high, or the per-cell tensor
   rotation fails) - go back to the eigenmode results (`../phase30_cell_eigenmode/`).
5. Also record |S11| at the feed (launch matching) for the three runs: a change of S11 between ACTIVE and CONTROL = the cells
   load the feed; the accepted-power normalisation must be used for every gain.
Epsilon_r caution: run ACTIVE at eps_r 10.2 and, as a sensitivity run, 11.2 (Rogers 'design Dk') - the model's X_bar/beta moves by
+45 ohm / +0.27 k0 for +1 in eps_r (`results/phase30/refine/eps11p2_sensitivity.json`), so a measured eps_r is needed before a build.

## STEP package (eps_r 11.2, 10.05 GHz; supersedes the AEDT scripts, removed because they did not run on the user machine)
Manual setup, predictions and the table to copy back: `HFSS_TESTBOARD_SETUP.md`. Port table: `ports_and_excitations.csv`. Geometry check: `step_check.json`. Generators: `scripts/research/phase30_testboard.py` (design; env P30_EPS=11.2 P30_F0=10.05e9; table results/phase30/library/uniform_table.json) and `phase30_testboard_step.py`. STEP files are gitignored; copies in C:\hmsa_sync\data\phase30_testboard\.
