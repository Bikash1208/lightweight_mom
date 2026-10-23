# Phase 20 R5 — RIS design mode on the R4 CP-2 library (plane wave first, then offset horn)

**Branch:** `phase20-ris` (R5 work may sit on `phase20-r5` created from it; heavy compute on `compute/<job>` branches per the infra brief). **Written:** 2026-10-23, from the user's answers to the R5 outline.
**Status:** brief only. Nothing in items 1-9 is started; **stop before item 1** and wait for the user. Item 0 (audit) was done while writing this brief (Appendix A).
Labels throughout: `MODEL_ONLY`, `HFSS_CALIBRATED` (baseline cell, inside the HFSS-measured domain only), `NOT_HFSS_VALIDATED`, `PROVISIONAL_DATASHEET`, `VIA_OMITTED_APPROXIMATION`, `EXTRAPOLATED_BEYOND_R3_TARGETS` (theta > 30), `NOT_AVAILABLE_FROM_CURRENT_DATA`.

## 0. Decisions (user, 2026-10-23) and inherited constraints
1. **Aperture:** 20 x 20 cells (300 mm at the 15 mm pitch). **HFSS check on a 10 x 10 sub-board (150 mm) first; the full board only if it agrees.**
2. **Illumination:** plane wave first, then the offset horn, both inside R5.
3. **Beam goal:** a single CP beam steered 0-50 deg in the phi = 0 and phi = 90 planes; report gain, SLL and the opposite-handed specular lobe. Multi-beam later (not R5).
4. **Polarisation:** CP first. The LP 1-bit interim mode is later (not R5). CP handedness: the incident wave is R in the propagation-relative basis; the wanted reflected wave is R (`r_RR`); `docs/theory/CONVENTIONS.md`.
5. **Calibration:** the current HFSS data (run 1 and 2: theta 0/30, phi 0, 8-12 GHz; `results/phase20_r3b/hfss_delta_table*.npz`) is enough to start. The user's coming **diode-tolerance cases** arrive as the same CSV format as before plus a **C_off column**. phi = 45 HFSS only if R5 shows sensitivity to it.
6. **Tolerance:** nominal codings with a **reported sensitivity band** (no tolerance-corner optimisation: the R3b lesson).
7. **Bandwidth:** report the **1 dB and 3 dB gain bandwidth**; starting goal ~10.3-10.7 GHz, not locked.
8. **HFSS board export is R5's last item**, predictions first.
9. Inherited: operating frequency **10.5 GHz**; cell = baseline CP-2 (octagon 4.601 mm, stubs 1.748 x 1.375 mm, pitch 15 mm, MADP-000907-14020), HFSS steps 98/92/102/68 deg, theta-30 opposite-handed -11.9 dB; **steering +-50 deg guaranteed over the band (grating-free to ~11.3 GHz), 60 deg only near band centre (<= 10.54 GHz)**;
   isolated-cell predictions are reliable for regular gradients, **not for random codings** (R3: order levels off 1-4 dB, up to ~10 dB): codings are verified at array level; library domain f 8-12 GHz, theta <= 45 (theta > 30 EXTRAPOLATED), phi 30 deg grid, higher-order refusal with the 3 % Wood guard; library mesh-correction residual 12 deg near 10.3-10.8 GHz.
10. **Compute rules (infra brief):** the Mac is coordinator-only. Anything above about a minute is a server job `jobs/<job>.md` (base commit, exact commands, workers x `OPENBLAS_NUM_THREADS=1`, outputs, stop condition, provenance record, manifest, checkpointing, `compute/<job>` branch). The Mac runs the design logic (seconds), targeted tests and the snapshot gate; the full fast suite runs on the server as the first step of the next job.
11. Repo rules: one commit per item, additive only, `make test-fast` green (on the server when the Mac may not run it; say so in the commit), bit-identical for every pre-existing result, honest labels, no fabricated HFSS data, no ML in the solver core. **Do not:** change CP-2 geometry or diode values, fit a surrogate and call it physics, build multi-beam or LP modes, call anything HFSS-validated without a user HFSS run.

## 1. Items

### R5-1 Calibrated library layer (`src/ris/calibrated_library.py`)
One object combining `StateLibrary` (corrected table) with the HFSS `DeltaTable`: `response(f, theta, phi, state)` returns the structure-frame 2x2 circular response, a **label** (`HFSS_CALIBRATED` inside theta <= 30, phi = 0 and 8-12 GHz; `HFSS_CALIBRATED_TRANSFERRED_UNCORRECTED_PHI` for phi != 0 (delta applied as a phi-independent structure-frame offset, flagged); `EXTRAPOLATED_BEYOND_R3_TARGETS` for theta > 30; `OutOfDomain` outside the grid / where a higher order propagates) and an uncertainty (library uncertainty + the 12 deg residual flag near 10.3-10.8 GHz).
* Delta from run 2 (0.05 GHz steps, plane correction) is the default; run 1 kept for comparison.
* **Tolerance data reader** for the user's C_off CSV (`read_hfss_reim` + the `Coff` column): per-state complex response vs C_off (0.0225 / 0.025 / 0.030 pF) at theta 0/30 -> a `ToleranceBand` object (interpolation of each state's response in C_off; MoM corner table `results/phase20_r4/tolerance_corners.json` as the fallback band, labelled MODEL_ONLY). Comparison script tolerance-HFSS vs MoM (prediction first: server job `r5-tolerance-mom`, 3 C_off x 4 states x theta 0/30 x 10.0-11.0 GHz at mesh 60/90, F = 2n).
* Acceptance: reproduces the HFSS points exactly at HFSS frequencies; labels and refusals tested; no change to `StateLibrary`.

### R5-2 Plane-wave forward model (`src/ris/forward.py`)
Reflected field of a finite coded aperture under plane-wave CP illumination (theta_i, phi_i; default normal incidence, optional 30 deg inside the calibrated domain).
* Cell grid 20 x 20 (parameterised); every cell has the same local incidence for a plane wave -> one library lookup per state at (f, theta_i, phi_i); full 2x2 response per cell (co-handed and opposite-handed parts).
* Aperture field = incident phase x cell response; **cell-sampled aperture** (piecewise-constant 15 mm cells: element factor sinc x array factor, so grating/quantisation lobes and the visible-region edge are inherent); obliquity factor as in `pattern/aperture_pattern.py` (documented MODEL choice; Appendix A) with a flag to switch it off.
* Outputs: far field over theta/phi (co-handed `E_R`, opposite-handed `E_L`), **gain** G = 4 pi U / P_inc with P_inc = incident power density x aperture area (includes diode loss, opposite-handed power and quantisation), directivity, aperture efficiency vs 4 pi A / lambda^2, pointing error, **SLL** (co-pol, in the phi = 0 and 90 cuts and the 2-D maximum), strongest quantisation / grating lobe, **opposite-handed specular lobe level re the main beam** and its angle, axial ratio on the beam.
* **Validation (tests, Mac):** uniform aperture = analytic sinc; continuous ideal gradient steers to the analytic angle with D = 4 pi A / lambda^2 (efficiency 1 within the obliquity convention); 2-bit ideal-ladder quantisation loss = 0.91 dB (`coded_synthesis.textbook_quantization_loss_db`); energy bookkeeping (|r_RR|^2 + |r_LR|^2 <= 1 per cell); **periodic gradient order amplitudes equal the isolated-cell prediction used in R3** (`results/phase20_r3/locality/`) to numerical precision; reciprocity of the pattern under phi -> phi + 180 for a symmetric coding.

### R5-3 Coding synthesis, plane wave (`src/ris/coding.py`)
* Required phase per cell for steering to (theta_s, phi_s) (phi_s = 0 and 90) from the given incidence; **state choice by the actual complex responses** (nearest in the complex plane including amplitude and the HFSS-measured step distortion), not the ideal 0/90/180/270 ladder.
* Refinement on the true pattern objective (the far field is linear in the cells: single-cell flips update O(directions)): maximise the gain at the target direction subject to the quantisation-lobe / SLL / opposite-handed constraints; deterministic seeds recorded.
* **Quantisation-lobe handling:** plain gradient vs random-offset / dithered phase vs period optimisation; report each (the plain 2-bit gradient has the strongest quantisation lobe ~ -9 dB).
* **Beam table:** theta_s = 0, 5, ..., 50 deg in phi_s = 0 and 90 at 10.5 GHz: gain, SLL, quantisation lobe, specular opposite-handed lobe, pointing error, scan loss; a 60 deg row near band centre marked as such (no band margin above 10.54 GHz).
* **Frequency / bandwidth:** fixed coding, f = 10.0-11.0 GHz (step 0.05): gain at the design direction and peak gain (beam squint), **1 dB and 3 dB gain bandwidth** (starting goal 10.3-10.7 GHz, not locked), the Floquet margin per frequency (+-50 deg guaranteed region).
* **Sensitivity band (item 6 of the answers):** the same coding evaluated with the C_off band (R5-1); report gain / pointing / SLL spread. No optimisation against corners.
* Report after this item (**stop point 1**): beam table, bandwidth, sensitivity band, comparison of the 2-bit design against the ideal and the quantisation theory.

### R5-4 Offset-horn illumination (`src/ris/horn.py`)
* Feed model: cos^q E-plane/H-plane pattern from an edge-taper input (default -10 dB at the aperture edge), phase centre, position (focal distance F and offset angle), CP handedness; spherical-wave incident phase and amplitude per cell; **local incidence angle and polarisation per cell**: the horn field is decomposed into the cell's TE/TM frame and the library's TE/TM 2x2 response (`StateLibrary.te_tm`) is applied; the reflected field is rebuilt in the circular basis.
* **Feasibility gate first (computed in this brief, to be implemented):** for the 20 x 20 aperture the local incidence range with the horn aimed at the centre is
  F = 240 mm: max theta_loc 40 / 47 / 50 / 52 deg for offset 0 / 15 / 20 / 25 deg; F = 270 mm: 37 / 44 / 47 / 50; F = 300 mm: 34 / 42 / 45 / 48 (fraction of cells above 30 deg: 6-59 %, above 45 deg: 0-17 %).
  Most of the aperture is therefore in the **EXTRAPOLATED** (theta > 30) region and the corners beyond 45 deg are outside the library for the usual 25 deg offset. The gate reports the coverage map and the allowed (F, offset) set for (i) fully inside the library (<= 45 deg), (ii) inside the HFSS-calibrated region (<= 30 deg, essentially a small offset only); **stop point 2:** the user decides between a smaller offset / larger F, a library extension to 50 deg (server job, ~12 h scale as R4), or accepting EXTRAPOLATED labels.
* Synthesis: the item-3 method with per-cell local responses and the horn path phase (focusing + steering); gain, spillover, illumination (taper) efficiency, SLL; same metrics and bandwidth as item 3.
* Validation: horn at large F reproduces the plane-wave result; spherical-wave path phase vs analytic; polarisation decomposition tests.

### R5-5 Array-level verification (server jobs; the R6 coupling check pulled in)
* Supercell MoM (`src/ris/supercell.py`, **F = 2n**) of (a) the **periodic gradient windows** matching the designed codings (n x 1 supercells, n <= 8, best rational approximants of the designed period) and (b) random-offset / dithered tiles (6 x 6 first) at theta = 0 and the plane-wave incidence; the order amplitudes and phases are compared with the isolated-cell prediction used by the design (gate 1 dB / 15 deg on the strongest orders, **reported, never tuned**). Costs (R3 at the minimum cutoff: gradient ~45 s, random 4 x 4 ~30 min per case; measured at F = 2n: the 4 x 1 gradient took 19-28 min per angle under load, i.e. ~25x; a 4 x 4 random case is estimated at several hours, a 6 x 6 tile more): server jobs `r5-supercell-gradients`, `r5-supercell-random`.
* Outcome per coding: `COUPLING_VERIFIED` (all gates pass), `COUPLING_CORRECTED` (design re-run with the supercell-measured order errors as a correction map; stated limits) or `COUPLING_FAILED_ISOLATED_CELL_MODEL` (reported as an error band on the pattern); no coding is called verified without these runs.
* Report (**stop point 3**): which codings pass, the error band on gain / SLL / specular lobe for those that do not.

### R5-6 RIS designer page (app)
"RIS designer" on the Phase 27 registry: family `ris_cp2_octagon_madp_rt5880_15`, frequency, aperture, steering (theta_s, phi_s), illumination (plane wave / horn parameters), results with labels and refusals (`NOT_AVAILABLE_FROM_CURRENT_DATA` outside the registry domain), sensitivity band, beam table, coupling status; heavy checks submitted as server jobs (the page shows the job request and ingests `hub` results; no heavy compute in the Mac app). Tests in the style of `tests/test_app_phase27_tensor_cp.py`.

### R5-7 HFSS export of a coded board - LAST ITEM, predictions first
* **10 x 10 sub-board (150 mm)** coded for a chosen steering (default 30 deg in phi = 0; the user may pick) with plane-wave incidence; one lumped RLC boundary per diode (400), states from the coding; AEDT script generator on the Phase 26/27 export pattern (`hfss_build_board.py`, layout CSV, README, `board.step` gitignored) - feasibility of the HFSS model size is checked with the user before generation.
* **Predictions committed first** (forward model: main-beam gain, pointing, SLL, specular opposite-handed lobe, bistatic pattern cuts; uncertainty = library + sensitivity band), a comparison script (`phase20_r5_hfss_compare.py`, comparison only), then the user runs HFSS. **Only if the 10 x 10 agrees** is the full 20 x 20 board generated (a follow-up decision, not part of this item).

### R5-8 Close-out
STATUS (capability row, open items), CHANGELOG, roadmap design-changes log, validation index (predictions and, when it exists, the HFSS comparison), registry note, archive this brief to `docs/history/briefs/`, stop and report (**stop point 4**).

## 2. Server jobs expected (each a `jobs/<job>.md`, all outputs with provenance and a manifest; Mac reviews and merges)
`r5-tolerance-mom` (R5-1), `r5-library-extension` (only if the user chooses it at stop point 2), `r5-supercell-gradients`, `r5-supercell-random` (R5-5), `r5-board-predictions` (R5-7, if the sub-board pattern grows beyond a minute), plus the fast suite as the first step of each. The equivalence check (`jobs/equivalence-check.md`) must have passed before any server result is used.

## 3. Acceptance and honesty rules
* Every number carries its label; HFSS_CALIBRATED only where the HFSS data exist; the opposite-handed term is always reported; extrapolated local angles never silently used.
* Predictions are written and committed before the user's HFSS run; comparison scripts report differences without calibrating.
* Seeds, grids and definitions of gain, SLL and bandwidth are written in the code and in the report; no result is reported without the check that produced it.

## 4. Main risks
HFSS ladder distortion (98/92/102/68 deg) and theta-30 opposite-handed -11.9 dB bound the pattern quality; the offset horn puts most cells at theta_loc > 30 deg (item 4 gate); phi != 0 has no HFSS calibration (phi = 45 run only if the sensitivity shows it matters); random-like codings may fail the isolated-cell model; the library residual (12 deg near 10.3-10.8 GHz) is of the order of the effects designed; 60 deg steering has no band margin above 10.54 GHz; supercell verification at F = 2n is expensive; a 400-boundary HFSS model may be impractical.

## Appendix A - item 0 audit (reuse inventory, done 2026-10-23)
* **Reuse unmodified:** `src/ris/state_library.py` (`StateLibrary.structure_circular / te_tm / guard / uncertainty`), `src/ris/hfss_calibration.py` (`read_hfss_reim`, `fit_reference_plane_mm`, `DeltaTable`), `src/ris/supercell.py` (`build_supercell`, `solve_orders`, `default_mesh(..., floquet_factor=2)`), `src/ris/polarization.py` (`circular_reflection`, `to_structure_frame`), `src/hologram/cell_family_registry.py` (`check_request`, RIS families),
  `src/phase8e/beam_reconstruction.far_field_pattern` (scalar aperture integral) / `src/pattern/aperture_pattern.compute_pattern` (cos(theta) obliquity convention, peak detection, directivity estimate) as the numerical baseline for the new vector CP forward model, `src/hologram/coded_synthesis.textbook_quantization_loss_db`, the metric helpers in `src/pattern/vector_aperture_pattern.py` (`circular_components`, `copol_component`, `strongest_secondary_beam`, `beam_local_peak`).
* **Not applicable to a plane-wave RIS:** `vector_aperture_pattern.compute_vector_pattern` / `cellwise_pattern.cell_contributions` (radial-current, feed-centred leaky-wave hologram model), `tensor_leaky_cp`, `leaky_forward`; `coded_synthesis` phase maps are leaky-wave (`beta_sw`) maps - only its quantisation helpers are reused. `src/app/routers/reflectarray.py` (library-cell selection with continuous phases) is a pattern for the page, not for the 4-state coding.
* **Pattern to copy:** `tensor_cp_service.py` / `tensor_cp_compare.py` / `src/export/tensor_cp_export.py` (design -> export -> HFSS comparison, request normalisation, jobs, honest refusals).
* **New code (all additive, `src/ris/`):** `calibrated_library.py`, `forward.py`, `coding.py`, `horn.py`; scripts under `scripts/research/phase20_r5_*.py`; tests with each.
