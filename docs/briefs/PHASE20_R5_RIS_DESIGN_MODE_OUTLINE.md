# Phase 20 R5 — RIS design mode: OUTLINE (not yet a brief; for the user's review, nothing started)

Branch: `phase20-ris` (or a new `phase20-r5` from it). Written 2026-10-23. **Outline only:** items, inputs, compute split and open questions; the full brief follows after the user's answers.

## 1. What R5 is
A design mode on the Phase 27 cell-family registry that turns the R4 CP-2 state library (+ the HFSS-calibrated delta) into a **2-bit coding map** for a finite RIS aperture, predicts the radiated pattern with honest labels, checks the coding at array level
(supercell), and exports a board for HFSS. Cell: the retained baseline CP-2 (R3b closed), registry family `ris_cp2_octagon_madp_rt5880_15` (CP) and `ris_lp_1bit_from_cp2_madp_rt5880_15` (LP interim, 1-bit codes of the same cell).

## 2. Decisions already taken (inputs, not to be re-opened)
* Operating frequency **10.5 GHz**; pitch 15 mm (lambda/pitch = 1.90). Library domain f 8-12 GHz, theta <= 45 (theta > 30 EXTRAPOLATED), phi 30 deg grid; higher Floquet orders refused (Wood guard 3 %); outside -> `NOT_AVAILABLE_FROM_CURRENT_DATA`.
* **Steering: +-50 deg guaranteed over the band** (grating-free to ~11.3 GHz at 15 mm); **60 deg only near band centre (<= 10.54 GHz)**; 2-bit quantisation lobes ~ -9 dB for a plain gradient (`results/phase20_r3b/grating_limits.json`).
* The **actual** state responses are used: HFSS steps 98/92/102/68 deg (not 90), theta-30 opposite-handed -11.9 dB (HFSS), amplitudes -0.6 dB at 10.5 GHz; the opposite-handed term carries no PB phase and goes to the specular lobe.
* Isolated-cell predictions are reliable for regular gradients, **not for random codings** at 15 mm pitch (R3: order levels off 1-4 dB, up to ~10 dB): R5 verifies codings at array level (supercell / coupling-aware).
* Labels: MODEL_ONLY; `HFSS_CALIBRATED` only for the baseline cell inside the HFSS-measured domain (theta 0/30, phi 0, run 1/2 data); everything else NOT_HFSS_VALIDATED. PROVISIONAL_DATASHEET diodes, VIA_OMITTED_APPROXIMATION.
* **Mac = coordinator only**: heavy compute (supercell verification, library extensions, large-array runs) are server jobs (`jobs/<job>.md`, OPENBLAS_NUM_THREADS=1 per worker, provenance on every result, checkpointed). The design logic itself must stay light enough for the Mac (O(cells x directions), seconds).

## 3. Items (proposed)
* **R5-0 audit and reuse plan** (Mac, small): what of `src/app/routers/reflectarray.py`, `src/pattern/cellwise_pattern.py` / `vector_aperture_pattern.py` (Phase 23: far field linear in each cell's contribution), `src/hologram/coded_synthesis.py` (Phase 15 E coded method), `auto_states.py`,
  `tensor_cp_service.py`/`tensor_cp_compare.py` (Phase 27 pattern: design -> export -> HFSS comparison) is reused unmodified; module layout `src/ris/design/`.
* **R5-1 calibrated library layer**: `StateLibrary` + `DeltaTable` as one object with a domain/label contract (HFSS_CALIBRATED inside theta <= 30, phi = 0 data; transferred uncorrected to phi != 0 and flagged; uncertainty = library + residual 12 deg near 10.3-10.8 GHz); diode-tolerance HFSS cases
  (the user's C_off = 0.0225/0.025/0.030 pF run) become an optional sensitivity band on the codings.
* **R5-2 illumination models**: plane wave (theta_i, phi_i, handedness) first; offset horn (feed pattern cos^q, F 240-300 mm, offset ~25 deg, edge taper, spillover) second. Every cell is evaluated at its **own local incidence angle and polarisation** (lookup (f, theta_loc, phi_loc)); cells whose local angle leaves the library domain are reported, not extrapolated silently.
* **R5-3 coding synthesis**: required reflection phase per cell (steering and, for the horn, path compensation) -> choose, per cell, the state whose **complex** response (amplitude and actual phase at that local angle) best matches the target; iterative refinement of the map on the true pattern objective (cheap because the far field is linear in the cells); quantisation-lobe suppression
  (random-offset / dithered phase, period optimisation); beam table over 0-50 deg (guaranteed) and up to 60 deg near band centre; frequency sweep 10.0-11.0 GHz; LP interim as 1-bit sign codes.
* **R5-4 forward model and report**: aperture field = illumination x per-cell complex response (co- and opposite-handed parts) -> far field; gain, pointing error, SLL, quantisation lobes, scan loss, CP axial ratio / XPD, **specular opposite-handed lobe level**, bandwidth; predicted-vs-ideal; honest labels on every number.
* **R5-5 array-level verification (R6 pulled in; server)**: supercell MoM (`src/ris/supercell.py`, F = 2n) of windows of the final coding (e.g. 4x4 / 6x6 tiles, gradient and the coded region's worst neighbourhoods), per-order level/phase errors against the isolated-cell prediction (gate 1 dB / 15 deg, reported never tuned), correction or rejection of bad neighbourhoods. Cost ~hours per tile at F = 2n: server jobs only.
* **R5-6 HFSS export and comparison**: a coded board (reduced size first, e.g. 8x8-10x10 cells) as AEDT script with per-cell lumped RLC states (Phase 26/27 export pattern), predictions committed first, comparison script, user runs HFSS. Full 20x20 = 1600 diode boundaries may not be practical in HFSS: size to be agreed.
* **R5-7 app page**: "RIS designer" on the registry (family, frequency, aperture, steering, illumination, bits) -> job -> results/labels -> export; heavy checks submitted as server jobs, results fetched via `hub`. Requests outside the domain return NOT_AVAILABLE_FROM_CURRENT_DATA.
* **R5-8 close-out**: STATUS, CHANGELOG, roadmap log, validation index entry; stop and report.

## 4. Heavy compute (server jobs) vs Mac
Server: supercell tiles (R5-5), any library extension (finer mesh-correction grid near 10.3-10.8 GHz, theta 15/45, phi refinement; each ~12 h scale as in R4), HFSS-board prediction runs, full-aperture pattern sweeps if they grow beyond seconds. Mac: design logic, small forward model, tests on targeted files, git/briefs/docs, review of server results.

## 5. Open questions for the user (answers shape the brief)
1. Aperture and cell count (300 x 300 mm = 20 x 20 cells at 15 mm?) and the HFSS board size we can actually solve.
2. First illumination: plane wave only, or the offset horn from the start (25 deg, F 240-300 mm, feed pattern/taper and handedness of the feed)?
3. Beam goals: single beam steering 0-50 deg (60 near centre), multi-beam, focusing? Any sidelobe/quantisation-lobe specification?
4. CP only first, with the LP 1-bit interim as a second mode, or both in the first release?
5. Calibration scope: is the existing HFSS data (theta 0/30, phi 0, 8-12 GHz) enough, or should a further HFSS run (phi = 45, theta 15/45, or a small coded supercell) be scheduled before R5-3?
6. Robustness to diode tolerance: nominal codings with a reported sensitivity band (proposed), or optimise the map against tolerance corners (the R3b lesson: the theta-30 / phi-45 corners are the MoM's least trustworthy condition)?
7. Bandwidth target around 10.5 GHz (the HFSS pass band is narrower than the library's: |r_RR| -1.5 dB at 10.0 and -1.7 dB at 11.0 GHz).
8. Is the HFSS export of a coded board part of R5 or a follow-on (R5b)?

## 6. Main risks
Actual ladder distortion (98/92/102/68 deg) and theta-30 opposite-handed -11.9 dB limit pattern quality; theta > 30 local angles (offset horn) are extrapolated; phi != 0 has no HFSS calibration; random-like codings fail the isolated-cell model (R3: 2 of 10); the library mesh-correction residual (12 deg) is comparable to the effects being designed;
steering to 60 deg has no band margin above 10.54 GHz; supercell verification at F = 2n is expensive.
