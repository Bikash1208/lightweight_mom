# Project Roadmap: `lightweight_mom`

Authoritative roadmap and project track from the project lead, paired with verified repository state.

---

## PARKED 2026-10-27 (Phase 28, branch `phase28-coded-holographic-v2`; the user switched projects, resume later)

**Architecture decided (design level, MODEL_ONLY except where stated):** side/sector-fed surface-wave hologram with a Yagi-type launcher (driven pin + reflector post + directors), SCALAR 1-bit, LINEAR polarisation, ONE MA4AGP907 (= MADP-000907-14020)
PIN diode per cell, individual control line per cell, 9 mm pitch, RT/duroid 5880 3.175 mm, 10.0 GHz, 150 mm hexagon (724 cells, 724 diodes). Diode model `MEASURED_VENDOR_SPARAM` (`src/ris/diodes_measured.py`, extraction
`scripts/research/phase28_diode_extract.py`; R 3.3 / 2.4 / 2.0 ohm at 5 / 10 / 15 mA, C 28.2 fF at 0 V, no series L visible). Lead cell: loop cell R4 (6.6 mm loop, 1.2 mm strip, diode in the middle of one side), **HFSS_CHECKED at ROT 0, normal incidence,
9 - 11 GHz** (one run: all sheet-admittance components within 4.3 % of the MoM, Delta B 4.22 mS against 4.05 mS; `hfss_results/phase28_loop_R4/NOTES.md`; calibration hook `src/ris/loop_cell_calibration.py`). Cells stay UNROTATED (the diode axis along x fixes a clean
x polarisation; model: single beam, SLL -14 to -21 dB, 21 - 23 dBi with a placeholder launch fraction; 1-bit sampling limits the scan toward the feed to about 30 deg). Earlier cell tracks (G, J, two-diode searches) are superseded; their
files stay as history.

**User note (2026-10-27 HFSS run):** the state files were swapped in the project: `OFF.csv` carries the R = 2.4 ohm state (lossy, isotropic) and `ON.csv` the C = 28.2 fF state (open, anisotropic). Which of the files or the boundaries was
swapped is still TO BE CONFIRMED by the user. (The first HFSS build used hand-set boundaries; `build_cell.py` is fixed since 9ea6897 and tested against an AEDT recording stub, but not yet re-run in AEDT.)

**State of the open items, and the next steps in order:**
1. **Rotated-load solver bug** (`transverse_conducts=False` on a rotated load gives ON = OFF). Cause found: the documented R3 limit (the 1e12 ohm transverse penalty locks the Cartesian rooftops on a rotated region). Done: a `UserWarning`
   in `check_load_geometry` (results bit-identical, fast suite 2900 passed), a minimal reproduction and a rotation-covariance physics test (conducting rotated loads: 2 - 7 % at 30 / 45 deg). NOT done: a real fix of blocked rotated loads
   (would need a basis aligned with the load or a projected constraint); not needed for the unrotated lead cell. Rotated-cell MoM results are good to 10 - 17 % only (staircase of the rotated geometry).
2. **Launcher design + HFSS package.** Model (coupled grounded posts through the TM0 wave, two feed calibrations) and package `hfss_exports/phase28_launcher/` exist with predictions (width 77 - 79 deg, S11 -8.7 / -15.9 dB, F/B 5 - 8 dB, predictions only).
   NOT done: any HFSS run; the model neglects near-field/space-wave mutual impedances and the loaded-surface edge. Next: user's HFSS run (bare slab first, then OFF-state loop surface), then re-tune.
3. **Phase 28 item 5:** app design mode (side-fed scalar 1-bit hologram synthesis, codebook per beam, forward-model predictions) and the full-board HFSS export (150 mm hexagon, 724 cells; smaller board first).
4. **Full-board HFSS check** (user): beam angle, gain, SLL, cross-pol against the forward model.
5. **Oblique / rotated-cell HFSS checks:** polarisation turned by 20 / 40 deg in the Floquet cell (the directionality check of `hfss_exports/phase28_loop_R4/`), the rotated cell, the measured diode reactance, package inductance.

**Scope:** design and solver/tool work verified in HFSS. NO fabrication, bias-network or FPGA work. Nothing here is merged to main, tagged or deleted; branch state is the last commit on `phase28-coded-holographic-v2`.

---

## 1. Planned Track (Authoritative)

1. **Deliverables:** (a) a monopole-fed coded (PIN-diode reconfigurable) holographic metasurface antenna at X-band with integrated monopole feed and FPGA beam steering (the proposal's "digital coded metasurface antenna"); (b) the HMSA Designer app with the tensor-CP flow and coded holographic design mode; (c) the MoM solver/tool itself, to be submitted later. *Dropped/parked:* reflectarray / plane-wave RIS designs with external feed, and transmitarrays.
2. **Done:** HMSA Phases 8–27 (app, multiport/coded, hexagon SLL optimiser, tensor-CP hexagon HFSS-validated at 10.5 GHz, 23.03 dBi RHCP); RIS R1 (loads, oblique), R2 (diode models), R3 (CP-2 cell; HFSS-checked at 10.5 GHz), R4 (state library; v2 full-grid library built on the server and passing the held-out gate, 2026-10-24); R5-1 to R5-3 (tolerance check, plane-wave forward model, coding synthesis); R5-5 (RIS designer page in the app); R5-6 infrastructure (RIS board exporter; 3×3 smoke package verified in AEDT).
3. **PARKED 2026-10-27 (see the PARKED section above):** Phase 28 (monopole-fed coded holographic metasurface antenna), branch `phase28-coded-holographic-v2`: **items 0 (re-plan), 1 (literature, verified DOIs only) and 2 (feasibility with the Phase 26 forward model: pitch, bit depth, required contrast in physical units, diode / shift-register / bias-line counts) done; at the item-2 stop point ■ awaiting the user's decision.** The first attempt (`review/gemini-*` on the hub) is kept as a record only: its item 0 was reused, its items 1-5 were reviewed and not used.
4. **Next, in order:**
   - Item 2 stop point approval (pitch, bit depth, polarization scope, diode count).
   - Item 3: Diode-loaded subwavelength cell candidates (dispersion/eigenmode with diode loads in each state, ±20% diode tolerance, HFSS eigenmode check package).
   - Item 4: State library on the server, held-out validation.
   - Item 5: Coded hologram synthesis per beam (codebook), forward-model predictions, app mode, HFSS board package, bias-layout sketch.
   - Fabrication (confirm flip-chip assembly with the lab); measurement.
5. **Open decisions:** **Operating frequency decided: 10.0 GHz (PI, 2026-10-26)**; the 10.54 GHz Phase 26 board and tables stay as history, the supercell table and the feasibility sweep are redone at 10.0 GHz; 1-bit vs 2-bit; polarization scope (LP vs CP); board A HFSS run (parked).
6. **Parked:** Plane-wave RIS R5 (R5-1/2/3/5/6 done, R5-4 stopped; `docs/briefs/_parked/PHASE20_R5_RIS_DESIGN_MODE.md`), CP-2 cell (MODEL_ONLY in registry); server jobs `r5-supercell-gradients`, `r5-supercell-random`, `r4-band-refine` (withdrawn); Butler sector scanner (Phase 25), unit-cell S-parameter tool (Phase 16), Phase 19 nylon lens E/F.
7. **Design Changes Log:**
   - **2026-10-24 (Phase 28 Decision):** Pivot to **monopole-fed coded holographic metasurface antenna** with integrated monopole feed and PIN-diode switched cells (digital coded metasurface antenna, FPGA-controlled). Reflectarray / plane-wave RIS designs with external feed (Phase 20 R5) and transmitarrays (Phase 19) are dropped. Carried over: Phase 26 monopole launch model (launch fraction 0.62 calibrated at 6.0 mm pin), leaky-wave forward model, tensor-CP synthesis, single-hexagon HFSS validation (150 mm hexagon, 6.4 mm lattice), Phase 20 R1/R2 diode circuit models and HFSS lumped-RLC representation. R5 parked; server jobs `r5-supercell-gradients`, `r5-supercell-random`, `r4-band-refine` withdrawn.
8. **Working rules:** briefs drive work; one commit per item; fast suite green; no silent physics changes; honest labels (`MODEL_ONLY`, `HFSS_CALIBRATED`, `NOT_AVAILABLE_FROM_CURRENT_DATA`); heavy compute only on the server via `jobs/` and `compute/<job>` branches; the Mac merges.

---

## 2. Verified Current State

Every item below cites ground truth verified on `phase20-ris` / `main` as of 2026-10-24.

| Track Item | Stated Status | Verified Evidence | Commit / Hash |
|---|---|---|---|
| **Phase 28: Monopole-fed Coded Holographic Antenna** | PARKED 2026-10-27 (was: in progress; items 0-2 done at 10.0 GHz, item-3 cell candidates done, at the item-3 stop point; MODEL_ONLY, nothing validated) | Brief `docs/briefs/PHASE28_CODED_HOLOGRAPHIC_ANTENNA.md`; branch `phase28-coded-holographic-v2`. | `dd8d81e` (brief on this branch) |
| **Deliverables (a): Coded Holographic Antenna** | In progress (feasibility study at 10.0 GHz done; converged candidate cells characterised, MODEL_ONLY, none HFSS-validated) | Integrated monopole launch + leaky-wave forward model (`src/hologram/leaky_forward.py`); diode models in `src/ris/diodes.py`. | Phase 28 active |
| **Deliverables (b): HMSA Designer App** | Tensor CP flow complete; RIS designer page built (parked) | App entry points `/tensor_cp.html` (`src/app/routers/tensor_cp.py`) and `/ris.html` (`src/app/routers/ris.py`, `src/ris/service.py`). | `4b87212` (Phase 27 close-out), `26b2433` (R5-5 RIS page) |
| **Deliverables (c): MoM Solver Tool** | in use; frozen core; the solver deliverable is not yet packaged | Periodic layered Green's function, rooftop basis, and impedance assembly in `src/solver.py`, `src/green_function.py`, `src/mom_matrix.py`, `src/sparameters.py`. | Tag `v1.1-handoff` |
| **HMSA Phases 8–27** | Done; single-hexagon tensor-CP board has "HFSS-validated prediction (one run, 10.5 GHz, one board)" (`docs/STATUS.md`) | End-to-end holographic synthesis, multi-port, hexagon SLL optimiser (`src/hologram/sll_optimizer.py`), tensor-CP single hexagon verified in HFSS (`hfss_results/hex150_tensor_cp_single/`). | `4b87212` (Phase 27 close-out), `58b864e` (Phase 26 single-hex close-out), `41f9ea0` (Board B export) |
| **RIS R1: Loads & Oblique** | MODEL_ONLY, NOT_HFSS_VALIDATED (`docs/STATUS.md`) | Oblique plane-wave excitation `compute_s_parameters_oblique` and lumped loads in `src/loads.py`, `tests/test_sparameters_oblique.py`. | `a04c2dd` (R1-1 oblique), `36b4d1a` (R1-2 loads), `3aa6b94` (R1 close-out) |
| **RIS R2: Diode Models** | MODEL_ONLY, PROVISIONAL_DATASHEET (`docs/STATUS.md`) | Skyworks SMP1345-079LF PIN diode models and bias sweeps in `src/ris/diodes.py`, `results/phase20_r2/README.md`. | `26551cb` (R2-1 diode models), `286e889` (R2 close-out) |
| **RIS R3: CP-2 Cell** | HFSS-checked at 10.5 GHz (`docs/STATUS.md`) | Baseline CP-2 octagon cell; HFSS check in `hfss_results/phase20_r3_cp_cell/NOTES.md` (`docs/STATUS.md`: "HFSS-checked at 10.5 GHz"). | `13463f5` (R3 close-out), `ba12f3b` (retained baseline CP-2) |
| **RIS R4: State Library v2** | GATE_PASS (`docs/STATUS.md`) | Full-grid library v2 (8942 points at F = 2n = 120) assembled from server chunks in `results/phase20_r4/cp2_library_v2_fullgrid.npz`. Held-out validation gate PASS (`compute/r4-fullgrid-validate`). | `7ddc843` (merge full grid), `1264bdc` (assemble v2), `37074da` / `f6688c9` (gate PASS) |
| **RIS R5 (R5-1 to R5-3, R5-5, R5-6)** | Parked (direction change 2026-10-24) | Tolerance MoM (`results/phase20_r5/TOLERANCE_GAP.md`), plane-wave forward model (`src/ris/forward.py`), coding synthesis (`src/ris/coding.py`), RIS page (`/ris.html`), board exporter infra (`src/export/ris_board_export.py`, 3×3 AEDT smoke verified). Brief moved to `docs/briefs/_parked/`. | `9d3deee`, `3f2cfe2`, `00528a4`, `26b2433`, `8033f90` |
| **RIS R5-4: Supercell Locality** | Parked / Stopped | Server jobs `jobs/r5-supercell-gradients.md` and `jobs/r5-supercell-random.md` marked WITHDRAWN 2026-10-24. | `7852dc2` |
| **Phase 26 Board A HFSS Run** | Parked | Board A export package in `hfss_exports/hex4_tensor_cp_static/`; simulation run parked per `docs/STATUS.md`. | `20ee0d2` (Board A export) |
| **Parked Ideas** | Parked | RIS R5 (`docs/briefs/_parked/PHASE20_R5_RIS_DESIGN_MODE.md`), Butler sector scanner (`docs/briefs/_parked/PHASE25_BUTLER_SECTOR_SCANNER.md`), FSS S-param tool (`docs/briefs/_parked/PHASE16_UNIT_CELL_SPARAM_TOOL.md`), nylon lens E/F (`docs/briefs/_parked/PHASE19_DIELECTRIC_TRANSMITARRAY_28GHZ.md`). | Commit `51ffbee` (Phase 25), `docs/STATUS.md` |
