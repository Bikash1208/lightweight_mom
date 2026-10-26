# Project Roadmap: `lightweight_mom`

Authoritative roadmap and project track from the project lead, paired with verified repository state.

---

## 1. Planned Track (Authoritative)

1. **Deliverables:** (a) a monopole-fed coded (PIN-diode reconfigurable) holographic metasurface antenna at X-band with integrated monopole feed and FPGA beam steering (the proposal's "digital coded metasurface antenna"); (b) the HMSA Designer app with the tensor-CP flow and coded holographic design mode; (c) the MoM solver/tool itself, to be submitted later. *Dropped/parked:* reflectarray / plane-wave RIS designs with external feed, and transmitarrays.
2. **Done:** HMSA Phases 8–27 (app, multiport/coded, hexagon SLL optimiser, tensor-CP hexagon HFSS-validated at 10.5 GHz, 23.03 dBi RHCP); RIS R1 (loads, oblique), R2 (diode models), R3 (CP-2 cell; HFSS-checked at 10.5 GHz), R4 (state library; v2 full-grid library built on the server and passing the held-out gate, 2026-10-24); R5-1 to R5-3 (tolerance check, plane-wave forward model, coding synthesis); R5-5 (RIS designer page in the app); R5-6 infrastructure (RIS board exporter; 3×3 smoke package verified in AEDT).
3. **In progress:** Phase 28 (monopole-fed coded holographic metasurface antenna), branch `phase28-coded-holographic-v2`: **item 0 (re-plan) done**; item 1 (literature, verified DOIs only) and item 2 (feasibility with the Phase 26 forward model: pitch, bit depth, required contrast in physical units, diode / shift-register / bias-line counts, stop point ■) in progress. The first attempt (`review/gemini-*` on the hub) is kept as a record only: its item 0 was reused, its items 1-5 were reviewed and not used.
4. **Next, in order:**
   - Item 2 stop point approval (pitch, bit depth, polarization scope, diode count).
   - Item 3: Diode-loaded subwavelength cell candidates (dispersion/eigenmode with diode loads in each state, ±20% diode tolerance, HFSS eigenmode check package).
   - Item 4: State library on the server, held-out validation.
   - Item 5: Coded hologram synthesis per beam (codebook), forward-model predictions, app mode, HFSS board package, bias-layout sketch.
   - Fabrication (confirm flip-chip assembly with the lab); measurement.
5. **Open decisions:** 10 vs 10.5 GHz with the PI/ISRO (both kept in feasibility sweep; 10.5 GHz baseline); 1-bit vs 2-bit; polarization scope (LP vs CP); board A HFSS run (parked).
6. **Parked:** Plane-wave RIS R5 (R5-1/2/3/5/6 done, R5-4 stopped; `docs/briefs/_parked/PHASE20_R5_RIS_DESIGN_MODE.md`), CP-2 cell (MODEL_ONLY in registry); server jobs `r5-supercell-gradients`, `r5-supercell-random`, `r4-band-refine` (withdrawn); Butler sector scanner (Phase 25), unit-cell S-parameter tool (Phase 16), Phase 19 nylon lens E/F.
7. **Design Changes Log:**
   - **2026-10-24 (Phase 28 Decision):** Pivot to **monopole-fed coded holographic metasurface antenna** with integrated monopole feed and PIN-diode switched cells (digital coded metasurface antenna, FPGA-controlled). Reflectarray / plane-wave RIS designs with external feed (Phase 20 R5) and transmitarrays (Phase 19) are dropped. Carried over: Phase 26 monopole launch model (launch fraction 0.62 calibrated at 6.0 mm pin), leaky-wave forward model, tensor-CP synthesis, single-hexagon HFSS validation (150 mm hexagon, 6.4 mm lattice), Phase 20 R1/R2 diode circuit models and HFSS lumped-RLC representation. R5 parked; server jobs `r5-supercell-gradients`, `r5-supercell-random`, `r4-band-refine` withdrawn.
8. **Working rules:** briefs drive work; one commit per item; fast suite green; no silent physics changes; honest labels (`MODEL_ONLY`, `HFSS_CALIBRATED`, `NOT_AVAILABLE_FROM_CURRENT_DATA`); heavy compute only on the server via `jobs/` and `compute/<job>` branches; the Mac merges.

---

## 2. Verified Current State

Every item below cites ground truth verified on `phase20-ris` / `main` as of 2026-10-24.

| Track Item | Stated Status | Verified Evidence | Commit / Hash |
|---|---|---|---|
| **Phase 28: Monopole-fed Coded Holographic Antenna** | In progress (item 0 done; items 1-2 in progress; nothing validated) | Brief `docs/briefs/PHASE28_CODED_HOLOGRAPHIC_ANTENNA.md`; branch `phase28-coded-holographic-v2`. | `fa22860` (brief) |
| **Deliverables (a): Coded Holographic Antenna** | In progress (feasibility study not yet done) | Integrated monopole launch + leaky-wave forward model (`src/hologram/leaky_forward.py`); diode models in `src/ris/diodes.py`. | Phase 28 active |
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
