# Project Roadmap: `lightweight_mom`

Authoritative roadmap and project track from the project lead, paired with verified repository state.

---

## 1. Planned Track (Authoritative)

1. **Deliverables:** (a) a 2-bit PIN-diode RIS at X-band with LP and CP as user-selectable families; (b) the HMSA Designer app with the tensor-CP flow and the RIS design mode; (c) the solver/tool itself, to be submitted later.
2. **Done:** HMSA Phases 8–27 (app, multiport/coded, hexagon SLL optimiser, tensor-CP hexagon HFSS-validated); RIS R1 (loads, oblique), R2 (diode models), R3 (CP-2 cell; HFSS-checked at 10.5 GHz), R4 (state library; v2 full-grid library built on the server and passing the held-out gate, 2026-10-24); R5-1 to R5-3 (tolerance check, plane-wave forward model, coding synthesis); R5-5 (RIS designer page in the app); R5-6 infrastructure (RIS board exporter; 3×3 smoke package awaiting an AEDT build check).
3. **In progress:** R5-4 supercell locality (server jobs r5-supercell-gradients, r5-supercell-random). After it: the 10×10 HFSS coded board package (R5-6) and the R5 close-out.
4. **Next, in order:** finish R5 (app RIS design mode, plane-wave illumination); HFSS check of a coded RIS board; stage-2 integrated antenna with its feed (decision pending: (a) monopole-fed reconfigurable holographic antenna — lead's preference, (b) reflectarray with a feed in front, (c) microstrip network); bias network in the RF model; fabrication (confirm flip-chip assembly with the lab); measurement.
5. **Open decisions:** 10 vs 10.5 GHz with the PI/ISRO (10.5 GHz accepted for now); 2-bit LP still needs its own cell (1-bit LP from CP-2 is the interim); stage-2 feed type; board A HFSS run (parked).
6. **Parked:** Butler sector scanner (Phase 25), unit-cell S-parameter tool (Phase 16), Phase 19 nylon lens E/F.
7. **Working rules:** briefs drive work; one commit per item; fast suite green; no silent physics changes; honest labels (`MODEL_ONLY`, `HFSS_CALIBRATED`, `NOT_AVAILABLE_FROM_CURRENT_DATA`); heavy compute only on the server via `jobs/` and `compute/<job>` branches; the Mac merges.

---

## 2. Verified Current State

Every item below cites ground truth verified from the live tip of `phase20-ris` (`f6688c9`).

| Track Item | Stated Status | Verified Evidence on `phase20-ris` | Commit / Hash |
|---|---|---|---|
| **Deliverables (a): 2-bit X-band RIS** | In progress (cell done, board in progress) | Baseline CP-2 octagon cell with MACOM MADP-000907-14020 diodes registered in `src/hologram/cell_family_registry.py`; diode RLC circuit models in `src/ris/diodes.py`. 1-bit LP interim code registered. | `f6688c9` (registry), `13463f5` (R3 registration), `1264bdc` (library v2) |
| **Deliverables (b): HMSA Designer App** | Feature complete for Tensor CP & RIS mode | App entry points `/tensor_cp.html` (`src/app/routers/tensor_cp.py`) and `/ris.html` (`src/app/routers/ris.py`, `src/ris/service.py`). | `4b87212` (Phase 27 close-out), `26b2433` (R5-5 RIS page) |
| **Deliverables (c): MoM Solver Tool** | Complete & frozen core | Periodic layered Green's function, rooftop basis, and impedance assembly in `src/solver.py`, `src/green_function.py`, `src/mom_matrix.py`, `src/sparameters.py`. | Tag `v1.1-handoff` |
| **HMSA Phases 8–27** | Done | End-to-end holographic synthesis, multi-port, hexagon SLL optimiser (`src/hologram/sll_optimizer.py`), tensor-CP single hexagon verified in HFSS (`hfss_results/hex150_tensor_cp_single/`). | `4b87212` (Phase 27), `58b864e` (Phase 26), `41f9ea0` (Board B export) |
| **RIS R1: Loads & Oblique** | Done | Oblique plane-wave excitation `compute_s_parameters_oblique` and lumped loads in `src/loads.py`, `tests/test_sparameters_oblique.py`. | `13f099f` |
| **RIS R2: Diode Models** | Done | Skyworks SMP1345-079LF PIN diode models and bias sweeps in `src/ris/diodes.py`, `results/phase20_r2/README.md`. | `9d645ce` |
| **RIS R3: CP-2 Cell** | Done | Baseline CP-2 octagon cell; HFSS check in `hfss_results/phase20_r3_cp_cell/NOTES.md` (`docs/STATUS.md:94`: "HFSS-checked at 10.5 GHz"). | `13463f5`, `ba12f3b` (retained baseline CP-2) |
| **RIS R4: State Library v2** | Done | Full-grid library v2 (8942 points at F = 2n = 120) assembled from server chunks in `results/phase20_r4/cp2_library_v2_fullgrid.npz`. Held-out validation gate PASS (`compute/r4-fullgrid-validate`). | `7ddc843` (merge full grid), `1264bdc` (assemble v2), `37074da` / `f6688c9` (gate PASS) |
| **RIS R5-1: Tolerance MoM** | Done | Server MoM tolerance predictions merged (`compute/r5-tolerance-mom`). Reviewed: mesh 90 reproduces HFSS sensitivity (`results/phase20_r5/TOLERANCE_GAP.md`). | `9d3deee` (merge), `14ac74d` (review) |
| **RIS R5-2: Forward Model** | Done | Plane-wave reflecting forward model in `src/ris/forward.py` tested against analytic references. | `3f2cfe2` |
| **RIS R5-3: Coding Synthesis** | Done | Coding synthesis (plain, dither, refined) and beam bandwidth tables in `src/ris/coding.py`, `results/phase20_r5/R5_3_REPORT.md`, `beam_table.json`. | `00528a4`, `0d46f83` (reference plane unified) |
| **RIS R5-5: RIS Designer Page** | Done | Web UI page `/ris.html` and backend router `/api/ris` operational on library v2 (`tests/test_app_phase20_r5_ris.py`). | `26b2433`, `f6688c9` |
| **RIS R5-6: Board Exporter Infra** | Done (Infra & 3×3 smoke test) | Board exporter `src/export/ris_board_export.py` with Load sheets, STEP, AEDT scripts; 3×3 smoke package in `hfss_exports/phase20_r5_ris_board_smoke3x3/`. | `977883b`, `e9e63dd` |
| **RIS R5-4: Supercell Locality** | In progress | Server jobs written and queued on Server 1: `jobs/r5-supercell-gradients.md` and `jobs/r5-supercell-random.md`. | `1a6b3d3` |
| **Phase 26 Board A HFSS Run** | Parked | Board A export package in `hfss_exports/hex4_tensor_cp_static/`; simulation run parked per `docs/STATUS.md:79`. | `58b864e` |
| **Parked Ideas** | Parked | Butler sector scanner (`docs/briefs/_parked/PHASE25_BUTLER_SECTOR_SCANNER.md`), FSS S-param tool (`docs/briefs/_parked/PHASE16_UNIT_CELL_SPARAM_TOOL.md`), nylon lens E/F (`docs/briefs/PHASE19_DIELECTRIC_TRANSMITARRAY_28GHZ.md`). | Commit `51ffbee` (Phase 25), `docs/STATUS.md:144` |
