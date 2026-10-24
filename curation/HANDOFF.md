# Project Handoff: `lightweight_mom` (Pass 5)

Authoritative onboarding, repository state, and handoff documentation prepared for the project lead, the parallel Claude session, and future engineers/AI agents.

---

## 1. Executive Summary & Active Branches

`lightweight_mom` (package `lwmom`) is a Python periodic Method-of-Moments solver for metasurface unit cells, built up into the **HMSA Designer** local web app (`http://127.0.0.1:8000`) for synthesizing holographic metasurface antennas and reconfigurable intelligent surfaces (RIS).

### Active Workspaces & Roles
- **Main Working Repository:** `/Users/biksmac/Documents/ISRO/HMSA/lightweight_mom` (macOS). The active branch `phase20-ris` is currently operated by the parallel Claude Code session. Ground rule: *Never touch, edit, push, or switch branches in that checkout.*
- **Curation Worktree:** `/Users/biksmac/Documents/ISRO/HMSA/lightweight_mom_curation` on branch `curation/2026-10-24` (branched from `1a6b3d3`). Holds all Pass 0–5 curation deliverables.
- **Compute Server Network:** Remote Windows servers (e.g. Server 1) perform all heavy computations via `jobs/<job>.md` on `compute/<job>` branches synced through the `hub` remote (`docs/briefs/INFRA_COMPUTE_SERVER_SYNC.md`). The local Mac acts strictly as coordinator and merges completed branches.

---

## 2. Planned Track & Current Status

### 2.1 The Three Deliverables
1. **(a) 2-bit PIN-Diode RIS at X-Band:** LP and CP as user-selectable families. Base CP-2 cell verified; 1-bit LP interim code registered; 10×10 board export infrastructure in progress.
2. **(b) HMSA Designer App:** Web app featuring the 3-step HMSA design wizard, Tensor CP flow (`/tensor_cp.html`), and RIS design mode (`/ris.html`).
3. **(c) MoM Solver Tool:** Periodic layered Green's function, rooftop basis functions, and matrix assembly core (`src/solver.py`, `src/sparameters.py`). Frozen core in use; solver deliverable to be packaged later.

### 2.2 Status Summary
- **HMSA Phases 8–27:** Completed and verified. Single-hexagon tensor-CP board has "HFSS-validated prediction (one run, 10.5 GHz, one board)" (`docs/STATUS.md:46`).
- **RIS Track (Phase 20):**
  - **R1 (Loads & Oblique):** Completed (`a04c2dd`, `36b4d1a`, `3aa6b94`). `MODEL_ONLY, NOT_HFSS_VALIDATED`.
  - **R2 (Diode Models):** Completed (`26551cb`, `286e889`). Skyworks SMP1345-079LF models in `src/ris/diodes.py`.
  - **R3 (CP-2 Cell):** Completed (`13463f5`, `ba12f3b`). Full-wave reference confirmed: "HFSS-checked at 10.5 GHz" (`hfss_results/phase20_r3_cp_cell/NOTES.md`, `docs/STATUS.md:94`).
  - **R4 (State Library v2):** Full-grid library (8942 points at $F = 2n = 120$) assembled from server chunks in `results/phase20_r4/cp2_library_v2_fullgrid.npz`. Validation gate: `GATE_PASS` (`37074da` / `f6688c9`).
  - **R5-1 (Tolerance Check):** Server MoM tolerance predictions reviewed (`results/phase20_r5/TOLERANCE_GAP.md`). MoM at mesh 90 reproduces HFSS $C_{\text{off}}$ ladder sensitivity ($\le 2.8^\circ$ vs HFSS $\le 2.2^\circ$).
  - **R5-2 (Forward Model) & R5-3 (Coding Synthesis):** Completed (`src/ris/forward.py`, `src/ris/coding.py`, `results/phase20_r5/R5_3_REPORT.md`).
  - **R5-5 (RIS Designer Page):** Built and operational (`src/ris/service.py`, `src/app/routers/ris.py`, `/ris.html`, commit `26b2433`).
  - **R5-6 (Board Exporter Infra):** Exporter infrastructure complete (`src/export/ris_board_export.py`, commit `977883b`); 3×3 smoke package built in AEDT without errors (`e9e63dd`).
  - **R5-4 (Supercell Locality):** In progress on Server 1 (`jobs/r5-supercell-gradients.md`, `jobs/r5-supercell-random.md`).
- **Parked Tracks:**
  - `_parked/PHASE19_DIELECTRIC_TRANSMITARRAY_28GHZ.md`: 28 GHz nylon lens (Parts A–D done; Parts E/F parked).
  - `_parked/PHASE25_BUTLER_SECTOR_SCANNER.md`: Phase 25 Butler sector scanner (paused).
  - `_parked/PHASE16_UNIT_CELL_SPARAM_TOOL.md`: Phase 16 unit-cell/FSS S-parameter tool.
  - Phase 26 Board A HFSS simulation run parked per `docs/STATUS.md:79`.

---

## 3. Cold-Reader Q&A

### Q1: What is this repository and what does it do?
> **Answer:** `lightweight_mom` (package `lwmom`) is a Python Method-of-Moments solver for periodic metasurface unit cells, built into a local web application ("HMSA Designer") that designs holographic/reflectarray metasurface antennas and reconfigurable intelligent surfaces (RIS) end to end.  
> **Source:** `docs/START_HERE.md:7-14`, `README.md:1-6`.

### Q2: What are the primary deliverables and project track?
> **Answer:** The project has three defined deliverables: (a) a 2-bit PIN-diode RIS at X-band with LP and CP families; (b) the HMSA Designer web app with tensor-CP and RIS design modes; (c) the standalone MoM solver deliverable.  
> **Source:** `docs/ROADMAP.md:7-13`.

### Q3: What is the current validation level of the core capabilities?
> **Answer:**
> - Unit cell solver (ungrounded): HFSS-validated off-resonance (Case G1 mean phase diff $-0.29^\circ$, $|S|$ diff $\sim 1.6 \times 10^{-5}$).
> - CP-2 RIS cell at 10.5 GHz: HFSS-checked at 10.5 GHz (`hfss_results/phase20_r3_cp_cell/NOTES.md`).
> - Single-hexagon tensor-CP board: HFSS-validated prediction (one run, 10.5 GHz, one board).
> - Oblique incidence & lumped loads: Model-only, verified analytically against Fresnel/Airy formulas ($10^{-10}$) and energy conservation.  
> **Source:** `docs/STATUS.md:29-94`, `docs/validation/README.md`.

### Q4: How do I run and test the codebase?
> **Answer:**
> - Setup: `python3 -m venv .venv && .venv/bin/pip install -e ".[dev,app,cad]"`
> - Start web app: `make app` (accessible at `http://127.0.0.1:8000`)
> - Fast test suite: `make test-fast` (`pytest -m "not slow"`, ~6–7 min, ~2900 tests)
> - Full test suite: `make test-full` (includes slow numerical solves, tens of minutes)  
> **Source:** `docs/START_HERE.md:60-67`, `CLAUDE.md:38-46`.

### Q5: What physical and coordinate conventions bite developers?
> **Answer:**
> 1. Time convention: $\exp(+j\omega t)$, with $\text{Im}(k_z) \le 0$ and $G = -Z$ (`docs/theory/CONVENTIONS.md`).
> 2. Reference planes: unit-cell library S-parameters reference 5 mm above the patch plane; array synthesis operates at the patch plane itself (`src/hologram/two_state_cells.py`).
> 3. Units: SI (meters, Hz) internally across all solver modules; mm and GHz only at UI and config boundaries.
> 4. Floquet cutoff: grounded cells with sub-wavelength lumped loads require $F = 2n$ summation cutoff (`results/phase20_r4/FLOQUET_AUDIT.md`).
> 5. Polarization & OAM sign rule: with a monopole feed, RHCP requires negative OAM charge and LHCP positive for vortex beams (`src/hologram/polarization_compensation.py`).  
> **Source:** `docs/START_HERE.md:70-85`, `docs/theory/CONVENTIONS.md`.

### Q6: How is heavy computation offloaded and synchronized?
> **Answer:** Heavy computational tasks (MoM sweeps, full-grid libraries, supercell solves) are never run on local Mac machines. Job requests (`jobs/<job>.md`) are dispatched to remote Windows servers via the `hub` remote on `compute/<job>` branches. The Mac acts strictly as coordinator and merges completed branches.  
> **Source:** `docs/briefs/INFRA_COMPUTE_SERVER_SYNC.md`, `jobs/README.md`, `CLAUDE.md:36`.

### Q7: Where are full-wave HFSS simulation data and benchmark references indexed?
> **Answer:** All external HFSS simulation exports, Touchstone S-parameter files, gain tables, and advisor run notes reside in `hfss_results/` and are indexed in `hfss_results/README.md`. Raw Touchstone files (`*.s4p`, `*.s36p`) are gitignored.  
> **Source:** `hfss_results/README.md`.

### Q8: Where are cell families and discrete-state libraries registered?
> **Answer:** The central cell registry is `src/hologram/cell_family_registry.py`. It registers cell geometry, substrate properties, state libraries, and validation regimes.  
> **Source:** `src/hologram/cell_family_registry.py`, `docs/START_HERE.md:32`.

---

## 4. Pending Actions for Next Session (Claude Session)

The curation branch `curation/2026-10-24` is fully prepared for merging into `phase20-ris`. The following handoff tasks should be executed by the Claude session:

1. **Merge Curation Branch:**
   ```bash
   git checkout phase20-ris
   git merge curation/2026-10-24 --no-ff -m "Merge branch 'curation/2026-10-24': repository audit, index updates, roadmap, and doc curation"
   ```
2. **Apply Protected File Proposals (`curation/PROPOSED_EDITS.md`):**
   - Apply the 9 relative link fixes in `docs/theory/GREEN_FUNCTION.md`, `MOM_FORMULATION.md`, and `S_PARAMETERS.md`.
   - Update `CLAUDE.md` "Physical model assumptions" section.
   - Update `docs/STATUS.md` header timestamp and Capability Matrix R4 row.
   - Insert missing 2026-10-24 entries into `CHANGELOG.md`.
   - Update links to parked brief `docs/briefs/_parked/PHASE19_DIELECTRIC_TRANSMITARRAY_28GHZ.md` in `docs/STATUS.md:98` and `docs/history/PROJECT_HISTORY.md:600`.
3. **Run Server-Side Full Test Suite:**
   Run full verification suite on the compute server post-merge.
