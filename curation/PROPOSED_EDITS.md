# Proposed Edits for Protected Documentation

This document records exact proposed text replacements for protected files (`docs/theory/*`, `CLAUDE.md`, `docs/STATUS.md`, and `CHANGELOG.md`). In accordance with the curation ground rules, these protected files are not edited directly by the curation worktree and are prepared here for handoff to the Claude session.

Every proposal includes the target file, exact old text, exact replacement text, and cited repository evidence.

---

## 1. `docs/theory/` Relative Link Fixes (9 links)

### 1.1 `docs/theory/GREEN_FUNCTION.md` (Line 3)

- **Target File:** `docs/theory/GREEN_FUNCTION.md`
- **Location:** Line 3
- **Old Text:**
```markdown
Implementation: [src/green_function.py](src/green_function.py) (`LayeredGreenFunction`).
```
- **New Text:**
```markdown
Implementation: [`src/green_function.py`](../../src/green_function.py) (`LayeredGreenFunction`).
```
- **Evidence:** The relative link target `src/green_function.py` from within `docs/theory/` resolves to non-existent path `docs/theory/src/green_function.py`. The source module resides at the repository root `src/green_function.py` (verified on `phase20-ris`).

---

### 1.2 `docs/theory/MOM_FORMULATION.md` (Lines 3–6)

- **Target File:** `docs/theory/MOM_FORMULATION.md`
- **Location:** Lines 3–6
- **Old Text:**
```markdown
Implementation: [src/basis_functions.py](src/basis_functions.py) (mesh, rooftops, transforms),
[src/mom_matrix.py](src/mom_matrix.py) (Z, V), [src/excitation.py](src/excitation.py)
(unpatterned-stack field, reference planes), [src/solver.py](src/solver.py)
(`assemble_system`, `build_system`, `solve`).
```
- **New Text:**
```markdown
Implementation: [`src/basis_functions.py`](../../src/basis_functions.py) (mesh, rooftops, transforms),
[`src/mom_matrix.py`](../../src/mom_matrix.py) (Z, V), [`src/excitation.py`](../../src/excitation.py)
(unpatterned-stack field, reference planes), [`src/solver.py`](../../src/solver.py)
(`assemble_system`, `build_system`, `solve`).
```
- **Evidence:** 4 link targets (`src/basis_functions.py`, `src/mom_matrix.py`, `src/excitation.py`, `src/solver.py`) are missing leading `../../` to resolve to repository root `src/` modules (verified on `phase20-ris`).

---

### 1.3 `docs/theory/S_PARAMETERS.md` (Lines 3–6)

- **Target File:** `docs/theory/S_PARAMETERS.md`
- **Location:** Lines 3–6
- **Old Text:**
```markdown
Implementation: [src/sparameters.py](src/sparameters.py) (extraction, data model),
[src/green_function.py](src/green_function.py) (`sheet_fields`, `stack_fields`),
[src/excitation.py](src/excitation.py) (`ReferencePlanes`, `stack_tangential_field`, now with `port`),
[src/solver.py](src/solver.py) (`solve` with several right-hand sides).
```
- **New Text:**
```markdown
Implementation: [`src/sparameters.py`](../../src/sparameters.py) (extraction, data model),
[`src/green_function.py`](../../src/green_function.py) (`sheet_fields`, `stack_fields`),
[`src/excitation.py`](../../src/excitation.py) (`ReferencePlanes`, `stack_tangential_field`, now with `port`),
[`src/solver.py`](../../src/solver.py) (`solve` with several right-hand sides).
```
- **Evidence:** 4 link targets (`src/sparameters.py`, `src/green_function.py`, `src/excitation.py`, `src/solver.py`) are missing leading `../../` to navigate from `docs/theory/` to `src/` (verified on `phase20-ris`).

---

## 2. `CLAUDE.md` Physical Model Assumptions

- **Target File:** `CLAUDE.md`
- **Location:** "Physical model assumptions" section (lines 50–57)
- **Old Text:**
```markdown
## Physical model assumptions

Rectangular lattice, one cell; air | zero-thickness PEC pattern at z = 0 | one isotropic substrate
| air, or (opt-in, `backing: pec`) air | pattern | substrate | PEC ground — no superstrate either
way. Normal incidence only; only the (0,0) Floquet order is a port mode. One production
(ungrounded) library: 4.8 mm period, rectangular PEC patch, RT/Duroid 5880, 10–15 GHz; a grounded
8 mm design point is also recorded (`src/phase12c/design_point.py`). See `docs/STATUS.md` for what
is and isn't validated against HFSS at this model.
```
- **New Text:**
```markdown
## Physical model assumptions

Rectangular or hexagonal unit cell; layered air | zero-thickness PEC pattern at z = 0 | isotropic substrate
| PEC ground (or ungrounded / air backing) — no superstrate. Both normal and oblique incidence
(`compute_s_parameters_oblique`, arbitrary θ, φ); only the (0,0) Floquet order is an external port mode.
Lumped RLC loads and PIN-diode circuit models (Skyworks SMP1345-079LF, MACOM MADP-000907-14020) at
arbitrary orientations on the pattern metallization (`src/loads.py`, `src/ris/diodes.py`). Floquet-order
summation cutoff requires F = 2n (mesh size n, e.g. n=60 -> F=120, n=90 -> F=180) for grounded cells with
sub-wavelength lumped loads (`results/phase20_r4/FLOQUET_AUDIT.md`). Production libraries:
(1) ungrounded 4.8 mm rectangular patch (RT5880, 10–15 GHz); (2) 2-bit grounded CP-2 octagon cell
(15 mm period, 10.5 GHz, 8942 points in `results/phase20_r4/cp2_library_v2_fullgrid.npz`); (3) 1-bit LP
interim code from CP-2; (4) tensor-CP pure-circular polarization cells. See `docs/STATUS.md` and
`docs/ROADMAP.md` for validated regimes and limits.
```
- **Evidence:**
  - Oblique incidence implemented in Phase 20 R1 (`src/loads.py`, `src/sparameters.py`, commit `13f099f`).
  - Lumped loads and diode models added in Phase 20 R1–R3 (`src/loads.py`, `src/ris/diodes.py`, `src/polygon_geometry.py`).
  - Production CP-2 grounded octagon cell library v2 assembled in Phase 20 R4 (`results/phase20_r4/cp2_library_v2_fullgrid.npz`, commit `1264bdc`).
  - Floquet cutoff audit establishes that $F = 2n$ is required for grounded loaded cells (`results/phase20_r4/FLOQUET_AUDIT.md`).

---

## 3. `docs/STATUS.md` Header & Capability Matrix R4 Row

### 3.1 Header Update

- **Target File:** `docs/STATUS.md`
- **Location:** Line 3
- **Old Text:**
```markdown
**Last updated:** 2026-10-23, by Claude (first HFSS check of the CP-2 RIS cell recorded, RIS operating frequency 10.5 GHz, R3b plan written and awaiting approval; Phase 20 R4 closed out on branch `phase20-ris`: CP-2 state library v1 + Floquet-cutoff audit; **all RIS numbers are now at Floquet cutoff F = 2n, the R3 values are superseded**; waiting for the user's HFSS run of `hfss_exports/phase20_r3_cp_cell/`. Earlier: Phase 20 R3 closed out on branch `phase20-ris`: CP-2 octagon cell with MADP diodes, LP as 1-bit codes of it, both registered MODEL_ONLY, HFSS check package pending; R4 brief next. Earlier: Phase 27 on branch `phase27-app-tensor-cp`, from `phase26-leaky-tensor-cp`: the tensor-CP design flow is in the app, with a cell-family registry; next is Phase 20 R1. Phase 26 on branch `phase26-leaky-tensor-cp`, from `phase25-butler-sector`: single-hexagon tensor-CP board checked against HFSS; 4-hex board scoped; Phase 24 curation merged into `main`, tagged `v1.1-handoff`; the Phase 25 locality caveat below stands). **Rule: every phase updates this
file** (last release tag/commit, capability matrix, open next steps) as part of its own final commit.
```
- **New Text:**
```markdown
**Last updated:** 2026-10-24, by Claude (R4 gate closed: library v2 full-grid 41x10x24 at F = 2n = 120 assembled and passing held-out validation gate on both random and phi-draws; R5-1 MoM tolerance check merged and reviewed: mesh 90 reproduces HFSS C_off sensitivity; R5-2 plane-wave forward model and R5-3 coding synthesis completed; R5-5 RIS designer app page operational; R5-6 board export infrastructure and 3x3 smoke package completed; R5-4 supercell locality server jobs queued; first HFSS check of CP-2 cell at 10.5 GHz recorded in hfss_results/phase20_r3_cp_cell/NOTES.md). **Rule: every phase updates this
file** (last release tag/commit, capability matrix, open next steps) as part of its own final commit.
```
- **Evidence:** Captures progress achieved on 2026-10-24 (commits `1264bdc`, `14ac74d`, `26b2433`, `977883b`, `e9e63dd`, `f6688c9`).

---

### 3.2 Capability Matrix R4 Row Update

- **Target File:** `docs/STATUS.md`
- **Location:** Table row for "2-bit RIS cells and the supercell / locality tool"
- **Old Text:**
```markdown
| 2-bit RIS cells and the supercell / locality tool (Phase 20 R3, `src/ris/cells.py`, `src/ris/supercell.py`, `src/ris/polarization.py`, `src/polygon_geometry.py`, rotated loads in `src/loads.py`; registry families `ris_cp2_octagon_madp_rt5880_15`, `ris_lp_1bit_from_cp2_madp_rt5880_15`) | **MODEL_ONLY, NOT_HFSS_VALIDATED**; check package `hfss_exports/phase20_r3_cp_cell/` pending; diode values PROVISIONAL_DATASHEET (MADP L_s assumed 0) | CP-2 (period 15 mm, RT5880 1.575 mm, 10 GHz), **current numbers at F = 2n** (R3 minimum-cutoff values superseded, `results/phase20_r4/FLOQUET_AUDIT.md`): theta 0 ladder deviation 10.7 deg, tolerance step errors 17.4 / 20.8 (C_T +20 %, open) / 19.3 deg, amplitude -1.16 dB, opposite-handed -15.2 dB (meets); theta 30: opposite-handed **-10.0 dB (open)**, amplitude -1.35 dB; **all R3 targets met only at the 10.0 GHz grid point for theta <= 7.5 deg**; state library `results/phase20_r4/cp2_library_v1.npz` (interpolation 1.81 deg / 0.166 dB at the 95th pct, gate pass; mesh-correction residual 12 deg / 0.96 dB near 10.3-10.8 GHz; theta > 30 EXTRAPOLATED; Wood-guard refusal); R3 (superseded) read: ladder deviation 12.1 deg, opposite-handed -13.6 / -14.0 dB; model uncertainty ~10 deg; LP-1/LP-2/LP-3 2-bit cells miss (LP-3 nearest); locality: 4-cell gradient passes (0.1 dB, 3 deg), random 4x4 codings 2/10 (**not rerun at F = 2n**; gradient rerun at F = 2n: unchanged) (`results/phase20_r3/OPTIMISATION.md`, `LOCALITY.md`, `LP_INTERIM.md`) |
```
- **New Text:**
```markdown
| 2-bit RIS cells and the supercell / locality tool (Phase 20 R3/R4, `src/ris/cells.py`, `src/ris/supercell.py`, `src/ris/polarization.py`, `src/polygon_geometry.py`, rotated loads in `src/loads.py`; registry families `ris_cp2_octagon_madp_rt5880_15`, `ris_lp_1bit_from_cp2_madp_rt5880_15`) | **HFSS-checked at 10.5 GHz** (`hfss_results/phase20_r3_cp_cell/NOTES.md`); state library v2 full-grid GATE_PASS (held-out gate: random p95 0.31 deg / 0.012 dB; phi-draw p95 2.46 deg / 0.109 dB); diode values PROVISIONAL_DATASHEET (MADP L_s assumed 0) | CP-2 (period 15 mm, RT5880 1.575 mm, 10.5 GHz), **numbers at F = 2n = 120** (`results/phase20_r4/FLOQUET_AUDIT.md`): HFSS check confirms ladder progression (98/92/102/68 deg); state library v2 (`results/phase20_r4/cp2_library_v2_fullgrid.npz`, 8942 points, 41x10x24) assembled and default for R5 (`CalibratedLibrary.load(version="v2")`); tolerance check (`results/phase20_r5/TOLERANCE_GAP.md`) confirms mesh-90 MoM reproduces HFSS C_off ladder shifts (<= 2.8 deg vs HFSS <= 2.2 deg); locality tool (`src/ris/supercell.py`) with server jobs `r5-supercell-gradients` and `r5-supercell-random` in progress |
```
- **Evidence:**
  - Status changed to "HFSS-checked at 10.5 GHz" per `docs/STATUS.md:94` and `hfss_results/phase20_r3_cp_cell/NOTES.md`.
  - State library v2 full grid assembled in `1264bdc`, held-out gate passed in `37074da` / `f6688c9` (`compute/r4-fullgrid-validate`).
  - Tolerance gap report in `results/phase20_r5/TOLERANCE_GAP.md` (commit `14ac74d`).

---

## 4. `CHANGELOG.md` Missing 2026-10-24 Entries

- **Target File:** `CHANGELOG.md`
- **Location:** Below main preamble, immediately preceding `## Unreleased - Phase 20 R5 brief (branch phase20-ris)`
- **Insertion Text:**
```markdown
## Unreleased - Phase 20 R5 progress & R4 library v2 validation (2026-10-24, branch `phase20-ris`)
- **R4 Gate PASS - Library v2 held-out validation** (`compute/r4-fullgrid-validate` @ `ed76d73`, merged in `f6688c9`): R4 item 4 held-out validation gate passes on both draws: 60 random points p95 0.31 deg / 0.012 dB; 40 points in resonance phi-band (75–105 deg) p95 2.46 deg / 0.109 dB. Full-grid library v2 (`results/phase20_r4/cp2_library_v2_fullgrid.npz`, 8942 points, 41 x 10 x 24 at F = 2n = 120) assembled from server chunks (`1264bdc`) is now the validated R5 default.
- **R5-1 MoM tolerance check review** (`compute/r5-tolerance-mom`, commits `9d3deee`, `14ac74d`): MoM at mesh 90 with F = 2n reproduces HFSS C_off ladder sensitivity (steps shift <= 2.8 deg vs HFSS <= 2.2 deg; `results/phase20_r5/TOLERANCE_GAP.md`). Mesh 60 does not at 10.5 GHz; R4 mesh-60 tolerance table at 10.5 GHz is superseded.
- **R5-2 & R5-3: Forward model and coding synthesis** (commits `3f2cfe2`, `00528a4`, `0d46f83`): Plane-wave reflecting forward model (`src/ris/forward.py`); coding synthesis with plain, dither, and refined algorithms (`src/ris/coding.py`); beam table, 1 dB / 3 dB gain bandwidth, and C_off sensitivity report (`results/phase20_r5/R5_3_REPORT.md`); reference plane unified (`PLANE_CORRECTION.md`).
- **R5-5: RIS designer web page** (commit `26b2433`): App page `/ris.html` and REST backend `/api/ris/*` (`src/ris/service.py`, `src/app/routers/ris.py`) added to HMSA Designer "More tools"; verified with `tests/test_app_phase20_r5_ris.py`.
- **R5-6: Board exporter infrastructure** (commits `977883b`, `e9e63dd`): RIS board exporter (`src/export/ris_board_export.py`) generating Load sheets, STEP models, and AEDT Python scripts; 3x3 smoke package created in `hfss_exports/phase20_r5_ris_board_smoke3x3/`.
- **Phase 26 Board B export** (commit `41f9ea0`): Export package for 4-port tensor-CP hexagon Board B.
```
- **Evidence:** Reconstructs the 2026-10-24 commit log on `phase20-ris` (`1264bdc`, `7ddc843`, `9d3deee`, `14ac74d`, `26b2433`, `977883b`, `e9e63dd`, `41f9ea0`, `ed76d73`, `f6688c9`).

---

## 5. Updates for Moved Briefs in Protected Documentation (Pass 3)

### 5.1 `docs/STATUS.md` (Line 98)

- **Target File:** `docs/STATUS.md`
- **Location:** Line 98 (Open Next Steps, item 1)
- **Old Text:**
```markdown
1. **Decision, Phase 19**: is ±30.00° a hard specification for the nylon lens? The error bar puts the beam at ≈30.3° median in the
   primary model (0.4° typical pointing error). Parts E/F of `docs/briefs/PHASE19_DIELECTRIC_TRANSMITARRAY_28GHZ.md` wait on this.
```
- **New Text:**
```markdown
1. **Decision, Phase 19**: is ±30.00° a hard specification for the nylon lens? The error bar puts the beam at ≈30.3° median in the
   primary model (0.4° typical pointing error). Parts E/F of `docs/briefs/_parked/PHASE19_DIELECTRIC_TRANSMITARRAY_28GHZ.md` wait on this.
```
- **Evidence:** `PHASE19_DIELECTRIC_TRANSMITARRAY_28GHZ.md` moved to `docs/briefs/_parked/` in Pass 3.

---

### 5.2 `docs/history/PROJECT_HISTORY.md` (Line 600)

- **Target File:** `docs/history/PROJECT_HISTORY.md`
- **Location:** Line 600
- **Old Text:**
```markdown
comparison) wait on the user's ±30.00° decision and an HFSS run. Brief: `briefs/PHASE19_DIELECTRIC_TRANSMITARRAY_28GHZ.md`.
```
- **New Text:**
```markdown
comparison) wait on the user's ±30.00° decision and an HFSS run. Brief: `briefs/_parked/PHASE19_DIELECTRIC_TRANSMITARRAY_28GHZ.md`.
```
- **Evidence:** Brief moved to `docs/briefs/_parked/`.
