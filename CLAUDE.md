# CLAUDE.md

Guidance for Claude when working in this repository.

## What this is

`lightweight_mom` (package `lwmom`) — a lightweight Python Method-of-Moments solver for periodic
metasurface unit cells, built up into a local web app ("HMSA Designer") that designs holographic/
reflectarray metasurface antennas end to end. **Read `docs/START_HERE.md` and `docs/STATUS.md`
first** — they carry the purpose, design flow, repo map, conventions and current validation state
that used to live here. The full phase-by-phase narrative is in `docs/history/PROJECT_HISTORY.md`,
with one write-up per phase alongside it in `docs/history/`; active work items are in
`docs/briefs/` (a finished brief moves to `docs/history/briefs/`).

## Hard rules

- **Never move or delete a release tag** (`submission-2026-10-06`, `phase16-complete`,
  `pre-cleanup-2026-10-07`, ...). Tag new milestones instead.
- **Briefs drive work**: one brief per phase/task in `docs/briefs/`, one commit per logical item,
  fast suite (`make test-fast`) green before every commit.
- **No silent physics/solver changes.** New phases reuse earlier, frozen modules unmodified (wrap
  or duck-type) — see "Code conventions" below for the documented exceptions. If you must touch a
  frozen module, say so explicitly and verify bit-identical output for every pre-existing
  config/sweep/result.
- **Honest labels, never fabricated data.** Report `NOT_AVAILABLE_FROM_CURRENT_DATA` instead of
  inventing an HFSS/measured number. Gates return explicit named states, never a bare pass/fail;
  don't upgrade a status without the evidence.
- **Never commit** raw HFSS Touchstone exports (`*.s36p`, `*.s4p`) or reference-paper PDFs
  (`docs/references/*.pdf`) — both are gitignored on purpose; cite papers by DOI instead
  (`docs/references/README.md`).
- **Don't rename or move Python packages under `src/`** — the research-phase packages
  (`phase8b`…`phase13`) are imported across the app; see `docs/ARCHITECTURE.md`'s module map for
  who imports what. Document a concern instead of refactoring the package.
- Every phase updates `docs/STATUS.md` (last tag/commit, capability matrix, open next steps) as
  part of its own final commit.
- **Roadmap alignment**: Check `docs/ROADMAP.md` for project deliverables, completed phases, in-progress tasks, and priority sequence.
- **Remote compute server workflow**: Heavy compute (MoM sweeps, full-grid libraries, supercell solves) is **never run on the local Mac**. Dispatch jobs via `jobs/<job>.md` to remote Windows compute servers through the `hub` remote on `compute/<job>` branches, following `docs/briefs/INFRA_COMPUTE_SERVER_SYNC.md` and `jobs/README.md`. The Mac coordinates and merges.
- **External HFSS ground truth**: Returned HFSS simulations, Touchstone files, and benchmark notes are indexed in `hfss_results/README.md`. Never invent full-wave data.

## Commands

```bash
python3 -m venv .venv && .venv/bin/pip install -e ".[dev,app,cad]"
make app                                               # local web app, http://127.0.0.1:8000
make test-fast                                         # pytest -m "not slow", ~6-7 min
make test-full                                         # everything, tens of minutes
.venv/bin/python -m pytest tests/test_mom_matrix.py -q # a single file, from the project root
make lint / make format                                # ruff check / ruff check --fix + format
```

## Code conventions

- Import as `from src.xxx import ...`; `pytest.ini` sets `pythonpath = .`, run from the project root.
- Config is strict: no electrical value has a code default; unknown keys are errors.
- **Documented, additive exceptions** to "reuse frozen modules unmodified": Phase 12A (optional PEC
  ground plane), 12B (mesh alignment/duplicate audit), 12C (grounded supercells), 12D (general
  N-port Touchstone reader) and Phase 20 R1 (oblique incidence `compute_s_parameters_oblique`, lumped loads `UnitCell.loads`/`src/loads.py`, `CompositeGeometry`, load sheets in `cad_export`; R2 adds the optional `LumpedLoad.z_of_f`, default `None` = unchanged; R3 adds rotated loads (`direction` = angle), `PolyComposite` (`src/polygon_geometry.py`) and their STEP export, 0/90 deg and axis-aligned cases unchanged; the normal-incidence path only gains no-op hooks, θ=0 without loads is bit-identical to the pre-R1 solver, checked against `tests/data/r1_theta0_reference.json`), and the Phase 28 rotated-diode fix (side quest, commit `fcc6691`: `loads.default_transverse_conducts`; `diode_load` and `diodes_measured.load` default `transverse_conducts=None` and `hfss_chain` applies the same rule -> False for x/y/0/90 deg (unchanged), True for a rotated axis; **documented behaviour change**: a rotated diode load that does not pass `transverse_conducts` used to freeze to an open circuit in every state (ON == OFF) and now switches, and `check_load_geometry` emits a `UserWarning` for a rotated load with an explicit `transverse_conducts=False` (wanted: it flags the old silent failure; explicit values are still honoured). Bit-identical for every other case: x/y loads, explicit rotated loads, `hfss_chain` along x and the CP-2 cell, checked against `main` fingerprints in `tests/data/rotated_load_main_reference.json` by `tests/test_loads_rotated_regression.py`) each extend a core module directly, always additively, verified
  bit-identical for every pre-existing case (see each `docs/history/PHASE12*.md`). `src/phase13/
  convergence.py`'s Richardson extrapolation is the one permitted additive numerics option, app-only.
- Tests check mathematical/numerical properties (reciprocity, energy conservation, symmetry,
  convergence, known-answer cases), not just "it runs." Add tests with every change.
- No ML / fitted surrogates in the solver core or `src/phase11/` (a test greps the source for this).

## Physical model assumptions

Rectangular lattice (`period_x`, `period_y`), one cell; air | zero-thickness PEC pattern at z = 0 | one
isotropic substrate | air, or (opt-in, `backing: pec`) air | pattern | substrate | PEC ground — no
superstrate either way. Normal incidence is the original path; **oblique incidence** (arbitrary θ, φ,
TE/TM and full 2x2, `compute_s_parameters_oblique`) is an additive Phase 20 R1 extension. Only the
(0,0) Floquet order is an external port mode. **Lumped RLC loads** (`src/loads.py`, rotated loads and
`PolyComposite` from R3) and PIN-diode circuit models (`src/ris/diodes.py`: Skyworks SMP1345-079LF,
MACOM MADP-000907-14020) sit on the pattern metal. **Floquet cutoff F = 2n** (n = mesh size, e.g.
n = 60 -> F = 120) is required for grounded cells with sub-wavelength lumped loads
(`results/phase20_r4/FLOQUET_AUDIT.md`); all RIS numbers are at F = 2n. Production libraries: the
ungrounded 4.8 mm rectangular patch (RT/Duroid 5880, 10-15 GHz); a grounded 8 mm design point
(`src/phase12c/design_point.py`); the 2-bit grounded CP-2 octagon cell (15 mm, 10.5 GHz, state
library v2, 8942 points, `results/phase20_r4/cp2_library_v2_fullgrid.npz`). See `docs/STATUS.md` and
`docs/ROADMAP.md` for what is and isn't validated against HFSS at this model.
