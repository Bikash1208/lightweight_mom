# Start here

For a human picking this up, or an AI agent resuming work with no prior context.

## Purpose, in 5 lines

`lightweight_mom` (package `lwmom`) is a Python Method-of-Moments solver for periodic metasurface
unit cells, built up into a local web app ("HMSA Designer") that designs holographic/reflectarray
metasurface antennas end to end. You upload a monopole feed's measured surface-wave field, pick
beam directions and a board, and the app synthesizes a hologram, predicts its radiation pattern,
and exports a layout ready for HFSS. Every output is labelled with how far it's actually been
validated (unit-cell-vs-HFSS, cross-checked, model-only) — see `docs/STATUS.md` for the real numbers.
This is an ISRO-internal research project, not a finished product.

## The design flow

```
monopole .fld (Step 1)  ->  beams & board (Step 2)  ->  pattern & export (Step 3)  ->  HFSS
   reference wave            hologram synthesis          far-field prediction          verify
   (uploaded, phase +         (object beams, board/      (vector model: co/cross-       in a
   optional magnitude)        substrate, design          pol, SLL, XPD, pointing/       real
                              method, ports, coded        level-matching reports)       solver
                              quantisation)                + CSV/JSON/DXF/STEP export
```

Run it: `make app`, then walk the 3-step wizard at `http://127.0.0.1:8000`. One worked example
lives at `inputs/reference_fields/efl_phase_EzZ_example.fld` (see `inputs/README.md`).

## Repo map

| Folder | What's in it | Touch when... |
|---|---|---|
| `src/` (flat modules) | Core periodic MoM solver: geometry, basis functions, Green's function, MoM matrix, S-parameters, unit-cell database | changing solver physics (rare — see `CLAUDE.md`'s reuse rule) |
| `src/dispersion/`, `src/fields/`, `src/hologram/`, `src/pattern/`, `src/reflectarray/`, `src/export/`, `src/verify/`, `src/hfss_compare/` | HMSA design layer: surface-wave dispersion, `.fld` import, hologram synthesis, far-field pattern, CAD/CSV export, full-array verification | adding a design capability |
| `src/hologram/ports.py`, `auto_states.py`, `sll_optimizer.py`, `sll_board.py`; `src/pattern/cellwise_pattern.py`; `src/export/exact_board.py`, `aedt_script.py` | Phases 21-23: multi-port layouts incl. polygon/hexagon ports (`hex4`), automatic coded-cell states, the cell-wise far-field model and the sidelobe optimiser ("Reduce sidelobes", `/api/sidelobes/optimize`), exact STEP + AEDT script for boards | multi-port, coded or sidelobe work (the Phase 23 rectangle gate: `scripts/research/phase23_rect_snapshot.py`) |
| `src/dielectric_cell/`, `src/transmitarray/` | Phase 19 (model-only): RCWA unit cell for 3D-printed dielectric (nylon) cells and the horn-fed transmitarray lens design/error-bar code; not wired into the app | the nylon lens (`docs/briefs/PHASE19_DIELECTRIC_TRANSMITARRAY_28GHZ.md`) |
| `src/phase8b/ … src/phase13/` | Frozen research-phase packages, still imported by the app (`docs/ARCHITECTURE.md` has the who-imports-what table) | **don't rename or move these** — document, don't refactor |
| `src/app/` | FastAPI backend + static HTML/JS frontend | adding/changing a page or endpoint |
| `tests/` | One `test_<module>.py` per module; `test_app_*.py` for the app | every change needs a test |
| `docs/` | `START_HERE.md` (this file), `STATUS.md`, `ARCHITECTURE.md`, `USER_GUIDE.md`, `theory/` (conventions), `validation/` (numbers + where), `history/` (phase write-ups + completed briefs), `briefs/` (active work items) | see `docs/README.md` for the full index |
| `examples/` | User-facing scripts (`make_paper_board_26ghz.py`, `build_unit_cell_library.py`, ...) | `scripts/research/` has the one-off research drivers instead |
| `configs/`, `inputs/` | Config-file presets; the example reference field | |
| `scripts/` | `research/` (per-phase drivers: `run_phase*.py`, `phase19_*.py`, `phase23_*.py`), `adhoc/` (rescued one-offs: `coded4.py`, `square4/`, FSS tools) — see each folder's README | reproducing a published number |
| `results/`, `hfss_exports/` | Generated outputs and HFSS export packages (each with its own `README.md` indexing status) | read `results/README.md` / `hfss_exports/README.md` before assuming something is stale |

## Where things are in the repo (depth 1)

`src/` code · `tests/` · `docs/` (this file, `STATUS.md`, `ARCHITECTURE.md`, `USER_GUIDE.md`, `DEVELOPING.md`, `theory/`, `validation/`,
`history/`, `briefs/`, `submission/`) · `examples/` · `scripts/` · `configs/`, `inputs/` · `results/` · `hfss_exports/` ·
`CHANGELOG.md`, `CLAUDE.md`, `Makefile` (`make clean`, `make distclean`, `make clean-app-outputs`: `docs/DEVELOPING.md` "Cleaning").
Active work: `docs/briefs/README.md` lists what is active and what is parked.

## How to run / test

```bash
python3 --version  # must be >= 3.11 -- substitute a 3.11+ interpreter below otherwise
python3 -m venv .venv && .venv/bin/pip install -e ".[dev,app,cad]"
make app          # start the web app at http://127.0.0.1:8000 -- RESTART after any backend change
make test-fast    # pytest -m "not slow", ~6-7 min (2550 tests), run before every commit
make test-full    # everything, including real MoM solves, tens of minutes
```

## Conventions that bite

- **`exp(+jωt)` time convention**, SI units internally (mm/GHz only in configs/UI). `Im(kz) <= 0`,
  `G = -Z`. Source of truth: `docs/theory/CONVENTIONS.md` — if code disagrees with it, that's a
  bug: fix one and add a test.
- **Two different reference planes**: the auto-built unit-cell libraries store S11 at a plane
  **5 mm above the patch**; a real board's reflection happens **at the patch plane itself**. Mixing
  these up silently shifts every phase by a frequency-dependent amount — see
  `src/hologram/two_state_cells.py`'s own module docstring.
- **mm vs m**: the UI and configs use mm/GHz; every solver/hologram function takes SI (m/Hz).
  Conversions happen at the app boundary, not inside `src/`.
- **IEEE right-hand-rule circular polarization naming**, and a project-specific **OAM sign rule**:
  with a monopole feed, RHCP needs a *negative* OAM charge and LHCP a *positive* one for a clean
  vortex (the opposite sign gives a cross-polar pencil, not a vortex) — see
  `src/hologram/polarization_compensation.py:oam_sign_mismatch`.
- **The app server must be restarted after any backend (`src/`) change** — there is no
  autoreload. The frontend is static HTML/JS served as-is; a browser refresh picks that up.

## The working method

- Work is driven by **briefs** in `docs/briefs/` (one markdown file per phase/task); a finished
  brief moves to `docs/history/briefs/`.
- **One commit per logical item**, fast suite (`make test-fast`) green before each commit.
- **Never move or delete a release tag** (`submission-2026-10-06`, `phase16-complete`, ...) —
  tag new milestones instead.
- Verify a UI change by actually **running the app in a real browser** (or headless), not just by
  reading the code — several real bugs in this project were only caught that way
  (`docs/DEVELOPING.md` "App tests specifically").
- Never fabricate data (HFSS/measured numbers): report `NOT_AVAILABLE_FROM_CURRENT_DATA` instead.
  Never commit raw HFSS Touchstone (`*.s36p`/`*.s4p`) exports or reference-paper PDFs — both are
  gitignored on purpose (`docs/DEVELOPING.md`).

## If you are an AI agent

Read this file, then `docs/STATUS.md` (what's validated, what's open right now), then the active
brief(s) in `docs/briefs/` (anything other than `_parked/` and `README.md`). `docs/ARCHITECTURE.md`
and `docs/USER_GUIDE.md` are the next two stops if you're touching the app specifically.
