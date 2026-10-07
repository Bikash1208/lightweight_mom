# Phase 17 brief — Handoff-ready repo: onboarding docs + professional file layout

Goal: a new person, an AI agent, or a fresh session can open this repo cold and, within minutes, know what the
project is, what works (and how well it is validated), how to run it, where everything lives, and exactly what to do
next — and the tree contains only what belongs there. **No physics/solver behaviour changes in this phase.**

## 0. Safety rules (do these first, follow throughout)
1. Tag the current state `pre-cleanup-2026-10-07` and work on branch `phase17-cleanup`; merge to main only when the
   fast suite is green and the fresh-clone check (§5) passes. Do not move `submission-2026-10-06`.
2. **Tracked research history is archived with `git mv`, never deleted.** Deleting tracked files is allowed only for
   true junk (duplicates, stray copies) and must be listed in the final report.
3. **Untracked/ignored data: list first, delete only after the user approves the list** (it may be irreplaceable —
   the raw HFSS `.s36p` exports were already lost once). Exception: caches (`__pycache__`, `.pytest_cache`,
   `.ruff_cache`, `.DS_Store`, `*.egg-info`) may be removed freely.
4. Do not rename or move Python packages under `src/` (phase packages import each other: e.g. `phase9` is imported
   from 24 other modules, `phase8b/8d_rev/8e/9/12d/13` are used by the app). Document instead (§2).
5. Commit per section; fast suite green before each commit.

## 1. Onboarding documents (the core deliverable)
Create/rewrite, short and current (history belongs in docs/history):
- **`README.md`** (front door, ≤ 1 screen + quickstart): what HMSA Designer / lightweight_mom is (one paragraph),
  screenshot, install, `make app`, `make test`, links to START_HERE, USER_GUIDE, STATUS, ARCHITECTURE.
- **`docs/START_HERE.md`** — for humans *and* AI agents resuming work. Sections: purpose in 5 lines; the design
  flow (monopole `.fld` → beams & board → pattern & export → HFSS); repo map (table: folder → what's in it → touch?);
  how to run/test; conventions that bite (exp(+jωt), reference planes (patch plane vs 5 mm), mm vs m, IEEE
  RHCP/LHCP labels and the OAM sign rule, the server must be restarted after backend changes); the working method
  (briefs in `docs/briefs/`, one commit per item, fast suite before commit, never move release tags, live browser
  check); "If you are an AI agent: read this, then STATUS.md, then the active brief".
- **`docs/STATUS.md`** — single source of truth for *now*: last release tag + commit; capability matrix with
  validation level per feature (HFSS-validated / independently cross-checked / model-only / not supported) —
  e.g. unit cell vs HFSS 0.13–0.52° (G1), ~4° at resonance (G4); vortex↔pencil physics; multi-port (model-only,
  HFSS pending); 1-bit coded (cell pair mesh-sensitive, HFSS pending); polarization efficiency ≈ 50 % (scalar-cell
  physics); known limitations (normal incidence only, single layer, PEC zero-thickness metal, no transmitarray, no
  finite-array MoM, 2-bit multi-port not wired); **open next steps** (HFSS checks of `hfss_exports/three_port_25deg`
  and `four_port_1bit_gap0/28`, unit-cell check of the 10.54 GHz 1-bit pair); parked ideas (link `docs/briefs/_parked/`:
  unit-cell/FSS tool with plane-defined S; also transmitarray, oblique incidence, tensor cells, 2-bit multi-port).
  Add a "Last updated: <date>, by <who>" line and a rule that every phase updates it.
- **`CLAUDE.md`** — rewrite to ≤ ~80 lines: what the repo is (3 lines), "read docs/START_HERE.md and docs/STATUS.md
  first", hard rules (tags, briefs, tests, no silent physics changes, honest labels, don't commit PDFs/.s36p),
  commands, where things are. Move the current long history text into `docs/history/PROJECT_HISTORY.md`.
- **`docs/ARCHITECTURE.md`** — update with a module map in three layers: (a) core periodic MoM solver
  (geometry, basis_functions, green_function, mom_matrix, solver, sparameters, …); (b) HMSA design (hologram,
  dispersion, pattern, fields, export, reflectarray, polarization); (c) app (src/app); plus (d) research-phase
  packages (`phase8b … phase13`) — frozen, still imported, with one line each on what they are and who imports them.
  Include the data-flow diagram of the app's 3 steps.
- **`docs/USER_GUIDE.md`** — make sure it matches today's UI (stepper navigation, Ports/defaults, level matching,
  pointing correction, polarization metrics, Coded = experimental) and add a "Taking a design to HFSS" section.
- **`docs/validation/README.md`** — one page indexing every validation with numbers and where the data is:
  Phase 12 HFSS unit-cell comparisons; independent free-space far-field cross-check of vortex/pencil and multi-port
  pointing; the square-loop FSS reference (periodic + finite-panel plane-S, `hfss_exports/fss_square_loop_10mm`);
  the Mishra 26 GHz reproduction.
- **`CHANGELOG.md`** — bring up to date (phases 13–16 + this one), Keep-a-Changelog style.

## 2. Docs layout
- `docs/history/` ← `PROJECT_HISTORY.md`, all `docs/phases/*.md` (git mv; keep names), and **completed** briefs
  (`docs/briefs/PHASE12*…PHASE16_FOLLOWUPS.md`) into `docs/history/briefs/`. Fix every inbound link (grep).
- `docs/briefs/` keeps only active briefs + `_parked/` + a `README.md` explaining the brief workflow.
- `docs/theory/` stays. `docs/screenshots/` stays (refresh the README screenshot if the UI changed).
- `docs/references/`: keep `.pdf` ignored; add `README.md` listing each paper with full citation + DOI and what it
  is used for (Ghosh OJAP 2026, Mishra TVT 2026, the Tier-1–4 list from CLAUDE.md).
- Add `docs/README.md` (index of the docs tree).

## 3. Code-adjacent files
- `examples/`: keep user-facing scripts (`run_sweep.py`, `build_unit_cell_library.py`, `inspect_unit_cell_library.py`,
  `make_submission_figures.py`, `make_paper_board_26ghz.py`, `print_config.py`, `sweep_config.yaml`) and git mv the
  research drivers (`run_phase8a…run_phase13_part0.py`, `export_*`) to `scripts/research/`; fix paths/imports and
  any Makefile/CI/test references; add `examples/README.md` and `scripts/README.md` (one line per script).
- `results/`: tracked research outputs stay (move under `results/research/` only if every reference is updated —
  otherwise leave in place and add `results/README.md` indexing each folder: which phase produced it, which script
  regenerates it, whether tests read it). `results/unit_cell_library/` heavy parts are already ignored — document the
  regenerate command.
- `results/submission_work/` (untracked scratch from Oct 5–6): **list contents for the user.** Proposed: move the
  reusable scripts into `scripts/adhoc/` with a README stating status (`fss/plane_s.py` + `fss/mom_save.py` +
  `fss/fssgeom.py` = finite-panel measurement-plane S for the square-loop FSS, parked; `three_port.py`, `coded4.py`,
  `four_port.py` = generators of the hfss_exports packages); delete the rest (logs, probes, `_tmp_example.fld`,
  `_moved_example.fld`, snapshots) after approval.
- `hfss_exports/`: add `hfss_exports/README.md` indexing every package (what, frequency, how generated, status:
  awaiting HFSS / validated). For the untracked packages (`three_port_25deg`, `four_port_1bit_gap0`,
  `four_port_1bit_gap28`, `fss_square_loop_10mm`): commit README/CSV/JSON/PNG; keep `*.step` out of git
  (add `hfss_exports/**/*.step` ignore EXCEPT the already-tracked phase12/two_beam STEP files) and record the
  regenerating script for each.
- Root: confirm `requirements*.txt` vs `pyproject.toml` are consistent (one source of truth; keep the lock file);
  `Makefile` targets documented in README (`make app`, `make test`, `make test-slow`, `make lint`, …); remove stray
  root files if any; `.DS_Store` stays ignored.
- `inputs/`: add `README.md` (what `efl_phase_EzZ_example.fld` is: the monopole reference field at 10.54 GHz).

## 4. Hygiene
- Dead code check: list modules under `src/` that nothing imports and no test covers (import graph); **report, do not
  delete** — the user decides.
- `ruff` clean, no stray `print` debugging, no TODOs without an owner/brief reference.
- Make sure the app shows its version/commit somewhere (footer "HMSA Designer · <tag or short sha>") so screenshots
  and bug reports are traceable.
- `make app`: print a hint "restart after pulling backend changes".

## 5. Acceptance checks
1. **Fresh clone:** clone into a temp dir, create a venv from the documented steps only, install, run the fast suite,
   start the app headless and load `/` and `/openapi.json`. Fix any doc step that fails.
2. **Cold-agent test:** spawn a subagent with NO conversation context, give it only the repo, and ask it to answer
   from the docs alone: (a) what does the app do; (b) how validated is the unit-cell solver and against what; (c) why
   do monopole-fed boards radiate a vortex and how is it fixed; (d) how to design and export a 2-port board;
   (e) what are the open next steps; (f) where are the HFSS packages awaiting simulation; (g) what must you never do
   (tags, PDFs, .s36p…). Fix the docs until every answer is correct; include the Q&A in the final report.
3. `git status` clean on the branch; tree listing (depth 2) in the report; list of everything moved/deleted.
4. Tag the merged result `v1.0-handoff` (separate from the submission tag).
