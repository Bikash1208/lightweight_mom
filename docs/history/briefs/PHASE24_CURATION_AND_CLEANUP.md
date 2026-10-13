# Phase 24 — second curation and cleanup pass (no physics/solver changes)

**Status:** active brief, written 2026-10-13. Read `docs/START_HERE.md`, `docs/STATUS.md`, `CLAUDE.md`, and the
completed first pass `docs/briefs/PHASE17_HANDOFF_AND_CLEANUP.md` (this phase follows the same method). Phases 18–23
were done in several parallel sessions; the repo needs consolidating so that `main` is again the single truth and a cold
reader can pick it up.

## Safety rules (follow before touching anything)
1. **Tag first:** `git tag pre-cleanup-2026-10-13 phase23-hex-sll` (and keep every existing tag; never move or delete a tag).
2. **No behaviour changes.** No edits to physics, solver, synthesis, pattern or export logic. Moves/renames must keep
   imports working; the fast suite and the Phase 23 rectangle-snapshot gate
   (`scripts/research/phase23_rect_snapshot.py`) must stay green/identical before and after.
3. **Nothing is deleted without the user's approval.** Produce a list first (path, size, why, how to regenerate), stop,
   and wait. This includes branches, untracked folders and caches outside the usual Python ones.
4. One commit per logical item; fast suite green before each.

## State found on 2026-10-13 (verify, don't trust)
- `main` = Phase 17 merge (`910a82d`). Current work tip: `phase23-hex-sll`, which already contains
  `phase19-transmitarray`, `phase21-multiport-layouts-coded` and `phase22-auto-cell-states` (33 commits ahead of main).
- `phase18-four-port-300mm`: 1 commit (`a41d6db`, 300 mm square 4-port boards + `scripts/adhoc/coded4.py` fine-grid
  metrics, cited by Phase 23 as the x/y-cut SLL definition) — **not** contained in phase23. Expect conflicts in
  `scripts/adhoc/coded4.py`, `CHANGELOG.md`, `docs/STATUS.md`, `hfss_exports/README.md`.
- `phase18-one-bit-cell-check` = main (empty); `phase17-cleanup` already merged.
- `docs/briefs/` still holds `PHASE17_HANDOFF_AND_CLEANUP.md` (finished — should be in history) and
  `PHASE19_DIELECTRIC_TRANSMITARRAY_28GHZ.md` (Parts E/F waiting on a user decision). Untracked: `PHASE20_RIS_ROADMAP.md`,
  `PHASE20_R1_SOLVER_LOADS_OBLIQUE.md` (RIS work, not started).
- Untracked export folders from other sessions: `hfss_exports/four_port_1bit_square_{gap0,gap28,300mm_gap0,300mm_gap28}/`.
  Gitignored leftover: `hfss_exports/one_bit_pairs_10p54ghz/` (4 STEP files from an abandoned single-cell check).
- Sizes: `.venv` 1.5 GB, `results/` 568 MB, `hfss_exports/` 243 MB, `docs/` 59 MB, `.git` 64 MB, tracked total ~98 MB.
  Caches: 29 `__pycache__` dirs, `.pytest_cache`, `.ruff_cache`, `lwmom.egg-info`; gitignored-but-heavy:
  `results/unit_cell_library/{geometry,database,_solve}`, `scripts/adhoc/fss/coef/` (~1 h to regenerate).
- Largest tracked files: STEP files in `hfss_exports/two_beam_20deg_phi0_180/` and `docs/submission/` (21, 13, 6.8,
  4.8 MB) and ~1 MB `layout.dxf` per export package.

## Work items

**1. Consolidate branches.**
- Merge `phase18-four-port-300mm` into `phase23-hex-sll` (resolve conflicts keeping both sides' content; the 300 mm
  square packages and Phase 23's citation of `coded4.py` must both survive). Run the fast suite.
- Then merge the result into `main` (`--no-ff`, message listing phases 18–23). Report the merge commit.
- Propose (don't execute until approved) deleting the fully-merged branches: `phase17-cleanup`,
  `phase18-one-bit-cell-check`, `phase18-four-port-300mm`, `phase19-transmitarray`, `phase21-multiport-layouts-coded`,
  `phase22-auto-cell-states`, `phase23-hex-sll`. Milestone tags preserve them: add `phase19-partD`, `phase21-complete`,
  `phase22-complete`, `phase23-complete` at the right commits.

**2. Briefs.**
- Move `PHASE17_HANDOFF_AND_CLEANUP.md` to `docs/history/briefs/`.
- `PHASE19_DIELECTRIC_TRANSMITARRAY_28GHZ.md`: leave active only if Parts E/F are really still open; otherwise move it
  and record the open decision (±30.00° hard spec?) in STATUS "Open next steps".
- Commit the two `PHASE20_*` RIS briefs under `docs/briefs/_parked/` (RIS is paused, not started) and list them in
  STATUS "Parked ideas".

**3. Untracked and leftover outputs** — list, then wait for approval:
- `hfss_exports/four_port_1bit_square_*` (4 folders): find which session/script made them (`scripts/adhoc/square4/`),
  whether Phase 21/22 reproduce them, and propose: commit (README/CSV/JSON/PNG, STEP stays ignored) and index, or delete
  as superseded.
- `hfss_exports/one_bit_pairs_10p54ghz/` and `scripts/adhoc/one_bit_pair_cells.py`: abandoned single-cell check of the
  10.54 GHz 1-bit pairs (STATUS open step 2). Propose: keep the script (documented, parked) and delete the STEP leftovers,
  or delete both.
- Any other untracked/ignored file that no README explains.

**4. Caches and regenerable data** — list sizes, then (after approval) delete: all `__pycache__`, `.pytest_cache`,
`.ruff_cache`, `lwmom.egg-info`. Do **not** delete `.venv`, `results/unit_cell_library/*` or `scripts/adhoc/fss/coef/`
without an explicit yes (they are slow to rebuild); report the regenerate command and time for each instead.
Add `make clean` (Python caches only) and `make distclean` (prints, doesn't delete, the heavy regenerable dirs).

**5. Repo growth policy (going forward only — no history rewrite).**
- New export packages: commit README, design report JSON, PNGs and the small CSV; keep `*.step`, `*.stl` **and**
  `layout.dxf`/`layout.json` over 1 MB out of git (`.gitignore` with the same "already-tracked files stay" exception
  Phase 17 used), each package README naming its regenerate command. Don't untrack anything already tracked.
- Report tracked size before/after.

**6. Docs refresh** (the cold-start set must describe phases 18–23):
- `docs/STATUS.md`: one consolidated capability matrix (Phase 19 nylon RCWA/lens, Phase 21 multi-port layouts/exact
  export, Phase 22 auto cell states, Phase 23 polygon ports + SLL optimiser), last tag, open next steps (HFSS runs
  pending: hex150 boards, square boards, nylon lens; Phase 19 decision; RIS parked), known limitations.
- `docs/START_HERE.md` repo map (new packages: `src/dielectric_cell`, `src/transmitarray`, `hologram/auto_states.py`,
  the SLL optimiser module, polygon ports), `docs/ARCHITECTURE.md` module map + who-imports-what, `docs/USER_GUIDE.md`
  check against the running app, `CHANGELOG.md` (one entry per phase 18–23), `docs/history/PROJECT_HISTORY.md`
  (add 18–23), `hfss_exports/README.md` and `results/README.md` (index every folder: what, which script, status).
- `docs/briefs/README.md` reflects what's active/parked.

**7. Hygiene (report, fix only zero-risk items).** `ruff check` count vs Phase 17's ~95 (don't bulk-fix); dead-code
scan as in Phase 17; no stray debug prints/TODOs; fast-suite runtime (still < 5 min?) and test count; any test that is
platform-sensitive (e.g. `tests/test_phase.py::test_branch_cut_negative_real_axis_is_plus_180_for_both_signed_zeros`
failed on Linux Python 3.12 but passes on the Mac) — mark/fix only if clearly a test bug.

**8. Acceptance.**
- Fresh clone of `main` into a temp dir, venv from the documented steps only, fast suite, app starts, `/` and
  `/openapi.json` load; the multi-port hex4 board and "Reduce sidelobes" run end to end in a headless browser.
- Cold-agent test (subagent with no context, repo only): (a) what the app does; (b) validation level of the unit-cell
  solver and of each Phase 19–23 capability; (c) how to design and export a 4-port hexagonal board with sidelobe
  optimisation; (d) which HFSS packages await simulation and in what order; (e) what's parked (RIS, Phase 19 E/F);
  (f) what must never be done. Fix docs until all answers are right; include the Q&A in the report.
- `git status` clean on `main`; tag `v1.1-handoff`; report the tree (depth 2), everything moved/deleted, sizes before/after.

## Do not
- Change any computed result, rewrite git history, force-push, move/delete tags, or delete anything unapproved.
- Start RIS (Phase 20) or Phase 19 E/F work.
