# Infra — server 1 implements, hub shares, Mac views (working model from 2026-10-24)

**Status:** active infrastructure brief. First written 2026-10-23 for a Mac-coordinator / Windows-compute split; **rewritten 2026-10-24** after the user's decision that server 1 is now the only working session
(code, tests, compute) and the Mac session is retired. Not a physics phase: no solver/physics change is allowed under this brief except portability fixes, each additive and bit-identical for every existing case.

## Roles
- **Server 1 (implementer, the only working session):** Windows Server 2019, 48 logical CPUs, 255.6 GB RAM, RTX 2000 (unused by the solver), Tailscale `100.88.61.96`, user `Bikash`, repo `C:\lightweight_mom`,
  Python 3.13 (Anaconda) venv `.venv`. Writes code, tests and briefs, runs all compute, commits on phase branches (currently `phase28-coded-holographic`), updates `docs/STATUS.md` and `CHANGELOG.md`, pushes to `hub`.
  Merges to `main` and release tags are made here too, after the fast suite and only when the user's stated condition holds.
- **Hub (shared repo):** SMB share `hmsa_sync` on server 1. Server path `C:\hmsa_sync`; `hmsa_sync\lightweight_mom.git` is the bare repo, remote name **`hub`**; `hmsa_sync\data\<job>\` holds large gitignored artefacts
  (libraries `.npz`, caches, task files, STEP, HFSS exports) with a manifest committed to git (rule 4).
- **Mac (viewer):** reads the hub (`/Volumes/hmsa_sync`, Finder ⌘K → `smb://100.88.61.96/hmsa_sync`, remote `hub` = `/Volumes/hmsa_sync/lightweight_mom.git`) to read reports and results. It does not commit or push. If the mount is missing
  the viewer remounts; nothing is ever pushed to another path or remote instead. The advisor reviews the server's reports through the user.
- Link: Tailscale. SSH is not available; SMB (port 445) is.

## Rules (all `CLAUDE.md` rules stay in force)
1. **Branches.** Phase work is committed on the phase branch (one commit per logical item, fast suite before each commit, stop at every ■ in a brief and report). Never rewrite history, never move or delete a tag,
   no deletion of branches/files/results without the user's approval. Heavy compute may use `compute/<job>` branches pushed to `hub`; they are merged or left as records. A job file `jobs/<job>.md` is **optional** now
   (the server may run its own compute without one); when one is written it states base commit, exact commands, workers, expected outputs and stop condition, and the server runs exactly that.
2. **Provenance on every result** (JSON sidecar or header; `scripts/research/provenance.py`): machine name, OS, git commit, dirty flag, Python and numpy versions, BLAS vendor, thread settings, worker count, run parameters, start/end time.
3. **Data manifest.** Large outputs live in `hmsa_sync\data\<job>\`; a manifest in git (`results/<phase>/<job>_manifest.json`: file, size, sha256, producing commit, machine) points to them. Text files written on Windows
   are normalised to LF by git (`.gitattributes`), so manifests of text files also give the LF-normalised sha256. Never commit large files.
4. **Fast suite before each commit:** `pytest -m "not slow" -n 24 -q` with `OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=1` (about 17-26 min; without the thread variables it was 4.7x slower). A docs-only commit does not
   change test outcomes; its message says so and refers to the last full run. Known failing tests on this machine (see "Known Windows-server test differences") are reported, not hidden.
5. **Threads and memory.** Workers × BLAS threads ≤ 48 including everything else running. Check free RAM before a job and use measured peaks, not estimates (below). A memory-guarded launcher
   (`scripts/infra/gradients_guarded.py` on `compute/r5-supercell-gradients`) is the pattern: start a task only if free RAM ≥ 30 GB + expected peak.
6. **Checkpointing:** every evaluation/task cached atomically (one file per task), resumable after a reboot or RDP disconnect; a long job stops new task starts at a stated wall-clock cap and lets running tasks finish
   (tasks do not checkpoint mid-run) until a final stop time.
7. **Honest labels:** server-authored results say so; nothing is called validated beyond its evidence.

## Measured facts (server 1, used for planning)
- Equivalence (2026-10-23, `compute/equivalence-check`): CP-2 cell, mesh 60, F = 2n, 4 states, θ = 0/30°, 10.5 GHz: Mac vs server **INTERCHANGEABLE**, max relative difference 4.3e-10 at 1 and 8 threads;
  server 1 thread vs 8 threads 1.1e-10.
- Cost of one 4-state task at F = 2n, 1 thread: mesh 60 89 s (0.43 GB peak), mesh 90 692 s (1.91 GB); mesh 120 at the rank-minimum cutoff 5.8 GB (at F = 2n not measured, extrapolated ~3500 s). With 40 workers each task runs
  about 1.5-1.8x slower than alone (memory bandwidth). Throughput for the R4 full grid: 8942 mesh-60/90 tasks in 8.52 h with 40 workers.
- Supercell MoM (R5-4, stopped): matrix order N = 1510 × nx × ny at 48 cells per period (exact); worker peak ≈ 1.28 × 16 N² bytes (4×4 tile: N = 24,160, 11.92 GB); cost grows steeply (3×1 25 min, 4×1 1.1 h with 3 threads).
- Phase 23 snapshot gate: reproducible on this machine (three runs at 1 thread give identical hashes, also identical at 8 threads). The Mac reference `rect_snapshot_after_C.json` stored hashes only, so the cross-machine
  difference of that synthesis cannot be measured from it; the gate is now numeric and server-based (`scripts/research/phase23_rect_numbers.py`, reference `results/phase23/rect_snapshot_server1.json`).

## Windows portability (state 2026-10-24)
- Every text write in `src/` outside the frozen `phase8b`-`phase13` packages uses `encoding="utf-8"` and `newline="\n"` (CSV writers keep the csv module's CRLF and gain the encoding); a scan test enforces it.
- `.gitattributes` (`* text=auto eol=lf`, CSV/HFSS txt `-text`, binaries marked); the Windows system git has `core.autocrlf=true`, which the attributes override for text files.
- `src/diagnostics.py` guards the Unix-only `resource` import (memory `NOT_AVAILABLE` on Windows).
- Commit identity on this machine is given per command (`git -c user.name=... -c user.email=...`); no git config is changed.

## Known Windows-server test differences (last-digit floating point, not logic)
Four fast-suite tests fail on this machine at `38e677f` because their references were produced on the Mac's CPU/BLAS and compare exactly:
`test_sparameters_oblique::test_normal_incidence_results_are_bit_identical_to_the_phase20_base_commit`, `test_tensor_cp_design::test_single_hexagon_cells_match_committed_layout_bit_identically`,
`test_tensor_cp_service::test_app_export_of_the_single_hexagon_matches_the_committed_phase26_folder`, `test_app_phase27_tensor_cp::test_single_hexagon_design_export_compare_end_to_end`.
Measured size (the service export's `predictions.json` against the committed reference): 27 of 69 numbers differ, max relative 1.4e-15. They stay as known items (report each run) until the user decides between
server-generated references and tolerance-based comparisons.
