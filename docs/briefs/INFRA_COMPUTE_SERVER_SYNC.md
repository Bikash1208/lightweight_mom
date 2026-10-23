# Infra — Mac coordinator + Windows compute servers (sync protocol)

**Status:** active infrastructure brief, written 2026-10-23. Applies to every session working on this repo from now on.
Not a physics phase: no solver/physics change is allowed under this brief except portability fixes listed below,
each additive and bit-identical on the Mac.

## Machines
- **Mac (coordinator):** `~/Documents/ISRO/HMSA/lightweight_mom`. Owns `phase20-ris` (and any phase branch), briefs,
  `docs/STATUS.md`, `CHANGELOG.md`, merges, close-outs.
- **Server 1 (compute):** Windows, 48 cores, RTX 2000 (unused by the solver), Tailscale `100.88.61.96`, user `Bikash`,
  repo `C:\lightweight_mom`, Python 3.13 (Anaconda) venv `.venv`. Server 2: same pattern when added.
- Link: Tailscale. SSH is not available; SMB (port 445) is.

## Shared hub (SMB share `hmsa_sync` on server 1)
- Server path `C:\hmsa_sync`; Mac mount `/Volumes/hmsa_sync` (Finder ⌘K → `smb://100.88.61.96/hmsa_sync`).
- `hmsa_sync/lightweight_mom.git` — bare git repo, remote name **`hub`** on both machines
  (Mac: `/Volumes/hmsa_sync/lightweight_mom.git`; server: `C:\hmsa_sync\lightweight_mom.git`).
- `hmsa_sync/data/<job>/` — large gitignored artefacts (libraries `.npz`, caches, STEP, HFSS exports).

## Rules
1. **Branches.** The server never commits to `phase20-ris` or any phase branch. Compute work goes on
   `compute/<job>` branches (e.g. `compute/r3b-rescore`, `compute/r4b-library`), created from the commit the Mac
   names in the job request, pushed to `hub`. The Mac fetches `hub`, reviews, and merges or cherry-picks.
2. **Job requests** come from the Mac session as a short file `jobs/<job>.md` committed on `phase20-ris` and pushed
   to `hub`: base commit, exact command(s), worker count, expected outputs, stop condition. The server runs exactly
   that and reports.
3. **Provenance on every result** (JSON sidecar or header): machine name, OS, git commit, dirty flag, Python and
   numpy versions, BLAS vendor (`numpy.show_config()`), thread settings, worker count, run parameters, start/end time.
4. **Data manifest.** Large outputs live in `hmsa_sync/data/<job>/`; a manifest committed to git
   (`results/<phase>/<job>_manifest.json`: file, size, sha256, producing commit, machine) points to them. Never commit
   the large files.
5. **Cross-machine equivalence (once, before any server result is used).** Run the same small fixed case on both
   machines (e.g. one CP-2 cell at mesh 60, F = 2n, 4 states, θ = 0/30°, 10.5 GHz; plus the Phase 23 snapshot gate)
   and record max |Δphase| and |Δ dB|. Expect last-digit differences (different CPU/BLAS); results are interchangeable
   if the differences are far below the model uncertainty (≤ 1e-6 relative), otherwise stop and report.
   Bit-identity gates stay **Mac-only**; on the server, report the difference size.
6. **Threads.** On the server set `OMP_NUM_THREADS` / `MKL_NUM_THREADS` / `OPENBLAS_NUM_THREADS` so that
   workers × threads ≤ cores; **every server job sets `OPENBLAS_NUM_THREADS=1` per worker (measured: `-n 24` without it was 4.7x slower)**; memory limits workers (mesh 120 ≈ 9 GB per worker) — check RAM before launching.
7. **Checkpointing** as in R4/R3b: every evaluation cached, resumable after a reboot/RDP disconnect.

## Mac push rule (added 2026-10-23)
**Before every push to `hub`** the Mac session checks that `/Volumes/hmsa_sync` is mounted and that `/Volumes/hmsa_sync/lightweight_mom.git` exists (`test -d /Volumes/hmsa_sync/lightweight_mom.git`).
If either is missing (the SMB mount drops after sleep / network changes), **stop and tell the user to remount** (Finder ⌘K → `smb://100.88.61.96/hmsa_sync`); never push to another path, a local copy or a different remote instead.
Same check before any fetch from `hub`.

## Portability items (Mac session implements; additive; bit-identical on the Mac)
- `src/diagnostics.py`: guarded import of the Unix-only `resource` module (memory reported as `NOT_AVAILABLE` on
  Windows), with a test.
- Any other Windows-only failure found by the server's first fast-suite run: report first, fix in the same way.
- Server fast-suite report 2026-10-23 (Windows, 48 cores, 255.6 GB, OpenBLAS 0.3.34): 13 failed / 85 errors = 93 from `resource` (fixed above), 2 last-digit BLAS differences (~1e-15, expected, not fixed) and 3 Windows-only test issues, fixed portably without physics change: `test_app_jsonutil::test_path_becomes_string` (compares with `str(path)`), `test_phase8b_conventions::test_coordinate_system_matches_conventions_md_text` (reads UTF-8 explicitly), `test_tensor_cp_service` export-vs-reference comparison (the gitignored reference `layout.json` / `layout.dxf` may be absent: the export must then still produce the file).

## First steps
1. Server: create the share and bare repo (done by the user), push nothing yet.
2. Mac: add remote `hub`, `git push hub phase20-ris` (and other active branches).
3. Server: `git remote add hub C:\hmsa_sync\lightweight_mom.git`, fetch, check out the pushed `phase20-ris`.
4. Both: run the equivalence case (rule 5); the Mac records the result in `docs/STATUS.md`.
