# WITHDRAWN 2026-10-24 (direction change): r5-supercell-gradients (Phase 20 R5-4a, server 1)

**Do not run this job.** Withdrawn 2026-10-24 per project direction change to monopole-fed coded holographic antenna (`docs/briefs/PHASE28_CODED_HOLOGRAPHIC_ANTENNA.md`). Reflectarray / plane-wave RIS designs with external feed are dropped. Partial results (3 of 11 tasks finished) stay on the hub branch `compute/r5-supercell-gradients` with their logs. Previous job specifications retained below for record.

**From:** Mac session (coordinator). **Branch on the server:** `compute/r5-supercell-gradients` from the base commit. Never commit to `phase20-ris`. **Base commit:** `7d4cc209271ae05a9fe43074c92518d19f00a11e` on `phase20-ris` (fetch from `hub`).

## START CONDITION
`jobs/equivalence-check.md` passed (it did: INTERCHANGEABLE) **and `jobs/r5-tolerance-mom.md` has finished** (queued behind it; the two supercell jobs may then run at the same time: this one uses at most 11 workers, `r5-supercell-random` at most 9, leaving the rest of the 48 cores free).

## Purpose
Array-level check of the designed gradient codings (the local-response assumption of the R5-3 forward model): supercell MoM of periodic windows, baseline CP-2, MADP diodes, f = 10.5 GHz, normal incidence, **Floquet cutoff F = 2n**, 48 mesh cells per 15 mm cell, every propagating order |m|,|n| <= 3 compared with the isolated-cell prediction (the one the design uses). Gate 1 dB / 15 deg on the strongest orders: reported, never tuned.
11 tasks (one supercell geometry each, no shared work): x-gradients 8x1 (13.8 deg), 7x1 (15.8), 6x1 (18.5), 5x1 (22.4), 4x1 (28.4), 7x1 q=2 (33.0), 3x1 (39.4), 8x1 q=3 (45.5), 5x1 q=2 (49.6); y-gradients 1x4 (28.4) and 1x5 q=2 (49.6) (the cell is not 4-fold symmetric, so phi_s = 90 is a separate case).

## Commands (PowerShell)
```
cd C:\lightweight_mom ; git fetch hub ; git checkout -b compute/r5-supercell-gradients 7d4cc209271ae05a9fe43074c92518d19f00a11e ; .venv\Scripts\activate
$env:OPENBLAS_NUM_THREADS="1"; $env:OMP_NUM_THREADS="1"; $env:MKL_NUM_THREADS="1"
python -m pytest -m "not slow" -q > results\phase20_r5\server1_fast_suite_r5_4a.txt     (first step; report failures other than the known BLAS last-digit ones; do not patch)
python scripts\research\phase20_r5_supercell.py gradients results\phase20_r5\supercell_gradients --workers 11 --smoke      (smoke: one task at n_per_cell 8, minutes; check it prints "done")
python scripts\research\phase20_r5_supercell.py gradients results\phase20_r5\supercell_gradients --workers 11
```
(create `results\phase20_r5` first; delete the smoke output before the real run). 1 BLAS thread per worker (infra rule 6).

## Cost, memory, checkpoint
Measured earlier (Mac, loaded): 4x1 at F = 2n, two load sets: 19-28 min. Cost grows steeply with n (assembly ~ unknowns^2 x F^2): the 7x1 and 8x1 tasks are estimated at **several hours each** (about 4-10 h on one thread); memory per task ~ unknowns^2 x 16 B, up to ~2-3 GB for 8x1. Each finished task writes `<name>.json` (with `seconds_total` and `peak_memory_kb`); a rerun skips finished tasks, so an interruption loses only the running tasks. **Push the finished task JSONs to `compute/r5-supercell-gradients` on `hub` whenever a task completes** (small files), so the Mac can start the analysis on partial results.

## Outputs
`results/phase20_r5/supercell_gradients/<task>.json` (per-order measured vs predicted level/phase for r_RR and r_LR, gate flag, provenance, timing, peak memory), `server1_fast_suite_r5_4a.txt`, manifest `results/phase20_r5/supercell_gradients_manifest.json` (file, size, sha256, producing commit, machine).

## Stop conditions
Smoke fails or the fast suite shows a new failure -> stop and report. Any task error -> stop that task, let the others finish, report (do not patch code). Wall clock beyond **12 h** with tasks still running -> stop and ask. Peak memory above 60 GB per task -> stop and report.

## What the Mac does with the result
Fetches `hub`, builds the per-coding coupling status (`COUPLING_VERIFIED` / `COUPLING_CORRECTED` / `COUPLING_FAILED_ISOLATED_CELL_MODEL`) and the error band on gain / SLL / specular lobe; stop point 2.
