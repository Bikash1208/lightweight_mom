# WITHDRAWN 2026-10-24 (direction change): r5-supercell-random (Phase 20 R5-4b, server 1)

**Do not run this job.** Withdrawn 2026-10-24 per project direction change to monopole-fed coded holographic antenna (`docs/briefs/PHASE28_CODED_HOLOGRAPHIC_ANTENNA.md`). Reflectarray / plane-wave RIS designs with external feed are dropped. Partial results (no task finished) stay on the hub branch `compute/r5-supercell-random` with their logs. Previous job specifications retained below for record.

**From:** Mac session (coordinator). **Branch on the server:** `compute/r5-supercell-random` from the base commit. Never commit to `phase20-ris`. **Base commit:** `7d4cc209271ae05a9fe43074c92518d19f00a11e` on `phase20-ris` (fetch from `hub`).

## START CONDITION
`jobs/equivalence-check.md` passed and `jobs/r5-tolerance-mom.md` finished (queued behind it). May run together with `r5-supercell-gradients` (this job uses at most 9 workers).

## Purpose
Same array-level check for the **non-periodic codings**: 4 x 4 tiles at f = 10.5 GHz, normal incidence, F = 2n, 48 cells per period, all propagating orders vs the isolated-cell prediction, gate 1 dB / 15 deg on the strongest orders (reported, never tuned). Codings (33): the five R3 random 4x4 codings repeated at F = 2n (R3 ran them at the minimum cutoff only), and **dithered** designed codings at 30 / 40 / 50 deg (phi_s = 0), 30 / 50 deg (phi_s = 90), plus the plain codings at 30 / 50 deg for contrast; for each design four 4x4 tiles are cut from the designed 20x20 array (corner, centre, far corner, edge). The dither seed is 0. A tile is repeated periodically by the supercell solver, so this tests the neighbour coupling inside and across tile borders, not the finite-array edge.
9 tasks (up to 4 codings per task share one impedance-matrix assembly): `tiles4x4_0` ... `tiles4x4_8`.

## Commands (PowerShell)
```
cd C:\lightweight_mom ; git fetch hub ; git checkout -b compute/r5-supercell-random 7d4cc209271ae05a9fe43074c92518d19f00a11e ; .venv\Scripts\activate
$env:OPENBLAS_NUM_THREADS="1"; $env:OMP_NUM_THREADS="1"; $env:MKL_NUM_THREADS="1"
python scripts\research\phase20_r5_supercell.py random results\phase20_r5\supercell_random --workers 9 --smoke      (smoke at n_per_cell 8; delete its output afterwards)
python scripts\research\phase20_r5_supercell.py random results\phase20_r5\supercell_random --workers 9
```
(the fast suite is run by `r5-supercell-gradients`; if this job starts alone, run it first as in that file). 1 BLAS thread per worker.

## Cost, memory, checkpoint
A 4x4 supercell at 48 cells per period has about 17k metal unknowns (dense matrix ~ 4.6 GB per task, LU in place; allow ~10 GB per worker). Measured at the minimum cutoff on the Mac: 5 random codings in one task 33 min; at F = 2n the assembly is much larger: **estimate 3-10 h per task** (not measured). Each task writes `tiles4x4_<k>.json` when done; reruns skip finished tasks; **push each finished task JSON to `compute/r5-supercell-random` on `hub`**. Report the first finished task's `seconds_total` and `peak_memory_kb` to the Mac as soon as it exists (the Mac then decides about a 6x6 tile run; **no 6x6 is part of this job**).

## Outputs
`results/phase20_r5/supercell_random/tiles4x4_<k>.json`, manifest `results/phase20_r5/supercell_random_manifest.json` (file, size, sha256, producing commit, machine).

## Stop conditions
Smoke fails -> stop and report. Task error -> stop that task, finish the others, report. A single task above **12 h** or peak memory above 60 GB -> stop and ask.

## What the Mac does with the result
As for the gradients: coupling status per coding; the dithered-vs-plain comparison decides whether dither can be the default from ~30 deg (the user's wish) with an array-level error band.
