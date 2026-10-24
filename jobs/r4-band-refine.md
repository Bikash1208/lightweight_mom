# Job request: r4-band-refine (library v2.1, resonance band 10.4-10.9 GHz, server 1)

**From:** Mac session (coordinator). **Branch on the server:** `compute/r4-band-refine` from the base commit. Never commit to `phase20-ris`. **Base commit:** `68814ae2aad573f67345cd7a4e6ce93da47a400a` on `phase20-ris` (fetch from `hub`).

## START CONDITION
Run **after both supercell jobs have finished** (`r5-supercell-gradients`, `r5-supercell-random`): this job uses up to 40 workers and the supercell tasks must not be starved or evicted. The v2 library stays the default until the Mac has assembled and reviewed v2.1.

## Why
Library v2 (0.1 GHz nodes) passes the R4 gate (p95), but every large error sits in the resonance band 10.4-10.9 GHz: the draw-B maximum is 19.3 deg at 10.75 GHz (theta 12.7, phi 101) and 5 deg at 10.45 GHz (`results/phase20_r4/README.md`). The band is where the mesh correction is largest and changes fastest (mesh 60 -> 90 moves r_RR by 37-47 deg at 10.5 GHz, 6 deg at 10.0 GHz, 0.5 deg at 11.0 GHz), so both the table nodes and the correction nodes (0.5 GHz apart) are too coarse there.

## Parts
**A. Frequency nodes (mandatory).** Five new frequencies, the ones in 10.4-10.9 GHz that are not on the 0.1 GHz grid: **10.45, 10.55, 10.65, 10.75, 10.85 GHz**, full theta 0-45 step 5 / phi 0-345 step 15 (217 points each), four states, **mesh 60, F = 2n**: 1085 tasks. Plus the **mesh-90 correction points** at those five frequencies (the R4 correction set: (theta, phi) = (0,0), (30,0), (30,90), (45,0), (45,90)): 25 tasks. (The mesh-60 correction points are grid points, so they need no extra tasks.)
**B. Validation points (mandatory).** 30 extra held-out points drawn only inside 10.4-10.9 GHz (all theta 0-45, all phi, seed 20261025, same draw rules as before), four states each, mesh 60, F = 2n, with the direct-solve references saved: 30 tasks. The Mac then evaluates v2 and v2.1 on draw A, draw B (their saved reference files are already in the repo) and these 30 points with `phase20_r4_validate.py --use-refs` (no new solves), and reports the gate plus the maximum error inside the band.
**C. Mesh-120 correction points (recommended, decision requested).** The library's correction is the three-mesh Richardson map (60/90/120); the Mac holds mesh-120 data only at the old nodes (10.0, 10.5, 11.0 ...). Without mesh 120 at the new frequencies the v2.1 correction in the band would have to be a mesh-90-only Richardson (p = 2), which differs from the three-mesh map by up to 24 deg (p95 16 deg) at the old nodes (`results/phase20_r4/fullgrid_review.json`, `mesh90_only_correction_residual`). 25 tasks, ~3500 s each, ~9 GB each. **Not part of the default run; run it only if the user confirms (it is the honest way to build v2.1).**

## Commands (PowerShell)
```
cd C:\lightweight_mom ; git fetch hub ; git checkout -b compute/r4-band-refine 68814ae2aad573f67345cd7a4e6ce93da47a400a ; .venv\Scripts\activate
$env:OPENBLAS_NUM_THREADS="1"; $env:OMP_NUM_THREADS="1"; $env:MKL_NUM_THREADS="1"
python scripts\research\phase20_r4_band_refine_run.py run C:\hmsa_sync\data\r4-band-refine --workers 40 --hours 3
python scripts\research\phase20_r4_band_refine_run.py consolidate C:\hmsa_sync\data\r4-band-refine
$env:PHASE20_R4_LIB="results\phase20_r4\cp2_library_v2_fullgrid.npz"; $env:PHASE20_R4_WORKERS="40"
python scripts\research\phase20_r4_validate.py a 30 results\phase20_r4\validation_v2_band.json --f-range 10.4,10.9 --seed 20261025 --save-refs results\phase20_r4\validation_v2_band_refs.npz
```
Part C only on the user's yes, after A and B are finished and pushed (so the 9-GB tasks do not compete with them):
```
python scripts\research\phase20_r4_band_refine_run.py run C:\hmsa_sync\data\r4-band-refine --workers 16 --hours 4 --mesh120
python scripts\research\phase20_r4_band_refine_run.py consolidate C:\hmsa_sync\data\r4-band-refine
```
(the second run skips every finished task and does the 25 mesh-120 ones.)

## Estimate (from the measured server timings: mesh 60 89 s alone / ~137 s with 40 workers, mesh 90 692 s alone, mesh 120 ~3500 s extrapolated)
* Part A: 1085 mesh-60 tasks x ~137 s / 40 = **~1.0 h**; the 25 mesh-90 tasks (~700-1000 s each) run first in the same pool and overlap with it: **total A about 1.2-1.5 h**. 1110 tasks.
* Part B: 30 tasks x ~90 s / 40 = **~3-5 min** (the driver run is finished first).
* Part C (if confirmed): 25 tasks x ~3500 s / 16 workers = 2 waves = **~2 h**, 144 GB RAM at 16 x 9 GB.
* Total with C: ~3.5-4 h wall; without C ~1.5 h.

## Resources, hard caps and stop conditions
* A and B: 40 workers x 1 BLAS thread; memory: mesh 90 ~2 GB, mesh 60 ~0.5 GB per worker, < 40 GB in total. **Hard cap 3 h** (`--hours 3`; the driver does not start a task that cannot finish before the stop).
* C: at most **16 workers x ~9 GB = 144 GB**; hard cap 4 h; stop and report if free memory falls below 30 GB.
* Any task error: the driver prints it and continues with the others; report the errors, do not patch code. A rerun resumes (per-task files).

## Outputs (small copies in git on the compute branch, data on the share)
`C:\hmsa_sync\data\r4-band-refine\` (`tasks\`, `chunks\grid_f*.npz`, `corr60_f*.npz`, `corr90_f*.npz`, (`corr120_f*.npz`), `MANIFEST.json` with sha256, `provenance_start/end.json`); in git: the manifest, provenance files, run log, `validation_v2_band.json` and `validation_v2_band_refs.npz`.

## What the Mac does next
Fetch, verify sha256, assemble **v2.1** with `phase20_r4_assemble.py` (grid = the 41 + 5 frequencies via `--f-list`, correction frequencies via `--corr-f-list`), evaluate v2 and v2.1 on the same points (draws A, B, the 30 band points) and report the gate and the in-band maximum error; v2 stays the default until then and v2.1 is not made the default without the user's yes.
