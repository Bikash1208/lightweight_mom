# Job request: r4-band-refine (library v2.1, resonance band 10.4-10.9 GHz, server 1) - revised: Part C approved and extended

**From:** Mac session (coordinator). **Branch on the server:** `compute/r4-band-refine` from the base commit. Never commit to `phase20-ris`. **Base commit:** `dc2e664dc19ae8b36d25eec7bc209df62e731493` on `phase20-ris` (fetch from `hub`).

## START CONDITION
Run **after both supercell jobs have finished** (`r5-supercell-gradients`, `r5-supercell-random`). v2 stays the default until the Mac has assembled and reviewed v2.1.

## Why
Library v2 (0.1 GHz nodes) passes the R4 gate (p95), but every large error sits in the resonance band 10.4-10.9 GHz (draw-B maximum 19.3 deg at 10.75 GHz, theta 12.7, phi 101; 5 deg at 10.45 GHz). The mesh correction there changes from 6 deg (10.0 GHz) to 37-47 deg (10.5 GHz) to 0.5 deg (11.0 GHz) and is sampled every 0.5 GHz, so the **correction map is the likely error source**, together with the 0.1 GHz table nodes. Both get 0.05 GHz nodes in 10.4-10.9 GHz.

## Stages (order fixed: a cap leaves a complete mesh-90 set)
| stage | content | tasks | mesh |
|---|---|---|---|
| **A grid** | table nodes at 10.45, 10.55, 10.65, 10.75, 10.85 GHz (the 0.05 GHz steps not on the 0.1 GHz grid), theta 0-45 step 5, phi 0-345 step 15, four states, F = 2n | 5 x 217 = 1085 | 60 |
| **B band references** | 30 held-out points drawn only inside 10.4-10.9 GHz (seed 20261025), four states, direct-solve references saved | 30 | 60 |
| **C1 corr90** | correction nodes at **every 0.05 GHz step in 10.4-10.9 GHz that has none yet**: 10.40, 10.45, 10.55, 10.60, 10.65, 10.70, 10.75, 10.80, 10.85, 10.90 GHz (10.50 exists from R4), the 5 R4 correction points each: (theta, phi) = (0,0), (30,0), (30,90), (45,0), (45,90) | 10 x 5 = 50 | 90 |
| **C2 corr120** | the same 50 points | 50 | 120 |

The mesh-60 correction values at the nodes are table points (the new grid tasks for 10.45 ... 10.85, the existing v2 grid chunks for 10.40/10.60/10.70/10.80/10.90): no extra mesh-60 tasks.

## Commands (PowerShell; the whole-job cap is written to `job_start.txt` by the first run, so `--hours 8` means 8 h for the entire job)
```
cd C:\lightweight_mom ; git fetch hub ; git checkout -b compute/r4-band-refine dc2e664dc19ae8b36d25eec7bc209df62e731493 ; .venv\Scripts\activate
$env:OPENBLAS_NUM_THREADS="1"; $env:OMP_NUM_THREADS="1"; $env:MKL_NUM_THREADS="1"
$D="C:\hmsa_sync\data\r4-band-refine"
python scripts\research\phase20_r4_band_refine_run.py run $D --stage grid --workers 40 --hours 8                      # A
$env:PHASE20_R4_LIB="results\phase20_r4\cp2_library_v2_fullgrid.npz"; $env:PHASE20_R4_WORKERS="40"
python scripts\research\phase20_r4_validate.py a 30 results\phase20_r4\validation_v2_band.json --f-range 10.4,10.9 --seed 20261025 --save-refs results\phase20_r4\validation_v2_band_refs.npz     # B
python scripts\research\phase20_r4_band_refine_run.py run $D --stage corr90 --workers 40 --hours 8                    # C1
python scripts\research\phase20_r4_band_refine_run.py run $D --stage corr120 --workers 16 --hours 8                   # C2 (16 workers x ~9 GB = 144 GB)
python scripts\research\phase20_r4_band_refine_run.py consolidate $D
```
Push the small files (validation JSON and references, manifest, provenance, run log) to `compute/r4-band-refine` after B and again after consolidation; the data stay on the share. `python ... status $D` shows the per-stage progress.

## Estimate (measured server timings: mesh 60 89 s alone / ~137 s with 40 workers, mesh 90 692 s alone, mesh 120 ~3500 s extrapolated, up to ~4500 s under load)
* **A:** 1085 x ~137 s / 40 = **~1.0-1.3 h**.
* **B:** 30 x ~90 s / 40 = **~5 min**.
* **C1:** 50 tasks, 40 workers, ~700-1000 s each = 2 waves = **~25-35 min**; ~2 GB per worker.
* **C2:** 50 tasks, 16 workers = 4 waves (16, 16, 16, 2) x ~3500-4500 s = **~3.9-5 h**, 144 GB RAM. (Not requested, but 17 workers would make it 3 waves, ~2.9-3.8 h, 153 GB; the file says 16.)
* **Total ~5.6-6.9 h; cap 8 h for the whole job.** The driver does not start a mesh-120 task that cannot finish before the cap (5400 s margin), a mesh-90 one with less than 2400 s left, a mesh-60 one with less than 300 s left.

## Resources, caps, stop conditions
* A, B, C1: 40 workers x 1 BLAS thread, < 40 GB in total. C2: **16 workers, free-memory floor 30 GB** (stop and report if free memory falls below it).
* **Hard cap 8 h for the whole job** (checkpointed and resumable: per-task files, a rerun skips finished tasks). If the cap is hit in C2, the mesh-90 set is complete and the finished mesh-120 nodes are kept; report which nodes are complete.
* A task error is printed and the others continue; report the errors, do not patch code.

## Outputs
`C:\hmsa_sync\data\r4-band-refine\` (`tasks\`, `chunks\grid_f*.npz`, `corr90_f*.npz`, `corr120_f*.npz`, `MANIFEST.json` with sha256 and the list of complete chunks, `provenance_*.json`, `job_start.txt`); in git: manifest, provenance, run log, `validation_v2_band.json`, `validation_v2_band_refs.npz`.

## What the Mac does next (no new solves)
Fetch, verify sha256, copy the 5 new grid chunks next to the 41 v2 ones, assemble **v2.1** with `phase20_r4_assemble.py` (`--f-list` = the 46 table frequencies, `--corr-f-list` = the 19 correction frequencies = the 9 R4 ones + the 10 new; mesh-60 correction chunks are cut from the grid chunks) and, for the attribution, **v2.1-oldmap** (the 46 table frequencies with the OLD 9-node correction map). Evaluate v2, v2.1-oldmap and v2.1 on the same points from saved references (`phase20_r4_validate.py --use-refs`): draw A, draw B and the 30 band points; report the gate and the maximum error inside the band, and which part of the in-band error comes from the grid (v2 -> v2.1-oldmap) and which from the correction map (v2.1-oldmap -> v2.1); in particular whether the 19.3 deg point of draw B was grid or correction. v2.1 is not made the default without the user's yes.
