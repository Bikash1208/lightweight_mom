# Job request: r4-fullgrid-validate (R4 item 4 gate for library v2, server 1)

**From:** Mac session (coordinator). **Branch on the server:** `compute/r4-fullgrid-validate` from the base commit. Never commit to `phase20-ris`. **Base commit:** `7d4cc209271ae05a9fe43074c92518d19f00a11e` (fetch from `hub`).

## Purpose
The R4 item 4 gate for the full-grid library `results/phase20_r4/cp2_library_v2_fullgrid.npz` (in the base commit): >= 50 held-out off-grid points solved directly (mesh 60, F = 2n, four states), interpolation error of the RAW table (`c_raw`) against them, 95th percentile <= 3 deg / 0.2 dB inside theta <= 30 deg. The Mac cannot give this number itself (the direct solves are heavy, the 2026-10-23 reference solutions were not stored), so it is a server job. Two draws, both saving their references:
* **A** the original draw: 60 points, seed 20261021, phi 0-360 (the same points as the v1 gate, so v1 and v2 are compared on identical points);
* **B** a phi = 90 focus: 40 points, seed 20261024, phi 75-105 deg, the plane the user asked about.

## Commands (PowerShell)
```
cd C:\lightweight_mom ; git fetch hub ; git checkout -b compute/r4-fullgrid-validate 7d4cc209271ae05a9fe43074c92518d19f00a11e ; .venv\Scripts\activate
$env:OPENBLAS_NUM_THREADS="1"; $env:OMP_NUM_THREADS="1"; $env:MKL_NUM_THREADS="1"
$env:PHASE20_R4_LIB="results\phase20_r4\cp2_library_v2_fullgrid.npz"; $env:PHASE20_R4_WORKERS="40"
python scripts\research\phase20_r4_validate.py a 60 results\phase20_r4\validation_v2_a.json --save-refs results\phase20_r4\validation_v2_a_refs.npz
python scripts\research\phase20_r4_validate.py a 40 results\phase20_r4\validation_v2_a_phi90.json --phi-range 75,105 --seed 20261024 --save-refs results\phase20_r4\validation_v2_a_phi90_refs.npz
```
Also run the same two commands with `$env:PHASE20_R4_LIB="results\phase20_r4\cp2_library_v1.npz"` (outputs `validation_v1_*`), so the old library is evaluated on the same points.

## Cost
100 points x 4 states, about 90 s each at mesh 60, F = 2n (one worker): 100 x 90 s / 40 = about 4 min per run, 4 runs about 20 min; ~0.5 GB per worker.

## Outputs (commit on the compute branch, push to `hub`)
The four result JSONs and four reference NPZs (small), manifest `results/phase20_r4/validation_v2_manifest.json` (file, size, sha256, producing commit, machine).

## Stop conditions
Any error -> stop and report, do not patch. Wall clock > 1 h -> stop and ask.
