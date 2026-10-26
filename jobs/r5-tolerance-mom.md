# Job request: r5-tolerance-mom (Phase 20 R5-1, server 1)

**From:** Mac session (coordinator). **Branch on the server:** `compute/r5-tolerance-mom` from the base commit. Never commit to `phase20-ris`. **Base commit:** `9d645cee187b15e3853aa7a270b04330f498d745` on `phase20-ris` (fetch from `hub`).

## START CONDITION
The cross-machine equivalence check (`jobs/equivalence-check.md`) has **passed** (all comparisons `INTERCHANGEABLE`, recorded by the Mac in `docs/STATUS.md`). If not, do not start.

## Purpose
MoM prediction of the diode-tolerance cases the user runs in HFSS, written **before** the HFSS comparison: baseline CP-2, MADP diode OFF capacitance C_off = 0.0225 / 0.025 / 0.030 pF, four one-hot states, theta = 0 and 30 deg, phi = 0,
f = 9.5-11.5 GHz in 0.05 GHz steps, meshes 60 and 90, Floquet cutoff F = 2n. 164 independent tasks (41 f x 2 theta x 2 meshes), 12 load sets each.

## Commands (PowerShell)
```
cd C:\lightweight_mom ; git fetch hub ; git checkout -b compute/r5-tolerance-mom 9d645cee187b15e3853aa7a270b04330f498d745 ; .venv\Scripts\activate
$env:OPENBLAS_NUM_THREADS="1"; $env:OMP_NUM_THREADS="1"; $env:MKL_NUM_THREADS="1"
python -m pytest -m "not slow" -q > results\phase20_r5\server1_fast_suite_r5_1.txt          (first step; report failures other than the 2 known BLAS last-digit diffs; do not patch)
python scripts\research\phase20_r5_tolerance_mom.py results\phase20_r5\tolerance_mom.npz --workers 40
```
(create `results\phase20_r5` first). 40 workers x 1 thread (infra rule 6); the mesh-90 tasks need ~3-4 GB each. Expected wall clock: ~30-60 min. The script is not checkpointed per task: if it is interrupted, rerun it (cost ~1 h).

## Outputs (commit on `compute/r5-tolerance-mom`, push to `hub`; the NPZ is ~1-2 MB)
`results/phase20_r5/tolerance_mom.npz` (array `c_structure[mesh, C_off, f, theta, state, 2, 2]`, provenance string inside: machine, OS, commit, dirty flag, versions, BLAS, threads, workers, start/end), `server1_fast_suite_r5_1.txt`, and a manifest `results/phase20_r5/tolerance_mom_manifest.json` (file, size, sha256, producing commit, machine).

## Stop conditions
Equivalence not passed -> do not start. Any failure -> stop and report; do not patch code. Wall clock > 3 h -> stop and ask.

## What the Mac does with the result
Fetches `hub`, then compares with the user's HFSS tolerance CSV (`scripts/research/phase20_r5_tolerance_compare.py`, comparison only) when it arrives.
