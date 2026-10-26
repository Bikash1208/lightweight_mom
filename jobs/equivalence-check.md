# Job request: equivalence-check (infra brief rule 5)

**From:** Mac session (coordinator). **To:** server 1 (Windows, 48 cores). **Branch to create on the server:** `compute/equivalence-check` from the base commit below. **Do not** commit to `phase20-ris`.

**Base commit:** `13f099f0aad699d01c53a0de08468aa7d4f50a66` on `phase20-ris` (fetch it from `hub`).

## Purpose
Show once, before any server result is used, that the same small fixed case gives interchangeable numbers on the Mac and on the server (CPU / BLAS differences only). The case: the baseline CP-2 cell (`results/phase20_r3/opt/cp2_madp_final.json`),
MADP nominal diodes, **mesh 60, Floquet cutoff F = 2n = 120**, the four one-hot states, f = 10.5 GHz, (theta, phi) = (0, 0) and (30, 0). The Mac reference is `results/infra/equivalence_mac.json` (in the base commit).

## Setup (first time on the server)
```
cd C:\lightweight_mom
git remote add hub C:\hmsa_sync\lightweight_mom.git      (skip if present)
git fetch hub
git checkout -b compute/equivalence-check 13f099f0aad699d01c53a0de08468aa7d4f50a66
.venv\Scripts\activate
```

## Commands (PowerShell; run in this order)
1. **Fast suite once** (portability check; report, do not fix): `python -m pytest -m "not slow" -q -x --maxfail=50 > results\infra\server1_fast_suite.txt`
   Expected: green apart from Windows-only failures. Report the list of failing tests with the first traceback line of each; **do not change code**.
2. **Equivalence case, 1 thread** (workers = 1):
   `$env:OMP_NUM_THREADS="1"; $env:MKL_NUM_THREADS="1"; $env:OPENBLAS_NUM_THREADS="1"; python scripts\research\equivalence_case.py results\infra\equivalence_server1_t1.json`
3. **Equivalence case, 8 threads** (workers = 1; shows whether the thread count changes the digits):
   `$env:OMP_NUM_THREADS="8"; $env:MKL_NUM_THREADS="8"; $env:OPENBLAS_NUM_THREADS="8"; python scripts\research\equivalence_case.py results\infra\equivalence_server1_t8.json`
4. **Compare with the Mac reference** (and the two server runs with each other):
   `python scripts\research\equivalence_compare.py results\infra\equivalence_mac.json results\infra\equivalence_server1_t1.json > results\infra\compare_mac_vs_t1.json`
   `python scripts\research\equivalence_compare.py results\infra\equivalence_mac.json results\infra\equivalence_server1_t8.json > results\infra\compare_mac_vs_t8.json`
   `python scripts\research\equivalence_compare.py results\infra\equivalence_server1_t1.json results\infra\equivalence_server1_t8.json > results\infra\compare_t1_vs_t8.json`
5. **Phase 23 snapshot (informational only; the bit-identity gate stays Mac-only):** `python scripts\research\phase23_rect_snapshot.py results\infra\snapshot_server1.json`, then report whether
   `fc /b results\infra\snapshot_server1.json results\phase23\rect_snapshot_after_C.json` says identical or different (hashes: a difference is expected to be possible and is not a failure by itself).

## Expected outputs (commit on `compute/equivalence-check`, push to `hub`)
`results/infra/equivalence_server1_t1.json`, `equivalence_server1_t8.json`, the three `compare_*.json`, `server1_fast_suite.txt`, `snapshot_server1.json` (small files only), each case file carrying its provenance record
(machine, OS, commit, dirty flag, Python / numpy / scipy, BLAS, threads, workers, start / end). Also report in the final message: the wall-clock of step 2 and step 3 (how much faster 8 threads are) and the machine RAM.

## Stop conditions
* `max_rel_complex_diff` > 1e-6 in any comparison involving the Mac file: **stop, push what you have, report** (verdict `DIFFERENCE_ABOVE_TOLERANCE_STOP_AND_REPORT`).
* Any failure that prevents steps 2-4: stop and report; do not patch the repo.
* Nothing else is to be run for this job; no optimisation, no library builds.

## What the Mac does with the result
Records the differences in `docs/STATUS.md`; if all comparisons are `INTERCHANGEABLE`, server results (with provenance) may be used for heavy jobs; the Mac keeps all bit-identity gates.
