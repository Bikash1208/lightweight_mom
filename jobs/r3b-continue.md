# Job request: r3b-continue (Phase 20 R3b, server 1)

**From:** Mac session (coordinator; no heavy compute on the Mac any more). **Branch to create on the server:** `compute/r3b-continue` from the base commit below. Never commit to `phase20-ris`.
**Base commit:** `75ef08c85c69192381fb08c79acaa1e35e496cfa` on `phase20-ris` (fetch from `hub`). It contains every input: the optimiser and its **249-evaluation cache** (`results/phase20_r3b/opt_cache.jsonl`), the run-1 and run-2 HFSS delta tables,
the run-2 HFSS data, the F = 2n state library, the re-score / verification / package scripts. A copy of the cache and tables with a manifest is in `C:\hmsa_sync\data\r3b\` (sha256 in `MANIFEST.json`; the manifest in git is `results/phase20_r3b/r3b_manifest.json`).

## START CONDITION
The cross-machine equivalence check (`jobs/equivalence-check.md`) has **passed** (all three comparisons `INTERCHANGEABLE`, relative difference <= 1e-6, recorded by the Mac in `docs/STATUS.md`). If it has not, do not start; report.

## Environment (every step)
```
cd C:\lightweight_mom ; git fetch hub ; git checkout -b compute/r3b-continue 75ef08c85c69192381fb08c79acaa1e35e496cfa ; .venv\Scripts\activate
$env:OPENBLAS_NUM_THREADS="1"; $env:OMP_NUM_THREADS="1"; $env:MKL_NUM_THREADS="1"
python -c "import json,sys; sys.path.insert(0,'scripts/research'); import provenance as p; json.dump(p.provenance({'job':'r3b-continue'}, workers=40), open('results/phase20_r3b/r3b_continue_provenance_start.json','w'), indent=1)"
```
40 workers x `OPENBLAS_NUM_THREADS=1` (infra rule 6). Memory: mesh 120 is ~9 GB per worker (<= 20 such workers at once = 180 GB of 255.6 GB); steps 1-2 are small.

## Steps (in order; each is resumable: rerun the same command after an interruption)
0. **Fast suite once** (the Mac no longer runs it): `python -m pytest -m "not slow" -q > results\phase20_r3b\server1_fast_suite.txt`. Expected: green except the 2 known last-digit BLAS diffs. Anything else: report before continuing.
1. **Resume the three starts** from the cache (the Nelder-Mead replay of the 249 cached evaluations is instantaneous; the starts had completed baseline 95, diag_long 86, diag_longer_wide 68 evaluations of maxfev 120):
   `python scripts\research\phase20_r3b_opt.py many 120 40 baseline diag_long diag_longer_wide`
   (three threads, ONE shared pool of 40 single-thread workers; the 8 condition solves of every evaluation of every start run concurrently, up to 24 at a time). Expected: ~1-1.5 h (one evaluation round ~70 s).
   Outputs: `results/phase20_r3b/opt_<start>.json`, `opt_cache.jsonl` (appended, fsync per evaluation). Final report line: evaluations per start and best J per start.
2. **Re-score the finalists with the run-2 delta** (top 5 of each start by the run-1 objective): `$env:R3B_POOL="40"; python scripts\research\phase20_r3b_rescore.py 5`
   Outputs: `rescore.json` (both rankings: run-1 delta used in the optimisation and run-2 delta), `cp2_r3b_final.json` (winner on the run-2 objective), `hfss_delta_table_run2.npz` (rebuilt). **Report both rankings** in the final message.
   **STOP CONDITION A:** if the winner is the unchanged baseline cell (`line_diag` / `width_diag` equal to the axis values within 0.05 mm) or its run-2 objective is within 2 % of the baseline cell's, **stop after this step and report** (do not generate an HFSS package for an unchanged cell); also report the contribution split (nominal / corners / edges / Q) of the winner.
3. **Mesh 60 / 90 / 120 verification of the winner** (F = 2n at every mesh): `python scripts\research\phase20_r3b_verify.py results\phase20_r3b\cp2_r3b_final.json --workers 40`
   Outputs: `verification.json` (provenance inside), `VERIFICATION.md`. Expected ~1-2 h (mesh 120, 20 tasks). Report the table and the Richardson uncertainty.
4. **HFSS package for the new cell**: `python scripts\research\phase20_r3b_cp_package.py --final results\phase20_r3b\cp2_r3b_final.json --workers 40` -> `hfss_exports\phase20_r3b_cp_cell\` (predictions at F = 2n, tolerance cases at 10.5 GHz, README, STEP).
   `cell.step` is gitignored: copy it to `C:\hmsa_sync\data\r3b\phase20_r3b_cp_cell\` (with the rest of the folder) and list it in the manifest.
5. **Finish**: `python -c "..."` as in the environment block with `p.finish(...)` into `r3b_continue_provenance_end.json`; copy `opt_cache.jsonl` and the new JSON/MD files to `C:\hmsa_sync\data\r3b\`; update that folder's `MANIFEST.json`
   and write `results/phase20_r3b/r3b_continue_manifest.json` (file, size, sha256, producing commit, machine) into git. Commit the small files on `compute/r3b-continue` and **push to `hub`**.

## Expected outputs
`results/phase20_r3b/`: `opt_baseline.json`, `opt_diag_long.json`, `opt_diag_longer_wide.json`, extended `opt_cache.jsonl`, `rescore.json`, `cp2_r3b_final.json`, `verification.json`, `VERIFICATION.md`, `server1_fast_suite.txt`, provenance files, manifest;
`hfss_exports/phase20_r3b_cp_cell/` (CSV, README, step_parts.json; STEP via hmsa_sync); `C:\hmsa_sync\data\r3b\`.

## Stop conditions
* Equivalence check not passed -> do not start.
* Any step fails, or the fast suite shows failures other than the 2 known BLAS last-digit diffs -> stop and report; do not patch code.
* STOP CONDITION A (step 2) above.
* Memory: if a mesh-120 worker is killed (out of memory), rerun step 3 with `--workers 16`; if it fails again, stop and report.
* A step that would exceed ~12 h of wall clock in total: stop and ask (the estimate for steps 1-4 is ~4-6 h).

## What the Mac does with the result
Fetches `hub`, reviews `rescore.json` / `VERIFICATION.md`, records the findings in STATUS and the roadmap log, merges or cherry-picks the compute branch, and hands the package `hfss_exports/phase20_r3b_cp_cell/` to the user for the second HFSS run.

## Note from the Mac (what the cache shows before you start)
After 249 evaluations the three starts have made **little progress**: baseline start J 27.4 -> 26.5, diag_long best 49.6, diag_longer_wide best 74.4. The objective is dominated by the tolerance corners (21.6 of 26.5, mostly the theta-30 / phi-45 corners)
and the diagonal-stub variations make it worse, not better. Resuming to maxfev 120 is cheap (~1 h) and is run as asked; expect STOP CONDITION A to trigger. Do **not** change the objective on the server.
