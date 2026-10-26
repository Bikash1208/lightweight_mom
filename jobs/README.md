# `jobs/` — Server Compute Job Requests and Records

Job request specifications and execution records for compute-intensive tasks offloaded to remote Windows compute servers. The Mac functions strictly as coordinator; heavy MoM sweeps, full-grid libraries, and supercell verifications run on the servers via a shared Git remote (`hub`).

Governed by the infrastructure protocol in `docs/briefs/INFRA_COMPUTE_SERVER_SYNC.md`.

---

## 1. Job Lifecycle

```
Mac (Coordinator)             Git Hub (remote 'hub')          Windows Server 1 (AEM-LAB1)
─────────────────             ─────────────────────          ───────────────────────────
1. Writes jobs/<job>.md
   (base commit, tasks)  ───>  git push hub phase20-ris ───>  git fetch hub
                                                              git checkout -b compute/<job> <base>
                                                              Runs job script (multi-worker)
                                                              Records provenance & manifest
                         <───  git push hub compute/<job> <──
2. Fetches compute/<job>
   Verifies hashes / log
   Merges into phase20-ris
```

### Working Rules
- **Mac writes, never runs heavy compute:** The Mac coordinates branches, reviews output hashes, and merges.
- **Server creates branch:** Jobs run on `compute/<job>` branched from the specified base commit. Servers **never** commit directly to `phase20-ris`.
- **Environment standards:** Server tasks run with single-threaded BLAS per worker (`OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`) and write checkpointed `.npz` chunks accompanied by JSON manifests and provenance logs.
- **Job files stay in place:** Job files are permanent records and remain in `jobs/` even after completion or withdrawal.

---

## 2. Jobs Status Matrix

| Job File | Description | Target Machine / Base Commit | Branch on Hub (`git branch -a`) | Status | Notes |
|---|---|---|---|---|---|
| `equivalence-check.md` | Initial cross-machine numerical equivalence check on CP-2 baseline cell (mesh 60, F = 2n = 120, 10.5 GHz). | Server 1<br>Base: `13f099f` | `remotes/hub/compute/equivalence-check` (`7b9a41c`) | **PASSED** (`INTERCHANGEABLE (max rel. 4.3e-10)`) | Quoting `docs/STATUS.md`: "cross-machine equivalence INTERCHANGEABLE (max rel. diff 4.3e-10, 1 and 8 threads, 2026-10-23)". Recorded in Mac `docs/STATUS.md` and commit `f5fda93`. |
| `r3b-continue.md` | Proposed continuation of R3b optimiser with 40 workers and 3 Nelder-Mead starts. | Server 1<br>Base: `c11abee` | *(None)* | **WITHDRAWN** | Withdrawn by user decision (2026-10-23, commit `ba12f3b`). Baseline CP-2 cell retained; see `results/phase20_r3b/OUTCOME.md`. |
| `r4-full-grid-server.md` | Full-grid state library build: 8942 tasks across 41 frequencies, 10 polar angles, and 24 azimuth angles. | Server 1 (48 cores)<br>Base: `7852dc2` | `remotes/hub/compute/r4-full-grid` (`6990d7e`) | **DONE & MERGED** | All 8942 tasks completed in 8.52 h (0 errors); merged into `phase20-ris` at commit `7ddc843`. Assembled into Library v2 (`1264bdc`). |
| `r4-fullgrid-validate.md` | Held-out validation gate for Library v2 (50+ random parameter draws vs direct MoM solves). | Server 1<br>Base: `7d4cc20` | `remotes/hub/compute/r4-fullgrid-validate` (`ed76d73`) | **DONE** | Completed on Server 1; manifest and provenance committed on branch `compute/r4-fullgrid-validate` at `ed76d73`. |
| `r4-mesh120-corrections.md` | Independent server re-solve of 45 mesh-120 correction points at F = 240. | Server 1<br>Base: `7d4cc20` | *(None)* | **STANDBY (DO NOT RUN)** | User confirmed: archival standby only, not to be run. The 45 mesh-120 points already exist in Mac chunks (`corr120_f*.npz`) and on the hub. |
| `r5-supercell-gradients.md` | Array-level supercell MoM verification of designed phase-gradient codings (11 tasks, x/y gradients, F = 2n). | Server 1 (≤11 workers)<br>Base: `7d4cc20` | `compute/r5-supercell-gradients` *(pending)* | **QUEUED** | Queued behind `r5-tolerance-mom`. Tests local-periodicity assumption of the R5-3 plane-wave forward model. |
| `r5-supercell-random.md` | Array-level supercell MoM verification of non-periodic and dithered 4×4 tiles (9 tasks, 33 codings, F = 2n). | Server 1 (≤9 workers)<br>Base: `7d4cc20` | `compute/r5-supercell-random` *(pending)* | **QUEUED** | Queued behind `r5-tolerance-mom`. May run concurrently with `r5-supercell-gradients`. |
| `r5-tolerance-mom.md` | MoM prediction of diode capacitance tolerance cases ($C_{\text{off}} = 0.0225, 0.025, 0.030\text{ pF}$, meshes 60 & 90, F = 2n, 164 tasks). | Server 1<br>Base: `9d645ce` | `remotes/hub/compute/r5-tolerance-mom` (`7eb06a8`) | **DONE & MERGED** | Merged into `phase20-ris` at commit `9d3deee`. Review (`14ac74d`) confirmed MoM mesh 90 reproduces HFSS sensitivity; mesh 60 does not. |
