# Job record: r4-full-grid-server (server-authored, pending Mac review)

**Status: NOT RUN. Stopped by the 6 h stop condition (estimate 10-12 h).** Authorised by the user (Mac offline, 2026-10-23 night): add grid options to `scripts/research/phase20_r4_build_library.py` (defaults unchanged), then build the full approved grid with checkpointing, 40 workers x `OPENBLAS_NUM_THREADS=1`, outputs in `C:\hmsa_sync\data\r4-full-grid\` with a manifest; stop and report if the estimate exceeds 6 h.

**Branch:** `compute/r4-full-grid`. **Base commit:** `7852dc2` (`hub/phase20-ris` at the time). **Machine:** AEM-LAB1 (Windows Server 2019, 48 logical CPUs, 255.6 GB, Python 3.13.5, numpy 2.5.4, OpenBLAS 0.3.34).

## The only code change (outside `src/`; nothing in `src/` touched)
`scripts/research/phase20_r4_build_library.py`: new `configure_grid(f_step, theta_step, phi_step, f_list, theta_list, phi_list)` and `parse_grid_options(argv)`; `__main__` accepts optional `--f-step X --theta-step X --phi-step X`, `--f-list a,b,c --theta-list ... --phi-list ...` and `--dry-run` (prints the grid and task count, runs nothing). The old positional arguments (`grid|corr60|corr90|corr120|status`, worker count, `reverse|only:...`) are parsed exactly as before. With no option nothing changes; the module-level `GRID_*` constants are untouched, so `phase20_r4_assemble.py` (which imports them) still sees the R4 grid.
Added `tests/test_phase20_r4_grid_options.py` (4 tests: defaults equal the old constants; `configure_grid()` and the default steps change nothing; the full grid is 41 x 217 points; the option parser keeps positional arguments). `results/infra/time_f2n.py` is the timing probe below.

**Dry-run, defaults vs the old constants** (`python scripts\research\phase20_r4_build_library.py grid --dry-run`): `GRID_F` = 8.0..12.0 step 0.25 (17), `GRID_THETA` = 0..45 step 7.5 (7), `GRID_PHI` = 0..330 step 30 (12), 73 points per f, 1241 tasks; identical (string-equal) with explicit `--f-step 0.25 --theta-step 7.5 --phi-step 30`, and equal bit for bit to the literal old expressions `[8.0+0.25*i ...]`, `[7.5*i ...]`, `[30.0*i ...]`.
**Full approved grid** (`--f-step 0.1 --theta-step 5 --phi-step 15`): 41 f (8.0..12.0), theta 0..45 step 5 (10), phi 0..345 step 15 (24); theta 0 = the phi 0 point only -> 1 + 9 x 24 = 217 points per f, **8897 mesh-60 tasks**.

## Measured cost on this server (one task = all four states, F = 2n, 10.5 GHz, 1 thread, `results/infra/time_f2n.py`)
| mesh | F | theta, phi | seconds | peak RAM |
|---|---|---|---|---|
| 60 | 120 | 30, 0 | 89.5 | 0.43 GB |
| 60 | 120 | 45, 45 | 87.8 | 0.43 GB |
| 90 | 180 | 30, 0 | 692.1 | 1.91 GB |
| 120 | 240 | not measured here | ~3500 (extrapolated with R4's 120/90 ratio 616/121 = 5.1) | ~9 GB (R4 notes) |

(The earlier "mesh 60 = 8.9 s per evaluation" figures in the chat were at the solver's rank-minimum Floquet cutoff, not F = 2n; they do not apply to the library. R4's "12 s at F = 120" was with full multithreaded BLAS.)

## Estimate for what was asked
- grid, 8897 tasks x ~89 s = 792k core-s; correction set as coded (`CORR_F` 9 f x 5 points = 45 tasks per mesh): mesh 60 4k core-s, mesh 90 31k, mesh 120 ~158k (extrapolated). Total ~985k core-s.
- 40 workers, **no contention: 6.8 h**; with the contention measured earlier (40 single-thread workers ran each task ~1.8x slower than a lone process, at the minimum cutoff; not measured at F = 2n): **~10-12 h**. The mesh-120 correction tasks need ~9 GB each, so at most ~20 at once (180 GB).
=> **above the 6 h stop condition; nothing was launched** (no chunks, `C:\hmsa_sync\data\r4-full-grid\` not created, no manifest).

## What it would take (for the Mac to decide)
Options: (a) accept ~10-12 h; (b) coarser steps, e.g. `--f-step 0.2 --theta-step 5 --phi-step 15` (4440 -> ~half the grid cost) or `--phi-step 30`; (c) split over two servers. Planned commands if approved (chunks are written per frequency, atomically, and skipped when valid, so a rerun resumes; the script prints one line per finished frequency):
```
cd C:\lightweight_mom ; git checkout compute/r4-full-grid
$env:OPENBLAS_NUM_THREADS="1"; $env:OMP_NUM_THREADS="1"; $env:MKL_NUM_THREADS="1"
cmd /c mklink /J results\phase20_r4\chunks C:\hmsa_sync\data\r4-full-grid\chunks      (chunks dir is gitignored; the junction puts the outputs on the share without a code change)
python scripts\research\phase20_r4_build_library.py grid 40 --f-step 0.1 --theta-step 5 --phi-step 15
python scripts\research\phase20_r4_build_library.py corr60 40
python scripts\research\phase20_r4_build_library.py corr90 40
python scripts\research\phase20_r4_build_library.py corr120 20        (the correction set uses the existing CORR_F / CORR_POINTS; unchanged)
```
Caveats: `phase20_r4_assemble.py` imports the default `GRID_*`, so assembling a full-grid library needs the same options added there (not done); chunk files `grid_f<f>.npz` carry their `points`, so full-grid chunks never validate as R4 reduced-grid chunks (and vice versa) if both sit in one directory.
