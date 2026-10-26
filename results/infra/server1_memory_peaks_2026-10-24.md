# Server 1 measured memory (peak working set per process) — server-authored, pending Mac review

Method: `Process.PeakWorkingSet64` / `WorkingSet64` of every `python.exe` on AEM-LAB1 (Windows' own counters, same as Task Manager's "peak working set"), read at 2026-10-24 12:33 (job started 12:14:39), plus single-task measurements made earlier with `K32GetProcessMemoryInfo`. 255.6 GB RAM; whole machine in use at 12:33: 24.5 GB.

## Running jobs (a snapshot, not final peaks)
| Job | Process | Working set now | Peak so far | Note |
|---|---|---|---|---|
| `r5-supercell-random` (n_per_cell 48, F = 2n, 9 workers, 19 min CPU each) | 1 worker | **11.92 GB** | **11.92 GB** | this worker has allocated its dense matrix |
| | 8 workers | 0.16-0.21 GB | 0.26 GB | still in the assembly phase before their matrix is allocated; their peak will rise (expect about the same ~12 GB as the first worker, i.e. up to ~108 GB for 9) |
| `r5-supercell-gradients` | pytest (fast suite, first step of the job) | 0.78 GB | 5.1 GB | the supercell tasks have not started yet; **no task peak exists for gradients** |

So the peaks are **not** all <= 5 GB: the first random worker is at 11.9 GB (a little above the job file's ~10 GB), and the other eight random workers had not reached their peak. Because of this the gradients real run was NOT changed to 4 threads per worker (the condition was peaks <= 5 GB, and for gradients there is no measurement yet).

## Earlier single-task measurements (`results/infra/time_f2n.py`, one task = four states, 1 thread)
| Task | Floquet | Peak working set | Time |
|---|---|---|---|
| mesh 60 | F = 120 (2n) | 0.43 GB | 89 s |
| mesh 90 | F = 180 (2n) | 1.91 GB | 692 s |
| mesh 120 | F = 61 (rank minimum, not 2n) | 5.81 GB | 492 s |
| mesh 120 | F = 240 (2n) | **NOT MEASURED** (extrapolating the mesh-60 -> 90 growth, ~6 GB; R4's notes say ~9 GB) | ~3500 s (extrapolated, not measured) |

For `r4-band-refine` worker counts: mesh 90 needs ~2 GB per worker (40 workers = ~76 GB); mesh 120 at F = 2n has no measured peak yet, so the 9 GB figure is an estimate. One mesh-120 F = 2n task on one core would measure it (about an hour); say if you want that run.
