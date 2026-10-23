# Overnight summary, 2026-10-23 -> 24 (server 1, AEM-LAB1) — server-authored, pending Mac review

Rules kept: nothing committed to `phase20-ris`; no `src/` change; everything on `compute/*` branches pushed to `hub`. Commits use a one-off identity (`compute-server-1 (AEM-LAB1)`), no git config changed.

| # | Item | Result | Branch @ commit |
|---|---|---|---|
| 1 | `jobs/r5-tolerance-mom.md` | **Done.** 164 tasks (41 f x 2 theta x 2 meshes, C_off 0.0225/0.025/0.030 pF, F = 2n), 40 workers x 1 thread, 73.7 min. `results/phase20_r5/tolerance_mom.npz` (`c_structure` shape (2,3,41,2,4,2,2), all finite, provenance inside) + manifest. Fast suite at the base: 9 failed (2 known BLAS last-digit diffs + 7 Windows items, same as below), not patched. | `compute/r5-tolerance-mom` @ `7eb06a8` |
| 2 | Windows tests after `.gitattributes` (00528a4) | **Report written:** `results/infra/server1_windows_tests.md`. 8 of the 9 previously failing tests still fail: the byte-compare tests now fail because the *generated* files get CRLF on Windows (`Path.write_text` without `newline`), the three `board_export` tests fail with the job error `SOLVE_FAILED: 'charmap' codec can't encode character 'θ'` (writes without `encoding="utf-8"` in `src/export/`), 2 BLAS last-digit diffs, 1 `test_diagnostics` assertion that expects a memory value (None on Windows by design). `test_diagnostics::...maxrss...` now passes. | `compute/windows-tests` @ `8f75f9e` |
| 3 | R4 full grid | **Grid options added (script only, defaults unchanged, 4 tests); the run was NOT started: estimate 10-12 h > 6 h stop condition.** At F = 2n one mesh-60 task = 89 s single thread (not 12 s as in R4 notes), mesh 90 = 692 s, mesh 120 ~3500 s (extrapolated); 8897 grid tasks + the unchanged correction set ~ 985k core-s = 6.8 h at perfect 40-worker scaling, 10-12 h realistically. Details, planned commands and options (coarser steps, second server) in `jobs/r4-full-grid-server.md`. | `compute/r4-full-grid` (this commit) |

## Decisions for the Mac in the morning
1. R4 full grid: accept ~10-12 h, or choose a coarser grid (e.g. `--f-step 0.2`), or split over two servers.
2. Windows portability of `src/export/*` text writes (`encoding="utf-8", newline="\n"`; output-identical on the Mac) and a platform condition for the one `test_diagnostics` memory assertion — see `server1_windows_tests.md`.
3. The server is idle and ready for `r5-supercell-gradients` / `r5-supercell-random` (job requests are on the hub; not started: the instructions for tonight did not include them).
