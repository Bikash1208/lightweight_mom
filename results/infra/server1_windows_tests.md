# Server 1 (Windows) test findings after the `.gitattributes` re-checkout

**Server-authored, pending Mac review.** Machine AEM-LAB1 (Windows Server 2019, 48 logical CPUs, 255.6 GB), Python 3.13.5, numpy 2.5.4, OpenBLAS 0.3.34, `OPENBLAS_NUM_THREADS=1`.
Branch `compute/windows-tests` = `hub/phase20-ris` at `7852dc2` (contains `.gitattributes` from 00528a4). Re-checkout done as instructed: `git rm --cached -r . ; git reset --hard hub/phase20-ris` on a clean tree.
No code was changed. `core.autocrlf=true` (system gitconfig) is still set on the server; with the new attributes (`* text=auto eol=lf`) checked-out text files are LF (`git ls-files --eol`: `hfss_build_board.py`, `predictions.json`, `r1_theta0_reference.json` are `i/lf w/lf`; `layout.csv` is `-text`, CRLF as committed).

## 1. Byte-compare tests, rerun (the 9 tests that failed before)

Result: **8 failed, 5 passed** in the 9-test set (13 test items; `tests/test_diagnostics.py` counted whole). The attributes did **not** make the byte-compare tests pass; the earlier hypothesis (CRLF in the checked-out *reference* files) was only half of it. With LF references the remaining difference is on the **generated** side:

| Test | Now | Cause |
|---|---|---|
| `test_tensor_cp_service::test_app_export_of_the_single_hexagon_matches_the_committed_phase26_folder` | FAIL | generated `hfss_build_board.py` has CRLF (`b'...\r\n'`), reference LF: the exporter writes text with `Path.write_text(...)` (no `newline="\n"`), so Windows translates `\n` to `\r\n` |
| `test_app_phase27_tensor_cp::test_single_hexagon_design_export_compare_end_to_end` | FAIL | same `\r` vs `\n` difference (`At index 23 diff: b'\r' != b'\n'`) |
| `test_sparameters_oblique::test_normal_incidence_results_are_bit_identical_to_the_phase20_base_commit` | FAIL | BLAS last-digit difference (expected, ~1e-15), not line endings |
| `test_tensor_cp_design::test_single_hexagon_cells_match_committed_layout_bit_identically` | FAIL | BLAS last-digit difference (30.757828525351485 vs ...282), not line endings |
| `test_diagnostics::test_ru_maxrss_kb_normalizes_macos_bytes_but_not_linux_kilobytes` | **PASS now** | fixed by the Mac's guarded `resource` import |
| `test_diagnostics::test_instrumented_solve_returns_result_and_diagnostics` | FAIL | asserts a memory value `is not None`; on Windows it is `None` (the intended `NOT_AVAILABLE`), so the test needs a platform skip/condition |
| 3 x `board_export` tests (see 2) | FAIL | cp1252 encode error, below |

## 2. One `board_export` test with the job's own error text

`tests/test_app_phase21_multiport.py::test_export_board_endpoint_feeds_files_pins_and_png` (same for `..._writes_exact_step_script_csv_and_readme` and `test_app_phase23_hex4::test_sidelobe_endpoint_end_to_end`). pytest hides the job dict (`'progress': 0.0, ...`); the probe `results/infra/board_export_error_probe.py` (a scratch test that replays the same API calls and prints the job JSON, run as `python -m pytest -c pytest.ini --rootdir . -s <probe>`) shows:

```
"kind": "board_export", "status": "failed", "progress_message": "Writing the exact STEP...",
"error": {"code": "SOLVE_FAILED",
          "message": "'charmap' codec can't encode character 'θ' in position 816: character maps to <undefined>"}
```

Root cause: text files written with the locale encoding (cp1252 on Windows) while the content contains non-ASCII characters (here `θ`). Candidate sites (read-only grep, not changed): `src/export/exact_board.py` lines 141, 269, 295, 334 (`open(path, "w")` / `Path.write_text(...)` without `encoding="utf-8"`), `src/export/array_export.py` 157, 200, 300, 311, `src/export/tensor_cp_export.py` (many `write_text`). The same family (no explicit `encoding="utf-8"` and no `newline="\n"`) explains the CRLF byte-compare failures in section 1. A portable fix would pass `encoding="utf-8", newline="\n"` at those sites; it is output-identical on the Mac. Not done on the server (rule: no code changes).

---

## Update 2026-10-24 (server-authored, pending Mac review)

Threads: `OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=1`. Checkout: detached `phase20-ris` commits named below (hub). No code changed.

### A. The 9 previously failing tests, rerun at `1a6b3d3` (with the `src/export` fix) — `server1_rerun9_1a6b3d3.log`
**7 failed, 6 passed** (not just the 2 expected BLAS diffs). `test_diagnostics` now passes entirely.

| Test | Cause now |
|---|---|
| `test_sparameters_oblique::..._bit_identical_to_the_phase20_base_commit` | BLAS last digit (expected) |
| `test_tensor_cp_design::..._match_committed_layout_bit_identically` | BLAS last digit (expected) |
| `test_tensor_cp_service::test_app_export_of_the_single_hexagon_matches_the_committed_phase26_folder` | line endings fixed; remaining difference is numeric (`predictions.json` differs at one digit, `b'8' != b'7'` at index 488): last-digit BLAS effect in the predicted numbers, but it is a byte-compare so it fails |
| `test_app_phase27_tensor_cp::test_single_hexagon_design_export_compare_end_to_end` | numeric last digits in `aperture_pattern` / `realized_gain_rhcp_dbi` (e.g. 28.47791131215888 vs 28.477911312158874): BLAS |
| 3 x `board_export` (`test_app_phase21_multiport` x2, `test_app_phase23_hex4::test_sidelobe_endpoint_end_to_end`) | **still failing, same error as before**: job `SOLVE_FAILED: 'charmap' codec can't encode character '\u03b8' in position 816` (`server1_board_export_job_error_1a6b3d3.log`). The fix did not reach `src/app/routers/export.py` lines 297 and 303: `script_path.write_text(generate_aedt_script(...))` and `readme_path.write_text(build_board_readme(...))` have no `encoding="utf-8"` (the README text from `src/export/board_readme.py` contains `θ`). `src/export/exact_board.py` and `tensor_cp_export.py` are fixed. |

### B. Fast suite on `26b2433`, `pytest -m "not slow" -n 24` — `server1_fast_suite_26b2433_n24.log`
**7 failed, 2757 passed, 129 skipped in 16 min 53 s** (`-n 24` with 1 BLAS thread per worker; the same 7 tests as section A, nothing else, no failures that only occur under `-n`):
`test_sparameters_oblique::test_normal_incidence_results_are_bit_identical_to_the_phase20_base_commit`, `test_tensor_cp_design::test_single_hexagon_cells_match_committed_layout_bit_identically`, `test_tensor_cp_service::test_app_export_of_the_single_hexagon_matches_the_committed_phase26_folder`, `test_app_phase27_tensor_cp::test_single_hexagon_design_export_compare_end_to_end`, `test_app_phase21_multiport::test_export_board_endpoint_writes_exact_step_script_csv_and_readme`, `test_app_phase21_multiport::test_export_board_endpoint_feeds_files_pins_and_png`, `test_app_phase23_hex4::test_sidelobe_endpoint_end_to_end`.
Expected: the 2 BLAS diffs. Not expected: the other 5. (`-n 24` is fast with `OPENBLAS_NUM_THREADS=1`: 17 min against 30-45 min serial; the 4.7x slowdown earlier was the default BLAS threads.)

### C. For the Mac
1. `src/app/routers/export.py:297,303`: add `encoding="utf-8", newline="\n"` (the three `board_export` tests then reach the next step; unverified beyond that).
2. The two numeric byte-compare tests (`test_tensor_cp_service` folder compare, `test_app_phase27` end-to-end) compare floating-point output exactly; on a different CPU/BLAS they cannot pass bit-for-bit. Consider a tolerance on the Windows server, or mark them Mac-only like the other bit-identity gates (the infra brief already says bit-identity gates stay Mac-only).
