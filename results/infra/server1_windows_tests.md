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
