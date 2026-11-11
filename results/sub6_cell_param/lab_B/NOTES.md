# Lab B run notes

- Source commit: `d38113380c5338c8fe493370c975328e71279d3e`.
- Local AEDT returned `2025.2` (HFSS profile reports `2025.2.0`); PyAEDT `1.9.0`, Python `3.12.16`.
- VERIFY status: PASS against user targets and the v251 reference (0.5% frequency tolerance).

| VERIFY row | Lab f (GHz) | Lab Q | Reference f (GHz) | Solve start | Solve stop | Adaptive passes | Final solved tets | Raw eigenvalue JSON (local time) |
|---|---:|---:|---:|---|---|---:|---:|---|
| OFF along | 5.774606449 | 154.7665 | 5.77461 | 2026-11-10 23:28:26 | 2026-11-10 23:32:19 | 12 | 24,772 | `verify/VERIFY_OFF_along.json` (23:32:44) |
| OFF across | 6.203002303 | 193.6578 | 6.20300 | 2026-11-10 23:33:02 | 2026-11-10 23:36:19 | 10 | 25,877 | `verify/VERIFY_OFF_across.json` (23:36:34) |
| ON along | 4.866152529 | 130.4262 | 4.86615 | 2026-11-10 23:36:53 | 2026-11-10 23:43:27 | 15 | 54,945 | `verify/VERIFY_ON_along.json` (23:43:55) |
| ON across | 4.931428809 | 87.1711 | 4.93143 | 2026-11-10 23:44:13 | 2026-11-10 23:49:47 | 13 | 49,985 | `verify/VERIFY_ON_across.json` (23:50:01) |

Each `verify/VERIFY_*.prof` is the AEDT exported HFSS profile for that row and records the start, stop, adaptive passes, and tet count. Each corresponding JSON stores eight modes with `raw_frequency`, `raw_unit`, frequency and Q read by PyAEDT from `EigenCell : LastAdaptive`. AEDT's native eigenvalue data for the last row is in the local, uncommitted file `C:\lwmom_labB\hfss_exports\sub6_cell_param\lab_B_work\verify_project\CellParam_CAL_NOMESH.aedtresults\CellParam.results\DV234_SOL136_M0_V139.sd`, last modified 2026-11-10 23:49:46 (402 bytes). The runner calls `cleanup_solution()` between rows, so earlier native `.sd` files are not retained; the per-row raw exports and profiles are retained. No AEDT project or native result file is committed.

Provenance: `cell_param_pyaedt.py` calls `solve_row()` for each VERIFY row. `solve_row()` saves and runs AEDT, then calls `h.post.get_solution_data()` for Eigen Modes and Eigen Q from `EigenCell : LastAdaptive`; it writes the per-row JSON. `verify()` selects from those returned modes and writes `verify.json`. Only after that does `check_verify_vs_reference.py` read `verify_reference/verify_reference.json` to compare. The reference is not a source of lab result values.

AEDT v252 validation warned `Boundary 'Boundary': unexpected unit type on Resistance` on all three BAR90 builds. Readback in a separate new local AEDT 2025.2 instance, using a copy of the lab-built project, showed `BAR90_02LRH` stores `Resistance='Rs'`, `Reactance='Xs'`; OFF variables are `Rs='22.23ohm'`, `Xs='-335.4ohm'`; ON variables are `Rs='1.6ohm'`, `Xs='29.2ohm'`. The absorbing floor stores `Resistance='377ohm'`, `Reactance='0ohm'`. These are the intended values, so the warning is treated as v252 unit labelling under the user's 2026-11-10 instruction. See `boundary_readback.json`.

All three BAR90 projects built with status `BUILT`, floor `IMPEDANCE_377_OHM_PER_SQUARE`, and no mesh operations. Before QUEUE_B, measured CPU busy was 12.9% of 48 logical CPUs (~41.8 free). Initial cap: 3 solver processes x 8 cores, with at least 8 logical CPUs reserved; the launcher remeasures before each 20-row batch. Estimate from the documented 3.5–5 min per row: about 8–11 hours plus batch start-up, subject to shared-machine load.

Runner fixes on this lab-b branch:
1. Create `lab_B_work` before writing `lab_B_machine.json`; detect the explicitly requested v252 install.
2. Batch QUEUE_B in 20-row groups for progress pushes, retaining original CSV stems for resume and read-back checks.
3. Reserve logical CPU headroom and remeasure before each batch and retry; preserve cumulative parallel reports.
4. Copy per-row raw VERIFY exports, AEDT profiles, and build summaries into the result package for provenance.

## Final QUEUE_B result (2026-11-11 08:09 IST)

- Completed all 336 planned rows: `SCREEN_2.csv` 144/144, `SCREEN_6_short_nostub.csv` 96/96, and `SCREEN_6_short_trace.csv` 96/96. FAILED: 0; failure flags: none.
- Independent final audit checked all 336 expected run IDs, eight modes per row, project variants, build floors, and `readback_si` against intended variables. Result: PASS.
- Result batches were pushed about every 20 rows to `compute/sub6-cell-param-lab-b` under `results/sub6_cell_param/lab_B/`. Actual layouts: full 3x8, 4x8, or 5x8 according to free CPU; no-stub and trace 5x8. CPU was remeasured before each batch.
- A Codex server restart interrupted the 101–120 batch after five private row results. No queue Python or AEDT process from that batch remained. The five rows passed the same readback check and were recovered into the owned output directory before resume. The remaining 15 rows were solved on resume. No row was counted twice in the final set.
