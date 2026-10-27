# Phase 28 at 10.0 GHz (PI operating frequency, 2026-10-26)

All results here are `MODEL_ONLY` / `NOT_HFSS_VALIDATED`. The 10.54 GHz files in `results/phase28/` are history and unchanged.

| item | files |
|---|---|
| Supercell table at 10.0 GHz (same method as Phase 26: 4 slots x N = 3, 4, 5, 6, 8 (tier 1) + 2, 10, 12, 16 (tier 2), mesh 24, 6 mm diameter, 6.4 mm pitch) | `task_*.json` (36 checkpoints), `manifest.json`, `supercell_table_*.json` (assembled, sha256 in `feasibility_results.json` `table.files`), `orchestrator*.log` |
| Item-2 sweep rerun, same code (`PHASE28_FREQ_GHZ=10.0`, `PHASE28_TABLE_GLOB`) | `feasibility_results.json` (+ `.provenance.json`), `FEASIBILITY_REPORT.md` |
| 10.0 vs 10.54 GHz | `FREQUENCY_COMPARISON.md`, `frequency_comparison.json` |
| Item 3 cell candidates (stop point) | `CELL_CANDIDATES.md`, `cells/` (spec, rows cache, derived JSON per candidate), `search_summary.json`, `search/*.json.gz` (raw screening) |
| Fast suite log of the code commit | `fast_suite.log` |

Table coverage: all 36 tasks finished (no task skipped for the deadline). The table is a function of the modulation wavenumber K at four slot sizes; coverage of a design point is judged by the "far from every column" metric of the sweep, not by the range of N alone (K above 1.56 k0 is not tabulated).

Resume: `python scripts/research/phase28_supercell_f10.py run results/phase28/f10 --tier all --deadline <ISO time>`; the characterisations resume from `cells/*.rows.jsonl`.
