# results/phase28 - Phase 28 (monopole-fed coded holographic antenna), branch phase28-coded-holographic-v2

- `LITERATURE.md`, `crossref/` - item 1: nine papers, metadata from Crossref only (`scripts/research/phase28_literature.py`); content columns `NEEDS_FULL_TEXT`.
- `feasibility_results.json` (+ `.provenance.json`) - item 2: the sweep (`scripts/research/phase28_feasibility_sweep.py`, 76 cases, 1 BLAS thread, server 1) on the Phase 26 forward model
  (`src/hologram/coded_feasibility.py`); anchor to the committed Phase 26 single-hexagon prediction, quantisation / pitch / steering / aperture cases, switching-contrast derivation, control-hardware counts.
- `FEASIBILITY_REPORT.md` - generated from the JSON only (`scripts/research/phase28_feasibility_report.py`); `tests/test_phase28_feasibility_report.py` checks that every number in it is a value of the JSON.

Labels: MODEL_ONLY, IDEALISED_CELLS, NOT_HFSS_VALIDATED. 10.0 GHz: NOT_AVAILABLE_FROM_CURRENT_DATA (the supercell table exists at 10.54 GHz only). Item 2 ends at a stop point: no cell is designed before it is approved.
