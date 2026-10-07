# Ad hoc scripts

Rescued from an untracked scratch folder (`results/submission_work/`, Oct 5–6 2026) because their
own output is still referenced elsewhere in the repo. These are genuinely ad hoc — terse,
single-purpose, not covered by the project's usual test/docstring conventions — kept for
reproducibility, not as a supported CLI. Run from the project root.

| Script | Status | What it does |
|---|---|---|
| `fss_loop.py` | Parked | Periodic-MoM square-loop FSS unit cell (the `hfss_exports/fss_square_loop_10mm/square_loop_FSS_S11_S21.csv` generator). Usage: `python scripts/adhoc/fss_loop.py <mesh_n> <freq_ghz[,freq_ghz...]> [floquet_order]`. |
| `fss/fssgeom.py`, `fss/mom_save.py`, `fss/plane_s.py` | Parked (`docs/briefs/_parked/PHASE16_UNIT_CELL_SPARAM_TOOL.md`) | The finite-11x11-panel measurement-plane S-parameter technique for the same square-loop FSS (`hfss_exports/fss_square_loop_10mm/finite_panel_11x11_plane_S.csv`). `mom_save.py <mesh_n> <freq_ghz[,...]> [out_dir]` solves and caches per-cell current coefficients as `.npz` (not committed — regenerable, several minutes per frequency); `plane_s.py` reads those `.npz` files and computes the plane-defined S11/S21. `fssgeom.py` is the shared unit-cell geometry both import. |
| `fss/coef/` | Cache, not committed | Cached 0.1 mm-mesh MoM current coefficients for the square-loop FSS, one `.npz` per frequency (15 frequencies, 8.5–20 GHz; ~1 h total to regenerate with `mom_save.py`). `plane_s.py`'s own input. |
| `three_port.py` | Generator | Built `hfss_exports/three_port_25deg/`. `python scripts/adhoc/three_port.py <out_dir> <port_side_mm> <gap_mm>`. |
| `coded4.py` | Generator | Built `hfss_exports/four_port_1bit_gap0/` and `four_port_1bit_gap28/` (run twice, with `gap_mm=0` and `gap_mm=28`). `python scripts/adhoc/coded4.py <out_dir> <gap_mm> <polarization> <thickness_mm> [beta_mode]`. Loads `examples/make_paper_board_26ghz.py` as a module at runtime (`importlib`) — a real dependency, not a typo. |

All six were fixed in place to resolve the project root from their own file location
(`Path(__file__).resolve().parent...`) rather than relying on the current working directory or
the stray `results/submission_work/fss/src` symlink that `fss/*.py` originally depended on
(removed as junk — see the Phase 17 cleanup report). Everything else about them (argument
parsing, terse variable names, no tests) is unchanged from how they were actually used.
