# Examples

User-facing scripts — things you'd actually run to use this project, as opposed to the one-off
research drivers that built it (those live in `scripts/research/`, see `scripts/README.md`).
Run from the project root, e.g. `.venv/bin/python examples/print_config.py`.

| Script | What it does |
|---|---|
| `print_config.py` | Load `configs/config.yaml` and print the parsed unit cell (SI internally, mm/GHz for display). |
| `run_sweep.py` | Run the sweep defined in a config file and write the unit-cell database (resumable). |
| `build_unit_cell_library.py` | Build the Phase 7 production unit-cell library (18x18 L/W sweep, 201 frequencies). Long-running, resumable via `checkpoint.jsonl`. |
| `inspect_unit_cell_library.py` | Query an already-built unit-cell library: exact-geometry lookup, phase-candidate search, phase inversion. Read-only, never re-solves. |
| `make_paper_board_26ghz.py` | Reproduce Mishra et al.'s 26 GHz 1-bit coded 5-port board through the app's own synthesis path; writes `docs/submission/paper26/comparison.md`. |
| `make_submission_figures.py` | Generate the Phase 15 submission figures (vortex vs pencil, deliberate OAM vortex, multi-port board + STEP/HFSS export). |
| `sweep_config.yaml` | Example config consumed by `run_sweep.py`. |

See `docs/STATUS.md` for what each output is validated against, and `docs/validation/README.md`
for the full provenance of any number these scripts produce.
