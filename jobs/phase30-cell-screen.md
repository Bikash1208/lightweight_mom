# Job: phase30-cell-screen (Phase 30 item 2, MODEL_ONLY)

Runs on the server itself (this session is the server). Base: branch `phase30-new-cell`.

```
python scripts/research/phase30_cell_screen.py results/phase30/screen --mesh 48 --workers 44
```

Output `results/phase30/screen/cell_screen.json` (one row per cell; unresolved-gap rows flagged by `gap_mesh_cells < 1.5`).
Single-threaded BLAS per worker. ~800 cells x 3 MoM solves at mesh 48.
