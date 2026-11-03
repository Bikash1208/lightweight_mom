# Job: phase30-library (Phase 30 cell library, MODEL_ONLY)

Runs on the server itself. Base: branch `phase30-new-cell`. eps_r 11.2, f0 = 10.05 GHz, mesh 48, E1 family (squircle n=4 + slot, outline fixed, slot turns).

```
export P30_EPS=11.2 P30_F0=10.05e9
python scripts/research/phase30_library.py uniform results/phase30/library --workers 31   # 31 slots, uniform-cell table
python scripts/research/phase30_library.py lut     results/phase30/library --workers 44   # 31 slots x 12 a x 12 psi = 4464 tensors
python scripts/research/phase30_library.py heldout results/phase30/library --workers 44   # 60 off-grid points
```

Outputs: `uniform_table.json`, `angle_lut_v1.npz/json`, `heldout.json`. Single-threaded BLAS per worker, ~2 GB per worker.
