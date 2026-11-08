# EXTRAPOLATED - do not use as results

Partial mesh-96 refine run (10 geometries), stopped by decision after round 1 of the old node policy.
The old node rule (S.next_node: median of the OFF and ON betas) places nodes between the two states' betas, so every
alpha / beta / Delta derived from these checkpoints is a linear EXTRAPOLATION of Y(kt) (nearest node 0.03-0.05 k0 away).
Only checkpoints (ckpt/*.json: raw Y(kt) tensors at the nodes) exist; no island_rows.json was produced.
Re-run with ISL_NODEFIX=1 (see scripts/research/sub6_island_cell.py).
