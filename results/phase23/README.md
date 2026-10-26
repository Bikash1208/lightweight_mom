# Phase 23 results (model-only)

- `comparison_v2_held_out.json` -- **current** square vs hexagon vs hexagon + optimiser comparison, 4-port 1-bit, 10.54 GHz,
  RHCP beams at theta 30 deg, phi 45/135/225/315 deg. Every number is on the held-out 0.5 deg x 1 deg grid at 10.33 / 10.54 /
  10.75 GHz with the two cell states RE-SOLVED at each frequency; starting designs are Phase 16 pointing-corrected; absolute
  pointing is constrained to <= 0.3 deg at the design frequency; `*_robust` rows minimise the worst case over the three
  frequencies (and, for all-ports, over 20 seeded +-10 deg port-phase errors). Produced by
  `scripts/research/phase23_compare_v2.py`.
- `comparison_held_out.json`, `frequency_check_baselines.json` -- superseded first version (no pointing correction; the
  +-2 % rows held the cells' 10.54 GHz phases, which is wrong -- the square's phase moves ~13 deg per 0.2 GHz).
- `rect_snapshot_before.json` / `rect_snapshot_after_C.json` -- the bit-identical gate for rectangle layouts
  (`scripts/research/phase23_rect_snapshot.py`).
- `rect_snapshot_server1.json` / `.provenance.json` -- numeric reference made on server 1 (2026-10-26, 1 BLAS thread): every number of the six
  rectangle syntheses, for number-by-number comparison (`scripts/research/phase23_rect_numbers.py`; see `docs/DEVELOPING.md`).
