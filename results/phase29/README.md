# Phase 29 — seven-hexagon cluster (MODEL_ONLY, INTER_TILE_COUPLING_NOT_MODELLED, NOT_HFSS_VALIDATED)

Driver `scripts/research/phase29_seven_hex.py` (4.4 min), module `src/hologram/seven_hex_cluster.py`, data `seven_hex.json`
(per tile alone, all 12 pairs x 12 phases, every config). 10.54 GHz (supercell-table frequency; the HFSS run is quoted at 10.5).
Layout: vertex-at-+/-y hexagons sharing flat edges (gap 0; the Phase 26 board needed none for isolation), centres sqrt(3) R apart
at azimuths 0,60,..300 (259.8 mm for R = 150), one shared 6.4 mm lattice through the cluster centre, edge cells owned by one tile.
Gain referred to the accepted power (1 W per port, 2 W for a pair); anchored to the single-hexagon HFSS run only (OFFSET_DB).
Verification: the centre-origin tile at (30, 225) reproduces the committed Phase 26 cells to 1e-14 mm and 23.03 dBi (the anchor,
so this checks the pipeline, not the physics). "Merged" = exactly one RHCP lobe within 3 dB of the peak, off both tile beams.
