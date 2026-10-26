# Job request: r4-mesh120-corrections (OPTIONAL cross-check, server 1)

**From:** Mac session (coordinator). **Branch on the server:** `compute/r4-mesh120-corrections` from the base commit. Never commit to `phase20-ris`. **Base commit:** `7d4cc209271ae05a9fe43074c92518d19f00a11e`.

## Read first: these 45 points already exist
The mesh-120 correction points (9 frequencies x 5 (theta, phi) points, F = 2n, four states) were computed in R4 on the Mac; the chunks `corr120_f*.npz` are in `results/phase20_r4/chunks` (gitignored) and, since 2026-10-24, copied with a manifest to `C:\hmsa_sync\data\r4-mac-chunks\` (hub `data/r4-mac-chunks`). The full-grid library v2 uses them together with the server's mesh-90 chunks, which agree with the Mac's mesh-90 chunks to 1e-13. The "deferred mesh-120 points" of the server note are therefore already in the library; **this job is only an independent server re-computation (a second cross-machine confirmation of the most expensive mesh), not a prerequisite.** Run it only if you want that confirmation.

## Task
The 45 points of `CORR_F x CORR_POINTS` of `scripts/research/phase20_r4_build_library.py` at mesh 120, F = 240, all four states: `python scripts\research\phase20_r4_build_library.py corr120 <workers>` (writes `corr120_f*.npz` through the build script's chunk directory; point the chunk directory junction at `C:\hmsa_sync\data\r4-mesh120-server\chunks` as in `jobs/r4-full-grid-server.md`, do not overwrite the Mac's copy).

## Resources (memory-checked)
One task ~ 3500 s single thread (extrapolated from the 90/120 ratio), **~9 GB peak** (R4 notes; the server's own 6 GB estimate is the optimistic end). 255.6 GB RAM: **20 workers x 9 GB = 180 GB** leaves 75 GB for the OS and the other jobs; with the other server jobs running at the same time use **at most 12 workers** (108 GB). **Hard caps:** never more than 20 workers; stop and report if free memory falls below 30 GB at any time; total wall clock cap **5 h** (45 tasks / 20 workers x ~1 h = ~3 h, three waves).

## Outputs and comparison
The 9 chunk files + manifest (size, sha256, producing commit, machine); the Mac compares them with its own chunks (expected agreement ~1e-13 as for mesh 90).

## Stop conditions
Any task error -> stop that task and report; memory or wall-clock cap -> stop and report; do not patch code.
