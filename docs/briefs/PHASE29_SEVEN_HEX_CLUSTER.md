# Phase 29 — seven-hexagon HMSA cluster (centre broadside + six tilted outer beams, pair merging)

**Status:** active brief, written 2026-10-27. Branch `phase29-seven-hex` from `main`. The RIS track and Phase 28 stay parked.
Read `docs/START_HERE.md`, `docs/STATUS.md`, `CLAUDE.md` first. Labels of everything here: **MODEL_ONLY**, NOT_HFSS_VALIDATED.

## Concept (the user's design, built as specified)
- 7 hexagonal HMSA tiles side by side (centre + 6 around), each the Phase 26 tile (150 mm circumradius hexagon, 6.4 mm lattice,
  RT/duroid 5880 1.575 mm), each with its own centre monopole port -> 7 ports.
- Centre tile: broadside beam. Six outer tiles: beams 45 deg from broadside, azimuths 60 deg apart, each pointing outward along its own
  direction from the cluster centre.
- Tile flow: the HFSS-validated tensor-CP flow (`hex150_tensor_cp_single`, 23.03 dBi in HFSS) - scalar cells with a centre monopole split the beam.
- Frequency: the tile's validated frequency (supercell table: 10.54 GHz, the HFSS run is quoted at 10.5 GHz) unless the user says otherwise.
- Goal: when two neighbouring ports are excited together (centre + outer, outer + adjacent outer) the two beams should merge into one beam in between.

## Work (model only, <= 3 h compute, <= 40 lines of report, one cheap verification)
1. Design the 7 tiles with the existing Phase 26 flow (`src/hologram/tensor_cp_design.py`, unmodified); cluster layout (side-by-side hexagons,
   gaps as small as the ports allow). Each tile alone: gain, pointing, SLL, AR.
2. Pair excitation, forward model: superposition of the tiles' own fields at their real positions; inter-tile coupling NOT modelled (labelled).
   Centre+outer and outer+outer pairs, amplitude ratio 1:1, relative phase 0-360 deg in 30 deg steps. Report where the main peak(s) go,
   the combined gain, and whether a single in-between beam forms.
3. Limits: repeat step 2 with the outer tilt reduced (45, 30, 20, 15, 10 deg) and tile size (150 mm and one smaller): where do pairs merge?
4. Recommend, do not decide. Stop. The HFSS package for the 7-tile board comes after review.

## Do not
Modify frozen modules; claim HFSS-level numbers for tilts or sizes the HFSS run did not cover; fabricate coupling; start the HFSS package.
Gain stays anchored to the single-hexagon HFSS run only (`OFFSET_DB`), so every other number is an extrapolation of the model.

---
## Addendum 2026-10-27 - design point B (user decision GO)
ONE shared hexagonal aperture (size free, expected 80-120 mm), 7 monopole ports, 10.54 GHz, passive fixed hologram, polarisation free. Each port alone: a BROAD
beam (HPBW ~30-35 deg), port 1 broadside, ports 2-7 at ~30 deg tilt, azimuths 60 deg apart, neighbouring beams crossing near -3 dB; the 7 holograms superimposed
on the same cells; scanning = complex weights on the ports (2-3 active). Stage 1: design + per-port metrics + coupling estimate + weight table on theta 0-45 deg
(pass: pointing within +-3 deg and gain dip <= 3 dB for theta <= 30 deg; 30-45 band reported as is) + one cheap fresh-subagent mechanical verification.
Stage 2 (only if stage 1 does not fail): `hfss_exports/phase29_point_B/` (7-port AEDT script, predictions first, weight CSV, README with Edit Sources steps).
Code: `src/hologram/shared_aperture_multiport.py`, `scripts/research/phase29_pointB_sweep.py` (parallel, BLAS 1 per worker, memory-guarded), `phase29_pointB_stage1.py`.

## Addendum 2026-10-27 - option 2 (compact 7 tiles), user decision after point B stage 1 failed
Shared aperture dropped (pointing criterion not met). Build the compact 7-tile layout (R30 or R40, whichever meets the relaxed criterion with the better worst case), continuous scanning at low gain.
Relaxed criterion: pointing within +-8 deg worst case and gain dip <= 3 dB for theta <= 30 deg; 30-45 deg reported as is. Weights tuned per direction (<= 3 active ports, parallel, `phase29_compact7_tune.py`).
HFSS package `hfss_exports/phase29_compact7/` (`phase29_compact7_export.py`, `src/export/compact7_export.py`); one fresh-subagent verification; no merge, no tag.
