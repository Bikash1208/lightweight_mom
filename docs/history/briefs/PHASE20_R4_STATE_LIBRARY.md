# Phase 20 / R4 — CP-2 state library (frequency × angle × polarisation)

**Status:** active brief, written 2026-10-21. Branch `phase20-ris` (after R3 close-out `13463f5`). Read
`docs/briefs/PHASE20_RIS_ROADMAP.md` (incl. the R5 log entries), `docs/history/briefs/PHASE20_R3_TWO_BIT_CELLS.md`,
`results/phase20_r3/OPTIMISATION.md`, `LOCALITY.md`, `LP_INTERIM.md`. All repo rules apply (one commit per item,
`make test-fast` + snapshot gates green, physics tests, honest labels, additive only, no ML/fitted surrogates in the
solver core — interpolation of tabulated MoM data is allowed and must report its own error).

## Why
R5 (RIS design mode) needs, for every cell on the aperture, the reflection of each of the four CP-2 states at the local
frequency and incidence direction (θ, φ) — the full 2×2, from which the CP (2-bit) and LP (1-bit codes) responses
follow. R3 produced the cell at a few points; R4 turns it into a dense, registered, error-bounded library.

## Scope
Cell: `ris_cp2_octagon_madp_rt5880_15` exactly as registered (MADP-000907, L_s = 0 nominal). The LP 1-bit family is
derived from the same states (codes 1↔3 co-LP, 2↔4 cross-LP) — no separate solves.

## Work items (one commit each; stop points ■)

**0. Test-clock fix (infrastructure, not physics).** The machine's wall clock jumps (`time.time()` advanced 903 s in a
2 s sleep). Change job-waiting test helpers and any job-timeout logic to `time.monotonic()`. Tests:
the three affected tests pass in the full suite; a unit test that a wall-clock jump does not trip a timeout (monkeypatch).

**1. Symmetry reduction, proven before used.** The four one-hot states are 45° rotations of each other. Show
numerically (mesh 60, a few f/θ/φ points) that state k at azimuth φ equals state 0 at φ − 45°·k, in the rotated TE/TM
frame, to mesh-staircase tolerance (R3 saw ≈ ±3°). If it holds, compute only state 0 over φ ∈ [0°, 360°) — or the
smallest range the cell's own mirror symmetry allows (prove that too) — and generate the others by rotation. If it does
not hold well enough, report the error and compute all states.

**2. Grid and cost estimate. ■** Proposed grid: f = 8–12 GHz step 0.1 GHz; θ = 0, 5, …, 45° (θ > 30° labelled
`EXTRAPOLATED_BEYOND_R3_TARGETS`); φ step 15° over the reduced range; TE and TM. Before running, estimate wall time
and memory on this machine (36 GiB; mesh 120 ≈ 15.6 GB per angle is too heavy for the grid). Plan: grid at the coarse
mesh used in R3's search, plus a sparse set of corrected points (R3 extrapolation practice, first-order fallback) to
fit a **smooth mesh-correction map** (Δphase, Δ|r| vs f, θ), applied to the whole grid and reported with its residual.
**Stop and report the estimate if the full grid exceeds ~24 h**, with a reduced grid proposal.

**3. Build the library.** Store per (f, θ, φ, state): the 2×2 TE/TM reflection, plus derived CP (r_RR, r_LR, …) and LP
quantities, mesh level, correction applied, uncertainty. Format: one versioned file (HDF5 or npz + JSON manifest) under
`results/phase20_r4/`; large binaries follow the repo's size rules (gitignore if large, with a regeneration script).

**4. Interpolator with an honest error.** Interpolate in (f, θ, φ) on complex values (or magnitude/unwrapped phase —
justify). Validate on **held-out solved points not in the grid** (≥ 50 random points): report max and 95th-percentile
error in phase and dB. Gate: interpolation error ≤ 3° / 0.2 dB at the 95th percentile inside θ ≤ 30°; outside that,
report only.

**5. Library metrics report** (`results/phase20_r4/README.md`): per f and θ — 2-bit step errors, magnitude spread,
opposite-handed level, usable band (all R3 targets met), angular stability vs θ and φ; LP 1-bit codes' phase difference
and unwanted term; tolerance corners (±20 % C_T, L_s 0/0.1/0.2 nH) at a reduced grid (f = 9.5, 10, 10.5 GHz; θ = 0,
30°). Plots: phase/amplitude vs f per state; opposite-handed level maps vs (θ, φ).

**6. Registry update + HFSS hook.** Register the library against the family in the Phase 27 registry (domain = the
grid; outside it → `NOT_AVAILABLE_FROM_CURRENT_DATA`). Add an optional, labelled calibration hook so the user's CP-2
HFSS results (package `hfss_exports/phase20_r3_cp_cell/`) can later be compared and, if the user approves, applied as
an `HFSS_CALIBRATED` correction — the hook and a comparison script only; no calibration until data exists.

**7. ■ Close-out.** STATUS, CHANGELOG, roadmap log, archive this brief. Stop and report.

## Do not
- Change the CP-2 geometry or diode values (re-optimisation would be a new R3 round).
- Build R5 (coding map, array design, app UI) or any array-level coupling model — R5 owns that.
- Fit a surrogate model and call it physics; the library is tabulated MoM + interpolation with reported error.
