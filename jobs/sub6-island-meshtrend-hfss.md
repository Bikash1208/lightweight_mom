# S2 fixed island mesh trend + HFSS cell packages (2026-11-10)

User-authorized follow-up to c3178fe. Worktree C:/lwmom_retune;
branch compute/sub6-island-meshtrend-hfss. Original C:/lwmom_final artifacts
are preserved. MODEL_ONLY / NOT_HFSS_VALIDATED. Do not tune geometry.

## Mesh trend

Fixed BASELINE board row and G1. Nominal BAR50-02V OFF (0.11 pF || 4 kohm)
+ 0.6 nH, ON 3 ohm + 0.6 nH. Meshes 96/128 reused where present; missing
BASELINE 128 and both cells 160 computed, optional 192 launched after 160.
Mesh192 attempt failed allocating a 22.9-GiB complex matrix (39220 square),
before any 192 node completed. Record unavailable; no worker cap introduced.
All Floquet cutoffs 2n and unchanged Bloch nodes 1.017/1.097, both directions.
Report extrapolated beta explicitly; do not stop for ON mesh changes.
Physical gap 1.00 mm, nearest-cell mask: 8/11/14/17 cells at 96/128/160/192.
No solver src changes. Existing FFT assembler and exact small-load updates
as audited in the parent job; FFT row chunk reduced to 256 for scratch
efficiency, not a physical/numerical truncation. No worker cap: Pool()=64,
OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=1. Previous S2 Python.

Tool-helper refresh interrupted initial mesh160 attempt; completed BASELINE
128 nodes preserved, only missing nodes restarted as a hidden detached
process. Do not re-run completed checkpoints or switch branches mid-job.

## HFSS packages

Four: BASELINE/G1 x no-via/vias+FR4+stubs. User explicitly retained both
2.0 x 0.3 mm bottom stubs from ISL_vias2. Same via locations, 0.3-mm barrels,
0.6-mm bottom pads, revised 0.6-mm ground anti-pads, 0.8-mm FR4. Retain 30-mm
lower air extension with Radiation, PEC lid +150 mm; no near-stub PEC floor.
+90deg frame: ACROSS=ph_x, ALONG=ph_y. One shared geometry per OFF/ON;
one 1.00 x 0.35 mm Impedance sheet; OFF 5.424-j79.320, ON 1.05+j7.65 ohm/sq.
No AEDT scripts. Pre-export and re-import gates must pass.

Method found: scripts/research/sub6_final_cal.py and actual data in
results/sub6_final_board/hfss_data.json. Archived ISL_vias2 setup:
21fe444:hfss_exports/sub6_island_choke/HFSS_SETUP.md. Actual recorded
calibration used vias/FR4/stubs, not no-via. Eigenmode f(phi), Q -> beta
and alpha_x=omega/(2 Q vg). Historical ON slope came from a MoM branch
anchor; new packages require >=3 direct HFSS phase points to measure slope.

Do not equate uniform-cell Q loss alpha_x with tensor holographic leakage
alpha. No 3D MoM via prediction is available; mark unavailable instead of
substituting no-via numbers. Expected directional no-via sheet beta/loss
estimates are explicitly MODEL_ONLY and flagged when extrapolated.

## Reproduction

Use C:/lightweight_mom/.venv313/Scripts/python.exe from this worktree;
PYTHONDONTWRITEBYTECODE=1. Set three BLAS thread vars=1; RETUNE_FFT_CHUNK=256.

```
python scripts/research/sub6_island_meshtrend.py initial
python scripts/research/sub6_island_meshtrend.py 192
python scripts/research/sub6_island_meshtrend.py report
python scripts/research/sub6_island_cell_packages.py
python scripts/research/sub6_island_package_predictions.py
python scripts/research/sub6_island_finalize.py
```

CAD generator before prediction append; predictions include the one-page
setup's numeric reference table. Desktop copy includes only STEP/JSON/CSV/MD
and preserves any user-modified delivery files. Regeneration need not re-run
CAD once geometry gates/STEP SHA256 are verified; only numeric references
need updating after each mesh completes. Validate tests, commit and push
github+hub under the previously authorized ENV_ONLY exceptions.
