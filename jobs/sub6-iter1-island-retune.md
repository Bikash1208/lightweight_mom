# SUB6 iteration-1 island retune (S2, 2026-11-10)

Authorized directly by user. Base 1ea2d13; isolated C:/lwmom_retune;
branch compute/sub6-iter1-island-retune; push github + hub. Preserve the
original C:/lwmom_final checkout and do not access C:/lightweight_mom.

MODEL_ONLY / NOT_HFSS_VALIDATED. BAR50-02V only: OFF (0.11 pF || 4 kohm)
+ 0.6 nH; ON 3 ohm + 0.6 nH. OFF sensitivity all four corners of
C={0.105,0.12} pF and Rp={3,5} kohm. One diode, physical gap 1.00 mm.
Accepted discretisation: mesh 96 gap 8 cells; mesh 128 gap 11 cells.

Targets: alpha_OFF >= 2 Np/m; 10 log10(alpha_OFF/alpha_ON) >= 12 dB,
ideally >=14 dB. Report beta shift, top three mesh-96/128 results and
sensitivity; check resonance/stop-band clearance >=1.5 GHz from 5.8 GHz.
No board/HFSS export or alternative diode work in this job.

Use unchanged src/ solver; research-only geometry wrapper. Reproduce old
mesh-96/F128 baseline before search at F=2n. Pool() uses all logical CPUs,
with one BLAS thread per process. Checkpoint every completed node.
Initial one-variable grid varies slot length/width, island length, bottom
clearance (which sets island width), and diode x position along the island.
Do not call a finite sampled frequency check proof of absence of all modes.

Commands (isolated .venv313/Scripts/python.exe):

```
python -m pytest -m "not slow" -q -n logical --maxfail=1
python scripts/research/sub6_iter1_island_retune.py baseline
python scripts/research/sub6_iter1_island_retune.py search
python scripts/research/sub6_iter1_island_retune.py search2
# Select candidates from computed results into results/.../selected.json.
python scripts/research/sub6_iter1_island_retune.py refine
python scripts/research/sub6_iter1_island_retune.py frequency
```

Hard stop on suspicious results/test failures; report evidence without
inventing missing results. Original diode_sheets.csv difference is only
column order: sheet_z_thick moved immediately after sheet_centre_y;
all 240 rows identical when compared by column name.

## Actual execution and authorized exceptions

Use previous S2 interpreter C:/lightweight_mom/.venv313/Scripts/python.exe
from C:/lwmom_retune with PYTHONDONTWRITEBYTECODE=1, and
OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=1. RETUNE_FFT=1 reuses
the existing src.dispersion.leaky_supercell_mom FFT assembler; it was audited
against production sums. RETUNE_WOODBURY=1 uses the exact small-load update
identity, audited against direct full solves. RETUNE_OUTPUT_SUFFIX is empty
for normal stages; _wb and _transpose were used for baseline audit outputs.
The direct baseline audit used RETUNE_WOODBURY=0. No src/ changes.
Original production-sum runs were restarted before any nodes completed due
to severe paging; no worker cap was imposed. Subsequent runs completed.

21 one-variable and 24 combined geometry evaluations (40 unique geometries),
then top-three nominal/sensitivity at meshes 96/128. Search2 follows the
first screen's shorter island and narrower slot; all dimensions recorded.

User accepted original four ENV_ONLY failures under revised norm/hash
criterion; full suite found six matching CP2 hash cases, all covered by
the same measured criterion. Polygon passes under previous pytest; app
timeout passed isolated rerun. No regression tests/fixtures changed.
New research-only circuit check expectation corrected after user confirmed
the series-L model and OFF HFSS impedance correction.

Final outcome: HARD_STOP_MESH_NOT_CONVERGED, ON alpha changes 176/446/1672%.
No frequency screening or board/HFSS export launched. All raw values and
uncertainty flags retained in results/sub6_iter1_island_retune/REPORT.md
and mesh_report.json. Commit/push evidence to github and hub as authorized.
