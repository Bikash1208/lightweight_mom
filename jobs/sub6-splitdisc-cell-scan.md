# Job: sub6-splitdisc-cell-scan (SUB6_SPLITDISC_1BIT stage 1, MODEL_ONLY)

Brief: `docs/briefs/SUB6_SPLITDISC_1BIT.md` (read sections 1, 4, 5 first). Server session. Branch `compute/sub6-splitdisc-cell-scan`
from `sub6-splitdisc` (created from the head of `sidequest-hybrid-hms-ris`; report the hash).
Single-threaded BLAS per worker (`OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`), workers x threads <= cores,
check RAM before launching. Checkpoint every evaluation (resumable). Provenance sidecar per rule 3 of `INFRA_COMPUTE_SERVER_SYNC.md`.

## 0. First step
Fast suite (`make test-fast`); report failures before using any result.

## 1. Script to write (server): `scripts/research/sub6_splitdisc_cell_scan.py`
Reuse existing code only (no solver change): the loaded Bloch sheet tensor at kt = beta (`src/dispersion/tensor_sheet.bloch_sheet_tensor`,
as used by `sidequest_cell_search.bloch_loaded` / `sidequest_rotated_cell_bloch.py`), `src/loads.LumpedLoad` for the diode (rotated loads
via `direction` = angle, `transverse_conducts` default from `fcc6691`), `PolyComposite` / `Poly` (`src/polygon_geometry.py`) or the
existing slotted-circle class (`src/dispersion/slotted_cells.SlottedCircularPatch`, if a full-length slot splits the disc) for the geometry.
**State in the report which path was used and why.** Add a test file with the physics checks of section 4.

Geometry (square cell, period p, grounded substrate `backing: pec`, F = 2n Floquet cutoff):
- Variant S: two semicircles of radius r, flat edges facing across a straight gap g through the centre (gap axis along x at psi = 0).
- Variant T: radii r1 (half 1), r2 = q r1 (half 2), same gap.
- Diode: lumped load bridging the gap at the cell centre. **The gap and the load footprint are set by the diode package, not free knobs**
  (user decision 2026-11-07); the HFSS RLC sheet uses the same footprint:
  - MA4AGP907 (flip chip): g = 0.33 mm (circuit-pad spacing 0.013 in), sheet 0.33 mm across x 0.35 mm along the gap (chip width 0.343-0.368 mm).
    Source: MACOM MA4GP907 data sheet outline (chip 0.66-0.69 x 0.34-0.37 mm, pads 0.305 x 0.203 mm). **Taken from the MA4GP907 sheet; the
    MA4AGP907 outline is assumed identical (same 907 flip-chip family) - UNVERIFIED until the user confirms from the MA4AGP907 sheet.**
  - SC-79 parts (BAR50-02V primary, SMP1345-079LF fallback): g = 0.65 mm and sheet 0.65 x 0.35 mm from the Infineon SC79 footprint values
    (pads 0.35 mm, 1.35 mm overall): **ESTIMATE** until confirmed from the drawing.
  - `OFF_OPEN` is run at both gaps.

Grid (frequency 5.8 GHz):
- substrate RT5880 (eps_r 2.2, tan d 0.0009), h in {1.575, 3.175} mm; p = 10 mm.
- S: r from 2.0 mm to p/2 - 0.3 mm in 0.2 mm steps; g and footprint per diode as above (MA states at g 0.33, SC-79 states at g 0.65).
- T: r1 on the same grid, q in {0.6, 0.7, 0.8, 0.9}; BAR50 states only (g 0.65). MA4AGP907 is over budget (brief sec. 2b): its rows are a
  reference only (user-measured part), variant S only.
- States: `OFF_OPEN` (no load), `OFF_MA` (28.16 fF), `ON_MA` (2.44 ohm) [MA4AGP907 user-measured], `OFF_SMP` (0.19 pF + 0.7 nH series),
  `ON_SMP` (1.5 ohm + 0.7 nH) [SMP1345-079LF from `src/ris/diodes.py`, in-production fallback], and the **chosen diode (user, 2026-11-07)
  Infineon BAR50-02V**: `OFF_BAR50` (0.15 pF + 0.6 nH; 0 V, 1-1.8 GHz typ), `OFF_BAR50_HI` (0.20 pF + 0.6 nH; 0 V / 100 MHz, sensitivity),
  `ON_BAR50` (3.0 ohm + 0.6 nH; 10 mA typ), `ON_BAR50_MAX` (4.5 ohm + 0.6 nH; 10 mA max, loss sensitivity). Source: Infineon BAR50 series data sheet.
  `OFF_SMP_0V` (0.25 pF + 0.7 nH, ESTIMATE, 0 V sensitivity of the 1 V SMP value).
  SC-79 states (SMP, BAR50) use **g 0.65 mm, sheet 0.65 x 0.35 mm** (BAR50 footprint, ESTIMATE until read from the Infineon SC79 drawing).
  **BAR50-02V is the primary diode: run every grid point with it; SMP and MA rows at g 0.65 / 0.33 only for variant S.**
- Tensor at kt = beta along the gap axis AND across it (two runs), beta from the uniform lattice of that state (iterate to self-consistency as the
  existing scripts do; report the iteration count). psi = 0 for the grid.
- Mesh 64 for the grid (cell 0.156 mm: the 0.33 mm gap and 0.35 mm load are about 2 mesh cells; mesh 48 would give 1.6, too coarse).
  Convergence subset (mesh 64 / 96 / 128): 3 geometries x 3 states (`OFF_MA`, `ON_MA`, `OFF_OPEN`) x both h.
- Rotation subset: psi in {22.5, 45} deg for the same 3 geometries, states `OFF_MA` / `ON_MA`, mesh 64.
- Frequency subset: 5.6 / 6.0 GHz for the same 3 geometries (re-centring if the user chooses another band).
Choose the 3 convergence geometries as small / medium / largest r of variant S, g = 0.4, after the grid (report which).

## 2. Output per row (`results/sub6_splitdisc/cell_scan.json`, one row per geometry x state x kt-direction)
Y in the gap frame (complex, mS): Y_par, Y_perp, Y_off (off-diagonal); mean = (Y_par + Y_perp)/2; Delta = Y_par - Y_perp; beta/k0 (uniform lattice);
alpha_abs (ohmic + dielectric, Np/m); nearest |Y| extremum / resonance in 5.0-6.6 GHz if detected; mesh; gap mesh cells (flag `< 1.5`).
Derived per geometry: rho_MA = Delta_ON_MA / Delta_OFF_MA, rho_BAR50 (primary) and rho_SMP likewise (complex; report Re and Im); delta_mean = mean_ON - mean_OFF;
delta_beta/k0.

## 3. Gates (named states, never bare pass/fail)
- `SYMMETRY`: |Y_off| / |Delta| <= 0.02 at psi = 0 for every row -> `H1_SUPPORTED`; otherwise `H1_VIOLATED` (list rows). Rotation subset:
  tensor at psi equals R Y R^T within 3 % -> `ROTATION_COVARIANT` / `ROTATION_DEVIATES` (give the size).
- `CONVERGENCE`: beta/k0 change mesh 96 -> 128 <= 0.5 % and Delta change <= 5 % -> `MESH_CONVERGED` / `MESH_NOT_CONVERGED` (do not hide it).
- `RHO_SPAN` (report, do not decide): min / max Re(rho) over the grid per variant / h / diode; list the geometries with |rho| <= 0.2
  (B-cell candidates), 0.8 <= rho <= 1.2 (A-cell candidates) and Re(rho) <= -0.5 (sign-flip candidates); for each, delta_beta/k0 and alpha_abs ON.
- `RESONANCE`: any geometry/state with an extremum within 0.5 GHz of 5.8 GHz -> flag `NEAR_RESONANT`.

## 4. Tests to add (physics, not "it runs")
Reciprocity of the tensor (Y_xy = Y_yx within solver tolerance), mirror symmetry (variant S: Y_off ~ 0 at psi = 0), rotation covariance
psi = 0 vs 90 deg (exact swap of par/perp), `OFF_OPEN` equals a load with Z -> infinity limit, `ON_MA` -> PEC strap limit as R -> 0.

## 5. Report (`results/sub6_splitdisc/CELL_SCAN.md`, <= 2 pages) and stop
Gates first, then a table of the candidates, then plots: Delta_OFF, Delta_ON, rho and beta/k0 vs r (one panel per h and diode).
Labels: MODEL_ONLY, IDEAL_BIAS (diode = lumped load only, as the HFSS RLC sheet will be), NOT_HFSS_VALIDATED.
Board context for the report (no board work in this job): 5.8 GHz confirmed, beams at theta 25 deg, phi 0 / 180 deg. **Stop after the report.** No board design, no
architecture choice: the decision session decides stage 2.
