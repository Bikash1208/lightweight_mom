# phase30_cell_eigenmode - HFSS eigenmode check of the Phase 30 cell (MODEL_ONLY predictions, written BEFORE any HFSS run; scripts NOT run in AEDT by the author)

Cell: Rogers RO3010 (eps_r 10.2, tan d 0.0022) 0.635 mm on PEC ground, 3.2 mm square cell, full-fill super-ellipse n = 4 (major = minor = 3.05 mm,
gap 0.15 mm), centre slot along x, PEC patch (sheet). Same master/slave + PEC-lid construction as `hfss_exports/slotted_cell_eigenmode/`
(Phase 26; the boundary dialogs there were the part not verified in AEDT: check Primary/Secondary pairs and the phase delays).
Air column 60 mm with a PEC lid (field at the lid ~ e^(-kappa 60 mm) = 2e-5 for kappa = 0.83 k0: negligible; with 30 mm it would be 4e-3).

## Cells (picked from `results/phase30/refine/`)
| Script | Slot L x W mm | M_t (model) | Bloch phase per cell |
|---|---|---|---|
| `build_E1_strong_M0p64.py` | 1.89 x 0.58 | 0.64 | 52.72 deg |
| `build_E2_longthin_M0p62.py` | 2.29 x 0.10 | 0.62 | 51.86 deg |
| `build_E3_mid_M0p27.py` | 1.89 x 0.10 | 0.27 | 48.96 deg |
| `build_E4_weak_M0p12.py` | 1.53 x 0.10 | 0.12 | 47.69 deg |
| `build_S1_supercell_leaky.py` | E1 x 7 cells, slot axis rotating 0..154 deg in 25.7 deg steps, radiation boundary on top | -- | 6.74 deg per 22.4 mm supercell |

The phase per cell is beta(10.54 GHz) x 3.2 mm with the model's mean beta (E1: 1.302 k0). Each script sets `PHASE_X` (along x, slot-parallel); for the y direction
swap the two phases (`phase_x = 0`, `phase_y = PHASE_X`). Run each script, set the phases, Analyze; the bound mode is the lowest eigenmode whose
field decays away from the patch plane (TM-like, E_z maximum at the patch edges, nodal at the ground).

## Predictions: eigenfrequency of the bound mode at the stated Bloch phase (GHz)
| Cell | along x, mesh 48 | along x, mesh 64 | along y, mesh 48 | along y, mesh 64 |
|---|---|---|---|---|
| E1 | 10.90 | 11.14 | 10.26 | 10.43 |
| E2 | 10.91 | 11.16 | 10.25 | 10.53 |
| E3 | 10.72 | 10.94 | 10.38 | 10.64 |
| E4 | 10.62 | 10.84 | 10.46 | 10.69 |
Reading: the mean of the x and y values is ~10.58 GHz (mesh 48) = the design frequency 10.54 within the mesh error; the x/y split (0.64 GHz for E1) is the anisotropy.
The two meshes differ by 2 % (0.2 GHz): the model is NOT mesh-converged for these near-resonant cells (see `results/phase30/REFINE_REPORT.md`), so
treat the predictions as +-0.25 GHz and the x/y SPLIT (E1 0.64, E2 0.65, E3 0.34, E4 0.16 GHz at mesh 48; 0.71, 0.62, 0.31, 0.15 at mesh 64) as the more robust test.
Pass criterion proposed: HFSS f_x, f_y within 0.3 GHz of the mesh-64 values and the split within 30 % -> the Bloch-tensor dispersion is trustworthy for this family.
Run E1 at eps_r 10.2 and 11.2 (edit `eps` in the material block): the model's beta/k0 moves 1.30 -> 1.57 for +1 in eps_r, i.e. the eigenfrequency at fixed phase drops by roughly 0.8-1 GHz.

## S1 (leaky supercell): alpha
Radiation boundary on the air top (no lid), Bloch phase 6.74 deg over the 22.4 mm supercell (= the n = -1 harmonic radiating at broadside for beta = 1.30 k0). Model: alpha/k0 = 0.058 (first order) or
0.031 (Floquet harmonic balance, mesh 48), v_g = 0.176 c  ->  Im f = 0.109 GHz (Q 49) or 0.058 GHz (Q 91) (alpha = 2 pi Im f / v_g).  HFSS reports a complex eigenfrequency;
alpha_HFSS = 2 pi Im(f) / v_g with v_g = d omega / d beta from the dispersion of E1 (model: 0.176 c; check with two E1 eigen runs at neighbouring phases).
A Q between ~40 and ~100 supports the model; Q > 200 means alpha is much lower than predicted (the main risk to the whole design).
The supercell Bloch phase is near 0 deg, where the eigenmode solver returns the band-edge pair: use the one with the radiating (broadside, circularly polarised) field.
Unverified: the AssignRadiation face call and the master/slave pair on the 22.4 mm cell (same construction as the unit-cell scripts, not run).

## Files
`predictions.json` (all numbers above, all meshes, all roots with their TM fractions); generator `scripts/research/phase30_eigenmode_package.py`.
