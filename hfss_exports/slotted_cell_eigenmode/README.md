> **Status: optional package, NOT run.** By decision the project relies on the internal solver (validated against HFSS earlier); these scripts and predictions remain for anyone who wants the check.

# slotted_cell_eigenmode -- HFSS unit-cell eigenmode checks of the slotted-circle family (Phase 26 B2.5)

Purpose: ground truth for the model's dispersion of near-resonant cells (the Bloch sheet tensor at kt > k0,
`src/dispersion/tensor_sheet.py`), which has never been compared with HFSS. The predictions below were written
BEFORE any HFSS run; MODEL_ONLY, NOT_HFSS_VALIDATED.

## The five cells (6.4 mm cell, circular patch, RT/duroid 5880 eps_r 2.2, tan d 0.0009, 1.575 mm on PEC ground)

| Script | Patch | Slot (L x W) | Slot angle to x |
|---|---|---|---|
| `build_E1_D5p8_slot0.py` | D 5.8 mm | 5.0 x 1.6 | 0 deg |
| `build_E2_D5p8_slot45.py` | D 5.8 mm | 5.0 x 1.6 | 45 deg |
| `build_E3_D6p0_slot0.py` | D 6.0 mm | 5.0 x 1.6 | 0 deg |
| `build_E4_D6p0_slot45.py` | D 6.0 mm | 5.0 x 1.6 | 45 deg |
| `build_E5_D5p8_bare.py` | D 5.8 mm | none | -- |

Slot ends keep a 0.29 mm (D 5.8) / 0.39 mm (D 6.0) copper bridge to the circle edge. PEC lid 60 mm above the
substrate (field at the lid ~ e^(-kappa 60 mm) = 4e-2 of the surface value: a small, unquantified frequency perturbation;
move the lid to 90 mm to check).

## Bloch settings (per script, two runs)

beta = 1.04 k0(10.54 GHz) = 229.7 rad/m; cell 6.4 mm.

* **along x**: master/slave phase `phase_x = 84.24 deg`, `phase_y = 0 deg`;
* **diagonal**: `phase_x = phase_y = 59.57 deg`.

The scripts build geometry, PEC patch/ground/lid, the substrate and air column, the eigenmode setup (>= 8 GHz,
6 modes) and TRY to assign the two lattice pairs with `PHASE_X` / `PHASE_Y`; if your AEDT version rejects a
property name, assign Primary/Secondary pairs by hand (x faces: phase_x, y faces: phase_y). The script was
not run in AEDT by the author.

## Which eigenmode

The bound surface-wave mode is the LOWEST eigenmode (near 10.2-10.6 GHz) whose field is strongest at the patch
plane and decays upward (|E| at the lid << |E| at the patch; Ez large over the patch). On the slotted cells a
TE-like bound mode (Ez ~ 0, field circulating around the patch) lies 0.05-0.3 GHz above it: both are
listed in the table.

## Predicted eigenfrequencies [GHz] (mesh 48; spot check at mesh 64 in brackets)

| Cell | along x: TM-like | along x: TE-like | diagonal: TM-like | diagonal: TE-like |
|---|---|---|---|---|
| E1_D5p8_slot0 | 10.574 (10.572) | 10.698 | 10.272 (10.320) | -- |
| E2_D5p8_slot45 | 10.286 | -- | 10.559 | 10.812 |
| E3_D6p0_slot0 | 10.535 | 10.576 | 10.189 | -- |
| E4_D6p0_slot45 | 10.211 (10.173) | -- | 10.515 (10.503) | 10.716 |
| E5_D5p8_bare | 10.567 (10.566) | -- | 10.556 (10.556) | -- |

Uncertainty: the staircased rotated slot and the 0.3-0.4 mm bridges make the model mesh-sensitive at the
~0.05 GHz level (mesh 48 vs 64 differ by up to 0.05 GHz here); HFSS finite-thickness effects are absent
(zero-thickness PEC in the model). A difference larger than ~0.1 GHz is a finding about the model.
Please report for each run: the eigenfrequency of the bound mode and of the TE-like neighbour, the mode's
Q (finite for a lossy slab only) and, if easy, |Ez| along a vertical line over the patch centre.
