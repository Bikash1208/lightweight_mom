# HFSS validation set — grounded 8 mm reflectarray cell (Phase 12)

Generated from the project's own geometry code (`src/cad_export.py`, `src/phase12b/period_study.py`)
so the STEP files are the exact geometry the MoM solver uses. Units: **mm**.

## Stack (all files)

| Part (STEP name) | Region | Assign in HFSS |
|---|---|---|
| `Substrate` | z ∈ [−1.575, 0] | Rogers RT/Duroid 5880: εr = 2.2, tanδ = 0.0009 |
| `Metal` (patches) | z ∈ [0, 0.035] | **Perfect E** (or replace by a sheet at z = 0 with Perfect E) |
| `Ground` | z ∈ [−1.61, −1.575] | **Perfect E** (or Perfect E boundary on the substrate bottom face) |

The solver's metal is zero-thickness PEC. The 35 µm thickness exists only so the solids are valid CAD;
don't assign copper unless you deliberately want to test thickness/conductor loss.

Add an air box above the patch yourself (e.g. 15 mm), master/slave (linked) boundaries on the four
side walls, and a Floquet port on top.

## Cases

**Single cells** (`single_cells/`), 8 × 8 mm, square patch L = W, centred:

| ID | L = W (mm) | MoM S11 phase @ 12.5 GHz (ref 5 mm) | Role |
|---|---|---|---|
| G1 | 2.0 | −21.95° | flat part of the curve |
| G2 | 4.8 | −57.13° | knee |
| G3 | 5.6 | −115.22° | steep region |
| G4 | 6.0 | −168.74° (\|S11\| = 0.9957) | resonance — hardest case |
| G5 | 6.8 | +99.78° | past resonance |
| G6 | 7.6 | +59.48° | largest patch, 0.4 mm gap |

**3×3 supercells** (`supercells/`), 24 × 24 mm, 3 × 3 grid of 8 mm sub-cells (layouts in `manifest.json`):

| ID | Layout |
|---|---|
| S1_uniform_3x3 | all 5.6 mm — must reproduce G3 (checks the periodic setup) |
| S2_gradient_3x3 | columns (x) 5.2 / 5.6 / 6.0 mm, rows identical — smooth gradient |
| S3_contrast_3x3 | centre 5.6 mm, all eight neighbours 1.6 mm — strong contrast |

⚠ A 24 mm period equals λ0 at 12.49 GHz: above that, higher Floquet orders propagate. Include **all
propagating Floquet modes plus a few evanescent ones** in the port and report the (0,0) modes, or limit the
supercell sweep to 10–12.25 GHz. No MoM reference exists yet for the supercells — Phase 12C will produce it.

## Solver settings (to match the MoM reference)

- Normal incidence (θ = 0, φ = 0). Floquet port modes: keep TE00 and TM00 (singles); note which is x-polarised.
- **Phase reference plane: 5 mm above the patch (z = +5 mm).** Place the port there, or higher and de-embed
  to z = +5 mm. The reference CSV also gives the phase re-referenced to the patch plane (z = 0) if you
  prefer to de-embed all the way down.
- Frequency: 10–15 GHz, discrete, 0.25 GHz step (21 points).
- Convergence: adapt at 15 GHz, max ΔS ≤ 0.005.
- Time convention: HFSS uses exp(+jωt), same as the solver — no conjugation expected. If phases look
  mirrored (sign-flipped), flag it rather than conjugating silently.

## MoM reference

`single_cells/mom_reference_single_cells.csv`: 6 cases × 21 frequencies; S11_xx, S11_yy, S11_yx
(magnitude, phase in degrees on (−180, 180]) at the 5 mm plane (`ref5mm`) and at the patch plane
(`ref0mm`). Solver: grounded MoM, aligned 40 × 40 mesh (0.2 mm), Floquet order 40 — reproduces the
Phase 12B `period_sweep.csv` values exactly at the shared frequencies.

## What to send back

One CSV (or Touchstone) per case: frequency and complex S11 for both Floquet modes, including the
cross-polarised terms, with a header note saying which mode is x-polarised and where the reference plane is.
Put them in `hfss_exports/phase12_grounded_8mm/hfss_results/`.

## Missing raw Touchstone files

The raw `.s36p` exports under `hfss_results/` are not in git (`.gitignore`) and are no longer on disk,
so the tests that read them skip rather than run. To re-enable them, re-export the 36-port
Touchstone files from HFSS for each case listed in `manifest.json` and place them under
`hfss_results/` with the names the manifest gives.
