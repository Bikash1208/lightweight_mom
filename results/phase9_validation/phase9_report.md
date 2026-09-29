# Phase 9 report -- physical control-variable and local-periodicity validation

## Final status: `LOCAL_PERIODICITY_INVALIDATED`

- Gate A (control-variable validation process): `CLOSED`
- Gate B (local-periodicity validation process): `OPEN`
- S11 classification: `UNVALIDATED`
- S21 classification: `UNVALIDATED`
- Selected control variable: `None`
- Local periodicity: `INVALIDATED`
- Geometry inversion: `LOCKED`

## Mandatory questions (Sec. 41)

### 1. What exact electromagnetic architecture was validated?

This project's own MoM solver architecture (Phase 1-7, unmodified): PARTIALLY_REFLECTIVE_TRANSMISSIVE -- more transmissive than reflective (mean reflected power 0.1922, transmitted 0.8070). No HFSS architecture could be validated: UNKNOWN (no HFSS project exists in this repository).

### 2. What does S11 physically represent in that architecture?

specular (0,0) Floquet order returning through port 1 (top reference plane), referenced at the top plane (-z (toward the pattern) at z_top = +top_air_thickness (above the pattern)).

### 3. What does S21 physically represent?

specular (0,0) Floquet order continuing through port 2 (bottom reference plane, z_bottom = -(h + bottom_air_thickness) (below the substrate)).

### 4. Which quantities were compared against HFSS?

S11_xx and S21_xx of the 10 reference cells at 7 frequency points each (Sec. 7, 11) -- see `s11_comparison.csv`/`s21_comparison.csv`.

### 5. What are the complex, magnitude, and phase errors?

NOT AVAILABLE FROM CURRENT DATA: no HFSS reference file found at /Users/biksmac/Documents/ISRO/HMSA/lightweight_mom/results/hfss_reference/hfss_reference_cells.csv

### 6. Are those errors frequency-dependent?

NOT AVAILABLE FROM CURRENT DATA (no HFSS reference exists to test this).

### 7. Are they geometry-dependent?

NOT AVAILABLE FROM CURRENT DATA (no HFSS reference exists to test this).

### 8. Does S11 agree sufficiently with the physical architecture?

`UNVALIDATED` -- NOT_AVAILABLE_FROM_CURRENT_DATA: no HFSS reference file found (hfss_reference.py)

### 9. Does S21 agree sufficiently with the physical architecture?

`UNVALIDATED` -- NOT_AVAILABLE_FROM_CURRENT_DATA: no HFSS reference file found (hfss_reference.py)

### 10. What happens to the centre cell when neighboring cells change?

Its induced-current dipole moment measurably changes; see `supercell_results.csv`/`lpa_error_metrics.csv` for every tested contrast level: 3x3_weak: relative complex error 0.1591, phase error 8.46 deg; 3x3_moderate: relative complex error 0.0599, phase error 3.20 deg; 3x3_strong: relative complex error 0.997, phase error 52.49 deg; 3x3_symmetric: relative complex error 0.1307, phase error 5.28 deg; 5x5_moderate: relative complex error 0.09026, phase error 3.91 deg; 5x5_strong: relative complex error 0.7436, phase error 15.85 deg.

### 11. Does the local-periodic approximation remain valid?

`INVALIDATED` -- at least one physically-realistic, in-sweep-range contrast level produced a deviation at least 50% of the signal's own magnitude (a self-referential, not an externally-chosen, criterion; Sec. 24).

### 12. Over what geometry/frequency range is it valid?

Not established over any range: this run's classification is a single overall verdict (`INVALIDATED`), not a range-restricted one -- `PARTIALLY_VALIDATED` was not reached, so no restricted applicability range can be honestly stated (Sec. 23).

### 13. What assumptions remain unvalidated?

Whether S11xx or S21xx is the physically correct holographic control variable (no HFSS data exists to test either); whether the local-periodic approximation holds under realistic neighbor contrast (measured evidence above suggests it does not, for strong contrast); the reference/object-wave formulation questions Phase 8B-8E already left open.

### 14. Is geometry inversion scientifically unlocked?

**NO**.

### 15. Is final STEP generation scientifically unlocked?

**NO** (tied to the same gate as geometry inversion, Sec. 5, 30).

### 16. What is the exact next phase?

Geometry inversion remains blocked. Two independent next steps are required: (1) obtain real HFSS (or other independent full-wave) reference data for at least the 10 reference cells at the 7 validation frequencies, since none exists in this project and Rule 1 forbids fabricating it; (2) given the measured local-periodicity invalidation under strong neighbor contrast, develop a coupled-cell or heterogeneous-library capability (Sec. 34) rather than continuing to use the isolated periodic lookup table for holographic synthesis.

## Rule compliance

- No boundary condition, port definition, polarization convention, normalization factor, empirical correction, phase offset, amplitude normalization, S-parameter transformation, HFSS result, validation threshold, or literature equation was invented (Rule 1).
- No mathematics was modified to force S11 or S21 agreement, or to force local-periodicity agreement (Rule 2).
- Every comparison retained complex, magnitude and phase error separately (Rule 3).
- No amplitude error was normalized away (Rule 4).
- Geometry inversion status: `LOCKED` (Rule 5).
