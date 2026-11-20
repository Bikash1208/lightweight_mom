# v3.1 gate field diagnostic — stored HFSS fields only

## STATUS

COMPLETE. Reopened the three converged HFSS 2026 R1 projects in fresh, private AEDT sessions; exported stored fields at 5.68 GHz. **No solve or geometry operation was run.** The 16-cell gate does not materially attenuate the near-top field behind the blocked half.

## GATES

Each project read back `Setup1 : LastAdaptive`, 5.68 GHz, and 518 PEC metal objects. The complex vector E field exported at 1 mm radial steps, eight azimuths, and z=2.035 mm (2.000 mm above the top copper). All 3 × 8 × 151 requested samples are in `radial_profiles.csv`; only r=0 in the PEC feed pin is `NaN`. A first parser assertion treated those expected pin samples as an error; the parser was corrected and the same stored field was re-exported. No re-solve was needed.

## DATA

The table is **10 log10(mean |E|²_case / mean |E|²_BOTH)** over r=40–140 mm, using complex phasor E. Ez uses the last complex component pair of the exported `<Ex,Ey,Ez>` vector.

| phi, degrees | A-only total dB | A-only Ez dB | B-only total dB | B-only Ez dB |
|---:|---:|---:|---:|---:|
| 0 | -0.21 | -0.14 | +0.39 | +0.56 |
| 45 | -0.49 | -0.45 | +0.16 | +0.27 |
| 90 (x=0) | -0.28 | -0.24 | -0.28 | -0.19 |
| 135 | -0.47 | -0.40 | -0.74 | -0.69 |
| 180 | -0.86 | -0.80 | +0.21 | +0.17 |
| 225 | -0.53 | -0.43 | +0.06 | +0.10 |
| 270 (x=0) | -0.18 | -0.08 | -1.03 | -0.93 |
| 315 | +0.28 | +0.34 | -0.90 | -0.81 |

| State / half | Total E² ratio dB, r40–140 | Ez² ratio dB, r40–140 | Total E² ratio dB, r6–29 |
|---|---:|---:|---:|
| A-only, **blocked x<0** (135/180/225) | **-0.63** | -0.55 | -0.13 |
| A-only, open x>0 (315/0/45) | -0.15 | -0.09 | -0.01 |
| B-only, **blocked x>0** (315/0/45) | **-0.03** | +0.09 | +0.09 |
| B-only, open x<0 (135/180/225) | -0.13 | -0.12 | -0.19 |
| A-only, x=0 phi90 / phi270 | -0.28 / -0.18 | -0.24 / -0.08 | -1.67 / -1.12 |
| B-only, x=0 phi90 / phi270 | -0.28 / -1.03 | -0.19 / -0.93 | -0.36 / -0.35 |

The blocked-half mean power-density proxy falls by only 0.63 dB in A-only and 0.03 dB in B-only. Individual lines can rise because coherent fields interfere. Thus the gate has a local effect but **does not block the outward near-top field** in either state. This is a near-top E-field comparison, not a modal decomposition or a direct bound-wave Poynting-flux measurement.

Via cells, from the supplied `gate_for` column (i,j):

- A-only, B vias: **(-1,-2), (-2,-1), (-1,-1), (0,-1), (-2,0), (-2,1), (-1,1), (0,1)**.
- B-only, A vias: **(0,-2), (1,-2), (1,-1), (2,-1), (2,0), (2,1), (0,2), (1,2)**.

## FLAGS

All fields use the saved port excitation in each project; no source edits were made. The total phasor |E| is formed from all six real/imaginary vector components, rather than AEDT's phase-0 `Mag_E` scalar. The profile includes near and radiated fields at z=2.035 mm; it cannot isolate TM surface-wave power. Net Poynting flux was optional and was not calculated. No other AEDT session or project was accessed.

## FILES

`radial_profiles.csv` contains state, phi, radius, position, total |E|, |Ez|, and their squares. `field_ratios.json` has unrounded line/group ratios and via coordinates. Three `*_E_complex.fld` files and `*_export_status.json` files preserve the raw calculator exports and own-session readbacks.
