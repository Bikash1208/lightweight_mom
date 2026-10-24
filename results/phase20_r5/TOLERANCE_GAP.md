# R5-1: MoM vs HFSS diode-tolerance sensitivity (does F = 2n close the gap?)

Data: server run `jobs/r5-tolerance-mom` (`tolerance_mom.npz`, merged from `compute/r5-tolerance-mom` 7eb06a8, sha256 checked against its manifest; meshes 60 and 90, **F = 2n for both**, C_off 0.0225 / 0.025 / 0.030 pF, states 1-4, theta 0 / 30) against the user's HFSS tolerance run. Script `scripts/research/phase20_r5_tolerance_gap.py`, numbers in `tolerance_gap.json`. Differences are taken inside each data set (C_off minus nominal), so the plane de-embedding and the library mesh correction cancel. COMPARISON ONLY, MODEL_ONLY.

## Result: the cutoff does not close it; the mesh does
Largest change of any ladder step caused by the C_off change (deg), 10.5 GHz:

| theta | C_off | HFSS | MoM mesh 90, F = 2n | MoM mesh 60, F = 2n |
|---|---|---|---|---|
| 0 | 0.0225 | 1.0 | 1.7 | 31.1 |
| 0 | 0.030 | 1.9 | 2.4 | 31.9 |
| 30 | 0.0225 | 1.1 | 1.8 | 42.0 |
| 30 | 0.030 | 2.2 | 2.8 | 42.8 |

* **Mesh 90 reproduces the HFSS sensitivity**: step changes <= 2.8 deg (HFSS <= 2.2), and the per-state phase changes with C_off agree to ~1 deg (theta 0, C_off = 0.030: HFSS -5.5 / -6.1 / -4.3 / -3.6, MoM -5.6 / -5.3 / -4.1 / -3.2).
* **Mesh 60 does not, at 10.5 GHz only** (31-43 deg), although it is fine at 10.25 and 10.75 GHz (<= 3.4 deg). The mesh-60 cell has its resonance slightly displaced from the converged one, so 10.5 GHz sits on a steep slope and a small C_off change moves the phase tens of degrees. F = 2n is already used in both columns, so the cutoff is not the cause; the discretisation (diagonal stubs) is.
* The R4 corner table (`results/phase20_r4/tolerance_corners.json`) was computed at mesh 60 with the *nominal* mesh correction added; its 10.5 GHz entries (C_T +-20 %: 19-30 deg state-phase shifts) are this mesh-60 artefact and are **superseded** for 10.5 GHz by the mesh-90 numbers above. Its 9.5 and 10.0 GHz entries (4-8 deg) are not affected in this way.
* A note on wording: the "17-28 deg step errors" in the R3/R4 text are the *absolute deviation of the steps from 90 deg* (nominal and corners alike); the HFSS "<= 2-3 deg" is the *change with C_off*. The comparable quantities are the changes, listed above.

Consequence for R5: the sensitivity band built from the HFSS C_off run (R5-3) is supported by the converged MoM; the cell is robust to the OFF capacitance at 10.5 GHz (main-beam gain moves <= 0.1 dB).
