# CP-2 diode-tolerance cases in HFSS (user run 2026-10-23), advisor's quick read

File `hfss_cp2_tolerance_Coff_0225_025_030.csv`: same project and reference plane as run 2 (residual de-embed offset
-0.482 mm vs the patch top, corrected below as exp(-j2·k0·0.482 mm·cos θ) on top of run-2 data), parametric
Coff = 0.0225 / 0.025 / 0.030 pF (−10 % / nominal / +20 %; OFF capacitance of all OFF diodes), states 1-4,
θ = 0 and 30°, 8-12 GHz step 0.05 GHz. Coff = 0.025 agrees with run 2 to |ΔS| ≤ 0.11 (HFSS run-to-run mesh variation;
steps agree within ~1 deg).

## Result (circular terms, repo formulas)
| f GHz | θ | Coff pF | steps (deg) | max |step−90| | |r_RR| dB | opposite-handed dB |
|---|---|---|---|---|---|---|
| 10.5 | 0 | 0.0225 | 97.7 / 91.7 / 102.1 / 68.5 | 21.5 | −0.68…−0.60 | −15.0 |
| 10.5 | 0 | 0.025 | 97.4 / 92.6 / 102.5 / 67.5 | 22.5 | −0.63…−0.58 | −16.7 |
| 10.5 | 0 | 0.030 | 96.7 / 94.4 / 103.2 / 65.6 | 24.4 | −0.68…−0.50 | −17.6 |
| 10.5 | 30 | 0.0225 | 96.5 / 91.4 / 100.2 / 72.0 | 18.0 | −0.88…−0.57 | −10.6 |
| 10.5 | 30 | 0.025 | 96.1 / 92.4 / 100.6 / 70.8 | 19.2 | −0.77…−0.55 | −11.8 |
| 10.5 | 30 | 0.030 | 95.2 / 94.6 / 101.6 / 68.6 | 21.4 | −0.60…−0.53 | −15.1 |
(10.25 and 10.75 GHz in `advisor` output: same picture; band edges degrade opposite-handed to −7…−12 dB.)

## Reading
* **In HFSS the cell is robust to the diode OFF capacitance:** −10 % / +20 % moves individual steps by ≤ 2-3 deg and
  |r_RR| by ≤ 0.2 dB. The MoM predicted 17-28 deg step errors for the same corners — the model's tolerance sensitivity
  is largely a model artefact (same region where the MoM staircases the diagonal stubs). This is why the R3b objective,
  dominated by those corners, could not improve.
* The real, tolerance-independent issue is the **short 4→1 step (≈ 66-72 deg)** and the **θ = 30° opposite-handed
  level (≈ −11…−12 dB nominal)**. +20 % Coff happens to improve the opposite-handed level at 10.5 GHz.
* L_s (+0.1 nH) was not run in HFSS.
