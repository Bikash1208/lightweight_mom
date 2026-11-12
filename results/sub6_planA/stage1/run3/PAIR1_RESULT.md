# Stage 1 Run 3 — V pair

## STATUS

V OFF and V ON were solved in separate, attached AEDT 2025.2 sessions, one at a time. Both completed the specified four passes and were saved and extracted. Neither reached the 0.02 delta-S target; the numbers below are four-pass results, not converged results.

## GATES

| Gate | V OFF | V ON |
|---|---|---|
| Saved build and AEDT readback | PASS | PASS |
| Full design validation, one port | PASS | PASS |
| Private AEDT temp readback | `C:/lwmom_labB/aedt_temp/ANCHOR_V_OFF` | `C:/lwmom_labB/aedt_temp/ANCHOR_V_ON` |
| Setup after removing imported `Auto1` | `Setup1` only, 5.8 GHz, max 4 passes | `Setup1` only, 5.8 GHz, max 4 passes |
| HFSS CPU request/readback | 48/48 | 48/48 |
| FFD nonzero and present | Yes | Yes, after reopening saved solution for extraction |

## DATA

| Metric | V OFF | V ON |
|---|---:|---:|
| Adaptive passes | 4 | 4 |
| Final tetrahedra | 1,677,410 | 1,619,088 |
| Final delta-S | 0.083473 | 0.228100 |
| AEDT profile peak RAM | 53.9 GB | 55.2 GB |
| Solve time | 3,427.66 s | 3,519.89 s |
| S11 at 5.8 GHz | −0.01453 dB | −0.00823 dB |
| GainRHCP / GainLHCP at (25°, 0°) | 7.591 / 0.194 dB | 10.914 / 0.392 dB |
| GainRHCP / GainLHCP at (25°, 180°) | 16.566 / −1.606 dB | 12.437 / −2.176 dB |
| GainRHCP / GainLHCP at (24°, 180°) | 16.566 / −1.605 dB | 12.438 / −2.178 dB |
| Peak RHCP | 16.621 dB at (24°, −178°) | 12.887 dB at (28°, 2°) |
| Radiated-power fraction, theta > 90° | 0.13898 | 0.28377 |

## FLAGS

- AEDT reported `Adaptive Passes did not converge` for both states at max passes = 4. The controller's `analyze_setup` call returned success and the saved solutions were extractable. No further passes were added.
- V OFF's original postprocess attempt failed because the profile parser expected AEDT's internal format rather than its exported text format. Extraction was repaired and run from the saved solution, without another solve.
- V ON's first attached-session metadata export produced an all-zero FFD with no port power data. Direct element-pattern export and AEDT point-gain queries on the saved solution were nonzero; after reopening, the metadata export returned the correct `Port1:1` source and nonzero FFD. The final row uses that reopened export. No model or solution was changed.
- AEDT's high S11 values are reported as solved; the user will assess RF suitability. The 240-diode anchor remains above the eventual 200-diode goal.

## FILES

- `ANCHOR_V_OFF.json`, `ANCHOR_V_ON.json`: per-run metrics and provenance.
- `ANCHOR_V_OFF/`, `ANCHOR_V_ON/`: exported `.prof`, `.conv`, `.ffd`, metadata XML and Touchstone files.
- Private AEDT projects, diagnostic logs, temporary files and prior failed export artefacts remain local under `C:\lwmom_labB\stage1_run3_work` and `C:\lwmom_labB\aedt_temp`.
