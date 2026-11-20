# V3.1 B-ALL-DIODE-V2, BAR90 OFF

HFSS 2026.1, one attached Direct solve at 5.68 GHz. Source was the unchanged V2 archive `BOARD_R150_V31_B_ALL_DIODE_V2.aedtz`, SHA-256 `54203F4835F49CA5F14CE0B799B03A96EC5D5119AA7C35103C6A9CF2E83363C3`. Only design variables changed: `Zr=2.5ohm`, `Zx=12ohm`. The model retained 250 DIODE sheets. Setup1 requested 6–12 adaptive passes and maximum delta S 0.02; six passes completed normally.

| Pass | Max magnitude delta S |
|---:|---:|
| 1 | unavailable (first pass) |
| 2 | 0.15061 |
| 3 | 0.02112 |
| 4 | 0.00775 |
| 5 | 0.00322 |
| 6 | 0.00253 |

S11 = **−17.6434 dB**; accepted power fraction = **0.982795**. Solve wall time 1593.8 s, peak reported memory 56.2 GB, 41 cores granted.

| Far-field quantity | Gain (dBi) | Direction (deg) |
|---|---:|---|
| RHCP peak | 8.643 | (32, 0) |
| LHCP peak | 0.329 | (36, 178) |
| LHCP | −9.347 | (26, 180) |
| LHCP | −3.458 | (30, 180) |
| LHCP | 0.216 | (36, 180) |

First-null SLL and other metrics are in `metrics_summary.json`. Raw FFD and S1P are in `OFF_BAR90_extract/farfield/`; the solver profile is `OFF_BAR90_solve.prof`. No other state was solved in this run.
