# Supercell resume 4 (2026-11-19)

## STATUS

D3p6 bridged 15-degree HFSS row complete; 30 degrees is the next continuity check. B2p2 45-degree row skipped by decision; no additional anchor HFSS solve.

## GATES

The common-frequency ON alpha >= 1.2 Np/m and ON/OFF >= 12 gates are pending D3p6 and, later, the OFF twin. Neither 5.8 GHz nor a particular beam angle is imposed. No operating point is selected yet.

## DATA

Anchor diagnostic: HFSS alpha uses local adjacent-phase slopes and the existing phase-matched Q baseline. Model alpha is reevaluated at each point's own radiated-harmonic angle. The model values use the same own 85/95-degree lidded TM beta pairs as the prior screen.

| Supercell phase | f GHz | local interval | v_g Mm/s | harmonic angle | HFSS alpha Np/m | chain alpha Np/m | delta/abs(B) | Q_base |
|---:|---:|---|---:|---:|---:|---:|---:|---:|
| 15 | 5.752876 | 15-30 | 158.595 | 32.56 deg | 2.0340 | 1.2958 | 0.283 | 1939.74 |
| 30 | 5.681668 | 15-30 | 158.595 | 31.41 deg | 1.6140 | 1.6775 | 0.339 | 2013.42 |
| 30 | 5.681668 | 30-45 | 193.735 | 31.41 deg | 1.3212 | 1.6775 | 0.339 | 2013.42 |
| 45 | 5.594682 | 30-45 | 193.735 | 30.35 deg | 1.0110 | 2.0171 | 0.470 | 2087.10 |

At 45 degrees the folded unit-cell phase is 84.375 degrees, 0.625 degrees below the Q baseline's 85-degree endpoint. Clamping Q_base to the measured 85-degree mean changes it from 2087.10 to 2062.54 and changes locally extracted alpha from 1.01105 to 1.01053 Np/m (0.05%). The opposite frequency slope survives all requested checks. The lowest-frequency delta/abs(B)=0.470 exceeds the stated ~0.4 first-order validity limit; its y beta was also extrapolated. **The chain is a 5.8 GHz screen only; the supercell is the reference.** Full readbacks: `anchor_resume4_diag.csv` and `.json`.

B2p2 D3p0 45 degrees is skipped. Its completed 15/30-degree candidate track gives provisional HFSS alpha 0.272/0.188 Np/m; the beta-to-sheet chain at 5.8 GHz gives 0.336-0.353 Np/m across theta 25-40 degrees. Every value is far below the 1.2 Np/m ON gate. This skips an expensive third phase without claiming a complete band or a verified mode track.

D3p6 bridged 15-degree HFSS row: 12 modes returned, highest 6.258287 GHz; 15 adaptive passes, max 102,994 solved tets, 18m32s solution time, normal completion. Candidate modes near the two folded x branches are listed below. These are **unassigned** until frequency continuity and saved fields separate them.

| Mode | f GHz | Q |
|---:|---:|---:|
| 3 | 5.664117 | 102.6 |
| 4 | 5.751826 | 92.2 |
| 5 | 5.793895 | 576.0 |
| 7 | 5.898020 | 723.6 |
| 8 | 6.054786 | 185.1 |
| 9 | 6.146695 | 114.1 |

The folded unrotated x-branch estimates at phase 15 are 5.7449 GHz (x-low) and 5.9032 GHz (x-high). Both windows contain plausible modes; phase 30 is needed. The 12-mode window ends below 6.5 GHz, so the upper spectral span is incomplete.

## FLAGS

- D3p6 has two x branches. The model-only 5.8-GHz predictions are 0.262-0.311 Np/m (x-low) and 3.032-3.308 Np/m (x-high) across theta 25-40 degrees. Both require HFSS mode matching.
- D3p6 phase 15 used a fresh local AEDT 2025.2 session, 5.4 GHz minimum, 12 requested modes, a 3600 s exact-own-PID watchdog, and 60 s heartbeats. Mode ambiguity supports a phase-30 continuation; phase 45 remains a later decision.
- The OFF twin and its copper-area audit remain last. Modal leakage does not determine which mode a centre pin excites.

## FILES

`anchor_resume4_diag.csv`, `.json` and `scripts/research/sub6_anchor_resume4_diag.py` reproduce the no-HFSS diagnostic. D3p6 run artifacts will be listed after completion.
