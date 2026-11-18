# Supercell resume 4 (2026-11-19)

## STATUS

D3p6 bridged 15/30/45-degree HFSS rows complete. B2p2 45-degree row skipped by decision; no additional anchor HFSS solve. No OFF-twin solve was started because a common viable f_op has not been established.

## GATES

The D3p6 upper-branch candidate is below the ON alpha >= 1.2 Np/m gate at all measured points, but its 45-degree mode identity is uncertain; the gate is **not certified**. The common-frequency ON/OFF >= 12 gate remains untested because the OFF twin has not been run. Neither 5.8 GHz nor a particular beam angle is imposed. No operating point is selected yet.

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

D3p6 bridged 30-degree row: 12 modes through 6.228252 GHz, max 100,780 solved tets, 20m32s solution time, normal completion. Candidate surface-wave frequencies/Q are mode 3 5.616420/128.6, mode 4 5.727140/212.3, mode 5 5.793549/94.8, mode 7 5.936188/823.3, and mode 8 6.005657/223.1. The negative-fold x-low/x-high unrotated predictions are 5.6972/5.8306 GHz. Modes crowd both windows, and Q changes sharply between some plausible 15-to-30 continuations. Phase 45 and saved-field matching are needed before assigning the two branches or computing alpha from a slope.

D3p6 bridged 45-degree row: 12 modes through 6.232631 GHz, max 118,236 solved tets, 21m18s solution time, normal completion. All three 12-mode windows stop below 6.5 GHz; this is an explicit spectral-span limit of the requested narrower setup.

The table gives the two negative-fold x-branch **candidates**. v_g uses adjacent phase slopes at 15/45 and the 15-to-45 central slope at 30, with the 92.8 mm supercell period. Q_base is the phase-matched mean of own unrotated PML40 across/along Q. The last angle is the radiated -1 harmonic; beta/k0 is the unfolded bound wave. Every f/Q is raw HFSS 2025.2.

| Track | phi_s | Mode | f GHz | Q | Max tets | beta/k0 | v_g Mm/s | Q_base | alpha Np/m | harmonic angle |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| x-low mean candidate | 15 | 3 | 5.664117 | 102.6 | 102,994 | 1.1169 | 106.2 | 1794.2 | 1.539 | 33.13 deg |
| x-low mean candidate | 30 | 3 | 5.616420 | 128.6 | 100,780 | 1.1025 | 121.9 | 1835.4 | 1.047 | 31.82 deg |
| x-low mean candidate | 45 | 2 | 5.554680 | 183.4 | 118,236 | 1.0905 | 137.5 | 1876.6 | 0.624 | 30.59 deg |
| x-high candidate | 15 | 8 | 6.054786 | 185.1 | 102,994 | 1.0449 | 109.4 | 1524.9 | 0.825 | 30.75 deg |
| x-high candidate | 30 | 8 | 6.005657 | 223.1 | 100,780 | 1.0310 | 96.1 | 1470.8 | 0.747 | 29.54 deg |
| x-high candidate | 45 | 7 | 5.968528 | 317.2 | 118,236 | 1.0149 | 82.7 | 1416.6 | 0.555 | 28.27 deg |

Field-magnitude profiles sampled at 16 points across all eight cells have normalized adjacent-phase cosine similarity 0.994/0.995 for x-low modes 3->3->2 and 0.993/0.784 for x-high modes 8->8->7. The 45-degree x-high assignment weakens near another mode: choosing mode 8 there gives alpha 0.205 Np/m, but that mode better continues the positive-fold branch (30 mode 7 -> 45 mode 8, field cosine 0.985). An **alternative x-low candidate** modes 4->4->4 has alpha 3.380/0.953/0.839 Np/m and field cosine 0.948/0.993; it is retained in `d3p6_resume4_tracks.csv`, not silently dropped. The nearly flat modes 5->5->5 are a resonance/hybrid track; dividing by their near-zero slope gives a misleading 27.3 Np/m at phase 15, so they are not used for leakage gates.

For the x-low mean candidate, alpha(f) at 5.555/5.616/5.664 GHz is 0.624/1.047/1.539 Np/m. For x-high, alpha(f) at 5.969/6.006/6.055 GHz is 0.555/0.747/0.825 Np/m. **No sampled frequency interval establishes a single-moded tracked branch with alpha >= 1.2 Np/m**: x-low exceeds 1.2 only at one endpoint, x-high at none, and nearby modes remain. The upper-mode shortfall triggers the specified driven-strip flag; pin-launch coupling is unknown.

## FLAGS

- D3p6 has two x branches. The model-only 5.8-GHz predictions are 0.262-0.311 Np/m (x-low) and 3.032-3.308 Np/m (x-high) across theta 25-40 degrees. Both require HFSS mode matching.
- D3p6 phases 15 and 30 used separate fresh local AEDT 2025.2 sessions, 5.4 GHz minimum, 12 requested modes, a 3600 s exact-own-PID watchdog, and 60 s heartbeats. The crowded 30-degree spectrum supports a phase-45 continuation.
- D3p6 phase 45 used the same setup in a third fresh local AEDT session. The field comparison uses only sampled |E| magnitudes, not a full complex-mode overlap; this limits certainty near the 45-degree upper-mode crossing. No new OFF-twin or driven-strip run was started.
- The prior unrotated PML40 spurious-alpha check for D3p6 branch 2 was -0.0008 Np/m across and +0.0746 Np/m along, under the 0.1 limit but material beside the extracted 0.555-0.825 Np/m upper-branch range.
- The OFF twin and its copper-area audit remain last. Modal leakage does not determine which mode a centre pin excites.

## FILES

`anchor_resume4_diag.csv`, `.json` and `scripts/research/sub6_anchor_resume4_diag.py` reproduce the no-HFSS diagnostic. `d3p6_resume4_tracks.csv`, `.json` and `scripts/research/sub6_d3p6_resume4_analyze.py` reproduce the modal table and retain alternate assignments. Per-phase `SUPER_B2p2_D3p6_ph{15,30,45}_PML40_resume4_*` files contain builds, raw modes, profiles, field samples, validation, and watchdog records. Private AEDT projects remain outside Git.

## HANDOFF

- STATUS: D3p6 bridged 15/30/45 complete; OFF twin pending common f_op.
- CELL JSON: source spec at 242f3fb; own geometry/boundary readbacks saved per phase.
- beta/k0, alpha: both candidate tracks in the table; x-high stays below 1.2 Np/m.
- M_t, loss ratio: not determined by modal leakage alone.
- Evidence: HFSS 2025.2 modes and sampled fields; prior beta-to-sheet chain MODEL_ONLY.
- Next request: driven-strip excitation test for the upper mode or a decision on f_op before OFF twin.
