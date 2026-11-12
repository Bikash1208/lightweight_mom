# SCREEN_C40 comparison (AEDT 2025.1)

**SCREEN_C40 settings:** PEC lid and floor 40 mm, 12 modes from 4.5 GHz, MaxDeltaFreq 0.1 %, MinimumConvergedPasses 2, MaximumPasses 15, 4 cores per job, tetrahedra count and passes recorded per row. Read-back in every row JSON (`setup_readback`): MinimumPasses reads **3** (4 was requested; AEDT did not accept 4), the other settings read as set.

## Step 3 repeatability gate: PASS (but see flag)
LS5p5 OFF across 95 deg solved twice on two fresh builds (`repeat_A`, `repeat_B`):
| Config | Branch 1 A | Branch 1 B | difference | worst mode f / Q difference (12 modes) |
|---|---|---|---|---|
| LS5p5 | 5.76081 GHz | 5.76081 GHz | 0.0000 % | 0.0000 % / 0.000 % |
| LS5p5 + via inset 0.20 mm | 5.77590 GHz | 5.77590 GHz | 0.0000 % | 0.0000 % / 0.000 % |
Limit 0.3 %. The two copies are bit-identical, so this proves solver determinism, not independence from the mesh. LS5p5 reached the 15-pass maximum (29,416 tetrahedra, not converged to 0.1 % twice); the via-inset run stopped at 14 passes (25,249). In this standard the via inset moves branch 1 by only +0.26 %, against -2.3 % in SCREEN_BOX40: that earlier drop was screen noise.

## Step 4 (island configs only): CS0-CS3 NOT RUN
The centre-shorted slot configs CS0-CS3 are not defined in anything I have; no CS prompt reached this session, so only BASE, LS5p5 and LS5p5 + WS2p4 (all with the 0.20 mm disc-via inset, `_VI`) were run: OFF and ON x (along 80, along 88, across 88, across 100) = 24 rows, all COMPLETE. 6 rows hit the 15-pass maximum (listed in the JSON). Only two phases per direction exist, so almost every beta/k0 at 5.8 GHz is extrapolated (`EXTRAP`), BASE's across branch by about 30 deg: treat every number as indicative. See `ranking_c40.csv` and `targets_analysis_c40.json` (all pairings, T1-T5).

| Config | beta/k0 OFF along / across | beta/k0 ON along / across | \|OFF delta\| (T1, needs >= 0.14) | ON \|delta\| (T2, <= 0.02) | T3 | T4 | T5 OFF across: Q, slope GHz/deg, v_g/c, att dB/m |
|---|---|---|---|---|---|---|---|
| BASE | 1.078 / 1.446 (a) | 1.141 / 1.026 | 0.367 pass | 0.114 fail | pass | pass | 385, 0.0201, 0.280, 4.9 |
| BASE (OFF along candidate b) | 1.217 / 1.446 | 1.141 / 1.026 | 0.229 pass | 0.114 fail | pass | pass | same |
| LS5p5 | 1.144 / 1.193 | 1.017 / 1.330 (or ON along b 1.174) | 0.049 fail | 0.314 (0.156 with b) fail | pass | pass | 949, 0.0312, 0.434, 1.3 |
| LS5p5 + WS2p4 | 1.140 / 1.187 | 1.163 / 1.277 | 0.047 fail | 0.114 fail | pass | pass | 1015, 0.0310, 0.432, 1.2 |

No config passes T1-T4; BASE passes T1, T3 and T4 (3 of 4). BASE's OFF across branch is `EXTRAP` and its slope at 5.8 GHz is 0.020 GHz/deg, i.e. at the SW limit. Step 5 (full-settings winner confirmation) was NOT started: it needs the best of all 7 configs.

## Files
`rows/` (24 rows), `repeat_A/`, `repeat_B/`, `modes_c40.csv` (336 mode lines), `ranking_c40.csv`, `targets_analysis_c40.json`.
