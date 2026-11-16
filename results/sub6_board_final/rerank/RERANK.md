# Layout re-rank with the real cell (code only, MODEL_ONLY, no HFSS)
Engine `jobs1_iter1_study.py` via `board_sweep3.py` (area-uniform placement, rot180, k = 1, R = 180 mm cell area). Engine "off" = switched cells strong = real two-beam state (diodes ON); engine "on" = weak = real one-beam state (diodes OFF).
Anchor: alpha_cell = 2.16 Np/m at the strong state (cell-gate Delta_ON +0.1040), so alpha is unchanged. Weak-state amplitude ratio r = |Delta_OFF|/|Delta_ON| = 0.0273/0.1040 = 0.262 (also r = 0.335 old and r = 0.15 for comparison).

| N | cap | two-beam balance A-B dB | two-beam SLL dB | one-beam SLL, r = 0.335 | r = 0.262 | r = 0.15 |
|---|---|---|---|---|---|---|
| 100 | 0.5 | +9.3 | -18.3 | -18.8 | -18.8 | -19.0 |
| 100 | 2.0 | +15.7 | -20.1 | -19.9 | -20.0 | -20.0 |
| 200 | 0.5 | +1.0 | -16.1 | -11.4 | **-13.8** | -17.7 |
| 200 | 1.0 | +4.2 | -17.8 | -14.6 | -17.0 | -18.2 |
| 200 | 2.0 | +7.1 | -18.7 | -17.8 | -18.5 | -18.6 |

Only N = 200, cap 0.5 meets balance <= 2 dB (+1.0 dB). There, one-beam SLL is -13.8 dB at r = 0.262; -15 dB needs r ~ 0.23 (linear interpolation between r = 0.262 and 0.15), i.e. Delta_OFF ~ 0.024 at Delta_ON 0.104, about 12.8 dB contrast, not 15 dB. Gap to target: about 1.2 dB of SLL, or Delta_OFF 0.027 -> 0.024. Model offsets to HFSS (+10.2 dB OFF / +7.9 dB ON SLL) are not applied.
