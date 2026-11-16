# Disc-radius screen, island cell, Ls 5.5, Ws 2.4, SCREEN_C40B (MODEL_ONLY)
Same sheets (OFF_ELS10, ON_ELS), phases along 88/95/102, across 90/98/106, disc via inset 0.20, 36 rows. Screen against screen only (baseline = LS5p5_WS2p4_VI from the Ls sweep).

| disc r mm | D_OFF | D_ON | r = D_OFF/D_ON | contrast dB | b OFF along | Q ON along/across | RES/HYBRID 5.5-6.1 |
|---|---|---|---|---|---|---|---|
| 5.435 (baseline) | +0.0296 | +0.0799 | 0.370 | 8.6 | 1.155 | 455 / 242 | ON along HYBRID |
| 5.3 | +0.0258 | +0.0633 | 0.408 | 7.8 | 1.133 | 539 / 260 | none |
| 5.2 | +0.0176 | +0.0501 | 0.351 | 9.1 | 1.120 | 512 / 276 | none |
| 5.0 | +0.0121 | +0.0446 | 0.271 | 11.3 | 1.100 | 710 / 316 | none |

Shrinking the disc lowers D_OFF and D_ON together, D_ON falls below the 0.08 floor, and r improves only at 5.0 (0.370 -> 0.271). Not monotonic (5.3 is worse than baseline). No disc radius gives D_OFF ~ 0 with |D_ON| >= 0.08.

Board model check (`results/sub6_board_final/rerank/grid/`, R 180/198, N 150-250, caps 0.3/0.5/0.8, balance <= 2 dB): baseline cell best worst-case SLL -13.9 dB (R 198, N 200, cap 0.3, balance +1.6). Smaller discs (alpha scaled by (D_ON/0.0799)^2) reach balance <= 2 dB only at N 250, cap 0.3, R 180, with worst SLL -13.1 dB (r 0.271) or -14.6 dB (r 0.19, optimistic, limited by two-beam SLL). The disc lever does not beat the baseline cell.
