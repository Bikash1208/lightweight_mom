# Ls 6.0 mm island cell at full settings (PEC150, 20 modes) - MODEL_ONLY
Phases along 80/85/88/95, across 88/95/100/110; OFF_ELS10 and ON_ELS; 16 rows, none NOCONV. Delta = beta/k0 across - along. F2 hung and was restarted as F5 for OFF across 95 and 110.

| | along | across | Delta |
|---|---|---|---|
| OFF beta/k0 at 5.8 GHz | 1.1607 | 1.2161 (SLOW_SW) | +0.0554 |
| ON beta/k0 at 5.8 GHz | 1.1923 (EXTRAP) | 1.3223 (SLOW_SW) | +0.1300 |

Contrast 7.4 dB (default pairing), Q OFF 1562 / 804, ON 346 / 168. T1 0.055, T2 0.130, T3 ok, no RES/HYBRID hits in 5.5-6.1 GHz.
Alternative pairings: ON along branch b (beta/k0 1.0096, near the light line) with ON across a gives Delta_ON +0.313 and 15.0 dB, but that pairs two different branches and is not credible as one mode; ON across b is an extrapolated, near-vertical branch (Q ~ 8e5), not used.

Versus the screen (C40B): Delta_OFF +0.040 -> +0.055, contrast 10.6 -> 7.4 dB. Ls 6.0 is worse at full settings. Full-settings Ls 5.5 (cell gate): Delta_OFF +0.0273, Delta_ON +0.1040, 11.6 dB. No Ls reaches 15 dB or Delta_OFF ~ 0.
