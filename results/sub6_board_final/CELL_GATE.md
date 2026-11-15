# Cell gate, full settings (MODEL: HFSS Eigenmode, 150 mm lid and floor, normal convergence, 20 modes from 4.0 GHz)
Cell LS5p5 + WS2p4 (island, disc via inset 0.202). OFF_ELS10 = 21.6 - j297.5, ON_ELS = 2.8 + j14.6 Ohm/sq. 16/16 rows COMPLETE.
Branches tracked per TARGETS.md (Hungarian + Q penalty); hologram branch "a" in all four cases. beta/k0 = phi58/80.8 at 5.8 GHz (phi58 interpolated in the bracketing pair of phases).
Delta := beta/k0(across) - beta/k0(along) (the T1 quantity).

| state | along phi58 | along beta/k0 | across phi58 | across beta/k0 | Delta | flags |
|---|---|---|---|---|---|---|
| OFF_ELS10 | 93.55 | 1.1577 | 95.75 | 1.1850 | +0.0273 | none |
| ON_ELS | 94.86 | 1.1741 | 103.27 | 1.2781 | +0.1040 | ON across SLOW_SW (slope 0.012 GHz/deg, Q 231) |

|Delta_OFF| = 0.0273, |Delta_ON| = 0.1040, ratio = 3.8.
**GATE |Delta_ON| < 0.5 |Delta_OFF|: FAIL.** The ON state is anisotropic about 4x more than OFF, the opposite of the wanted behaviour. T1 (|Delta| >= 0.14) fails in OFF too.
Other branch rules: T2 (|ON across - ON along| <= 0.02) = 0.104 FAIL; T3 pass (SLOW_SW flag on ON across); T4 pass.
Second SW branch at 5.8 GHz is EXTRAP (not bracketed) in OFF across and ON across; it was not used.
Per-phase hologram-branch frequencies (GHz): OFF along 80/85/88/95 = 5.274/5.481/5.612/5.849; OFF across 88/95/100/110 = 5.544/5.779/5.921/6.164; ON along = 5.244/5.460/5.524/5.805; ON across = 5.452/5.642/5.760/5.882.
Rows: cell_gate/*.json; full analysis: cell_gate_analysis.json.
