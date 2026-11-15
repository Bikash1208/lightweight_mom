# Island cell slot-length sweep at SCREEN_C40B (MODEL_ONLY, HFSS eigenmode, 5.8 GHz)

Slot width 2.4, island Li = Ls - 2.0, gap 0.30, clearance 0.60, island via on the island centre-line, disc via inset 0.20.
Sheets OFF_ELS10 (21.6 - j297.5) and ON_ELS (2.8 + j14.6). Phases along 88/95/102, across 90/98/106 deg. 84 rows, none NOCONV.
Delta = beta/k0(across) - beta/k0(along). Contrast = 20 log10(|D_ON|/|D_OFF|). Q = hologram branch, OFF along/across, ON along/across.

| Ls | b OFF along | b OFF across | b ON along | b ON across | D_OFF | D_ON | contrast dB | Q OFF al/ac | Q ON al/ac | RES/HYBRID 5.5-6.1 GHz |
|---|---|---|---|---|---|---|---|---|---|---|
| 5.5 | 1.1551 | 1.1847 | 1.1780 | 1.2579 | +0.0296 | +0.0799 | 8.6 | 1577 / 926 | 455 / 242 | ON along HYBRID (6.10 GHz at 88 deg) |
| 6.0 | 1.1565 | 1.1968 | 1.1793 | 1.3151 EXTRAP | +0.0402 | +0.1358 | 10.6 | 1605 / 842 | 335 / 165 | none |
| 6.5 | 1.1618 | 1.2126 | 1.2157 | no branch | +0.0508 | n/a | n/a | 1579 / 747 | 167 / - | ON across HYBRID |
| 7.0 | 1.1646 | 1.2395 | 1.1039 | no branch | +0.0749 | n/a | n/a | 1511 / 658 | 181 / - | ON along and across HYBRID |
| 7.5 | 1.1742 | 1.2718 | 1.1314 | no branch | +0.0976 | n/a | n/a | 1308 / 561 | 315 / - | ON along HYBRID |
| 8.0 | 1.1940 | 1.3027 | 1.1455 | no branch | +0.1087 | n/a | n/a | 986 / 473 | 460 / - | none |
| 8.7 | 1.2537 SLOW_SW | 1.3959 EXTRAP | 1.1504 | 1.0142 EXTRAP | +0.1422 | -0.1362 (not trusted) | -0.4 | 541 / 388 | 560 / 137 | none |

## Zero crossings
- **D_OFF** is positive and rises monotonically with Ls (+0.030 to +0.142). It does not cross zero in 5.5-8.7 mm. A linear extrapolation of the 5.5-7.0 points reaches zero near Ls = 4.6 mm (Li = 2.6 mm), outside the swept range and not simulated. At Ls = 5.5, the nearest point, |D_ON| = 0.080.
- **D_ON**: no crossing. It is defined only at 5.5 (+0.080) and 6.0 (+0.136, EXTRAP). From 6.5 to 8.0 the ON across branch does not exist at 5.8 GHz (HYBRID, no SW candidate). At 8.7 the ON across point is an extrapolation from 6.08-6.6 GHz data with T3 failing, so its sign is not reliable. The sign change between 6.0 and 8.7 is a branch change, not a measured crossing.
- Target (D_OFF ~ 0, |D_ON| >= 0.08, contrast >= 15 dB): not reached. The best contrast is 10.6 dB at Ls = 6.0, where D_OFF is +0.040, not zero, and the ON across point is extrapolated.

Files: `ls_sweep_summary.json`, 84 row JSONs, this file. Analysis uses `targets_analysis.analyse` (T1-T5 rules) on the rows.
