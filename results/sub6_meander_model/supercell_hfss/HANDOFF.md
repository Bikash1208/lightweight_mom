# HFSS rotated-supercell restart handoff (2026-11-18)

**Latest status:** B2p2 30° completed; B2p2 45° timed out at the specified 1500 s limit. Work stopped. See `RESUME3_STOP.md` for the current result and model-only angle screen.

## STATUS
Paused for user-requested server restart. AEDT 2025.2, fresh local sessions only. B2p2 30° was stopped by exact recorded PIDs; B2p2 45° and D3p6 supercells did not start.

## GATES
- Lidded 85° calibration: anchor +0.164631%; B2p2 and D3p6 bridged <0.000001%; all within 0.5%.
- PML spurious α at 40 mm air on target branches, across/along (Np/m): anchor 0.0074/0.0181; B2p2 0.0155/0.0103; D3p6 0.0307/0.0071. All ≤0.1. A separate B2p2 6.30 GHz mode gives 0.119 and is flagged.
- Anchor α = 1.453 Np/m at 5.682 GHz, within the expected range there; **5.8 GHz anchor gate unproven**.
- B2p2/D3p6 leakage and ON/OFF ratio: pending.

## DATA
| Cell | φ_s | Mode | f GHz | Q | α Np/m | Max tets | Wall s |
|---|---:|---:|---:|---:|---:|---:|---:|
| Anchor | 15° | 4 | 5.752876 | 54.45 | 1.8311 | 62,624 | 1123.33 |
| Anchor | 30° | 4 | 5.681668 | 67.40 | 1.4530 | 73,052 | 1354.20 |
| Anchor | 45° | 3 | 5.594682 | 86.03 | 1.1119 | 46,101 | 931.36 |
| B2p2 | 15° | unresolved | 5.772/5.872/5.900 | ~310/292/312 | — | 97,833 | 1502.02 total |

Anchor |v_g| = 176.165 Mm/s from 15°–45°; Q_base at 30° = 2013.42 from own PML rows. A SaveFields 30° rerun reproduced mode 4 exactly. Its top-side |E| averages across eight discs, normalized to their mean: 1.19, 1.27, 1.06, 0.69, 0.97, 0.86, 0.84, 1.12. Nearby mode 3 concentrates in two cells (2.13 and 2.01), supporting the mode-4 surface-wave match.

## FLAGS
- Radiation top was rejected by AEDT Eigenmode. PML at 20 mm shifted anchor f by ~1.5%; 40 mm shifted it by −0.25%.
- AEDT caps Eigenmode setups at 20 modes. A 40-mode build failed before solve; fixed to 20 modes, minimum 5.3 GHz. Completed supercells reached above 6.5 GHz.
- Initial anchor runs did not save fields. One 30° rerun used SaveFields=True; future setups save fields.
- An external watcher supplies 60 s heartbeats and an exact-own-PID 1500 s timeout because PyAEDT's blocking solve stalls its thread heartbeat.
- Field export uses Eigenmode stored-energy sources 0/1. Failed exports left own sessions and locks; only exact logged PIDs and stale own locks were removed.
- Folded Q baseline interpolates 85°–95° PML rows; 45° endpoints extrapolate by 0.625°.
- B2p2 30° is ABORTED_USER_RESTART with no eigenvalue row. Open-twin area audit is prepared but not run.
- Supercell α gives modal leakage, not which mode a pin excites.

## FILES
This folder has PML checks, folded Q, per-row modes/builds/profiles/convergence/log/status files, field samples, and the exact-PID stop record. Private AEDT projects remain under hfss_work/rotated_supercell/projects/ and are not pushed.

## HANDOFF
- STATUS: PAUSED for server restart; no own AEDT solve active.
- CELL JSON: source spec at commit 242f3fb; local build readbacks here.
- β/k0: no new dispersion claim.
- α: anchor 1.831/1.453/1.112 Np/m at 5.753/5.682/5.595 GHz.
- M_t: not measured.
- Loss ratio: ON/OFF undetermined.
- Evidence: HFSS 2025.2 eigenmodes, PML baseline, saved-field samples.
- Next: resume B2p2 30°/45°, D3p6 15°/30°/45°, then open-twin geometry audit.

