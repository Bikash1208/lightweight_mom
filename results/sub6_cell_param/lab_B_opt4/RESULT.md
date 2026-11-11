# Option 4 decision result — corrected G-DISC

## STATUS

Corrected G-DISC 56/56 eigenmode rows solved on local AEDT 2025.2, zero FAILED rows. G-ISL-a, G-ISL-b, and G-ISL-noline are closed by the 2026-11-11 user decision; their earlier 168 rows remain for provenance, with no further solves or lid checks.

## GATES

- Both G-DISC native builds passed design validation. The disc via is at (-3.0, 0) mm, with 0.6 mm pad, 0.3 mm drill, 1.25 mm slot clearance and 2.135 mm rim clearance.
- Top and bottom 0.40 mm channels enclose the 0.10 mm island bus with 0.15 mm clearance per side. Both line pieces are united with ISLAND; neither is united with DISC.
- Wrapped DC graph without diode: `island(n) ↔ island(n±1)` through the periodic y faces; `disc(n) → ground` through its via. There is no island/bus to ground path without the diode. The full graph is recorded in each build JSON.
- Every accepted row matched intended design variables on AEDT readback; each row was eligible for one retry. No retry failures or FAILED rows occurred.

## DATA

The seven phases are 70, 75, 80, 85, 90, 95, and 100 degrees, for OFF/ON and along/across. The CSV reports all eight tracks per state and direction, at every phase, with Q, slope, class, vg/c, and every bracketed β/k₀ candidate. No Q≥1000 mode was excluded.

| G-DISC island length | OFF nearest RES/HYBRID | ON nearest RES/HYBRID | Bracketed β/k₀ range at 5.8 GHz, OFF along/across | Bracketed β/k₀ range, ON along/across |
|---|---|---|---|---|
| 4.348 mm | None among eight exported tracks | 5.3455 GHz, Q 431.1, distance 0.4545 GHz (HYBRID) | 0.9353–0.9980 / 0.8668–0.9983 | 0.8713–0.9982 / 0.8669–0.9998 |
| 2.500 mm | None among eight exported tracks | 5.2323 GHz, Q 469.5, distance 0.5677 GHz (HYBRID) | 0.9351–0.9979 / 0.8681–0.9983 | 0.8708–0.9981 / 0.8681–0.9996 |

## FLAGS

- The ON near-resonance distances are below the 1.0 GHz target if those HYBRID tracks are physical. No lid check was requested for corrected G-DISC, so box-mode rejection and definitive target decisions remain open.
- Multiple branches bracket 5.8 GHz on each axis. A single β/k₀, signed anisotropy, and ON/OFF anisotropy ratio are unresolved without branch identification. No extrapolation was used.
- The eight-mode eigenvalue export starts from the 4.5 GHz setup minimum; absence of a track in that export is not proof that no other mode exists.
- The initial G-DISC topology was rejected before any G-DISC row was pushed. Its local files were archived; corrected projects were rebuilt and validated.

## FILES

- `builds/G-DISC_L4.348_PEC150.json`, `builds/G-DISC_L2.500_PEC150.json`: corrected native build and wrapped DC graph.
- `rows/OPT4_G-DISC_*.json`: 56 AEDT result/readback records.
- `analysis/every_mode_track.csv`: all 448 corrected G-DISC mode points.
- `analysis/summary.json`: machine-readable brackets, nearest tracks, and unresolved metric flags.
