# B-ALL-DIODE — HFSS 2026 R1 build only

## STATUS

**COMPLETE; NO SOLVE.** Copied the converged B-ALL-SHORT AEDT project to a new project and changed only the 250 slotted B-cell centre strips. All 16 gate discs remain plain, with no vias.

## GATES

| Audit | Result |
|---|---:|
| Slotted B cells modified | 250 |
| Sheets in the one `DIODE` impedance boundary | 250 |
| Area per sheet | 0.150000 mm² |
| Metal objects, before / after | 518 / 518 (516 cell discs + ground/feed) |
| Total metal volume, before / after CAD-history purge | 5071.486800 / 5071.486800 mm³ |
| Design validation | PASS |
| External STEP references in saved project and archive | 0 / 0 |
| Results and STEP files in archive | 0 |

Each cell's 0.5-mm central strip has a 0.3-mm cut along current across its full 0.5-mm width. The diode sheet lies at the top copper face, z = 0.035 mm, and meets the two copper strip halves at the cut edges. The cut and sheet rotate together with the original cell. The copper object stays one connected solid per site; the removed volume readback is 0.00525 mm³ per cell. The other geometry, feed, absorber, radiation boundary, infinite sphere, Setup1 and mesh settings read back unchanged.

The single boundary `DIODE` reads back `Resistance = Zr`, `Reactance = Zx` with 250 sheet objects. Default design variables are `Zr = 3.333ohm`, `Zx = -232.667ohm` (ON-A, surface ohms per square). Alternate settings: ON-B 37.05 / −559.0; OFF-D 2.667 / +48.667 ohm/square.

## FILES

Local AEDT project: `C:\lwmom_meander_model\hfss_work\full_board_gate_v31_diode\BOARD_R150_V31_B_ALL_DIODE.aedt` (9,898,315 bytes). Published archive: `BOARD_R150_V31_B_ALL_DIODE.aedtz` (1,385,112 bytes, SHA-256 `ED532E98AE14329F1E06A01EEF30668C31B6502678350E2DAE5BB8ECB0FAB162`). Its geometry is embedded: CAD import history was purged on all 522 imported bodies, and the archive contains no STEP or results. `audit.json` has every modified cell's mapping, area and volume readback; `validation.log` and `validation_embedded.log` record validation.
