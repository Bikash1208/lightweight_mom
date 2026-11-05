# SIDEQUEST_4QUADRANT - concept record (no work started)

Status: **CONCEPT ONLY, MODEL_ONLY, nothing built or simulated for it.** Recorded 2026-11-05 for later; it waits for the HFSS check of the 4-diode (C1) cell and for the K2-default board (`hfss_exports/sidequest_board_K2_default/`).

## Concept
One centre probe, one 150 mm hexagon, the board split into four quadrants (sectors of the surface-wave field; each quadrant is a hologram with its own states). One control circuit (all diodes of the board are ON together or OFF together).

| quadrant | default (circuit OFF, unpowered) | circuit ON |
|---|---|---|
| top-right (Q1) | LP beam | CP beam |
| top-left (Q2) | silent | CP beam |
| bottom-right (Q4) | silent | CP beam |
| bottom-left (Q3) | flexible: silent -> CP, or LP -> CP, whichever works best | CP beam |

* When the circuit is ON the default (LP) beam of Q1 turns off.
* "Silent" = the quadrant radiates no beam (its surface-wave power is not converted: low-leakage cells, or power absorbed/guided to the edge).
* The CP beams of the quadrants can point to the same or to different directions (open: 4 separate beams vs. one combined 4-quadrant CP beam; the combined beam needs a common phase reference and gives the full aperture gain).

## Why it is not covered by the K2-default cell
CP needs two orthogonal current components with a 90 deg phase difference. With a scalar cell (both diodes of a cell switched together, `m_x = m_y`) the cell current is `J = m rho_hat E`, a single linear polarisation. **CP needs per-axis switching**, i.e. the C1 4-diode cell (x-axis diodes and y-axis diodes on separate control states) so that `m_x` and `m_y` can differ in sign/magnitude between the OFF and ON states. The HFSS diagnostics of the crossed two-diode cell (`results/sidequest_hybrid/CELL_INDEPENDENCE.md`: only both diodes together switch) are the open question; the C1 4-diode cell package (`hfss_exports/sidequest_rt5880_cell_4diode/`) is the check.

## What each quadrant needs (design rule per quadrant)
* Q1: state OFF = LP hologram, state ON = CP hologram at the same or a different direction: per-cell pair (`m_x`, `m_y`) OFF -> ON with one ON/OFF word; the LP default hologram must vanish (beam suppressed) in ON.
* Q2, Q4: OFF = silent (all cells in a non-radiating static level or a hologram that cancels), ON = CP hologram.
* Q3: free choice, decided by the quadrant that gives the better gain/diode count.
* The common line count stays one if all switched cells go OFF -> ON (family A, see `results/sidequest_hybrid/K2_DEFAULT.md`); cells that would need ON -> OFF are static compromises (isolation cost) or need family-B cells (ON raises Y).

## Questions to answer before building
1. Does the C1 4-diode cell switch each axis independently in HFSS? (blocking)
2. Quadrant isolation and edge scattering: surface-wave power that crosses the quadrant boundaries (the quadrants are not independent apertures: the feed field is shared).
3. Gain per beam: a quadrant is one quarter of the aperture (about -6 dB against the full board at equal efficiency) unless the CP beams are combined.
4. Is "silent" achievable by cell choice alone (static levels with the lowest leakage) or does it need an absorbing/guiding region?
5. S11 of the probe with all four quadrants loaded.

## Dependencies
K2-default board results (HFSS) -> C1 4-diode cell check -> then this concept. Nothing here is started; no runs are planned.
