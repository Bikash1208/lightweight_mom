# OFF twin critical-path resume (2026-11-19)

## STATUS

Geometry audit and eight-cell dry build passed. OFF-twin phi_s=15 supercell solve completed; phi_s=30 is next. No driven-strip test has started.

## GATES

- Geometry: PASS. AEDT 2025.2 read back two PEC halves per cell, each 43.539602 mm2, total 87.079204 mm2 against expected 43.5+43.5 / 87.08.
- OFF leakage: pending the 15/30-degree supercells and phase-matched unrotated Q baseline. Gate is alpha_OFF <= alpha_ON/12, approximately 0.13 Np/m, at the ON x-low 15-degree point 5.664 GHz.
- Driven strip: conditional on OFF passing.

## DATA

The old end-run formula stopped at x=5.230679 mm, while the rim at the slot's inner y edge is x=5.366563 mm. It left a 0.135884 mm copper bridge. The corrected cut extends to x=+/-5.700 mm, beyond the 5.6 mm disc. After subtract, AEDT listed one PEC multi-lump object of volume 3.047772 mm3. `separate_bodies` returned two PEC objects, 43.539602 mm2 each. The full eight-cell dry build read back eight such pairs, 16 separate copper halves, lattice X/Y, PML top, PEC bottom, RT5880 substrate, and 20 requested modes from 5.4 GHz. Design validation passed; no solve was run in the dry build.

OFF-twin phi_s=15: HFSS 2025.2 returned 20 modes, highest 6.696874 GHz; 20 adaptive passes, max 91,679 solved tets, 22m27s solution time, normal completion. Near the 5.664 GHz ON x-low point, OFF modes 2/3 are at 5.662908/5.663298 GHz with Q 8.2/3.0. Their low Q makes frequency-only matching unsafe; the phi_s=30 row and saved fields will determine which, if either, follows the ON surface-wave branch. All 20 raw f/Q values are in the per-row modes JSON.

## FLAGS

- The prior 0 mm2 result came from reading a stale/disconnected object handle instead of enumerating the copper bodies after the edge-opening subtract. The end-run geometry also needed the rim-reaching correction above.
- The OFF sheet is an open-gap twin bound; it has no diode impedance. Its Q baseline must be obtained from own 2025.2 unrotated PML40 rows before applying the alpha gate.
- No B2p2, anchor bracketing, beam-angle, or f_op selection is included.

## FILES

`D3p6_OPEN_area_audit_resume5.json` and its log contain the single-cell object/volume audit. `SUPER_D3p6_OPEN_ph15_PML40_resume5_dry_build.json` and status/log/validation contain the eight-cell dry-build readback. Private AEDT projects are not pushed.
