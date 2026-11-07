# Job: sub6-slotdisc-cell-scan (MODEL_ONLY)

Setup: new worktree C:\lwmom_slotdisc, branch compute/sub6-slotdisc-cell-scan from sub6-splitdisc. Full worker pool (RAM-limited only),
single-thread BLAS, checkpoints, provenance. Do not touch session 1's grid.

Cell (user decision 2026-11-07): the Phase 26 slotted circular patch (src/dispersion/slotted_cells.SlottedCircularPatch:
disc radius r, centred rectangular slot L x W, rotated with the cell) PLUS a rim cut: at the +u end of the slot, a notch of width
g_c = 0.65 mm through the remaining rim to the disc edge. BAR50-02V diode bridges the notch: RLC sheet 0.65 mm across x 0.35 mm
along the rim (radial), centred in the rim, current across the notch. The ON state is close to the closed-slot Phase 26 cell.
Rim thickness (r - L/2) >= 0.4 mm. If SlottedCircularPatch cannot take the notch without changing frozen code, build it with
PolyComposite / Poly and say so.

Step 0, anchor (do first, report before the grid): take the Phase 26 10.54 GHz cell that went to HFSS (its exact r, L, W, period,
from the Phase 26 records), scale all lengths by 10.54/5.8, closed slot, no notch, RT5880 3.175 mm. Report its tensor, Delta,
beta/k0 and the closed-form leakage alpha at 25 deg (leakage_alpha_first_order_tensor), next to the original 10.54 GHz values.
This is the "worked in HFSS" reference.

## Part 2 (continued)
Grid, 5.8 GHz, mesh 64, F = 2n:
- Substrates: RT5880 3.175 mm with p = 11.6 mm (scaled anchor); RO3006 (take eps_r / tan d from the repo's existing RO3006 work,
  state the values) at its standard thickness closest to 2.5 mm, with p = 10 and 11.6 mm.
- r: disc-to-neighbour gap p - 2r in {0.4, 0.6, 0.8, 1.0} mm. L / (2r) in {0.5, 0.6, 0.7, 0.8, 0.9} (rim rule permitting).
  W in {0.4, 0.8, 1.2, 1.6} mm.
- States: NOCUT (closed slot, no notch), OFF_OPEN, OFF_BAR50 (0.15 pF + 0.6 nH), OFF_BAR50_HI (0.20 pF), ON_BAR50 (3.0 ohm + 0.6 nH),
  ON_BAR50_MAX (4.5 ohm).
- Tensor at kt = beta along and across the slot axis (beta of that state, self-consistent), psi = 0.
Outputs per row: Delta (gap frame), antisymmetric part A, mean, beta/k0, rho = Delta_ON/Delta_OFF, ON - OFF change of A,
and the closed-form leakage alpha at 25 deg for OFF and ON (the key number). Flag rows near resonance (5.0-6.6 GHz).
Gates (named states): ANCHOR_REPRODUCED (scaled anchor alpha within 2x of the 10.54 GHz value scaled by 5.8/10.54) or not;
LEAKAGE_OK where alpha_leak >= 1 Np/m in the state that carries the beam, else LEAKAGE_LOW; MESH check 64/96 on the 3 best rows.
Pilot one row first: report the solve time and projected total; if the projection is over 24 h, stop and ask.
Report results/sub6_slotdisc/CELL_SCAN.md (gates first, best-leakage table, rho and A per row), commit, push to hub, STOP.

## Amendment 2026-11-07 (user)
- RO3006: single 1.28 mm laminate (the thickest standard single sheet; flag in the report that it is about half of the 2.5 mm asked for). NO 2-layer 2.56 mm stack.
- ADD Rogers TMM6, single laminate: eps_r 6.3 (design Dk), tan d 0.0023 (Rogers TMM6 page), thickness 2.54 mm (0.100 in), p = 10 mm only.
  Same r / L / W / state grid as RO3006.
- Anchor: hex150_tensor_cp_single majority cell (D 6.0, slot 4.8 x 1.6, p 6.4 -> scaled p 11.63 mm). Accepted.
- SlotDiscNotch (duck-typed geometry, in the script) accepted; it must reproduce SlottedCircularPatch exactly without a notch (kept as a test).

## Amendment 2 2026-11-07 (user: grid correction, W was not scaled from the anchor)
- W set {1.6, 2.4, 2.9, 3.6} mm (anchor slot 1.6 x 1.82 = 2.9 mm wide; 0.4 and 0.8 dropped). L/(2r) set {0.6, 0.7, 0.8, 0.9} (anchor 0.8).
- Gap p - 2r set {0.4, 0.6, 0.8, 1.0} plus 0.73 (anchor scaled). Rim r - L/2 >= 0.35 mm (was 0.4). Rows closest to the anchor first (W 2.9, L/2r 0.8, gap 0.73), then outward.
- Script deviation (flagged): rows whose slot-corner rim r - sqrt((L/2)^2 + (W/2)^2) is below 0.2 mm (about one mesh-64 cell) are skipped.
- Notch width: user wrote "notch 1.00 mm still apply"; the job text has g_c = 0.65 mm. Clarification requested; script option SUB6_NOTCH_MM (default 0.65).
