# Stage 2 diode-count ladder — plan only

No Stage 2 model or AEDT solve was run.

## Sources found

- `gh/compute/sub6-final-board-13pad:scripts/research/sub6_final_design.py` fixes the existing hexagonal N240 island-board design (`side=0.2 m`, `n=240`, placement `a`, taper `0.5`, exclusion enabled). It calls `sub6_board_model4.py`, which uses a layered-medium reference-wave march (`march_rays`/`cell_design`) and calibrated static and switched cell tables. The slot axis follows the desired current angle and radial reference-wave angle; the `flips` array adds 90° where the optimizer selects it.
- `sub6_final_layout.py` creates cell and via tables and 12 series strings of 20 switched cells. `sub6_final_13pad.py` and `sub6_final_13pad_export.py` reroute these to 12 drive pads plus a common return, with connectivity and clearance checks.
- `sub6_final_rf_std.py full_win` converts the existing `results/sub6_final_board/mesh96/fab/cell_table.csv` and `routing.json` into the N240 STEP with the 12 × 12 mm SMA window. It hardcodes the hexagonal outline, 0.8 mm FR4 bias dielectric, and via/barrel length.
- `gh/compute/sub6-jobs1-prep:scripts/research/jobs1_iter1_study.py` is the circular, partitioned A-static/B-switched **model** study. Its `Board.evaluate` separates A and B currents by a switched mask. It supports `k=1` and arbitrary N through `JOB_KS`/`JOB_NS`, and radii through `JOB_RADII`; `ROT180=1` swaps model beam azimuths to the HFSS frame. The checked-in `iter1_study.json` uses N40/50, k1/2/3/4/6, diameters 220/280/360 mm. I did not find a checked-in N80–150 result in the named branches.
- `gh/compute/sub6-iter1-island-retune:scripts/research/sub6_iter1_island_retune.py` retunes island-cell response, but does not export a circular board.

## Feasibility

The circular partitioned **design model** can be parameterized for N60/90/120/150/200 and k1. No existing script directly emits a circular, partitioned, routed RF STEP for those cases. The layout/router and STEP exporter need to be adapted to consume the circular study's mask, site positions, and slot axes; route the N switched island diodes with verified return/connectivity; replace the bias dielectric with RT5880 at the selected thickness; preserve the 12 × 12 mm SMA window; and re-run object, overlap, port, diode-sheet, and RF readback gates. This is development work, not a switch of an existing option.

## Inputs needed before Stage 2 board generation

1. Which circular diameters to generate for each N (all 220/280/360 mm, or a selected pairing)?
2. Which RT5880 bias thickness(es) to use after Stage 1, and whether to retain the 0.8 mm option.
3. Selection rule among candidate placement strategies and cap/taper values, or approval to choose by the stated OFF beam balance and ON B suppression targets.
4. Bias network constraints: maximum current per string, preferred number/placement of drive pads, connector placement, and whether the existing 13-pad/common-return scheme is acceptable.
5. Confirm that the N240 anchor's cell response may be transferred to a lower-count circular board, or supply/approve a recalibrated BAR90 cell table after Stage 1.

## FLAGS — choices not specified

The current sources default to a 25 mm central exclusion in the circular study, placement strategies `inner`, `outer`, `sector_axis`, `sector_perp`, `field_weighted`, `coherence`, static-cell cap grid 0.5–4.0, taper and feed-efficiency assumptions inherited from older models, and a 0.15 mm minimum copper clearance in the board router. None is yet selected for the requested ladder. The existing checked-in iteration-1 study is labelled model-only/not HFSS validated and originally used BAR50 and N40/50; BAR90 N60–200 predictions would need updated cell data. The N240 anchor exceeds the intended 200-diode maximum because it is a validation baseline.
