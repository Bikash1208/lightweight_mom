# Job: sub6-board-model (received 2026-11-08, part 1 of 2; part 2 pending)

Branch compute/sub6-board-model from compute/sub6-slotdisc-hfss-cell.

## PART 1/2
MODEL_ONLY, PENDING_HFSS_CELL_CHECK. Phase 26 tensor engine (march_rays/cell_design, unmodified, as in sub6-board-requirements), 5.8 GHz, RT5880 3.175, p 11.6, hexagon side 200 mm (also 150 mm), centre probe (launch fraction from src/dispersion/pin_launch.py at 5.8 GHz, state it). Beams A (25,0) and B (25,180), RHCP.
Cell tables from results/sub6_slotdisc/rows_RT5880.jsonl: STATIC cells = NOCUT rows (no diode; amplitude via L/2r, axis psi); SWITCHED cells = notch rows, OFF_BAR50 and ON_BAR50 tensors (both, with their own beta: OFF ~1.05, ON 1.11-1.14).

## PART 2/2
(not yet received; work not started)
