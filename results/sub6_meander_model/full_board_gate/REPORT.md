# Full-board gate validation — Step 1 feasibility stop

**STATUS:** STOPPED_STEP1. No AEDT project was built or solved; the requested three-state HFSS comparison is unavailable.

**GATES:** The prior R100 crop contained 84 active sites and peaked near 17 GB; its all-static solve took 346 s. The requested R150 aperture has 516 active sites, a 6.143× increase. A linear extrapolation is 104 GB and 35 min; a plausible direct-solver 1.5-power extrapolation is **259 GB and 88 min**. The host has 255.6 GiB physical RAM but only **176.1 GiB available** at the check, and each requested solve has a 60-min timeout. The superlinear estimate exceeds available RAM and the timeout, so a safe full-board fit is not established. The older 856-site final-board estimate of about 110 GB used a different cell/layout/mesh and cannot validate this new board's resource use.

**DATA / build blocker:** The referenced source on `gh/compute/sub6-cell-param-lab-b`, `results/sub6_board_final/hfss_build/`, is an R198 cell-area board with 856 sites, static plus island cells, theta 25°, rot180 frame, and an existing STEP. Its exporter `scripts/research/sub6_final_rf_std.py` reads that board's `mesh96/fab/cell_table.csv` for per-site slot dimensions and axes. The current R150 dither model defines 516 A/B site labels and ideal coherent phase weights, but no per-site physical rotation/phase table for identical anchor slots at 5.68 GHz and theta 30°. Simply replacing copper in the old STEP would simulate a different hologram. The requested BOTH-state beam gate would therefore be uninterpretable even if the mesh fit.

**FLAGS:** The model's A/RHCP and B/LHCP assignment is ideal array-factor input, not a calibrated anchor-cell CP response. The memory and time bounds are estimates, not a claim that HFSS measured 259 GB on this board. AEDT v261 executable is installed, but no AEDT process was launched. No other session or project was accessed, and no process was killed.

**FILES:** `FEASIBILITY.json` records the inputs and extrapolation. Source evidence: `gh/compute/sub6-cell-param-lab-b:results/sub6_board_final/crop_testbed/RESULT.md`, `.../hfss_build/RUN_README.md`, `.../hfss_build/manifest.json`, and `scripts/research/sub6_final_rf_std.py`; current branch `scripts/research/sub6_sector_dither_grid.py`.

**HANDOFF:** STOPPED_STEP1 | full R150 board | no HFSS result | physical anchor-cell hologram layout and a safe full-board memory plan are required before a three-state solve.
