# Job: sub6-slotdisc-hfss-cell (received 2026-11-08, parts 1 and 2 as sent)

Branch compute/sub6-slotdisc-hfss-cell from compute/sub6-slotdisc-cell-scan (d79a2cb).

## PART 1/2
Build an HFSS EIGENMODE package for two cells, reusing the conventions of scripts/research/phase26_eigenmode_package.py / sidequest_via_cell_package.py (no frozen code changes):
- CELL1 RT5880 3.175 mm, p 11.6, gap 0.73, L/2r 0.8, W 2.9, notch 1.00 mm (anchor + notch).
- CELL2 RT5880 3.175 mm, p 11.6, gap 1.0, L/2r 0.9, W 2.9, notch 1.00 mm.
States per cell: OFF, ON (diode = one lumped RLC sheet 1.00 x 0.35 mm across the notch, centred in the rim) and NOCUT (closed slot, no notch).
HFSS lumped RLC is parallel, so use the 5.8 GHz parallel equivalents of the series models; check values:
OFF (0.15pF||5k)+0.6nH -> R 3.88 kohm || C 0.170 pF; ON 3.0ohm+0.6nH -> R 163 ohm || L 0.612 nH. Recompute and state them.

## PART 2/2
Master/slave phase along the slot axis (x) and across it (y) separately, phase points chosen so f crosses 5.8 GHz (expected beta/k0 1.05-1.21 -> about 85-100 deg per cell), mode 1 (TM-like). Deliverables in hfss_exports/sub6_slotdisc_cell/: STEP per cell/state (or one parametrised project), HFSS_SETUP.md (step-by-step for the user), states.csv, predictions.json (MoM beta/k0 along/across per state at mesh 96 where available, else 64, with the phase-vs-f the MoM predicts), results_template.csv (phase, f, Q per run). Geometry check json (notch width, rim, sheet size, no overlaps). Commit, push to hub and github. Then STOP. RESULT format with the list of runs the user must do in HFSS.
