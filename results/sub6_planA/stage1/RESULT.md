# Stage 1 anchor — ABORTED

## STATUS

**ABORTED_CONTROLLER_LOST.** No RF row is accepted. All ten own AEDT process trees were stopped by their recorded exact PIDs after all ten Python controllers unexpectedly exited with code 1 while AEDT kept solving. The projects, build reports, and partial solver profiles are retained locally. No FFD, S11, or completed far-field metric exists.

## GATES

| Gate | Result |
|---|---|
| AEDT | PyAEDT reports 2025.2; own HFSS profile reports **HFSS Version 2025.2.0**. All instances were new local processes. |
| STEP import | 1,244 AEDT objects: 1,000 metal bodies corresponding to the STEP compound's 1,000 lumps, 240 D sheets, SUBSTRATE, SUBSTRATE_FR4, AIR, PORT. |
| Thickness edit | PASS for 0.508 mm trial: bias bottom −4.010 → −3.718 mm, bias top −3.175 mm fixed; all 480 via/pad bottoms −4.045 → −3.753 mm, top copper +0.035 mm fixed. Other requested thicknesses passed the same per-body readback. |
| Material/boundary/build | 10/10 built and read back. Top SUBSTRATE uses RT5880 (2.2/0.0009), bias layer vacuum or RT5880, AIR vacuum+radiation, metal PEC. One BAR90 impedance boundary covers 240 D sheets, 50 Ω PORT, 5.8 GHz setup, maximum four passes, Global sphere 2° over full theta/phi ranges. |
| AEDT validation | T5/OFF independently passed design validation with exactly one port; all ten solve workers reached `SOLVING` only after their own validation passed. |
| Solve/export | **0/10 complete; 0 FFD.** Controllers exited without final reports. The orphaned AEDT solves were stopped; no partial row was promoted to a result. |

## DATA

The furthest partial solve, T5/OFF, reached Adaptive Pass 3. Its initial meshing took 17m41s and logged 4,183,852 tetrahedra before coarsening; assembled pass meshes were 774,732 (Pass 1), 996,194 (Pass 2), and 1,280,176 (Pass 3). Completed convergence footnotes available before interruption were delta-S 0.229975 and 0.152284. These are **partial process data**, not RF results.

## FLAGS

- **Controller failure cause unproven.** All ten Python worker sessions exited code 1 without traceback or a final JSON update. Their own AEDT processes and HFSS solver children were still active. The coincident exit is suspicious; Stage 1 stopped under the user's hard-stop rule. No retry was attempted after this infrastructure failure.
- Ten concurrent matrix solves drove free RAM to about 28 GB and Windows swap usage to about 22 GB before the controller loss was discovered. Future scheduling must reserve more memory per row in addition to checking free CPU.
- The existing anchor contains 240 diodes, above the 200-diode goal; it was used solely as the specified switching baseline.
- `SUBSTRATE_FR4` remains the imported object name even when its assigned material is vacuum or RT5880. Its material and thickness were read back per build.
- The unspecified setup delta-S target was set to 0.02; this did not cause an accepted result. The 50 Ω port integration line was chosen radially on PORT from x=0.635 to 2.05 mm.
- Stage 2 is a plan only. No Stage 2 solve or board generation was run.

## FILES

- `preflight/preflight.json`, `edit_trial/edit_trial.json`, `via_trial/via_trial.json`: own AEDT import and edit evidence.
- `builds/`: ten per-run object, material, impedance, setup, and geometry readbacks.
- `ANCHOR_*/ANCHOR_*.json`: ten abort records with own PIDs and stop times.
- `partial_profiles/`: ten incomplete HFSS profiles, retained for audit only.
- `../STAGE2_PLAN.md`: requested diode-count ladder source review and inputs needed.

The private `.aedt` projects and `.aedtresults` remain under `C:\lwmom_labB\stage1_work` and are not committed.
