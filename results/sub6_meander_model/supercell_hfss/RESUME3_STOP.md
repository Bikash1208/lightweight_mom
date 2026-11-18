# Supercell resume: stopped at B2p2 45-degree timeout

## STATUS

STOPPED after the B2p2 D3p0 45-degree HFSS solve reached its specified 1500 s limit. No D3p6 solve or OFF-twin audit was started. Completed mode rows and model checks are preserved.

## GATES

- Lidded calibration and target-branch PML spurious-loss checks passed previously (see `HANDOFF.md`).
- The B2p2 45-degree solver timed out; its 15/30/45-degree branch, single-mode band, and ON alpha gate cannot be completed from valid rows.
- The common-frequency ON/OFF >=12 gate and joint operating frequency/beam angle remain undetermined. The 5.8 GHz frequency is not imposed as a gate.

## DATA

Anchor HFSS modal leakage, negative-fold branch. Beta/k0 is the unfolded bound-wave value; theta is the radiated first-harmonic angle. The group velocity uses the **92.8 mm supercell** slope.

| Cell | phase deg | mode | f GHz | beta/k0 | v_g Mm/s | Q | Q_base | alpha Np/m | harmonic theta deg |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Anchor | 15 | 4 | 5.752876 | 1.0997 | 176.165 | 54.45 | 1939.74 | 1.8311 | 32.56 |
| Anchor | 30 | 4 | 5.681668 | 1.0898 | 176.165 | 67.40 | 2013.42 | 1.4530 | 31.41 |
| Anchor | 45 | 3 | 5.594682 | 1.0827 | 176.165 | 86.03 | 2087.10 | 1.1119 | 30.35 |

B2p2 **candidate** negative-fold continuation, modes 5->5. Its provisional v_g uses only two phases and its mode identity is not confirmed by the missing 45-degree row.

| Cell | phase deg | mode | f GHz | beta/k0 | provisional v_g Mm/s | Q | Q_base | provisional alpha Np/m | harmonic theta deg |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| B2p2 D3p0 | 15 | 5 | 5.772310 | 1.0960 | 179.083 | 309.81 | 1860.01 | 0.2724 | 32.43 |
| B2p2 D3p0 | 30 | 5 | 5.691903 | 1.0878 | 179.083 | 414.39 | 1898.60 | 0.1884 | 31.35 |

The MODEL_ONLY anchor beta-to-sheet-to-alpha chain at 30 degrees gives 2.019/1.676/1.285 Np/m at 5.595/5.682/5.753 GHz versus HFSS 1.112/1.453/1.831 (model/HFSS 1.816/1.153/0.702). The old expected 1.158-1.480 Np/m range was calculated **at 5.8 GHz**, not 5.682 GHz. The frequency trend disagrees. Full theta=25/30/35/40 results for anchor, B2p2, and both D3p6 x branches are in `MODEL_CHECK.md` and `model_check.csv`.

## FLAGS

- B2p2 45-degree build, boundary readback, and design validation passed; no eigenvalue row was emitted. The independent watchdog recorded `TIMEOUT_1500S` and terminated only the own AEDT PID 18712 and its children 4008/4464/4504. The controller reported `RuntimeError: SOLVE_TIMEOUT`. PyAEDT then raised a secondary `NoneType.objectID` error while trying to release the already terminated desktop. Exact PID check found none of those processes alive afterward.
- The anchor Q_base phases are 88.125/86.25/84.375 degrees; the 45-degree baseline extrapolates 0.625 degrees below the 85-degree unit-cell endpoint. Its v_g is 176.165 Mm/s from the 15-to-45-degree supercell slope. These were checked first and do not explain the model's opposite frequency trend.
- Supercell alpha describes modal leakage; pin-launch coupling remains unmeasured. B2p2 provisional points cannot establish a single-mode alpha>=1.2 band. No operating frequency or beam angle is selected.

## FILES

- `MODEL_CHECK.md`, `model_check.csv`, `model_check.json`: model-only frequency and angle screening.
- `SUPER_B2p2_D3p0_ph30_PML40_resume1_*`: completed HFSS row and solver evidence.
- `SUPER_B2p2_D3p0_ph45_PML40_*`: validated build, timeout status, controller log, watchdog record; no modes JSON.
- Private AEDT projects remain under `hfss_work/rotated_supercell/projects/` and are not pushed.

## HANDOFF

- STATUS: STOPPED at 1500 s B2p2 45-degree solver timeout.
- CELL JSON: source spec at commit 242f3fb; own build readbacks saved.
- beta/k0: anchor and provisional B2p2 values above; D3p6 HFSS pending.
- alpha: anchor HFSS 1.112-1.831 Np/m over 5.595-5.753 GHz; B2p2 incomplete.
- M_t, loss ratio: not established by this modal run.
- Evidence: HFSS 2025.2 for completed modes; `MODEL_ONLY` for tensor-chain estimates.
- Next request: decision on the 1500 s limit or a reduced HFSS solve before any further supercell work.
