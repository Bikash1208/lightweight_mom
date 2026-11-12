# BAR90 island tuning — base build hard stop

## STATUS

STOPPED before any eigenmode solve. The user approved moving the island via's local x from −1.78 to −1.34 mm and putting its local y on each island's centre-line. That resolved the Round 1 via clearance conflict. The first base build attempt hit a script-side pre-solve guard error; the guard was fixed and the base build retried once. AEDT 2025.2 then rejected the specified bottom radiation boundary (`AssignRadiation`, `RAD_BOTTOM`). The required boundary could not be read back, so no substitute boundary was used and no solve was started.

## GATES

| Gate | Result |
|---|---|
| Via-pad geometry | PASS for all seven Round 1 configurations after the approved move; see `geometry_gate.csv` |
| G1, all model/variable/boundary readbacks | FAIL: AEDT rejected `RAD_BOTTOM` before design validation or solve |
| G2, clean SW branch | NOT RUN |
| G3, T1–T4 ranking | NOT RUN |

## DATA

The exported local via was (−1.78, −0.53) mm with pad radius 0.30 mm ([sub6_final_layout.py:29,35](../../../scripts/research/sub6_final_layout.py)). The approved local x = −1.34 mm leaves a minimum **0.11 mm** pad-to-island-end margin at Lᵢ = 3.5 mm. For every configuration, local y is `(c−g)/2`, the transverse island centre-line. The minimum transverse margin is 0.35 mm at g = 1.0, c = 0.6. No Round 1 configuration fails the analytic pad check after this move.

The base geometry was created in a fresh, local AEDT 2025.2 instance with private temp under `C:\lwmom_labB\aedt_temp\bar90_tune\BASE`. It used the specified disc, slot, island, two via barrels/pads, radius-0.50 anti-pads, 0.254 mm RT5880 bias layer, no stubs, 0.30 × 0.60 mm diode sheet, PEC lid at z = +150.035 mm, and air bottom 30 mm below bottom copper. The diode boundary, two lattice pairs, and PEC lid were assigned. At the bottom face, the PyAEDT call to assign `Radiation RAD_BOTTOM` raised `AEDTRuntimeError: Failed to create boundary Radiation RAD_BOTTOM`; AEDT's gRPC call was `AssignRadiation`. No design validation or eigenmode analysis was reached. The failed local project was not used for results.

| Attempt | Phase reached | Outcome |
|---|---|---|
| 1 | Fresh project opened, before geometry build | Script guard incorrectly treated AEDT's empty `.aedtresults` folder as a solved project. Logged; guard changed to check actual existing analysis sweeps. |
| 2 | Base geometry and boundary assignments | Bottom radiation assignment rejected by AEDT. Hard stop before `validate_full_design` and before solve. |

## FLAGS

1. The first failure was a script/API error before solve, handled under the allowed single retry. The second failure is an unresolved required-boundary error; changing to PEC, impedance, or a different air-box topology would change the requested model and was not attempted.
2. No `GetMessages` call was made. The own AEDT process tree from the failed retry was stopped by exact recorded PID and creation-time match. No other AEDT process or project was touched.
3. There are **zero completed eigenmode rows**. No sensitivity table, Round 2 ranking, winner check, or BAR90 tuning conclusion is available.
4. Local controller logs, status JSON, and 60-second heartbeat are under `C:\lwmom_labB\bar90_tune_work`; the failure summary is included in `build_abort.json`. No `.py` file was pushed.

## FILES

`results/sub6_cell_param/bar90_tune/RESULT.md`, `geometry_gate.csv`, `build_abort.json`.
