# SD gate eigenmode run — ABORTED by user

## STATUS

ABORTED_BY_USER on 2026-11-11. Gate concept closed. No additional gate solves authorized. Option 1b remains closed.

## GATES

All accepted rows passed exact intended-variable/AEDT readback and eight-mode export. The 40, 150 and 300 mm PEC-lid builds passed native validation; AEDT returned `2025.2`. The active queue and its one remaining worker/AEDT instance were stopped by exact PID after the user stop.

## DATA

150 planned unique rows: 148 accepted, 1 FAILED after both permitted 30-minute attempts, 1 ABORTED during its second attempt. Finished rows from the interrupted 300 mm batch were validated and retained. The 150 mm gate sweep has 96 accepted rows; the 150 mm low-phase supplement has 6; the 40 mm lid has 24; the 300 mm lid has 22 accepted rows.

## FLAGS

- `GATEEIG_SD_w2p9_ylines_PEC300_ON_across_180`: FAILED after two 30-minute timeouts. No third attempt was made.
- `GATEEIG_SD_w2p9_ylines_PEC300_OFF_across_180`: ABORTED on user stop after one 30-minute timeout; its second attempt was interrupted. No result was accepted.
- The 100° C-variant point was run once while satisfying both requested phase ranges. C variants were solved only at the 150 mm lid. No branch or performance conclusions are assigned to this closed concept.

## FILES

- `rows/GATEEIG_*.json`: 148 accepted result/readback records.
- `FAILED_ROWS.json`, `ABORTED_ROWS.json`, `row_manifest.json`: complete accounting of 150 planned rows.
- `status.json`: machine-readable ABORTED state.
- `parallel_reports/*.json`: batch CPU allocations, including the interrupted batch.
- `builds/*.json`: native build summaries for 40, 150 and 300 mm PEC lids.
