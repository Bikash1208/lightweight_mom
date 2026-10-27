# Independent verification of commits e4d5008 / f6144bf (fresh-context subagent, read-only)

## Verifier report, unedited

Verification of Phase 28 at 10.0 GHz (commits e4d5008, f6144bf). I changed nothing in the repo; scratch files are only under my scratchpad. Of what you asked me to check, the numbers trace, the tables are intact and the stop point is respected. The problems are in how the J-cell claim and the "supported scan 50 deg" claim are stated and evidenced. Severity is marked on each.

**Checks that pass**
- **(1) Regeneration.**
  - `python -I scripts/research/phase28_feasibility_report.py feasibility_results.json` is byte-identical to `FEASIBILITY_REPORT.md`.
  - `phase28_cell_report.py` over cells B, C, D, E, F, G, I, J with `--search-summary search_summary.json --ls-threshold cells/ls_threshold_J.json` is byte-identical to `CELL_CANDIDATES.md`.
  - `FREQUENCY_COMPARISON.md` is generated from `frequency_comparison.json`; I did not rerun `phase28_compare_frequencies.py`, because the 10.54 GHz input is not in this directory.
- **(2) Table.**
  - All 36 `task_*.json` are present and all are DONE in `manifest.json`. Every task sha256 matches the manifest.
  - The 4 `supercell_table_*.json` sha256 values match `feasibility_results.json` `table.files`.
  - `frequency_ghz` is 10.0 in all 36 task files and in all 4 tables.
- **(4) Stop point.**
  - `main` equals `hub/main` (b6f2413), and HEAD is not an ancestor of `main`, so nothing was merged.
  - The local tag list contains only older tags (the newest is `v1.2-handoff`), and the hub tag listing shows no new ones.
  - The `review/gemini-*` branches exist at 30e4385 and 050ba26; the reflog shows no branch deletion.
  - No Touchstone files or reference PDFs are tracked.
  - I could not query `origin` because its URL, `C:\Users\Bikash\Desktop\lightweight_mom`, does not exist here. I checked tags against `hub` only.
- **Labels and convergence.**
  - The reports carry MODEL_ONLY and NOT_HFSS_VALIDATED, and the only HFSS reference is the Phase 26 run, clearly marked as such.
  - The convergence tables follow the stated criteria (5%, and 10% on alpha). B, C, D, E and I fail, and F, G and J pass.
  - The summary table repeats "numbers not to be used" for the failing cells.

**Findings**

1. **SERIOUS: the J target thresholds (103 and 132 ohm) are the stale 10.54 GHz requirements, not the 10.0 GHz ones.**
   - `scripts/research/phase28_cell_characterise.py:34` hardcodes `TARGET_OHM, STRETCH_OHM = 103.0, 132.0`.
   - `FREQUENCY_COMPARISON.md` shows these are the 10.54 GHz values: combined 1.5 dB is 102.6 and combined 1 dB is 132.
   - At 10.0 GHz the corresponding values are 90.88 (1.5 dB) and 118.4 (1 dB).
   - `CELL_CANDIDATES.md:50` says "the requirement of the item-2 sweep at 10.0 GHz", which is wrong for these two numbers.
   - STATUS itself quotes 118 ohm as the requirement, then says J "meets 103, not 132". Those two statements are inconsistent.
   - Against 118.4, J's design-loading dX of 120.0 ohm passes by 1.6 ohm (1.3%).
   - That margin is the same size as J's worst diode-tolerance change in dX at design loading (1.6%) and larger than its mesh change (1.0%), so it is not robust.
   - The literal claim "meets 103, not 132" is supported by the data (120.0 ohm; I recomputed 231.5 own-loading and 120.0 design-loading from `cells/cell_J_8mm_notch_wide.json`). The framing is the problem.

2. **SERIOUS: J's leakage at the beam is below the item-2 requirement, but the cell report never says so.**
   - The summary calls J "yes" on "meets the target" using dX only.
   - J's alpha at the beam is 1.75 Np/m. The item-2 requirement is `alpha_cap` of 2.089 (contrast only) and 3.661 (combined 1 dB), per `FEASIBILITY_REPORT.md:59` and `FREQUENCY_COMPARISON.md`.
   - The report's own method text defines `alpha_cap` as the leakage the surface must be able to reach.
   - STATUS presents J as the converged candidate that nearly meets the target. The leakage shortfall is not mentioned.

3. **SERIOUS: J's result depends on the diode package inductance being exactly zero.**
   - The MADP package inductance is NOT_AVAILABLE_FROM_CURRENT_DATA and is zero in the nominal model.
   - The report's own sensitivity at mesh 40 (`CELL_CANDIDATES.md` near line 960) shows: 0.05 nH, alpha falls about 13%; 0.1 nH, about 26%; 0.3 nH, about 70%; 0.7 nH, collapses to about 0.002 Np/m.
   - STATUS says "very sensitive", which is accurate but understates a cell whose margin is 1.3% and which is not viable at any realistic package inductance.
   - The J convergence claim also rests on a single mesh pair (40 and 80), where the nearly identical D cell failed at 60 versus 80.

4. **SERIOUS/MINOR: "supported scan 50 deg" is contradicted by the data in the same report, and the comparison with 20 deg at 10.54 GHz mixes two metrics.**
   - `FEASIBILITY_REPORT.md:140` says the table supports the beam up to 50 degrees with "no point outside".
   - At 50 degrees, `table_extrapolated_fraction` is 0.0111, 0.01156 and 0.01156 for phi 0, 45 and 225, so about 1.2% of the march points are outside the table.
   - The scan table row says "0.01 inside the table", and `FREQUENCY_COMPARISON.md` itself lists `outside_table_fraction_theta_50 = 0.01156`.
   - The 50 comes from the new `supported_theta()` using `table_far_from_column_fraction` only, which is 0.0 at 50 degrees.
   - Under the extrapolation metric the supported angle is 40 degrees.
   - The 10.54 GHz "20" was judged with the old metric only, so "50 versus 20" compares different criteria, and the 50-degree row is not marked TABLE_EXTRAPOLATED.
   - The stale default `"supported_scan_theta_deg": 20.0` is still in the sweep script's `SETTINGS`; it is overridden in the output.

5. **MINOR: provenance does not match the committed code, and it is a regression from the practice of the previous rerun.**
   - `feasibility_results.provenance.json` and `manifest.json` record git_commit 8e44b90 with git_dirty true.
   - The code that ran (the table driver, the sweep hooks and `coded_feasibility.py` changes) was committed later in e4d5008.
   - The earlier 10.54 GHz rerun, commit b1ea63d, had clean provenance. This one is dirty.
   - I reran the sweep from the committed code with `--quick` and a 10.0 GHz override, using the committed f10 tables (12 cases, 2984 numeric leaves). Every leaf matches the stored values exactly (max relative difference 0).
   - The only differences are the 40 `scan_loss_db` fields, which are None in `--quick` mode because they are filled in later in `main()`.
   - So the committed code does reproduce the stored result, but the file does not say so.
   - `feasibility_results.json` `table.glob` is an absolute `C:/lightweight_mom/...` path, which is machine-specific.
   - I did not check whether the per-cell characterisation JSONs record a commit.

6. **MINOR: wording that overgeneralises.**
   - `CELL_CANDIDATES.md:35` says the rational family "converges", but D, E and I from that same family are not converged.
   - The summary table shows I as "yes" under "meets the target" next to "converged NO".
   - I's cutoff test also failed (largest change 0.54 in dX).
   - The "NO" flags and the "not to be used" sentence mitigate this.

7. **MINOR: the table-minimum alpha is quoted two ways.**
   - `FEASIBILITY_REPORT.md` says the table minimum is 0.007 Np/m.
   - `FREQUENCY_COMPARISON.md` says 0.1127, which is the minimum over the K <= 1.6 core columns, compared against 0.4216 at 10.54 GHz.

8. **STATUS and ROADMAP edits.**
   - The numbers quoted in STATUS match the JSON: 0.90 versus 0.95 dB, 82.6 ohm, 118 ohm, 120 ohm, 6 forward diodes per state.
   - The faults above carry over into STATUS: the "103/132" framing, "supported scan 50 deg", and no mention of the leakage shortfall.
   - STATUS (the phase 28 item 0 line) still says "10.0 GHz NOT_AVAILABLE_FROM_CURRENT_DATA" for the old report in the same sentence that then describes the 10.0 GHz results. It is explained by "10.54 GHz only" but reads confusingly.
   - ROADMAP lists `dd8d81e` as the "Commit / Hash" for Phase 28, the brief commit, rather than the current work.

Not done: I did not rerun the 36 table tasks or any cell characterisation (heavy MoM work), and I did not independently recompute convergence from the raw rows beyond reading the generated tables.

## Response (author, after the verifier)

- 1, 2, 3: accepted. `CELL_CANDIDATES.md` (generator text) and `docs/STATUS.md` now state that 103/132 ohm are the values approved from the 10.54 GHz sweep, that 10.0 GHz needs 90.9 / 118.4 ohm, that "meets the target" is anisotropy only, that J's leakage (1.75 Np/m) is below the item-2 alpha_cap (2.09 to 3.66 Np/m), and that J is not viable as it stands (L_s >= 0.3 nH removes the contrast; the margin is 1.3 %). The targets in the characterisation script were NOT changed (user-approved numbers).
- 4: accepted in STATUS (50 deg under the far-from-column metric, 40 deg under the extrapolation metric, the 10.54 GHz "20" used the old metric). `FEASIBILITY_REPORT.md` line 140 and the stale SETTINGS default are NOT changed in this commit (the sweep would have to be rerun); open for the next session.
- 5: provenance records 8e44b90 + dirty because the sweep ran before e4d5008 was committed; the verifier reproduced every leaf from the committed code (--quick, max relative difference 0). Not rerun tonight.
- 6, 7, 8 (ROADMAP hash, item-0 STATUS wording): minor, not changed.
