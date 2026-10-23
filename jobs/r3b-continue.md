# WITHDRAWN: r3b-continue (decision of the user, 2026-10-23)

**Do not run this job.** The R3b optimisation is not resumed; the earlier content of this file (resume of the three Nelder-Mead starts with 40 workers, run-2 re-score, mesh 60/90/120 verification,
HFSS package for a new cell) is withdrawn and replaced by this notice (the old text is in git history: commit `c11abee`).

Reason (user): the objective was dominated by the theta-30 / phi-45 tolerance corners - the MoM's least trustworthy condition (staircased diagonals, no HFSS calibration there) - and was insensitive to the
five parameters; splitting the diagonal stub from the axis stub made it worse (249 cached evaluations: baseline start 27.4 -> 26.5, diag_long 49.6, diag_longer_wide 74.4).

Outcome: **the baseline CP-2 cell is retained** (HFSS-checked at 10.5 GHz); the R4 state library stays valid for it. See `results/phase20_r3b/OUTCOME.md`.

**The server currently runs only `jobs/equivalence-check.md`.** A new job request will be written when there is something to run.
