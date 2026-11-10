# Initial-suite environment diagnosis (2026-11-10, S2)

STATUS: ENV_ONLY_ACCEPTED_BY_DECISION_SESSION. No src/ or tests/ changes
during this diagnosis. No retune computations, commit, or push.

## Environments

| | Previous S2 | New isolated retune |
|---|---|---|
| Interpreter | C:/lightweight_mom/.venv313/Scripts/python.exe | C:/lwmom_retune/.venv313/Scripts/python.exe |
| Python | 3.13.17 | 3.13.17 |
| NumPy | 2.5.3 | 2.5.4 |
| SciPy | 1.18.1 | 1.18.1 |
| pytest | 9.1.1 | 9.2.0 |

Previous interpreter identified explicitly in results/jobs1/STATUS.md.
No venv/conda directory was found at the top level of C:/lwmom_final;
no jobs1 provenance/manifest filename was found. Previous interpreter used
only as an existing runtime; work directory and outputs remain in retune.
PYTHONDONTWRITEBYTECODE=1; OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=1.
Neither environment was modified during this diagnosis.

## Four-test rerun in previous environment

3 failed, 1 passed in 22.00 seconds. Polygon test passes under pytest 9.1.1
and fails with nested-data TypeError under 9.2.0: ENV_ONLY for that failure.
All three requested S11 tests still fail SHA256 equality.
Log: previous_env_four_tests.log.

## Numerical comparisons

Captured six complete 2x2 S11 matrices for ideal and MADP kits in states
0000, 0100, 1100, using the existing regression generator with a research-only
wrapper that adds raw values to its fingerprint outputs. No solver changes.

For all six, previous-env vs new-env max absolute S11 difference = 0;
max relative S11 difference = 0 (identical complex entries).

| case | abs difference of Frobenius norms / stored reference norm |
|---|---:|
| ideal 0000 | 3.1489213877503463e-16 |
| ideal 0100 | 1.577008039616045e-16 |
| ideal 1100 | 9.476506348929280e-16 |
| MADP 0000 | 3.1551331050331960e-16 |
| MADP 0100 | 7.9774410351954930e-16 |
| MADP 1100 | 3.2219830309626387e-16 |

All six S11 fingerprints differ from the stored reference in both envs.
The initial maxfail run reported only three of these.

**Max |dS11| and max relative S11 difference against the stored reference:
NOT_AVAILABLE_FROM_CURRENT_DATA.** tests/data/rotated_load_main_reference.json
stores only shape, SHA256 and scalar Frobenius norm. Its generating script
explicitly saves fingerprints without arrays. The reference was introduced
in commit 5c446866992e320adc255ca08a9206d8b2a5b40b in that form. Repository
and C:/hmsa_sync filename searches found no matching raw reference arrays.

A small norm difference is not an upper bound on entrywise complex errors;
under the initial criterion no ENV_ONLY classification or <=1e-9 claim
was justified for the S11 failures. The revised acceptance below supersedes
that initial block without asserting an entrywise error bound.
The two environments reproduce the available evidence equally well.

## Required to resume under the user's decision rule

Obtain the original complex S11 entries corresponding to the stored hashes,
or an explicit revised numerical acceptance test from the decision session.
No alteration of the stored fixture or tests is proposed.

## Revised acceptance (decision session, 2026-11-10)

User explicitly accepts platform hash differences as ENV_ONLY when old/new
environments produce bit-identical results and every relative norm difference
from the stored norm is <=1e-12. Both conditions hold: elementwise difference
0 for all six S11 matrices, maximum relative norm difference 9.48e-16.
ENV_DIAGNOSIS_PASS. The polygon pytest API failure is also ENV_ONLY.
Use the previous S2 environment with one BLAS thread per worker. The four
known initial-suite ENV_ONLY failures do not block the job's commit; record
them in its message. Do not change tests or fixtures. Later improvement:
store complex reference entries, not only hashes. This criterion does not
claim an entrywise error bound against the unavailable original entries.

Files: env_previous.json, env_new.json (complete raw S11 matrices and versions),
env_comparison.json, previous_env_four_tests.log. Raw numbers are MODEL_ONLY /
NOT_HFSS_VALIDATED; they are diagnostics of the CP2 regression, not SUB6 results.
