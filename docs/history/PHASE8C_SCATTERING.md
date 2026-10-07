# Unit-cell scattering architecture & control-quantity validation (Phase 8C)

Determines, from the *real* Phase 7 scattering data (not a filename, not a doc citation alone),
which physical quantity — if any — can serve as the holographic control variable. Read-only with
respect to Phases 1-7; no GUI, no final MTS STEP, no aperture calculation.

```
results/unit_cell_library/_solve/database.json (Phase 7's own full linear-basis S-matrix,
        read-only, still present -- xx/xy/yx/yy for every family, not just the "_xx" the
        delivered library exports)
        |
        v
topology_audit.py --> power_balance.py --> scattering_analysis.py (S11xx AND S21xx, full library)
        |                                          |
        v                                          v
reciprocity_symmetry.py               phase_only_test.py, architecture_test.py (reuses Phase 8B)
        |                                          |
        +------------------+-----------------------+
                            v
                 compatibility_matrix.py (factual table, no ranking) --> final_status()
                            |
                            v
                 results/phase8c_scattering/
```

## Why this reads the internal Phase 7 checkpoint, not the delivered library

The delivered `results/unit_cell_library/database/unit_cell_library.json` only exports the
co-polarized `_xx` component per S-family (`unit_cell_library.py::_cell_doc`). Answering the
brief's Sec. 11 cross-polarization question honestly needs the full linear-basis matrix
(`xx`/`xy`/`yx`/`yy`), which Phase 7 *did* compute and which is still sitting, unmodified, in
`results/unit_cell_library/_solve/database.json` (`database.Database`'s own native format).
`phase8c/data_source.py` reads it read-only and cross-checks a sample of its own `_xx` values
against the delivered library on every run (`tests/test_phase8c_data_source.py`) — proof this is
genuinely the same solve, not stale data.

## The measured result

Run against the real 324-cell library at 12.5 GHz:

| | S11xx | S21xx |
|---|---|---|
| Physical interpretation | reflection coefficient, top plane | transmission coefficient, port1→port2 |
| Phase coverage | 24.84° | 43.85° |
| Largest circular gap | 335.16° | 316.15° |
| Magnitude range | [0.231, 0.831] | [0.556, 0.973] |
| Phase-magnitude correlation | −0.92 | +0.98 |
| Distinct complex states | 169 of 324 | 169 of 324 |
| Architecture test | `REFLECTION_ARCHITECTURE_UNVALIDATED` | `TRANSMISSION_ARCHITECTURE_UNVALIDATED` |
| Validation status | `INSUFFICIENT_EVIDENCE` | `INSUFFICIENT_EVIDENCE` |

**S21xx has *more* phase coverage than S11xx (43.85° vs. 24.84°). This is reported as a
measurement, not used as a selection criterion (Sec. 1's explicit rule) — neither candidate is
"chosen" on that basis, or on any basis, because the physical architecture decision (is the
intended MTS reflective or transmissive?) that would validate one of them remains unmade.**

Two genuinely new, real findings from this phase's numeric checks (not assumptions):

- **Reciprocity holds to ~1e-15** (max|S21xx−S12xx| = 1.14e-15 over all 324 cells) — the exact
  precision `MOM_FORMULATION.md`'s "Z is complex-symmetric at normal incidence" claim predicts,
  now numerically confirmed on real solved data, not just asserted.
- **L↔W mirror symmetry holds to ~4e-15** (S11xx(L,W) = S11yy(W,L) over 171 pairs) — a real,
  testable physical symmetry of the rectangular-patch/square-lattice geometry, confirmed.
- **Cross-polarization is negligible** (mean R_cross, T_cross < 1e-6 of incident power) — checked
  using the full linear-basis matrix, not assumed because a bare `S11xx` "looked" clean.
- **Power balance**: mean R≈0.19, T≈0.81 across the whole library (not just Phase 6's case01),
  mean(R+T)=0.9992, consistent with the substrate's small loss tangent (0.0009) and the confirmed
  absence of any higher-order propagating Floquet mode in 10–15 GHz (onset ≈62.46 GHz).

## Final scientific status: `CONTROL_VARIABLE_UNVALIDATED`

Per the fixed decision order (topology → scattering meaning → reflection/transmission
interpretation → Phase 8B compatibility → complex-response adequacy → library adequacy → final
status, brief Sec. 62): architecture is resolved (partially reflective/transmissive, no ground
plane — a real, measured fact); the fundamental Floquet mode alone propagates; reciprocity holds;
but neither S11xx nor S21xx reaches `VALIDATED` (both `INSUFFICIENT_EVIDENCE`), and Formulation A
remains `LIBRARY_CONTROL_VARIABLE_INCOMPATIBLE` (Phase 8B, unchanged). The result is
`CONTROL_VARIABLE_UNVALIDATED` — not because the library's phase coverage is small (that would be
the *different*, unjustified conclusion Sec. 51 explicitly warns against drawing), but because the
physical interpretation itself remains an open architecture decision this project has not made.

Per Phase 8C's own stop condition (Sec. 63): **do not proceed to geometry inversion.** The next
required step is an explicit architecture decision (reflective vs. transmissive MTS), not a
software workflow transition.

## Running it

```bash
.venv/bin/python examples/run_phase8c_validation.py     # results/phase8c_scattering/
```

## Verification

`tests/conftest.py` (session-scoped fixture for the ~470 MB/~20 s checkpoint load, shared across
all Phase 8C test files); `tests/test_phase8c_data_source.py`, `test_phase8c_topology_audit.py`,
`test_phase8c_power_balance.py`, `test_phase8c_scattering_analysis.py`, `test_phase8c_phase_only_test.py`,
`test_phase8c_architecture_test.py`, `test_phase8c_reciprocity_symmetry.py`,
`test_phase8c_reference_plane_sensitivity.py`, `test_phase8c_compatibility_matrix.py`,
`test_phase8c_run.py` — 48 tests, synthetic edge cases (circular near-wrap, duplicate-state
detection, phase/magnitude separation, non-reciprocal synthetic data) plus direct checks against
the real 324-cell library.
