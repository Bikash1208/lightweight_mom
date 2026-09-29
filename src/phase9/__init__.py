"""Phase 9: physical control-variable validation (9A) and local-periodicity validation (9B).

Both workstreams are validation, not geometry generation (Sec. 1, 36): no final MTS STEP is
produced, and geometry inversion stays locked (`gate.py`) unless both gates open.

## Phase 9A (`control_variable.py`, `hfss_reference.py`, `error_metrics.py`)

Compares the existing, unmodified solver's `S11_xx`/`S21_xx` (read from the real Phase 7
checkpoint, never re-solved) against HFSS reference data -- **which does not exist anywhere in
this repository** (grep-verified, same finding as every earlier phase). Per Rule 1 ("if
information is unavailable, explicitly report UNRESOLVED / NOT AVAILABLE FROM CURRENT DATA. Do
not silently assume"), the comparison infrastructure is fully built and tested against synthetic,
explicitly-labelled test data, but the production run honestly reports every HFSS-comparison
metric as unavailable and both `S11`/`S21` classifications as `UNVALIDATED` for that reason --
never fabricated.

## Phase 9B (`supercell_geometry.py`, `dipole_extraction.py`, `supercell_solver.py`,
`local_periodicity.py`)

A genuinely new, real capability: builds a heterogeneous N x N "supercell" as one larger-period
`UnitCell`-duck-typed object (multiple sub-patches, each its own L/W) and solves it through the
**exact same, unmodified** `solver.assemble_system`/`solve` MoM pipeline Phases 1-7 already use --
no coupled-cell solver is invented, no new EM formulation, no new Green's function. The center
sub-cell's induced-current dipole moment (`RooftopBasis.evaluate`, unmodified, reused directly) is
compared between the isolated (periodic single-cell) and embedded (heterogeneous supercell) cases
-- a real, computed local-periodicity error, not asserted.

## `gate.py`

Sec. 29-30's two independent gates and the geometry-inversion lock; Sec. 42's four-value final
status taxonomy.
"""
