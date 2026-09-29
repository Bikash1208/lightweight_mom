"""Part 3: import of the HFSS field-calculator export format used for HMSA reference fields
(`fld_parser.py`), on top of which `reference_field_grid.py` builds a
`phase8e.reference_field_input.ReferenceFieldGrid` (reused unmodified) and does the
validation/sanity-check work PHASE13_WEBAPP_PROMPT.md's "Reference-field import" section asks for
(grid-regularity, local-wavenumber estimate vs the analytic TM0 baseline, cells-per-fringe
warning)."""
