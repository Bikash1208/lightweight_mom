"""Phase 12A: optional PEC ground plane (grounded-substrate stack).

This package is the *orchestration/reporting* layer only. The actual physics (the grounded
closed forms, one-port S-parameters, database/CAD support) lives in the core modules it reuses
unmodified: ``src/green_function.py``, ``src/geometry.py``, ``src/config.py``,
``src/excitation.py``, ``src/sparameters.py``, ``src/database.py``, ``src/sweep.py``,
``src/runner.py``, ``src/discretization.py``, ``src/cad_export.py`` (each extended, per
PHASE12A_GROUND_PLANE_PROMPT.md's own instruction, with an optional, backward-compatible
``backing``/``pec`` code path rather than rewritten).

Scope (the brief's own words): "This phase adds a capability and measures what it gives; it does
not build a production library, does not re-run Phases 8-11, and does not unlock geometry
inversion." Geometry inversion and STEP export remain LOCKED unconditionally (``gate.py``).

No ML anywhere in this package (there is nothing here that would need it: every number is a
direct MoM solve or a closed-form/analytic check against it).
"""
