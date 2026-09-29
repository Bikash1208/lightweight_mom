"""Phase 12D -- resolve the local-periodicity step limit, and ingest the real HFSS data.

Orchestration/reporting only, per this project's own convention (CLAUDE.md "New phases reuse
earlier phases unmodified"). The real capability additions this phase makes live in core modules,
extended additively, never edited in place:

- ``phase9/hfss_reference.py``: a general N-port Touchstone (``.sNp``) reader
  (``read_touchstone``/``TouchstoneNetwork``), for the real HFSS Floquet-port export now on disk
  (``hfss_exports/phase12_grounded_8mm/``).
- ``sparameters.py`` is reused UNMODIFIED: ``scattering_from_fields`` (already accepts an
  arbitrary supplied mean-current array, not only a solved one) is the basis of this phase's
  dipole-to-reflection-coefficient conversion (``reflection_metric.py``, Task A) and its
  Floquet-order power decomposition (``floquet_power.py``, Task B) reuses
  ``green_function.LayeredGreenFunction.sheet_fields``/``stack_fields`` and
  ``basis_functions.RooftopBasis.current_coefficients`` unmodified.
- ``phase12c.grounded_supercell``/``phase12c.locality_study`` (Task A/B) and
  ``phase12c.design_point`` (the unchanged grounded 8 mm design point) are reused, not edited.

Modules here: ``reflection_metric.py`` (Task A's dipole<->S11 conversion), ``step_curve.py``
(Task A's finer-grid step-size study), ``floquet_power.py`` (Task B's per-order power capability),
``aperture_test.py`` (Task B's orchestration), ``hfss_ingest.py`` (Task C), ``gate.py``,
``specification.py``, ``plots.py``, ``report.py``, ``run.py``.
"""
