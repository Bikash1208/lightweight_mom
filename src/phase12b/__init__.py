"""Phase 12B: mesh alignment, duplicate-geometry audit, grounded period study.

Orchestration/reporting only. The actual new *capability* (mesh-alignment helpers,
duplicate-geometry detection) lives in ``src/discretization.py`` and ``src/sweep.py``
(``SweepPlan.mesh_duplicate_groups``/``n_mesh_distinct``), extended additively -- every
pre-Phase-12B sweep/config/result is unaffected and verified unchanged. No ML anywhere here
either: every number is a direct MoM solve, an analytic/closed-form check, or a read of
already-solved data (the Phase 7 checkpoint, Phase 12A's own results).
"""
