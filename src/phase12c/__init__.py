"""Phase 12C: does local periodicity hold for the grounded 8 mm reflectarray cell?

Orchestration/report package only (mirrors phase9/10/11/12a/12b): every new *capability* this
phase needs lives in a core module, extended additively, never rewritten --
``src/phase12c/grounded_supercell.py`` wraps (does not edit) ``phase9.supercell_geometry`` to add
``backing="pec"``; ``src/cad_export.py`` gains supercell STEP export; ``src/phase9/hfss_reference.py``
gains one-port (no S21) CSV ingestion. See PHASE12C_GROUNDED_LOCALITY.md.
"""
