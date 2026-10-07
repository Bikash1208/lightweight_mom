"""Part 3, app page 6 ("Export"): full-array STEP export, layout CSV/JSON, and an HFSS-ready
README, reusing `cad_export.export_supercell_step`/`phase12c.grounded_supercell.
build_grounded_supercell` (unmodified) for the geometry. DXF export of the patch layer
(PHASE13_WEBAPP_PROMPT.md's "if ezdxf is available") was never built -- `ezdxf` is not a project
dependency and no code here reads or writes DXF; a documented scope gap, not a silent one (see
docs/history/PHASE13_WEBAPP.md)."""
