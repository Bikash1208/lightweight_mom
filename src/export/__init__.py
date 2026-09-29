"""Part 3, app page 6 ("Export"): full-array STEP export, layout CSV/JSON, and an HFSS-ready
README, reusing `cad_export.export_supercell_step`/`phase12c.grounded_supercell.
build_grounded_supercell` (unmodified) for the geometry and `manifest.RunManifest` for provenance.
DXF export uses the optional `ezdxf` package if present (`errors.cad_unavailable` otherwise, same
pattern as the `cad` extra)."""
