"""Phase 8C: unit-cell scattering architecture and control-quantity validation.

Read-only with respect to Phases 1-7: no solver, mesh classifier, Green's function, or Phase 7
data is modified. Where Phase 8A/8B's delivered `unit_cell_library.json` only exports the
co-polarized "_xx" component per S-family, this package reads the full linear-basis 4x4 S-matrix
(all of xx/xy/yx/yy) from the still-present internal Phase 7 solve checkpoint
(`results/unit_cell_library/_solve/database.json`, `database.Database`'s own native format,
unmodified, produced by the same `runner.run_sweep` Phase 7 already ran) -- real, already-computed
data, never re-solved, never fabricated. See `PHASE8C_SCATTERING.md` for the narrative summary.
"""
