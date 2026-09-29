"""Phase 11: physics-based heterogeneous MoM / nonlocal coupling study.

No machine learning, no data-fitted surrogate, anywhere in this package (Sec. 2's non-negotiable
rule). Every number reported here comes from a real solve of

    Z I = V

through the exact same, unmodified Phase 1-3 MoM machinery Phases 1-10 already use
(`solver.assemble_system`/`solver.solve`), applied to composite (heterogeneous) geometries via
Phase 9/10's own duck-typed `SupercellUnitCell` -- reused directly, not reimplemented. Where Phase
10's own wrapper (`phase10.grid_solver.solve_grid`) does not expose something Phase 11 needs (a
configurable excitation polarization, Sec. 15), a NEW Phase-11-owned wrapper is added
(`reference_solver.py`) rather than modifying Phase 10's frozen code, following the same
reuse-without-modification discipline every earlier phase in this project has used.

Reuse map:
  ``phase9.supercell_geometry.SupercellUnitCell``/``build_supercell``, ``phase9.dipole_extraction.
    local_dipole_moment``, ``phase9.local_periodicity.NOISE_FLOOR_RELATIVE``/
    ``SUBSTANTIAL_DEVIATION_RELATIVE`` (Phase 11L's convergence tolerances -- declared, sourced,
    reused, never re-invented).
  ``phase10.environment_geometry`` (the validated domain, CENTER_CELLS_MM, weak/moderate/strong
    contrast helpers), ``phase10.grid_solver.solve_grid`` (default x-polarization studies),
    ``phase10.diagonal_study``/``neighbor_sensitivity`` (reused directly for Sec. 9's "reproduce or
    independently re-establish" instruction -- literally the same functions, same physics).
  ``phase8b.object_wave``/``holographic_modulation`` (Sec. 24's three-beam holographic target),
  ``phase8e.beam_reconstruction.far_field_pattern`` (Sec. 23's beam-level far-field integral).
  ``phase9.hfss_reference`` (Sec. 25-26's HFSS-availability gate -- the same, unchanged check).
"""
