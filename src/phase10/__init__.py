"""Phase 10: environment-aware electromagnetic response model.

Entry scientific status (from Phase 9): `LOCAL_PERIODICITY_INVALIDATED` -- the isolated
periodic-cell model can be wrong by an amount comparable to the signal's own magnitude under
realistic strong-contrast neighbor conditions. Phase 10 asks: how much neighborhood information
does the electromagnetic system actually require, and can a finite, validated environment-aware
model replace the isolated-cell lookup for holographic synthesis?

Reuse map (nothing in Phase 1-9 is modified):
  ``phase9.supercell_geometry.SupercellUnitCell``/``build_supercell`` -- the composite,
    duck-typed multi-patch unit cell, solved through the exact same unmodified
    ``solver.assemble_system``/``solver.solve`` pipeline Phases 1-3 already use.
  ``phase9.dipole_extraction.local_dipole_moment`` -- the center/local-cell response readout,
    reusing ``basis_functions.RooftopBasis.evaluate`` (Phase 2, unmodified).
  ``phase9.error_metrics`` -- complex/magnitude/wrapped-phase error primitives.
  ``phase9.supercell_solver.solve_isolated_cell`` -- the isolated-periodic-cell reference.
  ``phase8e.beam_reconstruction.far_field_pattern`` -- the standard planar-aperture radiation
    integral, reused unmodified for Sec. 28-29's beam-level validation.
  ``phase.wrapped_difference_deg``/``principal_phase_deg``, ``frequency.free_space_wavenumber``.

The "center-cell electromagnetic response" S_c used throughout this package is the same physical
observable Phase 9B already validated: the local induced-current dipole moment (a complex number,
units A.m) of one sub-cell's footprint, extracted from a solved heterogeneous-supercell MoM
system. This is NOT S11 or S21 (Phase 9A's control-variable question remains unresolved, per
Rule/Sec. 44 of that phase, for lack of any HFSS reference) -- it is the same well-defined,
standard EM quantity (the zeroth-order/DC spatial moment of the induced surface current) that
Phase 9B's own correctness proof (a homogeneous 3x3 supercell reproduces the isolated cell to
~1e-14) already established as numerically sound. Sec. 1's generic "S" is instantiated concretely
as this quantity throughout Phase 10, and that choice is stated explicitly everywhere it matters
so no reader mistakes it for a validated S-parameter.
"""
