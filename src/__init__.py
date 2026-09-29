"""Lightweight spectral-domain periodic MoM solver for metasurface unit cells.

Phase 1: project foundation (configuration, constants, materials, geometry,
lattice/Floquet bookkeeping, excitation).
Phase 2: periodic layered-medium spectral dyadic Green's function
(green_function.py, see GREEN_FUNCTION.md).
Phase 3: rooftop MoM discretization, excitation vector, impedance matrix and the
solve Z I = V (basis_functions.py, mom_matrix.py, solver.py, see MOM_FORMULATION.md).
Phase 4: fields at arbitrary z and S-parameter extraction at the reference planes
(green_function.py, sparameters.py, see S_PARAMETERS.md).
Phase 5: frequency/geometry sweeps, discretization validation, resumable checkpoints and the
unit-cell database with linear and circular polarization data (sweep.py, discretization.py,
runner.py, database.py, polarization.py, phase.py, see DATABASE.md).
Phase 6: STEP export of representative unit cells for manual external validation, built directly
from the same geometry the solver consumes (cad_export.py, validation_cases.py,
validation_export.py, see VALIDATION_EXPORT.md). Purely additive: no EM formulation, Green's
function, MoM matrix, solver or S-parameter code is touched, so this phase does not change any
number the solver produces (see __version__ below) and needs the optional
requirements-cad.txt dependencies only for its own module/tests.
Phase 12A: an optional PEC ground plane at z = -h (backing="pec"; default "air", every prior
config/call site unchanged and verified bit-identical -- see __version__ below), so the solver
can also model a grounded reflectarray/holographic-impedance-surface stack: one-port
S-parameters (no port 2), a one-port database record, and grounded CAD export
(green_function.py, geometry.py, config.py, excitation.py, sparameters.py, database.py,
sweep.py, runner.py, cad_export.py; see PHASE12A_GROUND_PLANE.md and GREEN_FUNCTION.md Sec. 12A).
NOT implemented yet: holographic synthesis, oblique-incidence ports.
"""

# Part of every cell identifier and checkpoint key: bump it whenever a change alters the numbers
# the solver produces, so stale checkpoints are never reused. Phase 6 (STEP export) does not
# change any solver number, so it does not bump this.
__version__ = "0.5.0"
