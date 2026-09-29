"""Phase 13 Part 2: solver diagnostics captured on every app-facing solve. Wraps
`sparameters.compute_s_parameters` (unmodified) -- measures what's cheaply observable from OUTSIDE
it (timing, mesh/Floquet bookkeeping already public on `RooftopBasis`/`MeshSettings`, and
post-hoc checks on the returned `SParameters`) rather than editing the solver core to add internal
instrumentation, per CLAUDE.md's "physics/numerics untouched" rule.
"""

import resource
import time
from dataclasses import dataclass

import numpy as np

from .basis_functions import RooftopBasis
from .discretization import geometry_edges_aligned
from .errors import solve_failed
from .green_function import LayeredGreenFunction
from .logging_setup import get_logger
from .sparameters import compute_s_parameters, surface_wave_proximity

_log = get_logger("diagnostics")

# `numpy.linalg.cond` on the dense Z matrix is O(n^3) -- "where cheap" (brief's own wording): only
# attempted below this many unknowns, else reported as skipped rather than silently stalling a
# real-time app request.
CONDITION_ESTIMATE_MAX_UNKNOWNS = 300


@dataclass(frozen=True)
class SolveDiagnostics:
    run_id: str | None
    cells_x: int
    cells_y: int
    dx_m: float
    dy_m: float
    mesh_aligned: bool
    floquet_order: int
    n_unknowns: int
    total_time_s: float
    peak_memory_kb: float | None
    condition_estimate: float | None
    condition_estimate_skipped_reason: str | None
    grating_orders_propagating: bool
    surface_wave_margin: dict | None       # `sparameters.surface_wave_proximity`'s own dict, or None (ungrounded)
    energy_conservation_residual: float | None   # |1 - power_out_per_incident| for x-pol, None if not computable
    reciprocity_residual: float | None            # |S21_xx - S12_xx|, None for a grounded (one-port) result


def _grating_orders_propagate(unit_cell, green, floquet_order):
    """Whether any Floquet harmonic other than (0,0) is propagating (a plain, public
    re-derivation of the same physics `sparameters._warn_if_grating_orders` checks, exposed here
    as a diagnostic value rather than only as a `UserWarning` side effect)."""
    m, n = unit_cell.lattice.mode_grid(floquet_order)
    kx, ky = unit_cell.lattice.transverse_wavenumbers(0.0, 0.0, m, n)
    kz = green.normal_wavenumbers(kx, ky)[0]
    propagating = kz.imag == 0.0
    i0 = int(floquet_order)
    return bool(np.sum(propagating) > 1) if propagating[i0, i0] else bool(np.any(propagating))


def instrumented_compute_s_parameters(unit_cell, frequency, mesh, planes, run_id=None, **kwargs):
    """Solves via `compute_s_parameters` (unmodified) and returns `(SParameters, SolveDiagnostics)`.
    Raises `errors.LwmomError` (code `SOLVE_FAILED`) instead of letting a low-level exception
    propagate unlabeled -- the one place this module doesn't just observe but also translates."""
    basis = RooftopBasis(unit_cell, mesh.cells_x, mesh.cells_y)
    period_x, period_y = unit_cell.lattice.period_x, unit_cell.lattice.period_y
    geometry = getattr(unit_cell, "geometry", None)
    aligned = True
    if geometry is not None and geometry.kind == "rectangle":
        aligned = geometry_edges_aligned(geometry, unit_cell.lattice, mesh)

    backing = getattr(unit_cell, "backing", "air")
    green = LayeredGreenFunction(unit_cell.substrate, frequency, backing=backing)

    mem_before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    t0 = time.perf_counter()
    try:
        result = compute_s_parameters(unit_cell, frequency, mesh, planes, **kwargs)
    except Exception as exc:  # noqa: BLE001 -- deliberately broad: this IS the translation boundary
        _log.error("solve failed: %s", exc, extra={"cells_x": mesh.cells_x, "cells_y": mesh.cells_y})
        raise solve_failed(exc) from exc
    total_time = time.perf_counter() - t0
    mem_after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_memory_kb = float(max(mem_after, mem_before))  # ru_maxrss is already a running peak (POSIX)

    condition_estimate, skip_reason = None, None
    if basis.n_unknowns <= CONDITION_ESTIMATE_MAX_UNKNOWNS:
        from .mom_matrix import assemble_impedance_matrix
        z = assemble_impedance_matrix(basis, green, 0.0, 0.0, mesh.floquet_order)
        condition_estimate = float(np.linalg.cond(z))
    else:
        skip_reason = f"n_unknowns={basis.n_unknowns} > {CONDITION_ESTIMATE_MAX_UNKNOWNS} (too expensive)"

    sw_margin = None
    if backing == "pec":
        sw_margin = surface_wave_proximity(unit_cell, green, mesh.floquet_order)

    energy_residual = None
    reciprocity_residual = None
    try:
        power = result.power_out_per_incident()
        energy_residual = float(abs(1.0 - power[0, 0]))
    except Exception:
        pass
    if not result.grounded:
        try:
            # Reciprocity (CONVENTIONS.md: reciprocal => S = S^T): S21_xx must equal S12_xx.
            reciprocity_residual = float(abs(result.component("S21_xx") - result.component("S12_xx")))
        except Exception:
            pass

    diagnostics = SolveDiagnostics(
        run_id=run_id, cells_x=mesh.cells_x, cells_y=mesh.cells_y,
        dx_m=period_x / mesh.cells_x, dy_m=period_y / mesh.cells_y, mesh_aligned=aligned,
        floquet_order=mesh.floquet_order, n_unknowns=basis.n_unknowns, total_time_s=total_time,
        peak_memory_kb=peak_memory_kb, condition_estimate=condition_estimate,
        condition_estimate_skipped_reason=skip_reason,
        grating_orders_propagating=_grating_orders_propagate(unit_cell, green, mesh.floquet_order),
        surface_wave_margin=sw_margin, energy_conservation_residual=energy_residual,
        reciprocity_residual=reciprocity_residual,
    )
    _log.debug("solve complete", extra={"diagnostics": diagnostics.__dict__})
    return result, diagnostics
