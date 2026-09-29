"""Task B: the aperture-level test -- periodic phase-gradient supercells (anomalous reflectors),
solved as a FULL supercell (not just the centre-cell dipole Phase 12C's own periodic-supercell case
used), decomposed into every propagating Floquet order of the SUPERCELL's own lattice
(``floquet_power.py``), and compared against the local-periodicity (lookup-table) prediction. This
decides Phase 13 (brief "Why" point 2), not Phase 9's dipole-moment proxy.

Grid construction reuses ``phase12c.locality_study.build_periodic_phase_gradient_supercell``
UNMODIFIED (same N-cells-per-2pi-cycle construction, cells picked from a real, freshly-solved
isolated-cell curve -- "the lookup table"). Floquet order for the full-supercell solve uses the
SAME reduced (rank-sufficient-plus-margin, not the 2x-total machine-precision) convention Phase
12C's own Task D used for its large, qualitative-effect spot-checks
(``total_cells // 2 + margin``) -- justified here for the identical reason: Task B's comparisons
(fraction of power in the intended order, spurious-order levels) are large, qualitative effects
where a small residual floquet-truncation uncertainty (independently bounded below, this module's
own convergence check) cannot flip the conclusion, and the full-supercell solve at N=6-12 cells is
already at this project's practical cost ceiling even at this reduced order (measured, Sec. below).
"""

import time
from dataclasses import dataclass

import numpy as np

from .. import solver as solver_mod
from ..basis_functions import MeshSettings
from ..excitation import PlaneWaveExcitation
from ..green_function import LayeredGreenFunction
from ..phase import wrapped_difference_deg
from ..phase9.dipole_extraction import local_dipole_moment
from ..phase12c import design_point as dp
from ..phase12c.grounded_supercell import build_grounded_supercell
from ..phase12c.locality_study import CurvePoint, build_periodic_phase_gradient_supercell
from .floquet_power import floquet_order_powers, intended_order_for_gradient
from .reflection_metric import equivalent_reflection_from_dipole_moment

APERTURE_FLOQUET_MARGIN = 8   # matches phase12c.spot_check's own Task-D convention


def aperture_floquet_order(total_cells):
    return total_cells // 2 + APERTURE_FLOQUET_MARGIN


def build_lookup_curve(mesh_cells_per_subcell, l_values_mm):
    """A fresh, real isolated-cell S11_xx curve at THIS module's own mesh -- "the isolated-cell
    lookup table" the brief asks cells to be picked from."""
    from ..basis_functions import RooftopBasis
    from ..sparameters import compute_s_parameters

    floquet = 2 * mesh_cells_per_subcell   # isolated single-cell solves are cheap; use the fully
                                            # validated convention (Phase 12C), not the reduced one.
    points = []
    for l_mm in l_values_mm:
        cell = dp.unit_cell_at(l_mm * 1e-3)
        mesh = MeshSettings(mesh_cells_per_subcell, mesh_cells_per_subcell, floquet)
        s = compute_s_parameters(cell, dp.DESIGN_FREQUENCY_HZ, mesh, dp.REFERENCE_PLANES,
                                  warn_grating=False, warn_surface_wave=False)
        n_unknowns = RooftopBasis(cell, mesh_cells_per_subcell, mesh_cells_per_subcell).n_unknowns
        points.append(CurvePoint(l_mm, s.magnitude("S11_xx"), s.phase_deg("S11_xx"), n_unknowns))
    return tuple(points)


@dataclass(frozen=True)
class PerCellPhaseError:
    j: int
    l_mm: float
    predicted_phase_deg: float     # lookup-table (isolated-cell) prediction
    actual_phase_deg: float        # MoM full-supercell solve, per-cell equivalent reflection
    phase_error_deg: float


@dataclass(frozen=True)
class ApertureTestResult:
    n_cells: int
    mesh_cells_per_subcell: int
    floquet_order: int
    grid_mm: tuple
    achieved_span_deg: float
    is_full_2pi: bool
    n_unknowns: int
    runtime_s: float
    intended_order: tuple
    fraction_intended: float
    fraction_specular: float
    total_propagating_fraction: float
    spurious_orders: tuple         # ((m, n), fraction), every OTHER propagating order, sorted desc
    intended_order_phase_deg: float
    predicted_intended_phase_deg: "float | None"   # generalized-Snell's-law prediction, if defined
    per_cell_phase_errors: tuple   # PerCellPhaseError, one per column


def run_aperture_test(n_cells, mesh_cells_per_subcell, l_values_mm):
    """Solves the FULL n_cells x n_cells periodic phase-gradient supercell (Sec. above), decomposes
    its reflected power into every Floquet order of ITS OWN lattice, and reports the design metric.
    Raises only on a genuine solver failure; a too-expensive size is the CALLER's job to have
    already ruled out by timing (this function always attempts the solve; ``run.py`` decides which
    sizes to call it with and reports infeasibility for the rest, per the brief's own instruction)."""
    t0 = time.perf_counter()
    curve = build_lookup_curve(mesh_cells_per_subcell, l_values_mm)
    grid_mm, achieved_span_deg, is_full_2pi, targets = build_periodic_phase_gradient_supercell(
        curve, n_cells, mesh_cells_per_subcell,
    )
    grid_m = [[(l * 1e-3, l * 1e-3) for (l, _) in row] for row in grid_mm]
    supercell = build_grounded_supercell(n_cells, dp.PERIOD_M, dp.PERIOD_M, dp.SUBSTRATE, dp.METAL, grid_m)

    total_cells = n_cells * mesh_cells_per_subcell
    floquet = aperture_floquet_order(total_cells)
    excitation = PlaneWaveExcitation(frequency=dp.DESIGN_FREQUENCY_HZ, polarization="x")
    mesh = MeshSettings(total_cells, total_cells, floquet)
    system = solver_mod.assemble_system(supercell, excitation, mesh)
    currents = solver_mod.solve(system.z, system.v)

    green = LayeredGreenFunction(dp.SUBSTRATE, dp.DESIGN_FREQUENCY_HZ, backing="pec")
    max_order = min(n_cells * 2, floquet)   # cover a few multiples of n_cells' own grating orders
    power = floquet_order_powers(system.basis, currents, green, supercell.lattice, max_order)

    intended = intended_order_for_gradient(n_cells, gradient_axis="x")
    fraction_intended = power.fraction_at(*intended)
    fraction_specular = power.fraction_at(0, 0)
    intended_phase = power.phase_deg_at(*intended)

    all_orders = [(int(m), int(n)) for m, n in zip(power.m.ravel(), power.n.ravel())]
    fracs = {(int(m), int(n)): float(f) for m, n, f in zip(power.m.ravel(), power.n.ravel(), power.fraction.ravel())}
    prop = {(int(m), int(n)): bool(p) for m, n, p in zip(power.m.ravel(), power.n.ravel(), power.propagating.ravel())}
    spurious = sorted(
        ((mn, f) for mn, f in fracs.items() if prop[mn] and mn not in ((0, 0), intended)),
        key=lambda kv: -kv[1],
    )

    # Per-cell (per-column) phase-error profile: MoM's actual local response vs the lookup-table
    # (isolated-cell) prediction, one dipole extraction per distinct column (Phase 9's own
    # local_dipole_moment, reused; matches Phase 12C's own periodic-supercell convention).
    per_cell = []
    for j in range(n_cells):
        l_mm = grid_mm[0][j][0]
        cx = -dp.PERIOD_M * n_cells / 2 + (j + 0.5) * dp.PERIOD_M
        dipole = local_dipole_moment(system.basis, currents, cx, 0.0, dp.PERIOD_M / 2, dp.PERIOD_M / 2)
        actual = equivalent_reflection_from_dipole_moment(
            dipole.px, dp.PERIOD_M, dp.SUBSTRATE, dp.DESIGN_FREQUENCY_HZ, dp.REFERENCE_PLANES,
        )
        predicted = next(c for c in curve if c.l_mm == l_mm)
        err = abs(wrapped_difference_deg(actual.phase_deg, predicted.s11_phase_deg))
        per_cell.append(PerCellPhaseError(j, l_mm, predicted.s11_phase_deg, actual.phase_deg, err))

    runtime = time.perf_counter() - t0
    return ApertureTestResult(
        n_cells=n_cells, mesh_cells_per_subcell=mesh_cells_per_subcell, floquet_order=floquet,
        grid_mm=tuple(tuple(row) for row in grid_mm), achieved_span_deg=achieved_span_deg,
        is_full_2pi=is_full_2pi, n_unknowns=system.basis.n_unknowns, runtime_s=runtime,
        intended_order=intended, fraction_intended=fraction_intended, fraction_specular=fraction_specular,
        total_propagating_fraction=power.total_propagating_fraction, spurious_orders=tuple(spurious),
        intended_order_phase_deg=intended_phase, predicted_intended_phase_deg=None,
        per_cell_phase_errors=tuple(per_cell),
    )
