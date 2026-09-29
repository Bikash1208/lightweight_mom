"""Phase 12B Task B: re-assess Phase 12A's grounded physics-sanity sweep and refinement study on
mesh-aligned discretizations, so staircasing (a discretization artefact) is no longer entangled
with genuine under-resolution.
"""

from dataclasses import dataclass

from .. import library_schema
from ..basis_functions import MeshSettings
from ..discretization import feature_edges_aligned, max_abs_difference, smallest_aligned_cells
from ..phase12a.sweep_study import (
    DESIGN_FREQUENCY_HZ,
    SWEEP_FREQUENCIES_HZ,
    _planes,
    grounded_library_cell,
)
from ..phase_extraction import phase_coverage
from ..sparameters import compute_s_parameters

PERIOD_M = library_schema.LibraryConfig.default().period_m
COARSE_L_VALUES_MM = library_schema.L_VALUES_MM  # 1.0-4.4mm, 0.2mm step
FINE_L_VALUES_MM = tuple(round(3.6 + 0.1 * i, 1) for i in range(11))  # 3.6-4.6mm, 0.1mm step
ALL_L_VALUES_MM = tuple(sorted(set(COARSE_L_VALUES_MM) | set(FINE_L_VALUES_MM)))

ALIGNED_STUDY_CELLS = (24, 48, 96)  # Sec. "Redo the grounded refinement study"
REPRESENTATIVE_LW_MM = 2.8  # C3, this project's own already-used center geometry

# Evidence-based choice (Phase 12B measurement, PHASE12B_MESH_PERIOD_STUDY.md): at the aligned
# 96-cell mesh, floquet_order=48 agrees with floquet_order=96 to ~0.06deg / ~2e-4 relative on
# S11_xx for the representative patch -- i.e. already well converged -- while costing ~3x less
# than floquet_order=96 (12.8s vs 43s measured per solve). Used for the full aligned sweep below;
# the dedicated floquet-convergence study (function below) verifies this choice on its own.
ALIGNED_SWEEP_CELLS = 96
ALIGNED_SWEEP_FLOQUET_ORDER = 48


def aligned_mesh_for_family(l_values_mm=ALL_L_VALUES_MM, period_m=PERIOD_M):
    """The smallest aligned cells_x=cells_y (both axes, since every patch here is square) for the
    given L=W family, or raises if none is found within discretization's own search bound (an
    honest failure, not a silent fallback to an unaligned mesh)."""
    features = [(l * 1e-3, 0.0) for l in l_values_mm]
    n = smallest_aligned_cells(period_m, features)
    if n is None:
        raise ValueError(f"no aligned mesh found for {l_values_mm!r} within the search bound")
    return n


@dataclass(frozen=True)
class AlignedRefinementPoint:
    cells: int
    floquet_order: int
    s11_xx_mag: float
    s11_xx_phase_deg: float


def aligned_mesh_refinement_study(
    length_width_m=REPRESENTATIVE_LW_MM * 1e-3,
    frequency_hz=DESIGN_FREQUENCY_HZ,
    cells_list=ALIGNED_STUDY_CELLS,
):
    """Sec. 'Redo the grounded refinement study for L=2.8mm on aligned meshes (24/48/96 cells)':
    solve at each *aligned* cells_x=cells_y with floquet_order = cells (Phase 7's own '~1x cells'
    convention, library_schema.py), and report successive |S11| differences -- the joint
    mesh+Floquet analogue of discretization.refinement_study, but every mesh here is verified
    aligned first (raises if not -- never silently solve on an unaligned mesh here)."""
    for n in cells_list:
        if not feature_edges_aligned(PERIOD_M, n, length_width_m, 0.0):
            raise ValueError(f"cells={n} is not aligned for length_width_m={length_width_m!r}")
    cell = grounded_library_cell(length_width_m)
    planes = _planes()
    points = []
    s11s = []
    for n in cells_list:
        mesh = MeshSettings(n, n, n)
        s = compute_s_parameters(cell, frequency_hz, mesh, planes, warn_grating=False, warn_surface_wave=False)
        points.append(AlignedRefinementPoint(n, n, s.magnitude("S11_xx"), s.phase_deg("S11_xx")))
        s11s.append(s.s11)
    diffs = [max_abs_difference(s11s[i + 1], s11s[i]) for i in range(len(s11s) - 1)]
    return tuple(points), tuple(diffs)


@dataclass(frozen=True)
class FloquetConvergencePoint:
    floquet_order: int
    s11_xx_mag: float
    s11_xx_phase_deg: float


def floquet_only_convergence(
    length_width_m=REPRESENTATIVE_LW_MM * 1e-3,
    frequency_hz=DESIGN_FREQUENCY_HZ,
    cells=ALIGNED_SWEEP_CELLS,
    floquet_orders=(48, 64, 96),
):
    """Sec. 'also check Floquet convergence separately': mesh fixed (aligned) at ``cells``, only
    ``floquet_order`` varies -- isolates Floquet-sum truncation from mesh/staircasing error."""
    cell = grounded_library_cell(length_width_m)
    planes = _planes()
    points, s11s = [], []
    for order in floquet_orders:
        mesh = MeshSettings(cells, cells, order)
        s = compute_s_parameters(cell, frequency_hz, mesh, planes, warn_grating=False, warn_surface_wave=False)
        points.append(FloquetConvergencePoint(order, s.magnitude("S11_xx"), s.phase_deg("S11_xx")))
        s11s.append(s.s11)
    diffs = [max_abs_difference(s11s[i + 1], s11s[i]) for i in range(len(s11s) - 1)]
    return tuple(points), tuple(diffs)


@dataclass(frozen=True)
class AlignedSweepPoint:
    length_width_mm: float
    frequency_hz: float
    s11_xx_mag: float
    s11_xx_phase_deg: float
    s11_yx_mag: float


def run_aligned_grounded_sweep(
    l_values_mm=ALL_L_VALUES_MM,
    frequencies_hz=SWEEP_FREQUENCIES_HZ,
    cells=ALIGNED_SWEEP_CELLS,
    floquet_order=ALIGNED_SWEEP_FLOQUET_ORDER,
):
    """Sec. 'Re-run the Phase 12A sweep... on the aligned mesh, plus finer L steps near
    resonance': the same substrate/period/backing as Phase 12A, on the smallest mesh that aligns
    the *whole* combined coarse+fine L family (verified, not assumed -- ``aligned_mesh_for_family``
    would raise otherwise)."""
    if aligned_mesh_for_family(l_values_mm) > cells:
        raise ValueError(
            f"cells={cells} does not align the requested L family; "
            f"aligned_mesh_for_family reports {aligned_mesh_for_family(l_values_mm)}"
        )
    planes = _planes()
    points = []
    for l_mm in l_values_mm:
        cell = grounded_library_cell(l_mm * 1e-3)
        for f in frequencies_hz:
            mesh = MeshSettings(cells, cells, floquet_order)
            s = compute_s_parameters(cell, f, mesh, planes, warn_grating=False, warn_surface_wave=False)
            points.append(
                AlignedSweepPoint(
                    l_mm,
                    f,
                    s.magnitude("S11_xx"),
                    s.phase_deg("S11_xx"),
                    s.magnitude("S11_yx"),
                )
            )
    return tuple(points)


def coverage_by_frequency(points):
    freqs = sorted({p.frequency_hz for p in points})
    out = {}
    for f in freqs:
        table = [{"phase_deg": p.s11_xx_phase_deg} for p in points if p.frequency_hz == f]
        out[f] = phase_coverage(table)
    return out
