"""Phase 12B Task D: grounded rectangular-PEC-patch period study (4.8 / 8.0 / 12.0 mm), same
substrate as Phase 12A, aligned meshes throughout.
"""

import math
from dataclasses import dataclass

from ..basis_functions import MeshSettings
from ..constants import C0
from ..discretization import (
    grating_onset_frequency,
    grating_onset_frequency_oblique,
    max_abs_difference,
    smallest_aligned_cells,
)
from ..excitation import ReferencePlanes
from ..geometry import Geometry, UnitCell
from ..materials import Metal, Substrate
from ..periodic_bc import RectangularLattice
from ..phase_extraction import phase_coverage
from ..sparameters import compute_s_parameters, surface_wave_proximity

PERIODS_MM = (4.8, 8.0, 12.0)
SUBSTRATE_EPS_R = 2.2
SUBSTRATE_TAN_DELTA = 0.0009
SUBSTRATE_THICKNESS_M = 1.575e-3
FREQUENCIES_HZ = (10.0e9, 11.25e9, 12.5e9, 13.75e9, 15.0e9)
DESIGN_FREQUENCY_HZ = 12.5e9
FEED_ANGLES_DEG = (0.0, 15.0, 30.0, 45.0)   # a documented set (brief: "report, don't assume a feed angle")

# One data point per period: (L step, lower/upper fraction of P), chosen (not guessed) to keep
# the aligned mesh and the point count computationally tractable while covering ~0.2P to ~0.94P
# (PHASE12B_MESH_PERIOD_STUDY.md records the measured per-period aligned mesh size and cost).
_COARSE_STEP_MM = {4.8: 0.2, 8.0: 0.4, 12.0: 0.4}
_LO_FRAC, _HI_FRAC = 0.2, 0.94


def _substrate():
    return Substrate("RT5880", SUBSTRATE_EPS_R, SUBSTRATE_TAN_DELTA, SUBSTRATE_THICKNESS_M)


def l_family_for_period(period_mm, step_mm=None, lo_frac=_LO_FRAC, hi_frac=_HI_FRAC):
    """Coarse L=W family for one period, centred at 0, step chosen per-period (module docstring)."""
    step_mm = step_mm if step_mm is not None else _COARSE_STEP_MM[period_mm]
    lo = round(period_mm * lo_frac / step_mm) * step_mm
    hi = round(period_mm * hi_frac / step_mm) * step_mm
    n = round((hi - lo) / step_mm)
    return tuple(round(lo + i * step_mm, 6) for i in range(n + 1))


def aligned_cells_for_period(period_mm, l_values_mm=None, max_cells=4000):
    l_values_mm = l_values_mm if l_values_mm is not None else l_family_for_period(period_mm)
    features = [(l * 1e-3, 0.0) for l in l_values_mm]
    n = smallest_aligned_cells(period_mm * 1e-3, features, max_cells=max_cells)
    if n is None:
        raise ValueError(f"no aligned mesh found for period={period_mm}mm within {max_cells} cells")
    return n


def _unit_cell(period_mm, length_width_mm):
    lattice = RectangularLattice(period_mm * 1e-3, period_mm * 1e-3)
    geometry = Geometry("rectangle", length_width_mm * 1e-3, length_width_mm * 1e-3)
    return UnitCell(lattice, geometry, _substrate(), Metal(), backing="pec")


def _planes():
    return ReferencePlanes(5.0e-3, None, SUBSTRATE_THICKNESS_M)


@dataclass(frozen=True)
class PeriodSweepPoint:
    period_mm: float
    length_width_mm: float
    frequency_hz: float
    s11_xx_mag: float
    s11_xx_phase_deg: float


def run_period_sweep(period_mm, l_values_mm=None, frequencies_hz=FREQUENCIES_HZ, cells=None):
    """Direct MoM solves (grounded), on the smallest aligned mesh for the requested L family (or
    ``cells``, if the caller already knows it -- must still be aligned, verified). Returns
    (points, cells_used)."""
    l_values_mm = l_values_mm if l_values_mm is not None else l_family_for_period(period_mm)
    n = cells if cells is not None else aligned_cells_for_period(period_mm, l_values_mm)
    features = [(l * 1e-3, 0.0) for l in l_values_mm]
    if smallest_aligned_cells(period_mm * 1e-3, features, max_cells=n) != n:
        raise ValueError(f"cells={n} does not align every requested L for period={period_mm}mm")
    planes = _planes()
    points = []
    for l_mm in l_values_mm:
        cell = _unit_cell(period_mm, l_mm)
        for f in frequencies_hz:
            mesh = MeshSettings(n, n, n)
            s = compute_s_parameters(cell, f, mesh, planes, warn_grating=False, warn_surface_wave=False)
            points.append(PeriodSweepPoint(period_mm, l_mm, f, s.magnitude("S11_xx"), s.phase_deg("S11_xx")))
    return tuple(points), n


def coverage_by_frequency(points):
    freqs = sorted({p.frequency_hz for p in points})
    out = {}
    for f in freqs:
        table = [{"phase_deg": p.s11_xx_phase_deg} for p in points if p.frequency_hz == f]
        out[f] = phase_coverage(table)
    return out


def resonant_length_mm(points, frequency_hz=DESIGN_FREQUENCY_HZ):
    """The L (of those actually solved, at ``frequency_hz``) with the smallest |S11_xx| -- the
    measured resonance location, not an assumed one; ties broken by the first occurrence."""
    at_f = [p for p in points if p.frequency_hz == frequency_hz]
    if not at_f:
        raise ValueError(f"no points at frequency_hz={frequency_hz!r}")
    return min(at_f, key=lambda p: p.s11_xx_mag).length_width_mm


@dataclass(frozen=True)
class GratingSurfaceWaveReport:
    period_mm: float
    grating_onset_normal_hz: float
    grating_onset_oblique_hz: dict          # feed_angle_deg -> onset Hz
    specular_only_across_band: bool         # at normal incidence, whole FREQUENCIES_HZ band
    surface_wave_margins: dict              # frequency_hz -> nearest_margin_over_k0 (min over any L)


def grating_and_surface_wave_report(
    period_mm, representative_length_width_mm, frequencies_hz=FREQUENCIES_HZ, feed_angles_deg=FEED_ANGLES_DEG,
):
    cell = _unit_cell(period_mm, representative_length_width_mm)
    onset_normal = grating_onset_frequency(cell)
    onset_oblique = {deg: grating_onset_frequency_oblique(cell, math.radians(deg)) for deg in feed_angles_deg}
    specular_ok = all(f < onset_normal for f in frequencies_hz)
    margins = {}
    n = aligned_cells_for_period(period_mm, (representative_length_width_mm,))
    for f in frequencies_hz:
        from ..green_function import LayeredGreenFunction

        green = LayeredGreenFunction(cell.substrate, f, backing="pec")
        r = surface_wave_proximity(cell, green, floquet_order=n)
        margins[f] = r["nearest_margin_over_k0"]
    return GratingSurfaceWaveReport(period_mm, onset_normal, onset_oblique, specular_ok, margins)


@dataclass(frozen=True)
class PeriodConvergencePoint:
    cells: int
    s11_xx_mag: float
    s11_xx_phase_deg: float


def convergence_for_near_resonant_patch(period_mm, length_width_mm, frequency_hz=DESIGN_FREQUENCY_HZ, factors=(1, 2)):
    """Aligned-mesh convergence (floquet_order = cells, Phase 7's own convention) for one patch,
    at successively finer *aligned* meshes (the base aligned N, then 2x, 4x -- doubling an aligned
    mesh's cell count keeps every original edge aligned, since the new mesh lines are a superset)."""
    base = aligned_cells_for_period(period_mm, (length_width_mm,))
    cell = _unit_cell(period_mm, length_width_mm)
    planes = _planes()
    points, s11s = [], []
    for factor in factors:
        n = base * factor
        mesh = MeshSettings(n, n, n)
        s = compute_s_parameters(cell, frequency_hz, mesh, planes, warn_grating=False, warn_surface_wave=False)
        points.append(PeriodConvergencePoint(n, s.magnitude("S11_xx"), s.phase_deg("S11_xx")))
        s11s.append(s.s11)
    diffs = [max_abs_difference(s11s[i + 1], s11s[i]) for i in range(len(s11s) - 1)]
    return tuple(points), tuple(diffs)
