"""Part 3D / Phase 14 Part A (research item, kept for the record -- NOT the production dispersion
source; see `sheet_model.py` for that): a direct MoM surface-wave-eigenmode search over the
periodic array's own Z(kx, ky), and why it does not deliver a patch-loading-dependent dispersion
curve at any mesh size this project can afford to run interactively.

**The physics attempted**: a surface-wave eigenmode of the periodic grounded unit cell is (in
principle) a solution of this project's own Z(kx, ky) I = 0 (source-free) system at a REAL
(kx, ky) with kt = sqrt(kx^2+ky^2) > k0 (bound, non-radiating -- `CONVENTIONS.md`'s `Im(kz) <= 0`
branch choice already fixes the right analytic continuation for kt > k0, no new sign convention
needed). `find_surface_wave_mode` below scans the largest singular value of the assembled Z(kt)
across `sparameters.surface_wave_interval` and refines the peak -- **and that peak is real and
well-validated, but it is always the BARE-SLAB Green's-function pole, not a loaded-array mode**.

**Why, and how this was found (Phase 14's own investigation)**: `Z` is assembled from
`green.floquet_dyads`, whose (m, n) = (0, 0) term is the layered medium's own spectral Green's dyad
`LayeredGreenFunction.spectral_dyad` evaluated *exactly* at (kx0, ky0) -- and that dyad has a REAL
POLE at the bare grounded slab's own TM0 propagation constant (`dispersion.analytic_tm0`'s closed
form), *independent of any patch pattern*, since the dyad is a property of the layered medium
alone. That pole enters `Z` additively, so `svd(Z(kt))[0]` spikes there regardless of loading. Phase
14 tested this explicitly at patch fill fractions 0.05-0.92 and two mesh resolutions (8x8/floquet 5
and, to rule out a resolution artifact, 16x16/floquet 10, n_unknowns=264): **the peak's location is
identical to `analytic_tm0.tm0_beta` to numerical precision at every loading tested, and the
smallest singular value of Z(kt) -- where a genuine loaded-array eigenvalue, det(Z)=0, would show
up -- stays completely flat across the whole Brillouin zone, at every loading, at both mesh sizes**.
There is no numerically resolvable secondary feature to find at accessible mesh sizes; the original
claim in this docstring (removed) that the peak "shifts away from the bare value as the patch loads
the surface harder" was an over-read of one, as it turned out non-representative, weak-loading data
point, corrected here. Locating a periodic loaded structure's own eigenmode immediately next to a
Green's-function pole is a known hard problem in the numerical-electromagnetics literature (usually
addressed with complex-plane contour-integral eigensolvers, e.g. Beyn's algorithm, not a real-axis
scalar scan) -- genuinely out of scope for this project to solve from scratch.

**What's still real and kept**: `find_surface_wave_mode` correctly and robustly finds the bare-slab
TM0 pole (validated to <0.001 rad/m against `analytic_tm0.tm0_beta`, `tests/test_dispersion_mom.py`)
-- useful as a MoM-based cross-check of the analytic bare-slab value on a real, discretized cell,
and as the honest starting point for a future, more capable eigensolver. It is NOT used to build a
per-geometry dispersion table for Formulation A; `sheet_model.py`'s homogenized-sheet-admittance
model does that instead, using a completely different (and validated) technique that never needs to
resolve a pole-adjacent eigenvalue at all. `square_patch_unit_cell` below is still shared with that
module (a plain geometry helper, unrelated to the failed part).
"""

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize_scalar

from ..basis_functions import MeshSettings, RooftopBasis
from ..errors import LwmomError, config_invalid, dispersion_mom_no_clear_mode
from ..geometry import Geometry, UnitCell
from ..green_function import LayeredGreenFunction
from ..labels import DISPERSION_FROM_MOM, EXPERIMENTAL
from ..mom_matrix import assemble_impedance_matrix
from ..sparameters import surface_wave_interval
from .hfss_table import DispersionSample, DispersionTable

SCAFFOLD_LABELS = (DISPERSION_FROM_MOM, EXPERIMENTAL)

DEFAULT_MESH = MeshSettings(cells_x=8, cells_y=8, floquet_order=5)
DEFAULT_N_SCAN = 500
DEFAULT_MIN_PROMINENCE = 3.0  # peak / scan-median; below this the "peak" is scan noise, not a pole


@dataclass(frozen=True)
class SurfaceWaveMode:
    """One found MoM surface-wave eigenmode."""

    beta_rad_per_m: float
    x_s_ohm: "float | None"
    frequency_hz: float
    direction_deg: float
    prominence: float
    labels: tuple = SCAFFOLD_LABELS


def _peak_singular_value(basis, green, kt, direction_deg, floquet_order):
    kx0 = kt * math.cos(math.radians(direction_deg))
    ky0 = kt * math.sin(math.radians(direction_deg))
    z = assemble_impedance_matrix(basis, green, kx0, ky0, floquet_order)
    return float(np.linalg.svd(z, compute_uv=False)[0])


def _locate_peak(kts, vals, min_prominence):
    """Scan-then-refine peak search shared by the real search and its unit tests (a synthetic
    ``vals`` array can exercise the "no clear mode" path deterministically, without relying on a
    numerically fragile physical configuration to happen to be ambiguous)."""
    kts = np.asarray(kts, dtype=float)
    vals = np.asarray(vals, dtype=float)
    imax = int(np.argmax(vals))
    background = float(np.median(vals))
    if imax in (0, len(kts) - 1) or background <= 0:
        raise dispersion_mom_no_clear_mode(f"scan maximum at a boundary (kt={kts[imax]:.6g}) -- no interior peak found")
    prominence = float(vals[imax]) / background
    if prominence < min_prominence:
        raise dispersion_mom_no_clear_mode(
            f"peak/median prominence {prominence:.2f}x is below the {min_prominence:.2f}x threshold "
            f"-- this scan has no clear resonance, not a fabricated one"
        )
    return kts[imax - 1], kts[imax + 1], prominence


def find_surface_wave_mode(
    unit_cell,
    frequency_hz,
    direction_deg=0.0,
    mesh=DEFAULT_MESH,
    n_scan=DEFAULT_N_SCAN,
    min_prominence=DEFAULT_MIN_PROMINENCE,
):
    """Search `sparameters.surface_wave_interval` for the real, bound TM surface-wave eigenmode of
    `unit_cell` (must have `backing == "pec"`) at `frequency_hz`, along `direction_deg` (0 = +x).

    Raises `LwmomError` (`DISPERSION_MOM_NO_CLEAR_MODE`) if the scan shows no clear resonance --
    never a fabricated number. Returns a `SurfaceWaveMode`.
    """
    if getattr(unit_cell, "backing", "air") != "pec":
        raise config_invalid(
            message="A surface-wave eigenmode search needs a grounded cell (backing='pec'); "
            f"got backing={getattr(unit_cell, 'backing', 'air')!r}"
        )
    green = LayeredGreenFunction(unit_cell.substrate, frequency_hz, backing="pec")
    lo, hi = surface_wave_interval(green)
    basis = RooftopBasis(unit_cell, mesh.cells_x, mesh.cells_y)

    kts = np.linspace(lo * (1.0 + 2e-4), hi * (1.0 - 2e-3), n_scan)
    vals = np.array([_peak_singular_value(basis, green, kt, direction_deg, mesh.floquet_order) for kt in kts])
    bracket_lo, bracket_hi, prominence = _locate_peak(kts, vals, min_prominence)

    res = minimize_scalar(
        lambda kt: -_peak_singular_value(basis, green, kt, direction_deg, mesh.floquet_order),
        bounds=(bracket_lo, bracket_hi),
        method="bounded",
        options={"xatol": 1e-4},
    )
    beta = float(res.x)

    cos_d, sin_d = math.cos(math.radians(direction_deg)), math.sin(math.radians(direction_deg))
    _z_te_bare, z_tm_bare = green.equivalent_impedances(np.array([beta * cos_d]), np.array([beta * sin_d]))
    y_bare = 1.0 / z_tm_bare[0]
    x_s = float(-1.0 / y_bare.imag) if abs(y_bare.imag) > 1e-30 else None

    return SurfaceWaveMode(
        beta_rad_per_m=beta,
        x_s_ohm=x_s,
        frequency_hz=float(frequency_hz),
        direction_deg=float(direction_deg),
        prominence=prominence,
    )


def build_dispersion_table_from_mom(
    unit_cell_at,
    geometry_params,
    frequency_hz,
    direction_deg=0.0,
    mesh=DEFAULT_MESH,
    source_label="mom_dispersion_search",
):
    """Sweep `geometry_params` (e.g. a square patch side length, in metres) through
    `find_surface_wave_mode`, building a `DispersionTable` the same way `hfss_table.read_dispersion_csv`
    would -- `unit_cell_at(geometry_param)` builds the grounded `UnitCell` for one sample (the
    caller owns lattice/substrate/backing choices; this function only drives the sweep and the
    eigenmode search). A sample that raises `DISPERSION_MOM_NO_CLEAR_MODE` is skipped, not fabricated;
    raises `LwmomError` (`DISPERSION_MOM_NO_CLEAR_MODE`) itself if every sample is skipped.
    """
    samples = []
    for gp in geometry_params:
        cell = unit_cell_at(gp)
        try:
            mode = find_surface_wave_mode(cell, frequency_hz, direction_deg=direction_deg, mesh=mesh)
        except LwmomError:
            continue
        samples.append(
            DispersionSample(geometry_param=float(gp), beta_rad_per_m=mode.beta_rad_per_m, x_s_ohm=mode.x_s_ohm)
        )
    if not samples:
        raise dispersion_mom_no_clear_mode(
            f"no clear surface-wave mode found for any of {len(list(geometry_params))} geometry samples"
        )
    samples.sort(key=lambda s: s.geometry_param)
    return DispersionTable(
        source_path=source_label,
        frequency_hz=float(frequency_hz),
        samples=tuple(samples),
        label=DISPERSION_FROM_MOM,
    )


def square_patch_unit_cell(period_m, substrate, metal, patch_side_m):
    """Convenience `unit_cell_at` for `build_dispersion_table_from_mom`: a centred square patch,
    grounded, on a square lattice -- the shape the app's auto-library already uses
    (`library_sweep.build_square_library_cells`)."""
    from ..periodic_bc import RectangularLattice

    lattice = RectangularLattice(period_m, period_m)
    geometry = Geometry("rectangle", patch_side_m, patch_side_m)
    return UnitCell(lattice, geometry, substrate, metal, backing="pec")
