"""Part 3D: a MoM surface-wave-eigenmode capability, real -- this module used to be a scaffold
only (see git history); PHASE13D_POLISH_AND_PHYSICS.md asks for the real thing since it's "the
most valuable physics item left."

**The physics**: a surface-wave eigenmode of the periodic grounded unit cell is a solution of this
project's own Z(kx, ky) I = 0 (source-free) system at a REAL (kx, ky) with kt = sqrt(kx^2+ky^2) > k0
(bound, non-radiating -- `CONVENTIONS.md`'s `Im(kz) <= 0` branch choice already fixes the right
analytic continuation for kt > k0, since `Im(kz)` there is already negative/decaying by the
existing branch rule; no new sign convention needed). At such a point `Z` is (numerically) singular.
`Z` is assembled from `green.floquet_dyads`, whose (m, n) = (0, 0) term is the layered medium's own
spectral Green's dyad `LayeredGreenFunction.spectral_dyad` evaluated *exactly* at (kx0, ky0) -- and
that dyad has a REAL POLE at the bare grounded slab's own TM0 propagation constant (a textbook
transverse-resonance result; `dispersion.analytic_tm0` is its closed form). That pole enters `Z`
additively through the (0, 0) Floquet term, so **the largest singular value of Z(kt) blows up right
at the surface-wave pole** -- confirmed numerically (see `tests/test_dispersion_mom.py`): scanning
`svd(Z(kt))[0]` (not the smallest singular value / a raw `det`, both far noisier for a matrix this
size) across `sparameters.surface_wave_interval` shows a sharp, well-separated peak whose location
matches `analytic_tm0.tm0_beta` to <0.001 rad/m for a lightly-loaded cell, and shifts away from the
bare value as the patch loads the surface harder (a real, physical loading effect, not fitting
noise). This is the found beta; `find_surface_wave_mode` reports it, plus an equivalent surface
reactance derived from it (below), plus a `prominence` (peak / scan-median) so a flat, ambiguous
scan is refused rather than silently trusted -- "no silent fallbacks" (`errors.py`).

**Surface reactance**: `DispersionTable` (from `hfss_table.py`, reused here) needs both
`beta_rad_per_m` and `x_s_ohm` per geometry sample (`hologram/hmsa_synthesis.py`'s Formulation A
looks both up). There is no separate "reactance solve" here -- X_s is *defined* from the found beta
via the same transverse-resonance circuit the bare-slab TM0 pole itself satisfies: the loaded cell's
pole condition is `Y_TM,bare(beta) + Y_s = 0` (a shunt sheet admittance `Y_s` added at the patch
plane, in parallel with the bare slab's own TM admittance, transverse-resonating at the LOADED
beta rather than the bare one) -- so `Y_s = -Y_TM,bare(beta_loaded)`, i.e. `X_s = -1/Im(Y_s)`, using
`green.equivalent_impedances` (no patch) at the *found, loaded* beta. This is the standard
equivalent-sheet-impedance picture (Sievenpiper 1999; Fong et al. 2010) applied self-consistently
to this project's own full-wave pole, not a separate approximation.

**Validation, per the original scaffold's requirement 3**: `tests/test_dispersion_mom.py` checks a
lightly-loaded (small-patch) grounded cell's found beta against `analytic_tm0.tm0_beta` before
trusting this on a real patterned cell -- required before this label may read anything but
`EXPERIMENTAL`; it keeps that label regardless (a young numerical method, one board/frequency
combination's worth of validation), per `labels.EXPERIMENTAL`'s own definition.
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
        raise dispersion_mom_no_clear_mode(
            f"scan maximum at a boundary (kt={kts[imax]:.6g}) -- no interior peak found"
        )
    prominence = float(vals[imax]) / background
    if prominence < min_prominence:
        raise dispersion_mom_no_clear_mode(
            f"peak/median prominence {prominence:.2f}x is below the {min_prominence:.2f}x threshold "
            f"-- this scan has no clear resonance, not a fabricated one"
        )
    return kts[imax - 1], kts[imax + 1], prominence


def find_surface_wave_mode(unit_cell, frequency_hz, direction_deg=0.0, mesh=DEFAULT_MESH, n_scan=DEFAULT_N_SCAN, min_prominence=DEFAULT_MIN_PROMINENCE):
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


def build_dispersion_table_from_mom(unit_cell_at, geometry_params, frequency_hz, direction_deg=0.0, mesh=DEFAULT_MESH, source_label="mom_dispersion_search"):
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
        samples.append(DispersionSample(geometry_param=float(gp), beta_rad_per_m=mode.beta_rad_per_m, x_s_ohm=mode.x_s_ohm))
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
