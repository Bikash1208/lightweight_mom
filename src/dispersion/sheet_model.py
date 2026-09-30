"""Phase 14 Part A: a homogenized patch-sheet-admittance surface-wave dispersion model
(Luukkonen et al., "Simple and Accurate Analytical Model of Planar Grids and High-Impedance
Surfaces Comprising Metal Strips or Patches," IEEE Trans. Antennas Propag., 2008) -- the working
replacement for `mom_scaffold.py`'s direct MoM eigenmode search, which this project's own
numerical experiments found could not resolve a patch-loading-dependent mode at any accessible
mesh size (see that module's own docstring for the finding, kept as a documented research item).

**Why this works where the direct eigenmode search didn't**: the eigenmode search tried to find a
shifted POLE of the periodic array's own Z(kx, ky) directly, which turned out to be swamped by the
bare-slab Green's-function pole at every mesh size tested. This module never searches for a pole
of a periodic-array matrix at all. Instead it (1) extracts the patch LAYER's own shunt sheet
admittance from a quantity this project's MoM already solves robustly and has always solved
correctly -- the NORMAL-INCIDENCE reflection coefficient S11 of the periodic patch array -- then
(2) substitutes that admittance into the closed-form transverse-resonance condition for a TM
surface wave bound to the loaded, grounded structure, a plain scalar equation with no periodic-
array matrix in it at all, and solves it by bisection. Validated against `analytic_tm0.tm0_beta`
for a small patch (beta should approach the bare-slab value as the patch shrinks) and found to give
a real, monotonically increasing beta(L) -- see `tests/test_dispersion_sheet_model.py`.

**The physics** (this project's own exp(+j*omega*t), Im(kz) <= 0 conventions throughout):

1. `sparameters.compute_s_parameters` at normal incidence (kx0 = ky0 = 0, the function's default),
   with `ReferencePlanes(top_air_thickness=0.0, ...)` so S11 is already referenced exactly at the
   patch plane z=0 -- no manual de-embedding needed (verified: for a lightly-loaded cell this
   reproduces the bare-slab-only reflection coefficient computed from `_y_slab` below to 4+ sig
   figs).
2. `Y_in = Y0*(1 - Gamma)/(1 + Gamma)`, `Gamma = S11_xx` (`Y0 = 1/eta0`, `Y_in` is the total
   admittance of the loaded stack as seen from air at z=0).
   `Y_slab = _y_slab(kt=0, ...)`: the admittance of the grounded dielectric ALONE (no patch),
   looking down from z=0 -- a shorted-transmission-line closed form for the TM mode,
   `-j*(omega*eps0*eps_r/kz)*cot(kz*h)`, `kz = sqrt(eps_r*k0**2 - kt**2)`. This is NOT the same
   quantity as `green_function.LayeredGreenFunction.equivalent_impedances`, which is already the
   PARALLEL combination of the up- and down-looking admittances (its own docstring says so) --
   using it here would double-count the air side. Verified independently: `_y_air(kt) + _y_slab
   (kt)` reproduces `1/Z_TM` from `equivalent_impedances` to ~1e-17 (floating-point noise) at every
   kt tested, confirming both closed forms and the branch convention agree with this project's
   already-tested Green's function.
   `Y_sheet = Y_in - Y_slab` (the patch layer's own contribution, assuming the sheet and the slab
   admittance combine in parallel at z=0 -- the standard thin-sheet equivalent-circuit picture,
   and assuming `Y_sheet` itself is (nearly) independent of incidence angle -- valid for a period
   much smaller than a wavelength, the same regime this project's holograms already operate in,
   but an approximation, not exact: `labels.DISPERSION_FROM_SHEET_MODEL` says so).
3. `beta` solves `Im[Y_sheet + _y_air(beta) + _y_slab(beta)] = 0` for real beta in
   `(k0, Re(k_substrate))` (`sparameters.surface_wave_interval`'s own interval; at `Y_sheet = 0`
   this is exactly the bare-slab TM0 transverse resonance). The real part of the same expression at
   the found beta is reported as a residual (small for a low-loss board; not forced to zero, since
   a real S11 can carry a small real part this simplified model does not separately track).
   `X_s = eta0*sqrt((beta/k0)**2 - 1)` (the standard TM surface-reactance closed form).
"""

import numpy as np
from scipy.optimize import brentq

from ..constants import EPS0, ETA0
from ..errors import LwmomError, config_invalid, dispersion_sheet_model_no_root
from ..excitation import ReferencePlanes
from ..frequency import angular_frequency, free_space_wavenumber
from ..labels import DISPERSION_FROM_SHEET_MODEL, EXPERIMENTAL
from ..sparameters import compute_s_parameters, surface_wave_interval
from ..green_function import LayeredGreenFunction
from .hfss_table import DispersionSample, DispersionTable
from .mom_scaffold import square_patch_unit_cell

SCAFFOLD_LABELS = (DISPERSION_FROM_SHEET_MODEL, EXPERIMENTAL)

Y0 = 1.0 / ETA0

# Above this beta*period/pi, higher Floquet harmonics/grating effects are close enough to onset
# that the angle-independent, isolated-sheet assumption this model relies on is increasingly
# suspect -- flagged, not blocked (the brief's own "flag beta*P approaching pi").
STOP_BAND_WARNING_THRESHOLD = 0.9


def _y_slab(beta, omega, eps_c, k0, thickness_m):
    kz = np.sqrt(eps_c * k0**2 - beta**2)
    return -1j * (omega * EPS0 * eps_c / kz) / np.tan(kz * thickness_m)


def _y_air(beta, omega, k0):
    alpha = np.sqrt(beta**2 - k0**2 + 0j)
    return 1j * omega * EPS0 / alpha


class SheetModelMode:
    """One found (beta, X_s) point, plus diagnostics. A plain class (not frozen dataclass) kept
    minimal on purpose -- this is an internal return type, not a public schema."""

    def __init__(self, beta_rad_per_m, x_s_ohm, y_sheet, residual, beta_p_over_pi, stop_band_warning):
        self.beta_rad_per_m = beta_rad_per_m
        self.x_s_ohm = x_s_ohm
        self.y_sheet = y_sheet
        self.residual = residual
        self.beta_p_over_pi = beta_p_over_pi
        self.stop_band_warning = stop_band_warning
        self.labels = SCAFFOLD_LABELS


def sheet_admittance_from_normal_incidence(unit_cell, frequency_hz, mesh):
    """`Y_sheet` (complex, [S]) for one patch geometry -- step 1-2 above. `unit_cell` must be
    grounded (`backing == "pec"`; the transverse-resonance equation below assumes a shorted
    substrate line)."""
    if getattr(unit_cell, "backing", "air") != "pec":
        raise config_invalid(
            message="the sheet-admittance model needs a grounded cell (backing='pec'); "
            f"got backing={getattr(unit_cell, 'backing', 'air')!r}"
        )
    planes = ReferencePlanes(
        top_air_thickness=0.0, substrate_thickness=unit_cell.substrate.thickness, bottom_air_thickness=None
    )
    sp = compute_s_parameters(unit_cell, frequency_hz, mesh, planes, warn_grating=False, warn_surface_wave=False)
    gamma = sp.s11[0, 0]
    y_in = Y0 * (1.0 - gamma) / (1.0 + gamma)

    omega = angular_frequency(frequency_hz)
    k0 = free_space_wavenumber(frequency_hz)
    eps_c = unit_cell.substrate.eps_complex
    y_slab0 = _y_slab(0.0, omega, eps_c, k0, unit_cell.substrate.thickness)
    return y_in - y_slab0


def find_beta_from_sheet_admittance(y_sheet, substrate, frequency_hz, period_m=None, n_scan=2000):
    """Step 3: solves the loaded transverse-resonance condition for beta. Raises `LwmomError`
    (`DISPERSION_SHEET_MODEL_NO_ROOT`) if the scan finds no sign change -- never a fabricated
    number. `period_m`, if given, is used only to flag (not block) a mode approaching the
    Brillouin-zone edge (`STOP_BAND_WARNING_THRESHOLD`)."""
    omega = angular_frequency(frequency_hz)
    k0 = free_space_wavenumber(frequency_hz)
    eps_c = substrate.eps_complex
    green = LayeredGreenFunction(substrate, frequency_hz, backing="pec")
    lo, hi = surface_wave_interval(green)

    def f_imag(beta):
        return (y_sheet + _y_air(beta, omega, k0) + _y_slab(beta, omega, eps_c, k0, substrate.thickness)).imag

    betas = np.linspace(lo * (1.0 + 1e-4), hi * (1.0 - 1e-4), n_scan)
    vals = np.array([f_imag(b) for b in betas])
    sign_changes = np.where(np.diff(np.sign(vals)) != 0)[0]
    if len(sign_changes) == 0:
        raise dispersion_sheet_model_no_root(
            f"no sign change in the transverse-resonance equation's imaginary part over "
            f"({betas[0]:.6g}, {betas[-1]:.6g}) rad/m"
        )
    i = int(sign_changes[0])
    beta = float(brentq(f_imag, betas[i], betas[i + 1]))
    residual = float(
        (y_sheet + _y_air(beta, omega, k0) + _y_slab(beta, omega, eps_c, k0, substrate.thickness)).real
    )
    x_s = float(ETA0 * np.sqrt((beta / k0) ** 2 - 1.0))

    beta_p_over_pi = None
    stop_band_warning = False
    if period_m is not None:
        beta_p_over_pi = float(beta * period_m / np.pi)
        stop_band_warning = beta_p_over_pi >= STOP_BAND_WARNING_THRESHOLD

    return SheetModelMode(
        beta_rad_per_m=beta,
        x_s_ohm=x_s,
        y_sheet=y_sheet,
        residual=residual,
        beta_p_over_pi=beta_p_over_pi,
        stop_band_warning=stop_band_warning,
    )


def find_beta_for_cell(unit_cell, frequency_hz, mesh, period_m=None):
    """Steps 1-3 combined: `Y_sheet` from `unit_cell`'s own normal-incidence S11, then `beta`/`X_s`
    from the transverse-resonance equation. `period_m` defaults to `unit_cell.lattice.period_x`
    (square lattices only need one)."""
    y_sheet = sheet_admittance_from_normal_incidence(unit_cell, frequency_hz, mesh)
    if period_m is None:
        period_m = unit_cell.lattice.period_x
    return find_beta_from_sheet_admittance(y_sheet, unit_cell.substrate, frequency_hz, period_m=period_m)


def build_dispersion_table_from_sheet_model(unit_cell_at, geometry_params, frequency_hz, mesh, period_m=None, source_label="sheet_model_dispersion_search"):
    """Sweeps `geometry_params` (e.g. a square patch side length, in metres) through
    `find_beta_for_cell`, building a `DispersionTable` the same shape `hfss_table.
    read_dispersion_csv` would. A sample that raises `DISPERSION_SHEET_MODEL_NO_ROOT` is skipped,
    not fabricated; raises the same error itself if every sample is skipped."""
    samples = []
    for gp in geometry_params:
        cell = unit_cell_at(gp)
        try:
            mode = find_beta_for_cell(cell, frequency_hz, mesh, period_m=period_m)
        except LwmomError:
            continue
        samples.append(DispersionSample(geometry_param=float(gp), beta_rad_per_m=mode.beta_rad_per_m, x_s_ohm=mode.x_s_ohm))
    if not samples:
        raise dispersion_sheet_model_no_root(
            f"no bound-mode root found for any of {len(list(geometry_params))} geometry samples"
        )
    samples.sort(key=lambda s: s.geometry_param)
    return DispersionTable(
        source_path=source_label,
        frequency_hz=float(frequency_hz),
        samples=tuple(samples),
        label=DISPERSION_FROM_SHEET_MODEL,
    )
