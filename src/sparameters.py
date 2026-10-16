"""S-parameter extraction: solved current -> fields at the reference planes -> S11/S21 (and S12/S22).

Full definitions, equations and conventions: S_PARAMETERS.md. Summary (time factor exp(+j w t)):

* Two ports, both in air and both carrying the same transverse axes x, y:
  port 1 = the plane z = z_top (above the pattern), port 2 = the plane z = z_bottom
  (below the substrate); both come from ``ReferencePlanes`` (explicit in config.yaml).
* Only the specular (0,0) plane-wave order is a port mode; normal incidence.
* Wave amplitudes are tangential-E amplitudes referenced to the plane:

      port 1:  incident  a1 exp(+jk0 (z - z_top)) [travels -z],  reflected  b1 exp(-jk0 (z - z_top))
      port 2:  incident  a2 exp(-jk0 (z - z_bot)) [travels +z],  reflected  b2 exp(+jk0 (z - z_bot))

  (b1, b2) = S (a1, a2), each S_ij a 2x2 matrix [received pol, incident pol] in the basis (x, y).
* The amplitudes are *measured*: the total tangential E and H of the (0,0) harmonic at each
  plane (unpatterned-stack field + sheet field) are split into up- and down-going parts,

      E_up = (E - eta0 z x H)/2,   E_down = (E + eta0 z x H)/2,

  so the incident amplitude used in S = b/a is what the fields actually contain.
"""

import math
import warnings
from dataclasses import dataclass, replace

import numpy as np

from .basis_functions import MeshSettings, RooftopBasis
from .constants import C0, EPS0, ETA0, MU0
from .excitation import PlaneWaveExcitation, ReferencePlanes, stack_tangential_field
from .geometry import UnitCell
from .green_function import LayeredGreenFunction
from .mom_matrix import assemble_excitation_vector, assemble_impedance_matrix
from .solver import solve

PORTS = ("top", "bottom")  # port 1, port 2
POLARIZATIONS = "xy"  # index 0, 1
_COMPONENT_INDEX = {"xx": (0, 0), "yx": (1, 0), "xy": (0, 1), "yy": (1, 1)}  # (received, incident)
_FAMILIES = ("S11", "S12", "S21", "S22")


def decompose_air_fields(e_t, h_t):
    """Split tangential (E, H) of a (0,0) wave in air into (E_up, E_down).

    Up-going (+z) waves have H = z x E / eta0, down-going (-z) waves H = -z x E / eta0, hence
    E_up = (E - eta0 z x H)/2 and E_down = (E + eta0 z x H)/2 with z x H = (-H_y, H_x).
    """
    e_t = np.asarray(e_t, dtype=complex)
    h_t = np.asarray(h_t, dtype=complex)
    zxh = np.stack([-h_t[..., 1], h_t[..., 0]], axis=-1)
    return 0.5 * (e_t - ETA0 * zxh), 0.5 * (e_t + ETA0 * zxh)


@dataclass(frozen=True, eq=False)
class SParameters:
    """Complex S-parameters of one unit cell at one frequency, 2 polarizations, and (for the
    default ungrounded stack) 2 ports.

    Each of ``s11, s12, s21, s22`` is a (2, 2) complex array indexed [received, incident] with
    index 0 = x, 1 = y, so e.g. ``s21[1, 0]`` is the y-polarized field transmitted to port 2
    for an x-polarized wave incident on port 1 (cross-polarized). Diagonal = co-polarized.

    ``grounded=True`` (Phase 12A, a grounded stack, ``backing="pec"``): there is only one port
    (top). ``s12``, ``s21`` and ``s22`` are then ``None`` -- **undefined, not zero** (S21 is
    physically zero for a lossless grounded stack, but representing "undefined" as an actual
    zero would make that a silent assumption rather than a measured result) -- and every method
    that would read them (``matrix``, ``component`` for a family other than S11, ...) raises
    ``ValueError`` instead of silently returning 0. ``incident_amplitudes`` then has shape
    (1, 2) and ``mean_currents`` shape (1, 2, 2) (port 1 only).
    """

    frequency: float
    planes: ReferencePlanes
    grounded: bool
    s11: np.ndarray
    s12: np.ndarray | None
    s21: np.ndarray | None
    s22: np.ndarray | None
    incident_amplitudes: np.ndarray  # (n_ports, 2 pols): measured a at the reference plane
    incoming_residual: float  # largest incoming-wave amplitude that must vanish, / |a|
    mean_currents: np.ndarray  # (n_ports, 2 pols, 2): J~_00 [A/m] (zeros: no metal)

    def __post_init__(self):
        undefined = (self.s12, self.s21, self.s22)
        if self.grounded:
            if any(v is not None for v in undefined):
                raise ValueError("a grounded (one-port) SParameters must have s12=s21=s22=None")
            if self.incident_amplitudes.shape[0] != 1 or self.mean_currents.shape[0] != 1:
                raise ValueError("a grounded SParameters must have exactly 1 port of data")
        else:
            if any(v is None for v in undefined):
                raise ValueError("an ungrounded (two-port) SParameters needs s12, s21 and s22")

    # ------------------------------------------------------------------ views
    @property
    def k0(self):
        return 2.0 * np.pi * self.frequency / C0

    def matrix(self):
        """4x4 matrix [[S11, S12], [S21, S22]] over (port1 x, port1 y, port2 x, port2 y).

        Raises ``ValueError`` for a grounded (one-port) result: there is no port 2, so no 4x4
        matrix exists -- use ``.s11`` directly.
        """
        if self.grounded:
            raise ValueError("no 4x4 S-matrix exists for a grounded (one-port) stack (no port 2); use .s11")
        return np.block([[self.s11, self.s12], [self.s21, self.s22]])

    def component(self, name):
        """One complex entry, e.g. ``component("S21_yx")`` = y received at port 2, x incident at port 1.

        Raises ``ValueError`` for a family other than S11 of a grounded result (undefined, not 0).
        """
        family, _, comp = name.partition("_")
        if family not in _FAMILIES or comp not in _COMPONENT_INDEX:
            raise ValueError(f"unknown component {name!r}; use <{'|'.join(_FAMILIES)}>_<{'|'.join(_COMPONENT_INDEX)}>")
        block = getattr(self, family.lower())
        if block is None:
            raise ValueError(
                f"{family} is undefined for a grounded (one-port) stack (no port 2 exists); only S11 is defined"
            )
        return complex(block[_COMPONENT_INDEX[comp]])

    def magnitude(self, name):
        return abs(self.component(name))

    def magnitude_db(self, name):
        return 20.0 * np.log10(max(abs(self.component(name)), 1e-300))

    def phase_deg(self, name):
        """Phase in degrees, in (-180, 180] (output boundary; internally radians)."""
        return float(np.degrees(np.angle(self.component(name))))

    def to_record(self):
        """Flat dict (SI frequency; magnitude linear, phase in degrees) for a unit-cell database.

        For a grounded (one-port) result only the S11 entries are included (S12/S21/S22 are
        undefined -- Sec. "S12/S21/S22 undefined" above).
        """
        rec = {"frequency_hz": self.frequency}
        families = ("S11",) if self.grounded else _FAMILIES
        for fam in families:
            for comp in _COMPONENT_INDEX:
                key = f"{fam}_{comp}"
                rec[f"{key}_mag"] = self.magnitude(key)
                rec[f"{key}_phase_deg"] = self.phase_deg(key)
        return rec

    def power_out_per_incident(self):
        """(n_ports, 2) array [port, pol]: sum over all output ports/pols of |S|^2 for that
        incident wave. 1 for a lossless cell whose only propagating order is (0,0); < 1 with
        loss. For a grounded (one-port) result this is the reflected power fraction (the *only*
        power balance that exists, since there is no port 2): (1, 2), from S11 alone."""
        if self.grounded:
            return np.sum(np.abs(self.s11) ** 2, axis=0).reshape(1, 2)
        m = self.matrix()
        return np.sum(np.abs(m) ** 2, axis=0).reshape(2, 2)

    # --------------------------------------------------------- reference planes
    def at_planes(self, new_planes):
        """Re-reference to other planes by pure (0,0) propagation in air (exact, lossless air).

        With D1 = z_top' - z_top and D2 = z_bottom' - z_bottom (moves outward when D1 > 0, D2 < 0):

            S11' = S11 exp(-2j k0 D1),   S22' = S22 exp(+2j k0 D2),
            S21' = S12' = S21 exp(+j k0 (D2 - D1)).

        For a grounded (one-port) result ``new_planes`` must also be grounded (no bottom plane to
        move to or from); only the S11/D1 part of the above applies.
        """
        if not np.isclose(new_planes.substrate_thickness, self.planes.substrate_thickness, rtol=1e-12):
            raise ValueError("substrate thickness differs between plane sets")
        if self.grounded != new_planes.grounded:
            raise ValueError("cannot re-reference between a grounded (one-port) and an ungrounded (two-port) plane set")
        d1 = new_planes.z_top - self.planes.z_top
        k0 = self.k0
        if self.grounded:
            return replace(
                self,
                planes=new_planes,
                s11=self.s11 * np.exp(-2j * k0 * d1),
                incident_amplitudes=self.incident_amplitudes * np.exp(1j * k0 * d1),
            )
        d2 = new_planes.z_bottom - self.planes.z_bottom
        return replace(
            self,
            planes=new_planes,
            s11=self.s11 * np.exp(-2j * k0 * d1),
            s22=self.s22 * np.exp(+2j * k0 * d2),
            s21=self.s21 * np.exp(+1j * k0 * (d2 - d1)),
            s12=self.s12 * np.exp(+1j * k0 * (d2 - d1)),
            incident_amplitudes=self.incident_amplitudes
            * np.array([[np.exp(1j * k0 * d1)] * 2, [np.exp(-1j * k0 * d2)] * 2]),
        )


def _check_planes(green, planes):
    if not np.isclose(planes.substrate_thickness, green.substrate.thickness, rtol=1e-12, atol=0.0):
        raise ValueError(
            f"reference planes assume substrate thickness {planes.substrate_thickness}, "
            f"stack has {green.substrate.thickness}"
        )
    grounded = green.backing == "pec"
    if grounded != planes.grounded:
        raise ValueError(
            f"stack backing={green.backing!r} but reference planes "
            f"{'have no bottom port' if planes.grounded else 'expect a bottom port'} "
            f"(planes.bottom_air_thickness={planes.bottom_air_thickness!r}); "
            f"a grounded stack needs planes with bottom_air_thickness=None and vice versa"
        )


def _scattering_from_fields_two_port(green, planes, mean_currents):
    z_planes = np.array([planes.z_top, planes.z_bottom])
    s_blocks = {name: np.zeros((2, 2), dtype=complex) for name in _FAMILIES}
    amplitudes = np.zeros((2, 2), dtype=complex)
    residual = 0.0

    for port_index, port in enumerate(PORTS):
        for pol in range(2):
            incident = np.eye(2)[pol]
            e_stack, h_stack = green.stack_fields(0.0, 0.0, incident, z_planes, port)
            e_sheet, h_sheet = green.sheet_fields(0.0, 0.0, mean_currents[port_index, pol], z_planes)
            up, down = decompose_air_fields(e_stack + e_sheet, h_stack + h_sheet)  # (plane, comp)

            if port == "top":  # incident wave = down-going at plane 1; nothing incoming at plane 2
                a = down[0, pol]
                incoming_should_vanish = np.concatenate([[down[0, 1 - pol]], up[1]])
            else:  # incident wave = up-going at plane 2; nothing incoming at plane 1
                a = up[1, pol]
                incoming_should_vanish = np.concatenate([[up[1, 1 - pol]], down[0]])
            amplitudes[port_index, pol] = a
            residual = max(residual, float(np.max(np.abs(incoming_should_vanish)) / abs(a)))
            b1, b2 = up[0], down[1]  # outgoing: +z wave at plane 1, -z wave at plane 2
            s_blocks[f"S1{port_index + 1}"][:, pol] = b1 / a
            s_blocks[f"S2{port_index + 1}"][:, pol] = b2 / a

    return SParameters(
        frequency=green.frequency,
        planes=planes,
        grounded=False,
        s11=s_blocks["S11"],
        s12=s_blocks["S12"],
        s21=s_blocks["S21"],
        s22=s_blocks["S22"],
        incident_amplitudes=amplitudes,
        incoming_residual=residual,
        mean_currents=mean_currents,
    )


def _scattering_from_fields_grounded(green, planes, mean_currents):
    """One-port version of ``_scattering_from_fields_two_port``: only port 1 (top) exists, so
    only S11 is computed (Phase 12A; S12/S21/S22 are undefined, ``SParameters.grounded=True``)."""
    z_top = planes.z_top
    s11 = np.zeros((2, 2), dtype=complex)
    amplitudes = np.zeros((1, 2), dtype=complex)
    residual = 0.0

    for pol in range(2):
        incident = np.eye(2)[pol]
        e_stack, h_stack = green.stack_fields(0.0, 0.0, incident, z_top, "top")
        e_sheet, h_sheet = green.sheet_fields(0.0, 0.0, mean_currents[0, pol], z_top)
        up, down = decompose_air_fields(e_stack + e_sheet, h_stack + h_sheet)  # (comp,)

        a = down[pol]  # incident wave = down-going at the top plane
        incoming_should_vanish = np.array([down[1 - pol]])  # nothing incoming in the other pol
        amplitudes[0, pol] = a
        residual = max(residual, float(np.max(np.abs(incoming_should_vanish)) / abs(a)))
        s11[:, pol] = up / a  # outgoing: +z wave at the top plane

    return SParameters(
        frequency=green.frequency,
        planes=planes,
        grounded=True,
        s11=s11,
        s12=None,
        s21=None,
        s22=None,
        incident_amplitudes=amplitudes,
        incoming_residual=residual,
        mean_currents=mean_currents,
    )


def scattering_from_fields(green: LayeredGreenFunction, planes: ReferencePlanes, mean_currents):
    """S-parameters from the (0,0) fields at the reference planes.

    ``mean_currents`` has shape (2 ports, 2 pols, 2) for the default ungrounded stack: the
    Floquet coefficient J~_00 of the solved current for each excitation (all zeros for a cell
    without metal). For a grounded stack (``green.backing == "pec"``) there is only one port, so
    the shape is (1, 2, 2) and the result's ``s12/s21/s22`` are ``None`` (Phase 12A).
    """
    _check_planes(green, planes)
    grounded = green.backing == "pec"
    mean_currents = np.asarray(mean_currents, dtype=complex)
    expected = (1, 2, 2) if grounded else (2, 2, 2)
    if mean_currents.shape != expected:
        raise ValueError(f"mean_currents must have shape {expected}, got {mean_currents.shape}")
    if grounded:
        return _scattering_from_fields_grounded(green, planes, mean_currents)
    return _scattering_from_fields_two_port(green, planes, mean_currents)


def stack_s_parameters(green: LayeredGreenFunction, planes: ReferencePlanes):
    """S-parameters of the bare dielectric slab (no metal): the Fresnel/Airy slab result
    (grounded: the bare grounded-slab reflection, S11 only)."""
    shape = (1, 2, 2) if green.backing == "pec" else (2, 2, 2)
    return scattering_from_fields(green, planes, np.zeros(shape))


def _warn_if_grating_orders(unit_cell, green, kx0=0.0, ky0=0.0):
    lat = unit_cell.lattice
    m, n = lat.mode_grid(3)
    kx, ky = lat.transverse_wavenumbers(kx0, ky0, m, n)
    kz = green.normal_wavenumbers(kx, ky)[0]
    propagating = (kz.imag == 0.0) & ~((m == 0) & (n == 0))
    if np.any(propagating):
        warnings.warn(
            f"{int(propagating.sum())} higher Floquet order(s) propagate in air at this frequency "
            "(period >= wavelength): S-parameters describe the specular (0,0) order only and "
            "energy conservation will not hold for them",
            UserWarning,
            stacklevel=3,
        )


# Relative (of k0) tolerance for the surface-wave-proximity diagnostic below: "near" the interval
# means within this fraction of k0 of either edge, not just strictly inside it. A policy choice
# (like the a-priori mesh/wavelength checks in discretization.py), not a theorem.
SURFACE_WAVE_TOLERANCE_RTOL = 0.02


def surface_wave_interval(green: LayeredGreenFunction):
    """(k0, Re(k_substrate)) [rad/m]: the transverse-wavenumber interval containing a grounded
    dielectric slab's TM0 surface-wave pole (Phase 12A). A grounded slab's TM0 surface wave has
    NO cutoff -- it exists for any h > 0 at every frequency, unlike the ungrounded slab's guided
    modes (GREEN_FUNCTION.md Sec. 7) -- and its propagation constant kt always lies in this
    interval (Collin, *Field Theory of Guided Waves*; Pozar, *Microwave Engineering*, grounded-
    slab surface-wave section). This is the *interval* the brief asks to check Floquet harmonics
    against, not the pole's exact location (which needs solving a transcendental eigenvalue
    equation) -- an a priori check, in the same spirit as ``_warn_if_grating_orders`` above.
    """
    return float(green.k0), float(green.k_substrate.real)


def surface_wave_proximity(unit_cell, green, floquet_order, tolerance_rtol=SURFACE_WAVE_TOLERANCE_RTOL):
    """Sec. "Physics" (Phase 12A): whether any Floquet harmonic used in the impedance-matrix
    assembly has kt inside, or within ``tolerance_rtol`` of k0 from the edges of,
    ``surface_wave_interval`` -- only meaningful for a grounded stack (``green.backing=="pec"``).

    Returns a dict (also usable as a per-record diagnostic, like ``specular_only``):
    ``k0``, ``k_substrate_re`` (the interval), ``any_in_band`` (strictly inside),
    ``any_near_band`` (inside or within tolerance), ``nearest_margin_over_k0`` (the closest
    harmonic's distance to the nearest edge of the interval, in units of k0; 0 exactly at an
    edge, negative would mean inside -- reported as the unsigned distance to the nearest edge,
    with ``any_in_band`` telling you which side), ``tolerance_rtol``.
    """
    lo, hi = surface_wave_interval(green)
    m, n = unit_cell.lattice.mode_grid(floquet_order)
    kx, ky = unit_cell.lattice.transverse_wavenumbers(0.0, 0.0, m, n)
    kt = np.hypot(kx, ky)
    inside = (kt >= lo) & (kt <= hi)
    tol = tolerance_rtol * lo
    near = inside | (np.abs(kt - lo) <= tol) | (np.abs(kt - hi) <= tol)
    margin = float(np.min(np.minimum(np.abs(kt - lo), np.abs(kt - hi)))) if kt.size else float("inf")
    return {
        "k0": lo,
        "k_substrate_re": hi,
        "any_in_band": bool(np.any(inside)),
        "any_near_band": bool(np.any(near)),
        "nearest_margin_over_k0": margin / lo,
        "tolerance_rtol": float(tolerance_rtol),
    }


def _warn_if_surface_wave_near(unit_cell, green, floquet_order):
    """Like ``_warn_if_grating_orders``: a UserWarning (not an error) when a Floquet harmonic
    used in this assembly falls inside, or near, the grounded slab's TM0 surface-wave interval."""
    r = surface_wave_proximity(unit_cell, green, floquet_order)
    if r["any_near_band"]:
        warnings.warn(
            f"a Floquet harmonic is {'inside' if r['any_in_band'] else 'within tolerance of'} the "
            f"grounded slab's TM0 surface-wave interval (k0={r['k0']:.6g}, "
            f"Re(k_substrate)={r['k_substrate_re']:.6g} rad/m; margin "
            f"{r['nearest_margin_over_k0']:.3g} x k0): S-parameters near a surface-wave resonance "
            f"may be numerically sensitive",
            UserWarning,
            stacklevel=3,
        )


def compute_s_parameters(
    unit_cell: UnitCell,
    frequency,
    mesh: MeshSettings,
    planes: ReferencePlanes,
    warn_grating=True,
    warn_surface_wave=True,
):
    """MoM solve for both polarizations (and, for the default ungrounded stack, both ports),
    then S-parameters at ``planes``.

    Normal incidence. Z is assembled once; the right-hand sides (port x polarization; only 2,
    top-port-only, for a grounded ``unit_cell.backing == "pec"``, Phase 12A) share one
    factorization. ``warn_grating=False`` suppresses the UserWarning about propagating higher
    Floquet orders (sweeps record that condition per frequency instead); ``warn_surface_wave=False``
    similarly suppresses the Phase 12A surface-wave-pole-proximity UserWarning (grounded stacks
    only; a no-op for an ungrounded one).
    """
    # getattr: duck-typed "UnitCell-like" objects that predate Phase 12A (e.g.
    # phase9.supercell_geometry.SupercellUnitCell) have no `.backing` attribute at all -- absent
    # means "air", exactly like every other pre-Phase-12A context/config (CONVENTIONS.md).
    backing = getattr(unit_cell, "backing", "air")
    green = LayeredGreenFunction(unit_cell.substrate, frequency, backing=backing)
    _check_planes(green, planes)
    if warn_grating:
        _warn_if_grating_orders(unit_cell, green)
    if warn_surface_wave and backing == "pec":
        _warn_if_surface_wave_near(unit_cell, green, mesh.floquet_order)
    basis = RooftopBasis(unit_cell, mesh.cells_x, mesh.cells_y)
    z_matrix = assemble_impedance_matrix(basis, green, 0.0, 0.0, mesh.floquet_order)

    ports = ("top",) if backing == "pec" else PORTS
    v = np.empty((basis.n_unknowns, 2 * len(ports)), dtype=complex)
    for port_index, port in enumerate(ports):
        for pol in range(2):
            exc = PlaneWaveExcitation(green.frequency, polarization=POLARIZATIONS[pol])
            e_stack = stack_tangential_field(green, exc, port)
            v[:, 2 * port_index + pol] = assemble_excitation_vector(basis, e_stack, 0.0, 0.0)
    currents = solve(z_matrix, v)

    mean = np.empty((len(ports), 2, 2), dtype=complex)
    for port_index in range(len(ports)):
        for pol in range(2):
            mean[port_index, pol] = basis.current_coefficients(currents[:, 2 * port_index + pol], 0.0, 0.0)
    return scattering_from_fields(green, planes, mean)


def s_parameters_from_config(config):
    """``compute_s_parameters`` from a validated SimulationConfig (normal incidence, x and y both); for ``theta != 0`` the
    oblique solver (TE/TM basis, ``ObliqueSParameters``)."""
    if config.excitation.theta != 0.0:  # Phase 20 R1: oblique incidence -> TE/TM basis (ObliqueSParameters)
        return compute_s_parameters_oblique(
            config.unit_cell, config.excitation.frequency, config.mesh, config.reference_planes, config.excitation.theta, config.excitation.phi
        )
    return compute_s_parameters(config.unit_cell, config.excitation.frequency, config.mesh, config.reference_planes)


# =============================================================================================================
# Phase 20 R1 item 1: oblique incidence, S-parameters in the TE / TM basis (additive; the normal-incidence
# functions above are unchanged and still used for theta = 0 in the x/y basis).
# =============================================================================================================
#
# Basis (docs/theory/S_PARAMETERS.md Sec. 8): for the incident wave (theta, phi), p = (cos phi, sin phi) is the in-plane
# direction of the incident transverse wavevector and q = z x p = (-sin phi, cos phi).
#   TE: E along q (perpendicular to the plane of incidence);  TM: tangential E along p (E in the plane of incidence).
# Index 0 = TE, 1 = TM.  At theta = 0, phi = 0: TE = y, TM = x.
# The (0,0) harmonic of each wave is split into up- and down-going parts per mode with the mode characteristic
# admittances  Y_TE = kz0/(w mu0),  Y_TM = w eps0/kz0  (kz0 = k0 cos theta):  V_up/dn = (V +/- I/Y)/2.  Wave amplitudes are
# tangential-E amplitudes; the returned S are POWER-wave normalised, S_ij = sqrt(Y_i/Y_j) V_i/V_j (equal to the
# E-amplitude ratio for co-pol terms; a TE<->TM term carries sqrt(Y_TM/Y_TE) = 1/cos theta), so that for a lossless
# structure the sum over output modes and ports of |S|^2 is 1 when only the (0,0) order propagates.
TE, TM = 0, 1
_OBLIQUE_MODES = ("TE", "TM")


def oblique_basis_vectors(phi):
    """`(q, p)` = (TE, TM) in-plane unit vectors for azimuth `phi`: q = (-sin phi, cos phi), p = (cos phi, sin phi)."""
    p = np.array([math.cos(phi), math.sin(phi)])
    q = np.array([-math.sin(phi), math.cos(phi)])
    return q, p


def _mode_admittances(frequency, theta):
    k0 = 2.0 * math.pi * frequency / C0
    kz0 = k0 * math.cos(theta)
    w = 2.0 * math.pi * frequency
    return kz0 / (w * MU0), w * EPS0 / kz0


def decompose_oblique_fields(e_t, h_t, phi, y_te, y_tm):
    """Split tangential (E, H) of a (0,0) wave in air at (theta, phi) into per-mode up- and down-going E amplitudes.

    Returns `(up, down)`, each `(..., 2)` ordered (TE, TM). E_t = V_TM p + V_TE q ... and H_t = -I_TE p + I_TM q
    (`LayeredGreenFunction._apply_modes`), so  V_TE = E.q, I_TE = -H.p, V_TM = E.p, I_TM = H.q;
    up = (V + I/Y)/2, down = (V - I/Y)/2.
    """
    q, p = oblique_basis_vectors(phi)
    e_t, h_t = np.asarray(e_t, dtype=complex), np.asarray(h_t, dtype=complex)
    v_te, i_te = e_t @ q, -(h_t @ p)
    v_tm, i_tm = e_t @ p, h_t @ q
    up = np.stack([0.5 * (v_te + i_te / y_te), 0.5 * (v_tm + i_tm / y_tm)], axis=-1)
    down = np.stack([0.5 * (v_te - i_te / y_te), 0.5 * (v_tm - i_tm / y_tm)], axis=-1)
    return up, down


@dataclass(frozen=True, eq=False)
class ObliqueSParameters:
    """S-parameters of one unit cell at one frequency and one incidence direction, in the TE / TM basis.

    Each of `s11, s12, s21, s22` is a (2, 2) complex array `[received, incident]`, index 0 = TE, 1 = TM, POWER-wave
    normalised (see the module block comment). Grounded stack: only `s11` (the others `None`, undefined not zero).
    `reflection` is the full 2x2 reflection (co- and cross-pol) at the top reference plane.
    """

    frequency: float
    theta: float  # [rad]
    phi: float  # [rad]
    planes: ReferencePlanes
    grounded: bool
    s11: np.ndarray
    s12: np.ndarray | None
    s21: np.ndarray | None
    s22: np.ndarray | None
    incident_amplitudes: np.ndarray  # (n_ports, 2 modes) tangential-E amplitude measured at the plane
    incoming_residual: float
    mean_currents: np.ndarray  # (n_ports, 2 modes, 2) J~_00 [A/m]
    higher_orders_propagate: bool

    @property
    def k0(self):
        return 2.0 * math.pi * self.frequency / C0

    @property
    def kz0(self):
        return self.k0 * math.cos(self.theta)

    @property
    def reflection(self):
        return self.s11

    def component(self, name):
        """`"S11_TE_TM"` = TE received, TM incident (family, received, incident)."""
        fam, rx, inc = name.split("_")
        block = getattr(self, fam.lower())
        if block is None:
            raise ValueError(f"{fam} is undefined for a grounded (one-port) stack")
        return complex(block[_OBLIQUE_MODES.index(rx), _OBLIQUE_MODES.index(inc)])

    def magnitude_db(self, name):
        return 20.0 * np.log10(max(abs(self.component(name)), 1e-300))

    def phase_deg(self, name):
        return float(np.degrees(np.angle(self.component(name))))

    def power_out_per_incident(self):
        """(n_ports, 2): sum over all output ports/modes of |S|^2 for each incident mode (1 for a lossless cell whose
        only propagating order is (0,0))."""
        if self.grounded:
            return np.sum(np.abs(self.s11) ** 2, axis=0).reshape(1, 2)
        m = np.block([[self.s11, self.s12], [self.s21, self.s22]])
        return np.sum(np.abs(m) ** 2, axis=0).reshape(2, 2)

    def to_record(self):
        rec = {"frequency_hz": self.frequency, "theta_deg": math.degrees(self.theta), "phi_deg": math.degrees(self.phi),
               "higher_orders_propagate": self.higher_orders_propagate}
        for fam in ("S11",) if self.grounded else _FAMILIES:
            for rx in _OBLIQUE_MODES:
                for inc in _OBLIQUE_MODES:
                    key = f"{fam}_{rx}_{inc}"
                    rec[f"{key}_mag"] = abs(self.component(key))
                    rec[f"{key}_phase_deg"] = self.phase_deg(key)
        return rec

    def at_planes(self, new_planes):
        """Re-reference by (0,0) propagation in air at the oblique angle: phases use kz0 = k0 cos(theta)."""
        if not np.isclose(new_planes.substrate_thickness, self.planes.substrate_thickness, rtol=1e-12):
            raise ValueError("substrate thickness differs between plane sets")
        if self.grounded != new_planes.grounded:
            raise ValueError("cannot re-reference between a grounded (one-port) and an ungrounded (two-port) plane set")
        d1 = new_planes.z_top - self.planes.z_top
        kz = self.kz0
        if self.grounded:
            return replace(self, planes=new_planes, s11=self.s11 * np.exp(-2j * kz * d1),
                           incident_amplitudes=self.incident_amplitudes * np.exp(1j * kz * d1))
        d2 = new_planes.z_bottom - self.planes.z_bottom
        return replace(
            self, planes=new_planes, s11=self.s11 * np.exp(-2j * kz * d1), s22=self.s22 * np.exp(+2j * kz * d2),
            s21=self.s21 * np.exp(+1j * kz * (d2 - d1)), s12=self.s12 * np.exp(+1j * kz * (d2 - d1)),
            incident_amplitudes=self.incident_amplitudes * np.array([[np.exp(1j * kz * d1)] * 2, [np.exp(-1j * kz * d2)] * 2]),
        )


def _oblique_scattering(green, planes, mean_currents, kx0, ky0, theta, phi, higher_orders):
    _check_planes(green, planes)
    grounded = green.backing == "pec"
    mean_currents = np.asarray(mean_currents, dtype=complex)
    q, p = oblique_basis_vectors(phi)
    e_inc = (q, p)  # mode 0 = TE, 1 = TM: unit tangential E of the incident wave
    y_te, y_tm = _mode_admittances(green.frequency, theta)
    y = np.array([y_te, y_tm])
    ports = ("top",) if grounded else PORTS
    z_planes = np.array([planes.z_top] if grounded else [planes.z_top, planes.z_bottom])
    blocks = {name: np.zeros((2, 2), dtype=complex) for name in _FAMILIES}
    amplitudes = np.zeros((len(ports), 2), dtype=complex)
    residual = 0.0
    for pi, port in enumerate(ports):
        for mode in range(2):
            e_stack, h_stack = green.stack_fields(kx0, ky0, e_inc[mode].astype(complex), z_planes, port)
            e_sheet, h_sheet = green.sheet_fields(kx0, ky0, mean_currents[pi, mode], z_planes)
            up, down = decompose_oblique_fields(e_stack + e_sheet, h_stack + h_sheet, phi, y_te, y_tm)  # (plane, mode)
            if port == "top":  # incident = down-going at plane 1; nothing incoming at plane 2 or in the other mode
                a = down[0, mode]
                incoming = [down[0, 1 - mode]] + ([] if grounded else list(up[1]))
            else:  # incident = up-going at plane 2; nothing incoming at plane 1 or in the other mode
                a = up[1, mode]
                incoming = [up[1, 1 - mode]] + list(down[0])
            amplitudes[pi, mode] = a
            residual = max(residual, float(np.max(np.abs(incoming)) / abs(a)))
            scale = np.sqrt(y / y[mode])  # power-wave normalisation  sqrt(Y_out / Y_in)
            blocks[f"S1{pi + 1}"][:, mode] = up[0] / a * scale
            if not grounded:
                blocks[f"S2{pi + 1}"][:, mode] = down[1] / a * scale
    return ObliqueSParameters(
        frequency=green.frequency, theta=theta, phi=phi, planes=planes, grounded=grounded, s11=blocks["S11"],
        s12=None if grounded else blocks["S12"], s21=None if grounded else blocks["S21"], s22=None if grounded else blocks["S22"],
        incident_amplitudes=amplitudes, incoming_residual=residual, mean_currents=mean_currents, higher_orders_propagate=higher_orders,
    )


def oblique_stack_s_parameters(green: LayeredGreenFunction, planes: ReferencePlanes, theta, phi=0.0):
    """TE/TM S-parameters of the bare slab (no metal) at (theta, phi): the Fresnel/Airy result."""
    k0 = green.k0
    kx0, ky0 = k0 * math.sin(theta) * math.cos(phi), k0 * math.sin(theta) * math.sin(phi)
    shape = (1, 2, 2) if green.backing == "pec" else (2, 2, 2)
    return _oblique_scattering(green, planes, np.zeros(shape), kx0, ky0, theta, phi, False)


def _higher_orders_at(unit_cell, green, kx0, ky0):
    lat = unit_cell.lattice
    m, n = lat.mode_grid(3)
    kx, ky = lat.transverse_wavenumbers(kx0, ky0, m, n)
    kz = green.normal_wavenumbers(kx, ky)[0]
    return bool(np.any((kz.imag == 0.0) & ~((m == 0) & (n == 0))))


def compute_s_parameters_oblique(
    unit_cell: UnitCell, frequency, mesh: MeshSettings, planes: ReferencePlanes, theta, phi=0.0, warn_grating=True,
    warn_surface_wave=True,
):
    """MoM solve at incidence `(theta, phi)` [rad], 0 <= theta < pi/2, for both modes (TE, TM) and (ungrounded) both ports;
    returns `ObliqueSParameters` (TE/TM basis, full 2x2). Z is assembled once at the incident transverse wavenumber
    `(kx0, ky0)`; the right-hand sides share one factorisation. At theta = 0 the TE/TM labels are x/y rotated by `phi`.

    `warn_grating` warns (UserWarning) when a higher Floquet order propagates at this angle; the same fact is recorded in
    `higher_orders_propagate` (and `to_record`) so sweeps carry it.
    """
    if not 0.0 <= theta < math.pi / 2:
        raise ValueError(f"theta must satisfy 0 <= theta < pi/2, got {theta}")
    backing = getattr(unit_cell, "backing", "air")
    green = LayeredGreenFunction(unit_cell.substrate, frequency, backing=backing)
    _check_planes(green, planes)
    k0 = green.k0
    kx0, ky0 = k0 * math.sin(theta) * math.cos(phi), k0 * math.sin(theta) * math.sin(phi)
    higher = _higher_orders_at(unit_cell, green, kx0, ky0)
    if warn_grating and higher:
        _warn_if_grating_orders(unit_cell, green, kx0, ky0)
    if warn_surface_wave and backing == "pec":
        _warn_if_surface_wave_near(unit_cell, green, mesh.floquet_order)
    basis = RooftopBasis(unit_cell, mesh.cells_x, mesh.cells_y)
    z_matrix = assemble_impedance_matrix(basis, green, kx0, ky0, mesh.floquet_order)
    q, p = oblique_basis_vectors(phi)
    e_inc = (q, p)
    ports = ("top",) if backing == "pec" else PORTS
    v = np.empty((basis.n_unknowns, 2 * len(ports)), dtype=complex)
    for pi, port in enumerate(ports):
        for mode in range(2):
            e_stack = green.stack_fields(kx0, ky0, e_inc[mode].astype(complex), 0.0, port)[0]
            v[:, 2 * pi + mode] = assemble_excitation_vector(basis, e_stack, kx0, ky0)
    currents = solve(z_matrix, v)
    mean = np.empty((len(ports), 2, 2), dtype=complex)
    for pi in range(len(ports)):
        for mode in range(2):
            mean[pi, mode] = basis.current_coefficients(currents[:, 2 * pi + mode], kx0, ky0)
    return _oblique_scattering(green, planes, mean, kx0, ky0, theta, phi, higher)


def sweep_s_parameters_oblique(unit_cell, frequencies, mesh, planes, theta, phi=0.0):
    """`compute_s_parameters_oblique` over frequencies; returns `(results, records)` where each record carries
    `higher_orders_propagate` for that frequency (warnings suppressed: the flag is the record)."""
    results = [compute_s_parameters_oblique(unit_cell, f, mesh, planes, theta, phi, warn_grating=False, warn_surface_wave=False) for f in frequencies]
    return results, [r.to_record() for r in results]
