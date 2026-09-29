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

import warnings
from dataclasses import dataclass, replace

import numpy as np

from .basis_functions import MeshSettings, RooftopBasis
from .constants import C0, ETA0
from .excitation import PlaneWaveExcitation, ReferencePlanes, stack_tangential_field
from .geometry import UnitCell
from .green_function import LayeredGreenFunction
from .mom_matrix import assemble_excitation_vector, assemble_impedance_matrix
from .solver import solve

PORTS = ("top", "bottom")       # port 1, port 2
POLARIZATIONS = "xy"            # index 0, 1
_COMPONENT_INDEX = {"xx": (0, 0), "yx": (1, 0), "xy": (0, 1), "yy": (1, 1)}   # (received, incident)
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
    incident_amplitudes: np.ndarray    # (n_ports, 2 pols): measured a at the reference plane
    incoming_residual: float           # largest incoming-wave amplitude that must vanish, / |a|
    mean_currents: np.ndarray          # (n_ports, 2 pols, 2): J~_00 [A/m] (zeros: no metal)

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
            raise ValueError(
                "no 4x4 S-matrix exists for a grounded (one-port) stack (no port 2); use .s11"
            )
        return np.block([[self.s11, self.s12], [self.s21, self.s22]])

    def component(self, name):
        """One complex entry, e.g. ``component("S21_yx")`` = y received at port 2, x incident at port 1.

        Raises ``ValueError`` for a family other than S11 of a grounded result (undefined, not 0).
        """
        family, _, comp = name.partition("_")
        if family not in _FAMILIES or comp not in _COMPONENT_INDEX:
            raise ValueError(
                f"unknown component {name!r}; use <{'|'.join(_FAMILIES)}>_<{'|'.join(_COMPONENT_INDEX)}>"
            )
        block = getattr(self, family.lower())
        if block is None:
            raise ValueError(
                f"{family} is undefined for a grounded (one-port) stack (no port 2 exists); "
                f"only S11 is defined"
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
            raise ValueError(
                "cannot re-reference between a grounded (one-port) and an ungrounded (two-port) "
                "plane set"
            )
        d1 = new_planes.z_top - self.planes.z_top
        k0 = self.k0
        if self.grounded:
            return replace(
                self, planes=new_planes,
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
            up, down = decompose_air_fields(e_stack + e_sheet, h_stack + h_sheet)   # (plane, comp)

            if port == "top":       # incident wave = down-going at plane 1; nothing incoming at plane 2
                a = down[0, pol]
                incoming_should_vanish = np.concatenate([[down[0, 1 - pol]], up[1]])
            else:                   # incident wave = up-going at plane 2; nothing incoming at plane 1
                a = up[1, pol]
                incoming_should_vanish = np.concatenate([[up[1, 1 - pol]], down[0]])
            amplitudes[port_index, pol] = a
            residual = max(residual, float(np.max(np.abs(incoming_should_vanish)) / abs(a)))
            b1, b2 = up[0], down[1]         # outgoing: +z wave at plane 1, -z wave at plane 2
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
        up, down = decompose_air_fields(e_stack + e_sheet, h_stack + h_sheet)     # (comp,)

        a = down[pol]                       # incident wave = down-going at the top plane
        incoming_should_vanish = np.array([down[1 - pol]])   # nothing incoming in the other pol
        amplitudes[0, pol] = a
        residual = max(residual, float(np.max(np.abs(incoming_should_vanish)) / abs(a)))
        s11[:, pol] = up / a                # outgoing: +z wave at the top plane

    return SParameters(
        frequency=green.frequency,
        planes=planes,
        grounded=True,
        s11=s11, s12=None, s21=None, s22=None,
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


def _warn_if_grating_orders(unit_cell, green):
    lat = unit_cell.lattice
    m, n = lat.mode_grid(3)
    kx, ky = lat.transverse_wavenumbers(0.0, 0.0, m, n)
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
    unit_cell: UnitCell, frequency, mesh: MeshSettings, planes: ReferencePlanes, warn_grating=True,
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
    """``compute_s_parameters`` from a validated SimulationConfig (normal incidence, x and y both)."""
    if config.excitation.theta != 0.0:
        raise ValueError("S-parameter extraction supports normal incidence only")
    return compute_s_parameters(
        config.unit_cell, config.excitation.frequency, config.mesh, config.reference_planes
    )
