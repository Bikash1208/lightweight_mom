"""Phase 15 Part B: a VECTOR aperture-to-far-field model, replacing `aperture_pattern.py`'s scalar
one for HMSA. Methodology source: `docs/briefs/PHASE15_POLARIZATION_VORTEX_MULTIFEED.md` Part B.

## Why

`aperture_pattern.compute_pattern` is a bare scalar 2-D Fourier transform of the aperture field
(`phase8e.beam_reconstruction.far_field_pattern`) -- it has no notion of polarization at all, so it
cannot see the vortex this phase's brief measured in HFSS: a monopole-fed HMSA's surface current is
radial, `J(x,y) = rho_hat(x,y) * J_s(x,y)`, and radiating a radial current necessarily carries an
orbital-angular-momentum term in its circularly-polarized far-field components.

## Method (reciprocity via this project's own layered Green's function, not a free-space dipole)

For a grounded stack (the HMSA board: air | patch/current plane | substrate | PEC ground,
`green_function.LayeredGreenFunction(substrate, frequency_hz, backing="pec")`), the tangential field
radiated by a current sheet AT the current plane, for one spectral (Floquet) component, is exactly
`G(kx, ky) @ J~(kx, ky)` (`LayeredGreenFunction.spectral_dyad`, already validated against an
independent field-matching reference solver, `test_green_function_fields.py`). By reciprocity, that
SAME spectral response, evaluated at the *continuous* (not periodic-lattice-quantized) point
`(kx, ky) = k0*sin(theta)*(cos(phi), sin(phi))`, is proportional to the far field radiated toward
`(theta, phi)` -- this is the standard spectral-domain/reciprocity far-field extraction technique for
a planar current on a layered medium (e.g. Pozar's microstrip-antenna spectral-domain papers; the
`cos(theta)` prefactor below is the standard stationary-phase Jacobian factor of that same
extraction, not an extra approximation on top of it).

Steps, per output `(theta, phi)`:
1. `J~x, J~y` = the SAME direct-summation DFT `phase8e.beam_reconstruction.far_field_pattern`
   already uses for a scalar aperture field, applied once each to the vector current's `Jx`, `Jy`
   (`_aperture_dft`, generalizing that summation, not reinventing it).
2. `(Ex, Ey) = G(kx, ky) @ (J~x, J~y)`, `G = green_function.spectral_dyad(kx, ky)`.
3. Apply `cos(theta)` (this project's existing obliquity factor, `aperture_pattern.py`'s own,
   cited-Balanis, PHASE13D-derived -- reused unchanged: the SAME Jacobian factor, now applied on
   top of the correct layered-medium spectral response instead of the old free-space assumption).
4. Project Cartesian -> spherical (Ludwig-3, standard antenna theory, e.g. Balanis "Antenna
   Theory," aperture/equivalent-current far-field projection):
       E_theta = Ex*cos(theta)*cos(phi) + Ey*cos(theta)*sin(phi)
       E_phi   = -Ex*sin(phi) + Ey*cos(phi)
5. Local circular components, Ludwig-3's boresight-referenced construction (theta_hat, phi_hat
   literally equal x_hat, y_hat at boresight, and vary smoothly off-axis -- exactly what makes this
   substitution valid everywhere in the pattern, not just at theta=0): reuse
   `polarization.amplitudes_to_circular` UNCHANGED, treating `(E_theta, E_phi)` as a direct
   substitute for `(E_x, E_y)` (`circular_components`, below) -- this is what keeps the R/L
   handedness labels consistent with this project's one fixed convention (`polarization.E_R`/`E_L`,
   `CONVENTIONS.md`'s Phase 5 section), rather than a hand-rolled formula that could silently swap
   them (an earlier version of this module did exactly that: it built `E_R = (E_theta+j*E_phi)/sqrt(2)`
   directly instead of calling `amplitudes_to_circular`, which -- because that function's own
   internal convention is `a_R = (V_x - j*V_y)/sqrt(2)`, the conjugate-basis decomposition formula,
   not the naive `V_x + j*V_y` guess -- had R and L exactly swapped; caught by checking a
   known-polarization beam's radiated field at exact boresight against `polarization.handedness`
   during development, not by inspection).

## Which circular label is "co-polar" for a compensated beam

`copol_label` maps a `Beam.desired_polarization` ("rhcp"/"lhcp") to whichever of
`circular_components`' outputs the compensated beam actually radiates into. Derived (not assumed):
`hologram.polarization_compensation`'s phi_ext is built so `conj(p_hat_b) . rho_hat(x,y)`'s phase is
cancelled, i.e. the aperture's own `conj(p_hat_b) . J(x,y)` ends up with the correctly-steered
phase -- and `conj(p_hat_b) . J` is EXACTLY `amplitudes_to_circular`'s own decomposition coefficient
of `J` along `p_hat_b` (`a_v(J) = conj(v) . J` for an orthonormal basis vector `v`). At boresight the
medium's spectral Green's function is isotropic (`spectral_dyad`'s own `kt=0` branch: `G = -Z_avg *
identity`), so `E_far = -Z_avg * J` exactly and the decomposition commutes: the far field's own
`p_hat_b`-labeled coefficient is what's dominant -- i.e. `copol_label` returns the SAME label as
`p_hat_b` itself, with NO inversion (confirmed numerically at boresight during development: a
"rhcp"-compensated broadside beam radiates >60dB more power into `amplitudes_to_circular`'s own
`a_L` output than `a_R`, matching `p_hat_b = E_L` for "rhcp", `hologram.polarization_compensation.
_circular_jones_vector`). Off-boresight the medium's TE/TM anisotropy (`Z_TE != Z_TM`) mixes the two
circular components somewhat -- this is why XPD is finite (tens of dB, not infinite) away from
boresight in this module's own tests, a real, documented limitation of a scalar/isotropic unit cell
(the brief's own "tensor/anisotropic cells are the next step"), not a bug.
"""

from dataclasses import dataclass

import numpy as np

from ..frequency import free_space_wavenumber
from ..green_function import LayeredGreenFunction
from ..phase8e.beam_reconstruction import detect_local_maxima
from ..polarization import amplitudes_to_circular, handedness

_SQRT_HALF = 2.0**-0.5


def aperture_current_from_scalar(x_m, y_m, scalar_current, feed_xy_m=(0.0, 0.0), leakage_alpha_per_m=0.0):
    """`J = rho_hat(x,y) * scalar_current(x,y) * exp(-alpha*rho)` -- Part B's VECTORIZATION of the
    scalar aperture current magnitude Part A/PHASE14 already compute (`hmsa_synthesis`'s
    `aperture_field`, already documented there as `~ J_s`), about `feed_xy_m` (Part C: the assigned
    feed for a multi-feed design). `leakage_alpha_per_m` (the brief's Part B: "an exponential
    leakage decay e^{-alpha*rho}, alpha from the paper's eq. 11-12 form") is an OPTIONAL extra
    amplitude taper representing a leaky-wave aperture's own power loss along its radial extent --
    default 0 (off, `exp(0)=1` identically, no change to any pre-existing caller) since this
    project's own dispersion model (`dispersion.sheet_model`) does not itself derive an alpha value
    (Fong 2010/the paper's eq. 11-12 do, from the unit cell's own loss/coupling, which this
    project's sheet-admittance model does not currently extract) -- a caller who has an independent
    alpha estimate may supply it, documented as exactly that, not a value this project computes.
    Returns `(Jx, Jy)`, same shape as `x_m`."""
    x_m, y_m = np.asarray(x_m, dtype=float), np.asarray(y_m, dtype=float)
    x0, y0 = feed_xy_m
    rho = np.hypot(x_m - x0, y_m - y0)
    phi_prime = np.arctan2(y_m - y0, x_m - x0)
    scalar_current = np.asarray(scalar_current, dtype=complex)
    if leakage_alpha_per_m:
        scalar_current = scalar_current * np.exp(-leakage_alpha_per_m * rho)
    return scalar_current * np.cos(phi_prime), scalar_current * np.sin(phi_prime)


def _aperture_dft(x_flat, y_flat, values_flat, kx, ky, dx_m, dy_m):
    """The exact DFT-summation algebra of `phase8e.beam_reconstruction.far_field_pattern`,
    factored out so it can be applied once per vector current component."""
    shape = kx.shape
    phase = np.outer(kx.ravel(), x_flat) + np.outer(ky.ravel(), y_flat)
    return ((np.exp(1j * phase) @ values_flat) * dx_m * dy_m).reshape(shape)


def _vector_far_field(theta_deg, phi_deg, x_m, y_m, Jx, Jy, frequency_hz, green_function, dx_m, dy_m):
    """`(E_theta, E_phi)` at arbitrary `(theta_deg, phi_deg)` points (any shape) -- the module
    docstring's steps 1-4, factored out so both the regular-grid pattern (`compute_vector_pattern`)
    and an arbitrary small ring of points (`measure_vortex_winding`) share one implementation."""
    k0 = free_space_wavenumber(frequency_hz)
    theta = np.radians(np.asarray(theta_deg, dtype=float))
    phi = np.radians(np.asarray(phi_deg, dtype=float))
    kx = k0 * np.sin(theta) * np.cos(phi)
    ky = k0 * np.sin(theta) * np.sin(phi)

    x_flat = np.asarray(x_m, dtype=float).ravel()
    y_flat = np.asarray(y_m, dtype=float).ravel()
    jx_tilde = _aperture_dft(x_flat, y_flat, np.asarray(Jx, dtype=complex).ravel(), kx, ky, dx_m, dy_m)
    jy_tilde = _aperture_dft(x_flat, y_flat, np.asarray(Jy, dtype=complex).ravel(), kx, ky, dx_m, dy_m)

    g = green_function.spectral_dyad(kx, ky)  # shape kx.shape + (2, 2)
    j_tilde = np.stack([jx_tilde, jy_tilde], axis=-1)[..., None]  # (..., 2, 1)
    e_spectral = np.squeeze(g @ j_tilde, axis=-1)  # (..., 2): (Ex, Ey) at the current plane

    obliquity = np.cos(theta)  # PHASE13D's own factor, reused -- see module docstring
    ex = e_spectral[..., 0] * obliquity
    ey = e_spectral[..., 1] * obliquity

    e_theta = ex * np.cos(theta) * np.cos(phi) + ey * np.cos(theta) * np.sin(phi)
    e_phi = -ex * np.sin(phi) + ey * np.cos(phi)
    return e_theta, e_phi


def circular_components(e_theta, e_phi):
    """Ludwig-3 local circular basis, via `polarization.amplitudes_to_circular` UNCHANGED (module
    docstring: `(E_theta, E_phi)` substituted for `(E_x, E_y)`). Returns `(E_R, E_L)`, this
    project's own fixed convention (`a_R = (V_x - j*V_y)/sqrt(2)`, `a_L = (V_x + j*V_y)/sqrt(2)`) --
    use `copol_label`/`copol_component` to map a `Beam.desired_polarization` to the physically-
    matching one of these two."""
    e_theta, e_phi = np.asarray(e_theta, dtype=complex), np.asarray(e_phi, dtype=complex)
    a = amplitudes_to_circular(np.stack([e_theta, e_phi], axis=-1))
    return a[..., 0], a[..., 1]


def copol_label(desired_polarization):
    """"rhcp"/"lhcp" -> "E_R" or "E_L", whichever `circular_components` output is co-polar with a
    beam compensated for that handedness (`hologram.polarization_compensation.polarization_jones_vector`
    picks `p_hat_b` via the SAME `handedness()` call used here) -- see module docstring ("Which
    circular label is co-polar") for the derivation: at boresight the medium is isotropic, so the
    far field's own `p_hat_b`-labeled decomposition coefficient is what's dominant, i.e. NO
    inversion relative to `p_hat_b`'s own label."""
    if desired_polarization not in ("rhcp", "lhcp"):
        raise ValueError(f"copol_label is only meaningful for 'rhcp'/'lhcp', got {desired_polarization!r}")
    wanted = "RHCP" if desired_polarization == "rhcp" else "LHCP"
    for circ_label, field_name in (("R", "E_R"), ("L", "E_L")):
        if handedness(circ_label, "+z") == wanted:
            return field_name
    raise AssertionError("unreachable: neither E_R nor E_L matches the requested handedness")


def copol_component(e_theta, e_phi, desired_polarization):
    """The complex field of whichever `circular_components` output is co-polar with
    `desired_polarization` ("rhcp"/"lhcp"); for "hp"/"vp"/"slant_p45"/"slant_m45" this is the
    corresponding LINEAR projection of `(E_theta, E_phi)` (the same Jones vectors
    `hologram.polarization_compensation._LINEAR_JONES_VECTORS` uses, applied in the local
    Ludwig-3 (theta_hat, phi_hat) frame instead of the global (x_hat, y_hat) one)."""
    e_theta, e_phi = np.asarray(e_theta, dtype=complex), np.asarray(e_phi, dtype=complex)
    if desired_polarization in ("rhcp", "lhcp"):
        e_r, e_l = circular_components(e_theta, e_phi)
        return e_r if copol_label(desired_polarization) == "E_R" else e_l
    linear = {
        "hp": (1.0, 0.0),
        "vp": (0.0, 1.0),
        "slant_p45": (_SQRT_HALF, _SQRT_HALF),
        "slant_m45": (_SQRT_HALF, -_SQRT_HALF),
    }
    if desired_polarization not in linear:
        raise ValueError(f"unknown desired_polarization {desired_polarization!r}")
    pt, pp = linear[desired_polarization]
    return np.conj(pt) * e_theta + np.conj(pp) * e_phi


@dataclass(frozen=True)
class VectorPatternResult:
    theta_deg: np.ndarray
    phi_deg: np.ndarray
    e_theta: np.ndarray
    e_phi: np.ndarray
    e_r: np.ndarray  # circular_components(e_theta, e_phi) -- see copol_label for which handedness
    e_l: np.ndarray
    power_db: np.ndarray  # |E_theta|^2 + |E_phi|^2, relative to the pattern's own peak, 0dB at peak
    main_beam_theta_deg: float
    main_beam_phi_deg: float
    sll_db: "float | None"
    directivity_estimate_dbi: float


def compute_vector_pattern(
    x_m, y_m, Jx, Jy, frequency_hz, green_function: LayeredGreenFunction, dx_m, dy_m,
    theta_max_deg=90.0, n_theta=91, n_phi=73,
):
    """Vector far field over a hemisphere, on an `n_theta` x `n_phi` grid -- the vector counterpart
    of `pattern.aperture_pattern.compute_pattern`, same grid/peak-normalization/peak-detection
    conventions (`detect_local_maxima`, reused unmodified) so results are directly comparable."""
    theta = np.linspace(0.0, theta_max_deg, n_theta)
    phi = np.linspace(0.0, 360.0, n_phi, endpoint=False)
    theta_grid, phi_grid = np.meshgrid(theta, phi, indexing="ij")

    e_theta, e_phi = _vector_far_field(theta_grid, phi_grid, x_m, y_m, Jx, Jy, frequency_hz, green_function, dx_m, dy_m)
    e_r, e_l = circular_components(e_theta, e_phi)

    power = np.abs(e_theta) ** 2 + np.abs(e_phi) ** 2
    peak = power.max()
    if peak <= 0.0:
        raise ValueError("compute_vector_pattern: aperture field radiates zero power (all-zero input?)")
    power_db = 10.0 * np.log10(power / peak)

    peaks = detect_local_maxima(power, theta_grid, phi_grid, min_relative_db=-40.0, phi_periodic=True)
    peaks = sorted(peaks, key=lambda p: -p.relative_power_db)
    main = peaks[0]
    sll_db = peaks[1].relative_power_db if len(peaks) > 1 else None

    directivity = _directivity_estimate(power, theta_grid)
    return VectorPatternResult(
        theta_deg=theta_grid, phi_deg=phi_grid, e_theta=e_theta, e_phi=e_phi, e_r=e_r, e_l=e_l,
        power_db=power_db, main_beam_theta_deg=main.theta_deg, main_beam_phi_deg=main.phi_deg,
        sll_db=sll_db, directivity_estimate_dbi=directivity,
    )


def _directivity_estimate(power, theta_grid_deg):
    """Same definition/limitations as `aperture_pattern._directivity_estimate` (hemisphere-only,
    no back-radiation notion) -- duplicated rather than imported because `aperture_pattern`'s is
    module-private; both compute the identical `4*pi*P_max / integral(P sin(theta) dtheta dphi)`."""
    theta_rad = np.radians(theta_grid_deg)
    dtheta = np.radians(theta_grid_deg[1, 0] - theta_grid_deg[0, 0]) if theta_grid_deg.shape[0] > 1 else 1.0
    dphi = np.radians(360.0 / theta_grid_deg.shape[1])
    solid_angle_weight = np.sin(theta_rad) * dtheta * dphi
    integral = float(np.sum(power * solid_angle_weight))
    if integral <= 0.0:
        return float("nan")
    return 10.0 * np.log10(4.0 * np.pi * power.max() / integral)


def aperture_efficiency_estimate(directivity_dbi, area_m2, frequency_hz):
    """Phase 15 Part C: `directivity / D_max`, `D_max = 4*pi*area/lambda^2` the standard maximum
    directivity of a uniformly-illuminated aperture of physical area `area_m2` (Balanis, "Antenna
    Theory," aperture antennas chapter) -- used to report a per-port aperture efficiency (Mishra's
    own brief: "Aperture efficiency per port"). Returns `(efficiency_fraction, max_directivity_dbi)`;
    `efficiency_fraction` can exceed 1 for a design whose OWN `_directivity_estimate` (hemisphere-
    only, no back-radiation notion -- see that function's docstring) overstates directivity, which
    is reported as-is rather than clamped (an honest symptom of that approximation, not hidden)."""
    from ..frequency import free_space_wavelength

    wavelength_m = free_space_wavelength(frequency_hz)
    max_directivity = 4.0 * np.pi * area_m2 / wavelength_m**2
    max_directivity_dbi = 10.0 * np.log10(max_directivity)
    efficiency = 10.0 ** ((directivity_dbi - max_directivity_dbi) / 10.0)
    return float(efficiency), float(max_directivity_dbi)


@dataclass(frozen=True)
class VortexMeasurement:
    null_depth_db: float  # field power AT the beam direction, relative to the pattern's own peak
    winding_number: float  # net phase winding (cycles) around the beam direction, in the co-pol component
    ring_alpha_deg: np.ndarray  # the sample ring's own parameter (0-360deg), for plotting
    ring_phase_deg: np.ndarray  # unwrapped phase (deg) of the co-pol component around the ring


def measure_vortex_winding(
    x_m, y_m, Jx, Jy, frequency_hz, green_function, dx_m, dy_m,
    theta_b_deg, phi_b_deg, desired_polarization, pattern_peak_power=None,
    ring_radius=0.03, n_ring_samples=72,
):
    """Makes "vortex vs pencil" a MEASUREMENT (the brief's own requirement), not an assumption:
    samples the beam's own co-pol component (`copol_component`) on a small ring around
    `(theta_b_deg, phi_deg)` in (u, v) = (sin(theta)cos(phi), sin(theta)sin(phi)) space (flat and
    non-singular even at boresight, unlike a theta/phi ring), unwraps its phase around the ring, and
    reports the net winding in CYCLES (an integer for a true vortex; near 0 for a pencil beam) plus
    the field's own power AT the beam direction relative to the whole pattern's peak (`pattern_peak_power`,
    from a prior `compute_vector_pattern` call on the same design -- if omitted, relative to this
    ring's own peak instead, which still correctly reports a deep null as very negative).

    `ring_radius` [dimensionless, in (u,v) space] default 0.03 (~1.7deg at boresight) is small
    enough to stay within a single main lobe for the grids this module's own tests and the
    brief's Design-1 reproduction use; a caller with a much broader or narrower beam should pass a
    radius scaled to its own beamwidth."""
    u_b = np.sin(np.radians(theta_b_deg)) * np.cos(np.radians(phi_b_deg))
    v_b = np.sin(np.radians(theta_b_deg)) * np.sin(np.radians(phi_b_deg))
    alpha = np.linspace(0.0, 2.0 * np.pi, n_ring_samples, endpoint=False)
    u = u_b + ring_radius * np.cos(alpha)
    v = v_b + ring_radius * np.sin(alpha)
    rho = np.hypot(u, v)
    if np.any(rho > 1.0):
        raise ValueError("measure_vortex_winding: ring_radius takes the sample ring past the visible hemisphere")
    theta_ring_deg = np.degrees(np.arcsin(rho))
    phi_ring_deg = np.degrees(np.arctan2(v, u))

    e_theta, e_phi = _vector_far_field(
        theta_ring_deg, phi_ring_deg, x_m, y_m, Jx, Jy, frequency_hz, green_function, dx_m, dy_m
    )
    field_ring = copol_component(e_theta, e_phi, desired_polarization)
    # standard discrete winding number: sum of wrapped consecutive phase differences around the
    # CLOSED ring (including the last-sample-back-to-first-sample step), / 2*pi -- exact for a
    # smoothly-varying phase sampled finely enough that no true turn exceeds pi between samples,
    # and (unlike np.unwrap) makes no assumption about where the ring "starts".
    phases = np.angle(field_ring)
    diffs = np.diff(np.concatenate([phases, phases[:1]]))
    diffs_wrapped = (diffs + np.pi) % (2.0 * np.pi) - np.pi
    winding = np.sum(diffs_wrapped) / (2.0 * np.pi)

    e_theta_b, e_phi_b = _vector_far_field(
        theta_b_deg, phi_b_deg, x_m, y_m, Jx, Jy, frequency_hz, green_function, dx_m, dy_m
    )
    field_b = copol_component(np.array([e_theta_b]), np.array([e_phi_b]), desired_polarization)[0]
    peak = pattern_peak_power if pattern_peak_power is not None else np.max(np.abs(field_ring) ** 2)
    null_depth_db = 10.0 * np.log10(np.abs(field_b) ** 2 / peak) if peak > 0.0 else float("-inf")
    return VortexMeasurement(
        null_depth_db=null_depth_db, winding_number=float(winding),
        ring_alpha_deg=np.degrees(alpha), ring_phase_deg=np.degrees(np.unwrap(phases)),
    )


def first_null_half_width_deg(aperture_extent_m, wavelength_m, theta0_deg=0.0):
    """Approximate first-null half-width of a uniformly illuminated aperture of extent `D`,
    `asin(1.22*lambda/D)` (Balanis), widened by `1/cos(theta0)` for a steered beam. Used only to
    exclude the main lobe's own shoulder from the search for a strongest SECONDARY beam."""
    ratio = min(1.0, 1.22 * wavelength_m / aperture_extent_m)
    return float(np.degrees(np.arcsin(ratio)) / max(np.cos(np.radians(theta0_deg)), 1e-3))


def _angular_distance_deg(theta1_deg, phi1_deg, theta2_deg, phi2_deg):
    t1, t2 = np.radians(theta1_deg), np.radians(theta2_deg)
    dphi = np.radians(phi1_deg - phi2_deg)
    cos_d = np.clip(np.cos(t1) * np.cos(t2) + np.sin(t1) * np.sin(t2) * np.cos(dphi), -1.0, 1.0)
    return float(np.degrees(np.arccos(cos_d)))


def strongest_secondary_beam(result: VectorPatternResult, exclusion_deg):
    """The strongest local maximum of `result`'s power lying OUTSIDE `exclusion_deg` of the main
    beam (angular distance on the sphere), as `(theta_deg, phi_deg, level_db_rel_main)`, or `None`
    if no such maximum exists above the -40 dB search floor. Reuses `aperture_pattern`'s
    `detect_local_maxima` unmodified, the same peak finder `compute_vector_pattern` uses."""
    power = np.abs(result.e_theta) ** 2 + np.abs(result.e_phi) ** 2
    peaks = detect_local_maxima(power, result.theta_deg, result.phi_deg, min_relative_db=-40.0, phi_periodic=True)
    outside = [
        p for p in peaks
        if _angular_distance_deg(p.theta_deg, p.phi_deg, result.main_beam_theta_deg, result.main_beam_phi_deg)
        > exclusion_deg
    ]
    if not outside:
        return None
    best = max(outside, key=lambda p: p.relative_power_db)
    return (float(best.theta_deg), float(best.phi_deg), float(best.relative_power_db))


def mirror_beam_level_db(result: VectorPatternResult):
    """Level (dB, relative to the main beam) at the MIRROR direction `(theta0, phi0+180deg)` on the
    pattern grid -- the direction a symmetric 1-bit pair would also radiate into if the quantization
    were sign-symmetric. Nearest grid sample; `phi` is taken modulo 360."""
    power = np.abs(result.e_theta) ** 2 + np.abs(result.e_phi) ** 2
    peak = float(power.max())
    theta_axis = result.theta_deg[:, 0]
    phi_axis = result.phi_deg[0, :]
    i = int(np.argmin(np.abs(theta_axis - result.main_beam_theta_deg)))
    mirror_phi = (result.main_beam_phi_deg + 180.0) % 360.0
    dphi = np.abs((phi_axis - mirror_phi + 180.0) % 360.0 - 180.0)
    j = int(np.argmin(dphi))
    return float(10.0 * np.log10(max(power[i, j], 1e-300) / peak)), float(theta_axis[i]), float(phi_axis[j])


ORTHOGONAL_POLARIZATION = {
    "rhcp": "lhcp", "lhcp": "rhcp", "hp": "vp", "vp": "hp", "slant_p45": "slant_m45", "slant_m45": "slant_p45",
}
MAIN_BEAM_REGION_DB = -10.0  # the co-pol contour that bounds the main-beam region
_SEARCH_RADIUS_DEG = 12.0


@dataclass(frozen=True)
class BeamPolarizationMetrics:
    point_xpd_db: "float | None"  # co over cross at the grid point nearest the request
    integrated_xpd_db: "float | None"  # co over cross power summed over the main-beam region
    polarization_efficiency: "float | None"  # co-pol power / total power over the hemisphere
    sll_copol_db: "float | None"  # highest co-pol level outside the main-beam region, rel. its peak
    cross_peak_db: "float | None"  # highest cross-pol level over the hemisphere, rel. the co-pol peak
    region_fraction: float  # fraction of the hemisphere grid inside the main-beam region


def _direction_mask(theta_grid_deg, phi_grid_deg, theta_deg, phi_deg, radius_deg):
    t0, p0 = np.radians(theta_deg), np.radians(phi_deg)
    t, p = np.radians(theta_grid_deg), np.radians(phi_grid_deg)
    cos_angle = np.cos(t) * np.cos(t0) + np.sin(t) * np.sin(t0) * np.cos(p - p0)
    return cos_angle >= np.cos(np.radians(radius_deg))


def _main_beam_region(co_power, theta_grid_deg, phi_grid_deg, theta_deg, phi_deg):
    """The connected -10 dB region of the co-pol power around the beam: seeded at the co-pol
    maximum within `_SEARCH_RADIUS_DEG` of the request, grown on the (theta, phi) grid with phi
    periodic (wrapped by one column on each side so a region across phi=0/360 stays whole)."""
    from scipy import ndimage

    seed_window = _direction_mask(theta_grid_deg, phi_grid_deg, theta_deg, phi_deg, _SEARCH_RADIUS_DEG)
    masked = np.where(seed_window, co_power, -np.inf)
    si, sj = np.unravel_index(int(np.argmax(masked)), masked.shape)
    peak = co_power[si, sj]
    above = co_power >= peak * 10.0 ** (MAIN_BEAM_REGION_DB / 10.0)
    padded = np.concatenate([above[:, -1:], above, above[:, :1]], axis=1)
    labels, _ = ndimage.label(padded, structure=np.ones((3, 3), dtype=int))
    region = labels == labels[si, sj + 1]
    return region[:, 1:-1], peak


def beam_polarization_metrics(pat: VectorPatternResult, desired_polarization, theta_deg, phi_deg):
    """Polarization figures of merit for one requested beam, from a computed vector pattern.

    `desired_polarization` is the beam's own label ("rhcp"/"lhcp"/"hp"/"vp"/"slant_p45"/
    "slant_m45"); its orthogonal partner is the cross-pol. Every power is over the grid's own
    hemisphere with the sin(theta) solid-angle weight. A scalar cell cannot fix both components of a
    monopole's radial current, so the cross-pol share is reported as measured, not assumed."""
    if desired_polarization not in ORTHOGONAL_POLARIZATION:
        raise ValueError(f"beam_polarization_metrics needs a copol label, got {desired_polarization!r}")
    cross_pol = ORTHOGONAL_POLARIZATION[desired_polarization]
    e_theta, e_phi = pat.e_theta, pat.e_phi
    co = np.abs(copol_component(e_theta, e_phi, desired_polarization)) ** 2
    cross = np.abs(copol_component(e_theta, e_phi, cross_pol)) ** 2
    total = np.abs(e_theta) ** 2 + np.abs(e_phi) ** 2

    theta_rad = np.radians(pat.theta_deg)
    dtheta = np.radians(pat.theta_deg[1, 0] - pat.theta_deg[0, 0])
    dphi = 2.0 * np.pi / pat.theta_deg.shape[1]
    weight = np.sin(theta_rad) * dtheta * dphi

    region, co_peak = _main_beam_region(co, pat.theta_deg, pat.phi_deg, theta_deg, phi_deg)
    outside = ~region
    i = int(np.argmin(np.abs(pat.theta_deg[:, 0] - theta_deg)))
    j = int(np.argmin(np.abs(((pat.phi_deg[0, :] - phi_deg + 180.0) % 360.0) - 180.0)))
    point_cross = cross[i, j]
    point_xpd = float(10.0 * np.log10(co[i, j] / point_cross)) if point_cross > 0.0 else None

    region_co = float(np.sum(co * weight * region))
    region_cross = float(np.sum(cross * weight * region))
    integrated_xpd = float(10.0 * np.log10(region_co / region_cross)) if region_cross > 0.0 else None

    total_power = float(np.sum(total * weight))
    efficiency = float(np.sum(co * weight) / total_power) if total_power > 0.0 else None

    outside_co = float(np.max(co[outside])) if np.any(outside) else 0.0
    sll = float(10.0 * np.log10(outside_co / co_peak)) if outside_co > 0.0 else None
    cross_peak = float(10.0 * np.log10(np.max(cross) / co_peak)) if co_peak > 0.0 else None

    return BeamPolarizationMetrics(
        point_xpd_db=point_xpd,
        integrated_xpd_db=integrated_xpd,
        polarization_efficiency=efficiency,
        sll_copol_db=sll,
        cross_peak_db=cross_peak,
        region_fraction=float(np.mean(region)),
    )
