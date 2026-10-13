"""Phase 15 Part E: the Coded (1-bit / N-bit) design method (Mishra, Barad, Sarkar, IEEE TVT 2026).

## The phase formula is ALREADY this project's own phi_ext/phase-difference machinery

Mishra's eq. 1, `psi_pq = -k0(x sin(th0)cos(ph0) + y sin(th0)sin(ph0)) + beta_sw*rho + chi*atan2(y,x)`
(`chi = +-1` for LH/RH CP -- the same spiral compensation as Part A; `chi = 0` with a half-plane
flip for LP), decomposes EXACTLY into quantities this phase already computes:

- the steering term is `arg(phase8b.object_wave.single_beam_object_field(...))` (identical formula,
  `-k0*(kx0*x+ky0*y)` via `transverse_wavenumber_deg`);
- `beta_sw*rho` is `-arg(E_ref)` for the standard monopole reference wave
  (`phase8b.reference_wave.monopole_reference_field`'s own `exp(-j*k_swg*rho)`), so adding it is
  exactly `required_phase = arg(E_obj) - arg(E_ref)` -- Formulation B's own phase-difference
  hologram (`hmsa_synthesis.formulation_b_control_map`);
- the CP spiral / LP half-plane-flip term is EXACTLY `arg(hologram.polarization_compensation.
  polarization_compensation_term(..., oam_charge=0))` (Part A's own phi_ext, with no extra OAM).

So `formulation_b_control_map(reference_field, [beam], frequency_hz, feed_xy_m, cp_phase_offset_deg)`
called on a SINGLE beam already returns `psi_m(x,y)` in degrees -- the brief's own "unify the code
so both methods use one phi_ext implementation" falls out for free; `per_beam_phase_maps_deg` below
just calls it once per beam instead of combining beams first.

## The combination rule IS a special case of one formula

Mishra eq. 3 (multi-beam): `f = sum_m w_m*cos(psi_m)`; eq. 2 (1-bit): `B = 1 if f>0 else 0`; the
brief's own N-bit generalization: quantize `arg(sum_m w_m*exp(j*psi_m))` into `2**n_bits` levels.
These are the SAME rule: `f = sum_m w_m*cos(psi_m) = Re(sum_m w_m*exp(j*psi_m))`, and quantizing an
angle into 2 levels {0, 180} by nearest-level is exactly `sign(Re(.))` -- `f > 0` decides this
exact way. `coded_phase_map_deg` therefore implements ONE quantization formula for any `n_bits`,
not two.
"""

from dataclasses import dataclass, replace

import numpy as np

from ..phase import principal_phase_deg, wrapped_difference_deg
from ..phase8e.reference_field_input import ReferenceFieldGrid
from .hmsa_synthesis import GeometryCell, map_control_to_geometry_library
from .hmsa_synthesis import formulation_b_control_map as _psi_m_control_map


def per_beam_phase_maps_deg(reference_field, object_beams, frequency_hz, feed_xy_m=(0.0, 0.0), cp_phase_offset_deg=0.0):
    """psi_m(x,y) [deg] for each beam, independently (module docstring) -- `feed_xy_m` may be one
    shared `(x0, y0)` or a per-beam sequence, exactly as `polarization_compensation_term` already
    accepts (Part C: each beam radiating about its own assigned feed)."""
    from .polarization_compensation import _feed_xy_for_beam

    maps = []
    for i, beam in enumerate(object_beams):
        feed_xy_b = _feed_xy_for_beam(feed_xy_m, i)
        maps.append(_psi_m_control_map(reference_field, [beam], frequency_hz, feed_xy_b, cp_phase_offset_deg))
    return maps


def coded_phase_map_deg(
    reference_field, object_beams, frequency_hz, n_bits=1, feed_xy_m=(0.0, 0.0), cp_phase_offset_deg=0.0, weights=None
):
    """Returns `(continuous_deg, quantized_deg, level_map)`: `continuous_deg` is
    `arg(sum_m w_m*exp(j*psi_m))`, the un-quantized reference used to measure quantization loss;
    `quantized_deg`/`level_map` snap it to the nearest of `2**n_bits` evenly-spaced levels starting
    at 0deg. `weights`: real per-beam weight `w_m` (default `abs(Beam.weight)`, matching the
    hologram path's own power-weight convention) -- Mishra's eq. 3 weights are real, not complex;
    any phase of a complex `Beam.weight` is DISCARDED here, not silently honoured."""
    if n_bits < 1:
        raise ValueError(f"n_bits must be >= 1, got {n_bits!r}")
    psi_maps_deg = per_beam_phase_maps_deg(reference_field, object_beams, frequency_hz, feed_xy_m, cp_phase_offset_deg)
    if weights is None:
        weights = [abs(b.weight) for b in object_beams]
    complex_sum = sum(w * np.exp(1j * np.radians(p)) for w, p in zip(weights, psi_maps_deg, strict=True))
    continuous_deg = np.degrees(np.angle(complex_sum))

    n_states = 2**n_bits
    step_deg = 360.0 / n_states
    level_map = (np.round((continuous_deg % 360.0) / step_deg).astype(int)) % n_states
    quantized_deg = level_map * step_deg
    return continuous_deg, quantized_deg, level_map


@dataclass(frozen=True)
class CodedDesign:
    x_m: "object"
    y_m: "object"
    n_bits: int
    continuous_phase_deg: "object"  # un-quantized reference (quantization-loss comparisons)
    quantized_phase_deg: "object"  # IDEAL target phase per cell (0, step, 2*step, ...)
    level_map: "object"  # int array, which of the 2**n_bits states each cell uses
    cells: "tuple[GeometryCell, ...]"
    aperture_field: "object"  # E_ref * exp(j*ACHIEVED phase) -- the real design's radiated field
    labels: tuple
    two_state_states: "tuple | None"  # PHASE15 Part E: the 2 TwoStateCell`s used, n_bits==1 only
    unreachable_fraction: float = 0.0  # every cell takes one of the two real states, so none is unreachable
    path_iterations: int = 0  # state-dependent beta iterations used (0 when beta was not path-integrated)
    path_converged: bool = True
    keepout_cells_changed: int = 0  # cells forced to state 0 under the feed (apply_feed_keepout_one_bit)


MAX_PATH_ITERATIONS = 10


def synthesize_coded_one_bit(
    reference_field,
    object_beams,
    frequency_hz,
    two_state_choice,
    period_m,
    feed_xy_m=(0.0, 0.0),
    cp_phase_offset_deg=0.0,
    weights=None,
    surface_wave_betas_rad_per_m=None,
):
    """n_bits=1: the quantized {0, 1} map maps DIRECTLY onto `two_state_choice`'s own two physical
    states (no further library search -- these ARE the two available states, by definition of
    1-bit coding). `GeometryCell.geometry_param` holds the chosen STATE INDEX (0.0 or 1.0, a float
    so `GeometryCell`'s existing `"float | None"` type is unchanged) -- the caller maps a state
    index back to `two_state_choice.state_0/state_1` (each carrying its own `.geometry` dict,
    rectangle or annular-slot) for export; `achieved_value` is that state's own MEASURED reflection
    phase (not the ideal 0/180deg target), `phase_error_deg` the gap between them.

    `surface_wave_betas_rad_per_m`: `(beta_state_0, beta_state_1)`. The surface-wave phase is then
    path-integrated along each radial ray from the feed, `integral beta(s) ds`, where `beta(s)` is
    the beta of the state chosen at that point -- the two states' betas differ, so the phase is
    not `beta_mean * rho`. The state map and the phase are solved together by iteration (each
    pass re-quantizes with the phase from the previous pass's states)."""
    from ..fields.reference_field_grid import resample_to_cell_grid

    cell_field = resample_to_cell_grid(reference_field, period_m)
    states = (two_state_choice.state_0, two_state_choice.state_1)
    path_iterations = 0
    path_converged = True
    if surface_wave_betas_rad_per_m is None:
        continuous_deg, quantized_deg, level_map = coded_phase_map_deg(
            cell_field,
            object_beams,
            frequency_hz,
            n_bits=1,
            feed_xy_m=feed_xy_m,
            cp_phase_offset_deg=cp_phase_offset_deg,
            weights=weights,
        )
    else:
        beta_states = tuple(float(b) for b in surface_wave_betas_rad_per_m)
        beta_map = np.full(cell_field.x_m.shape, 0.5 * (beta_states[0] + beta_states[1]))
        path_converged = False
        while path_iterations < MAX_PATH_ITERATIONS:
            path_iterations += 1
            phase_rad = path_integrated_phase_rad(cell_field.x_m, cell_field.y_m, feed_xy_m, beta_map, period_m / 2.0)
            psi_reference = _phase_reference_field(
                cell_field, phase_rad, f"PATH_INTEGRATED_BETA({beta_states[0]:.6g},{beta_states[1]:.6g}rad/m)"
            )
            continuous_deg, quantized_deg, level_map = coded_phase_map_deg(
                psi_reference,
                object_beams,
                frequency_hz,
                n_bits=1,
                feed_xy_m=feed_xy_m,
                cp_phase_offset_deg=cp_phase_offset_deg,
                weights=weights,
            )
            beta_next = np.where(level_map == 0, beta_states[0], beta_states[1])
            if np.array_equal(beta_next, beta_map):
                path_converged = True
                break
            beta_map = beta_next
    ny, nx = cell_field.x_m.shape
    cells = []
    achieved_deg = np.empty((ny, nx))
    for i in range(ny):
        for j in range(nx):
            level = int(level_map[i, j])
            state = states[level]
            achieved = float(principal_phase_deg(state.value))
            achieved_deg[i, j] = achieved
            target = float(quantized_deg[i, j])
            cells.append(
                GeometryCell(
                    i,
                    j,
                    float(cell_field.x_m[i, j]),
                    float(cell_field.y_m[i, j]),
                    control_value=target,
                    geometry_param=float(level),
                    achieved_value=achieved,
                    phase_error_deg=float(abs(wrapped_difference_deg(achieved, target))),
                    reachable=True,
                )
            )
    ref_complex = cell_field.Ex if cell_field.amplitude_available else np.exp(1j * np.radians(cell_field.phase()))
    aperture_field = ref_complex * np.exp(1j * np.radians(achieved_deg))
    labels = ("CODED_DESIGN_METHOD",)
    return CodedDesign(
        x_m=cell_field.x_m,
        y_m=cell_field.y_m,
        n_bits=1,
        continuous_phase_deg=continuous_deg,
        quantized_phase_deg=quantized_deg,
        level_map=level_map,
        cells=tuple(cells),
        aperture_field=aperture_field,
        labels=labels,
        two_state_states=states,
        path_iterations=path_iterations,
        path_converged=path_converged,
    )


FEED_KEEPOUT_FRACTION = 0.6  # a cell is "under the feed" when its centre is within 0.6 periods of it
# in x AND y: a feed at a cell centre touches 1 cell, on an edge 2 cells, on a corner 4 cells


def apply_feed_keepout_one_bit(design, feed_xy_m, period_m):
    """Phase 21 C2: force the cell(s) under the feed to state 0 (the empty cell -- a monopole's
    coax passes through the board there, so no patch can sit on it). Returns `(design, n_changed)`:
    `n_changed` counts only cells that were state 1 before. The radiated `aperture_field` of a
    changed cell is rebuilt with state 0's own measured phase (its reference-wave amplitude and
    phase are kept), so the pattern prediction sees the real, keep-out-corrected board."""
    if design.n_bits != 1 or design.two_state_states is None:
        raise ValueError("apply_feed_keepout_one_bit needs a 1-bit coded design with two_state_states")
    fx, fy = feed_xy_m
    limit = FEED_KEEPOUT_FRACTION * period_m
    state_0 = design.two_state_states[0]
    achieved_0 = float(principal_phase_deg(state_0.value))
    cells = list(design.cells)
    level_map = np.array(design.level_map, copy=True)
    aperture_field = np.array(design.aperture_field, copy=True)
    n_changed = 0
    for k, c in enumerate(cells):
        if abs(c.x_m - fx) >= limit or abs(c.y_m - fy) >= limit or int(c.geometry_param) == 0:
            continue
        n_changed += 1
        aperture_field[c.i, c.j] *= np.exp(1j * np.radians(achieved_0 - c.achieved_value))
        level_map[c.i, c.j] = 0
        cells[k] = replace(
            c,
            geometry_param=0.0,
            achieved_value=achieved_0,
            phase_error_deg=float(abs(wrapped_difference_deg(achieved_0, c.control_value))),
        )
    return (
        replace(
            design,
            cells=tuple(cells),
            level_map=level_map,
            aperture_field=aperture_field,
            keepout_cells_changed=n_changed,
        ),
        n_changed,
    )


def path_integrated_phase_rad(x_m, y_m, feed_xy_m, beta_map_rad_per_m, step_m):
    """`phi(x,y) = integral_0^rho beta(s) ds` along the straight ray from the feed to each cell
    centre, the integrand taken from the cell nearest each sample point (samples every `step_m`,
    exact segment lengths). A uniform beta gives exactly `beta*rho`."""
    x0, y0 = feed_xy_m
    xs, ys = x_m[0, :], y_m[:, 0]
    dx = float(xs[1] - xs[0]) if len(xs) > 1 else step_m
    dy = float(ys[1] - ys[0]) if len(ys) > 1 else step_m
    ny, nx = x_m.shape
    phase = np.zeros((ny, nx))
    for i in range(ny):
        for j in range(nx):
            vx, vy = x_m[i, j] - x0, y_m[i, j] - y0
            rho = float(np.hypot(vx, vy))
            if rho == 0.0:
                continue
            n = max(1, int(np.ceil(rho / step_m)))
            edges = np.linspace(0.0, rho, n + 1)
            mids = 0.5 * (edges[:-1] + edges[1:])
            sx = x0 + mids * vx / rho
            sy = y0 + mids * vy / rho
            ii = np.clip(np.round((sy - ys[0]) / dy).astype(int), 0, ny - 1)
            jj = np.clip(np.round((sx - xs[0]) / dx).astype(int), 0, nx - 1)
            phase[i, j] = float(np.sum(beta_map_rad_per_m[ii, jj] * np.diff(edges)))
    return phase


def _phase_reference_field(cell_field, phase_rad, source):
    field = np.exp(-1j * phase_rad)
    return ReferenceFieldGrid(
        x_m=cell_field.x_m,
        y_m=cell_field.y_m,
        Ex=field,
        Ey=None,
        amplitude_available=True,
        source=source,
        dx_m=cell_field.dx_m,
        dy_m=cell_field.dy_m,
        z_reference_m=0.0,
        phase_deg=np.degrees(np.angle(field)),
    )


def synthesize_coded_n_bit(
    reference_field,
    object_beams,
    frequency_hz,
    n_bits,
    library_candidates,
    period_m,
    feed_xy_m=(0.0, 0.0),
    cp_phase_offset_deg=0.0,
    weights=None,
):
    """n_bits > 1: the quantized target phase (`2**n_bits` discrete levels) is matched against a
    CONTINUOUS phase library (e.g. the square-patch one) via the SAME
    `hmsa_synthesis.map_control_to_geometry_library`/`nearest_response_match` the continuous
    Formulation B path already uses -- an N-bit coded design is reachability-reported exactly like
    a continuous one (`GeometryCell.geometry_param`/`geometry_param_w` are the matched L_mm/W_mm,
    as usual)."""
    from ..fields.reference_field_grid import resample_to_cell_grid

    cell_field = resample_to_cell_grid(reference_field, period_m)
    continuous_deg, quantized_deg, level_map = coded_phase_map_deg(
        cell_field,
        object_beams,
        frequency_hz,
        n_bits=n_bits,
        feed_xy_m=feed_xy_m,
        cp_phase_offset_deg=cp_phase_offset_deg,
        weights=weights,
    )
    cells = map_control_to_geometry_library(cell_field.x_m, cell_field.y_m, quantized_deg, library_candidates)
    achieved_deg = np.array([c.achieved_value if c.achieved_value is not None else 0.0 for c in cells]).reshape(
        cell_field.x_m.shape
    )
    ref_complex = cell_field.Ex if cell_field.amplitude_available else np.exp(1j * np.radians(cell_field.phase()))
    aperture_field = ref_complex * np.exp(1j * np.radians(achieved_deg))
    labels = ("CODED_DESIGN_METHOD",)
    return CodedDesign(
        x_m=cell_field.x_m,
        y_m=cell_field.y_m,
        n_bits=n_bits,
        continuous_phase_deg=continuous_deg,
        quantized_phase_deg=quantized_deg,
        level_map=level_map,
        cells=cells,
        aperture_field=aperture_field,
        labels=labels,
        two_state_states=None,
    )


def beta_reference_field(cell_field, feed_xy_m, beta_rad_per_m):
    """Uniform-beta reference `exp(-j*beta*rho)` about the feed: `path_integrated_phase_rad` with
    a single beta, kept for the one-beta case (a state pair whose betas are the same)."""
    x0, y0 = feed_xy_m
    rho = np.hypot(cell_field.x_m - x0, cell_field.y_m - y0)
    return _phase_reference_field(cell_field, beta_rad_per_m * rho, f"SURFACE_WAVE_BETA({beta_rad_per_m:.4g}rad/m)")


def textbook_quantization_loss_db(n_bits):
    """The standard main-beam loss of N-bit phase quantization, `-10*log10(sinc^2(pi/2^N))`
    (uniformly distributed continuous phase, nearest-level rounding): 3.92 dB for 1 bit, 0.91 dB
    for 2 bits. A sanity reference for `quantization_loss_db`, not a measured value."""
    x = np.pi / (2**n_bits)
    return float(-10.0 * np.log10((np.sin(x) / x) ** 2))


def aperture_coupling_loss_db(ideal_aperture_field, continuous_aperture_field):
    """The coherent (main-beam) fraction of the continuous aperture the quantized one retains:
    `-10*log10(|<E_q, E_c>|^2 / (|E_q|^2 |E_c|^2))`, with equal aperture amplitudes this is exactly
    the textbook quantity above. Independent of any pattern grid or directivity normalisation."""
    e_q = np.asarray(ideal_aperture_field, dtype=complex).ravel()
    e_c = np.asarray(continuous_aperture_field, dtype=complex).ravel()
    num = abs(np.vdot(e_c, e_q)) ** 2
    den = np.vdot(e_c, e_c).real * np.vdot(e_q, e_q).real
    if den <= 0.0 or num <= 0.0:
        return float("inf")
    return float(-10.0 * np.log10(num / den))


def continuous_reference_aperture_field(reference_field, continuous_phase_deg):
    """The hypothetical zero-quantization-error aperture field (`E_ref * exp(j*continuous_phase)`)
    -- compared against a `CodedDesign.aperture_field` (which uses the REAL quantized+achieved
    phase) to measure quantization loss (Part E's own "report quantization loss vs continuous")."""
    ref_complex = (
        reference_field.Ex if reference_field.amplitude_available else np.exp(1j * np.radians(reference_field.phase()))
    )
    return ref_complex * np.exp(1j * np.radians(continuous_phase_deg))
