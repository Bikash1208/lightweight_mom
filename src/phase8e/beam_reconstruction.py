"""Sec. 18-19: continuous and discrete aperture -> far-field reconstruction, BEFORE any geometry
inversion. The aperture-to-far-field relation used here is the standard planar-aperture radiation
integral (e.g. Balanis, 'Antenna Theory,' aperture-antenna far-field pattern as the 2-D spatial
Fourier transform of the aperture distribution) -- a cited, standard antenna-theory relation, not
an invented one; it is evaluated by direct summation (not FFT) so the requested (theta, phi) grid
need not be regular in k-space.
"""

from dataclasses import dataclass

import numpy as np

from ..frequency import free_space_wavenumber


def far_field_pattern(x_m, y_m, aperture_field, frequency_hz, theta_deg, phi_deg, dx_m, dy_m):
    """F(theta, phi) = sum_xy A(x,y) * exp(+j*(kx*x + ky*y)) * dx*dy -- the planar-aperture
    far-field radiation integral, direction cosines kx=k0 sin(theta)cos(phi), ky=k0 sin(theta)
    sin(phi) (this project's one transverse-wavenumber formula, `phase8b.object_wave.
    transverse_wavenumber_deg`'s underlying algebra, reused via direct broadcast here since a
    scalar per-angle evaluation is simpler than reconstructing Beam objects for a whole angle
    grid). Returns a complex array shaped like ``theta_deg``."""
    k0 = free_space_wavenumber(frequency_hz)
    theta = np.radians(np.asarray(theta_deg, dtype=float))
    phi = np.radians(np.asarray(phi_deg, dtype=float))
    kx = k0 * np.sin(theta) * np.cos(phi)
    ky = k0 * np.sin(theta) * np.sin(phi)

    x_flat = np.asarray(x_m, dtype=float).ravel()
    y_flat = np.asarray(y_m, dtype=float).ravel()
    a_flat = np.asarray(aperture_field, dtype=complex).ravel()

    shape = kx.shape
    phase = np.outer(kx.ravel(), x_flat) + np.outer(ky.ravel(), y_flat)   # (Nf, Na), kx_i*x_j + ky_i*y_j
    field = (np.exp(1j * phase) @ a_flat) * dx_m * dy_m
    return field.reshape(shape)


@dataclass(frozen=True)
class DetectedPeak:
    theta_deg: float
    phi_deg: float
    relative_power_db: float


def detect_local_maxima(power, theta_grid_deg, phi_grid_deg, min_relative_db=-20.0):
    """Simple grid-based local-maximum detector (greater than all 4 axis-neighbours), above
    ``min_relative_db`` of the global peak -- sufficient to separate a handful of intended beams
    from sidelobes on a modest angular grid; not a general antenna-pattern peak-search algorithm."""
    power = np.asarray(power, dtype=float)
    peak_power = power.max()
    peaks = []
    ny, nx = power.shape
    for i in range(ny):
        for j in range(nx):
            v = power[i, j]
            if v < peak_power * 10 ** (min_relative_db / 10.0):
                continue
            neighbors = [
                power[i - 1, j] if i > 0 else -np.inf, power[i + 1, j] if i + 1 < ny else -np.inf,
                power[i, j - 1] if j > 0 else -np.inf, power[i, j + 1] if j + 1 < nx else -np.inf,
            ]
            if v >= max(neighbors):
                peaks.append(DetectedPeak(
                    float(theta_grid_deg[i, j]), float(phi_grid_deg[i, j]),
                    float(10.0 * np.log10(v / peak_power)),
                ))
    return sorted(peaks, key=lambda p: -p.relative_power_db)


@dataclass(frozen=True)
class BeamReconstructionMetrics:
    n_intended_beams: int
    pointing_error_deg: tuple            # per intended beam, angular distance to nearest detected peak
    matched_peak_power_db: tuple         # per intended beam, that peak's relative power [dB]
    unwanted_peaks: tuple                # DetectedPeak objects not matched to any intended beam
    highest_unwanted_relative_db: float | None
    detected_peaks: tuple


def _angular_distance_deg(theta1, phi1, theta2, phi2):
    """Great-circle-style angular distance between two (theta, phi) directions via their
    Cartesian unit vectors -- exact for any theta/phi, unlike a naive coordinate difference."""
    def unit(theta, phi):
        t, p = np.radians(theta), np.radians(phi)
        return np.array([np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)])

    u1, u2 = unit(theta1, phi1), unit(theta2, phi2)
    cos_angle = np.clip(np.dot(u1, u2), -1.0, 1.0)
    return float(np.degrees(np.arccos(cos_angle)))


def evaluate_beam_reconstruction(far_field, theta_grid_deg, phi_grid_deg, intended_beams, min_relative_db=-20.0):
    """Sec. 18: does the (continuous or discrete) far field actually produce the requested
    beam(s)? Reports per-beam pointing error and matched-peak power, plus any unmatched ("unwanted
    dominant") peak, never averaged into one hidden number (Sec. 19's explicit prohibition)."""
    power = np.abs(far_field) ** 2
    peaks = detect_local_maxima(power, theta_grid_deg, phi_grid_deg, min_relative_db)
    matched_indices = set()
    pointing_errors, matched_db = [], []
    for beam in intended_beams:
        if not peaks:
            pointing_errors.append(float("nan"))
            matched_db.append(float("nan"))
            continue
        distances = [_angular_distance_deg(beam.theta_deg, beam.phi_deg, p.theta_deg, p.phi_deg) for p in peaks]
        best = int(np.argmin(distances))
        matched_indices.add(best)
        pointing_errors.append(distances[best])
        matched_db.append(peaks[best].relative_power_db)
    unwanted = [p for i, p in enumerate(peaks) if i not in matched_indices]
    return BeamReconstructionMetrics(
        n_intended_beams=len(intended_beams), pointing_error_deg=tuple(pointing_errors),
        matched_peak_power_db=tuple(matched_db), unwanted_peaks=tuple(unwanted),
        highest_unwanted_relative_db=(max(p.relative_power_db for p in unwanted) if unwanted else None),
        detected_peaks=tuple(peaks),
    )


@dataclass(frozen=True)
class ReconstructionComparison:
    complex_error_mean: float
    complex_error_max: float
    amplitude_error_mean: float
    phase_error_deg_mean: float


def compare_continuous_vs_discrete(continuous_far_field, discrete_far_field):
    """Sec. 19: complex/amplitude/phase error between the continuous-formulation far field and the
    discrete-library-reconstructed far field, over the whole angular grid -- reported as summary
    statistics here; per-angle detail belongs in the plots (Sec. 23), not hidden behind one number
    (Sec. 19's explicit prohibition on hiding failures behind an averaged phase error alone --
    amplitude and complex error are reported alongside it)."""
    c = np.asarray(continuous_far_field, dtype=complex)
    d = np.asarray(discrete_far_field, dtype=complex)
    complex_err = np.abs(c - d)
    amp_err = np.abs(np.abs(c) - np.abs(d))
    with np.errstate(invalid="ignore"):
        phase_err = np.abs(np.angle(c) - np.angle(d))
        phase_err = np.minimum(phase_err, 2 * np.pi - phase_err)
    return ReconstructionComparison(
        complex_error_mean=float(np.mean(complex_err)), complex_error_max=float(np.max(complex_err)),
        amplitude_error_mean=float(np.mean(amp_err)),
        phase_error_deg_mean=float(np.degrees(np.nanmean(phase_err))),
    )
