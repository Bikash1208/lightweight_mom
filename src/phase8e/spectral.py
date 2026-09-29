"""Sec. 7-8: 2-D spatial-Fourier-transform analysis of the reference field and the object wave.
Nothing here discards evanescent content -- it is reported as its own fraction of the spectral
power, never silently dropped (Sec. 7's explicit instruction).
"""

from dataclasses import dataclass

import numpy as np

from ..frequency import free_space_wavenumber
from ..periodic_bc import normal_wavenumber
from ..phase8b.object_wave import transverse_wavenumber_deg
from ..phase8d_rev.candidate_period import required_transverse_spatial_frequency_per_m


@dataclass(frozen=True)
class SpectralAnalysis:
    kx_max_per_m: float           # Nyquist-limited spectral grid extent from the sample spacing
    ky_max_per_m: float
    k0_per_m: float
    dominant_kx_per_m: float      # transverse wavenumber of the peak spectral-power bin
    dominant_ky_per_m: float
    k_t_dominant_per_m: float
    propagating_fraction: float   # power-weighted fraction of spectrum with kx^2+ky^2 <= k0^2
    evanescent_fraction: float
    local_periodic_compatible: bool   # True iff essentially all power is propagating (kt <= k0)
    note: str


def reference_field_spectrum(field, frequency_hz):
    """Sec. 7: 2-D FFT of ``field.Ex`` (the polarization spectrum additionally reports ``Ey`` if
    it carries nonzero power). Raises ValueError if the field has no grid spacing (e.g. a
    phase-only or irregular input) -- spectral analysis needs an amplitude-and-phase field on a
    regular grid, never approximated from less."""
    if field.Ex is None or field.dx_m is None or field.dy_m is None:
        raise ValueError(
            "reference_field_spectrum: needs a complex field (Ex) on a regular grid "
            "(dx_m/dy_m); phase-only or irregular inputs are not supported here"
        )
    ny, nx = field.Ex.shape
    spectrum_x = np.fft.fftshift(np.fft.fft2(field.Ex))
    kx = np.fft.fftshift(np.fft.fftfreq(nx, d=field.dx_m)) * 2.0 * np.pi
    ky = np.fft.fftshift(np.fft.fftfreq(ny, d=field.dy_m)) * 2.0 * np.pi
    kx_grid, ky_grid = np.meshgrid(kx, ky, indexing="xy")

    power_x = np.abs(spectrum_x) ** 2
    power = power_x.copy()
    if field.Ey is not None and np.any(np.abs(field.Ey) > 0.0):
        spectrum_y = np.fft.fftshift(np.fft.fft2(field.Ey))
        power = power + np.abs(spectrum_y) ** 2

    k0 = free_space_wavenumber(frequency_hz)
    kt2 = kx_grid**2 + ky_grid**2
    propagating_mask = kt2 <= k0**2
    total_power = power.sum()
    propagating_fraction = float(power[propagating_mask].sum() / total_power) if total_power > 0 else float("nan")
    evanescent_fraction = 1.0 - propagating_fraction

    peak_idx = np.unravel_index(np.argmax(power), power.shape)
    dominant_kx, dominant_ky = float(kx_grid[peak_idx]), float(ky_grid[peak_idx])
    k_t_dominant = float(np.hypot(dominant_kx, dominant_ky))

    return SpectralAnalysis(
        kx_max_per_m=float(np.max(np.abs(kx))), ky_max_per_m=float(np.max(np.abs(ky))), k0_per_m=float(k0),
        dominant_kx_per_m=dominant_kx, dominant_ky_per_m=dominant_ky, k_t_dominant_per_m=k_t_dominant,
        propagating_fraction=propagating_fraction, evanescent_fraction=evanescent_fraction,
        local_periodic_compatible=bool(propagating_fraction > 0.99),
        note=(
            f"{propagating_fraction*100:.4g}% of spectral power is propagating (kt <= k0); "
            f"{evanescent_fraction*100:.4g}% is evanescent, reported, not discarded (Sec. 7)."
        ),
    )


@dataclass(frozen=True)
class ObjectWaveSpectralPoint:
    theta_deg: float
    phi_deg: float
    kx_per_m: float
    ky_per_m: float
    kz_per_m: complex
    is_propagating: bool


@dataclass(frozen=True)
class ObjectWaveSpectrum:
    beams: tuple                  # ObjectWaveSpectralPoint per beam
    k_required_max_per_m: float | None   # Sec. 8: interference spatial frequency across beam pairs
    all_beams_propagating: bool


def object_wave_spectrum(beams, frequency_hz):
    """Sec. 8: (kx, ky, kz) for each beam, from the actual beam specification
    (`phase8b.object_wave.transverse_wavenumber_deg`, `periodic_bc.normal_wavenumber` for the
    physical-branch kz -- both reused, unmodified), plus the multi-beam interference spatial
    frequency (`phase8d_rev.candidate_period.required_transverse_spatial_frequency_per_m`, reused
    unmodified, needs >=2 beams)."""
    k0 = free_space_wavenumber(frequency_hz)
    points = []
    for beam in beams:
        kx, ky = transverse_wavenumber_deg(frequency_hz, beam.theta_deg, beam.phi_deg)
        kz = normal_wavenumber(k0, kx, ky)
        is_propagating = bool((kx**2 + ky**2) <= k0**2)
        points.append(ObjectWaveSpectralPoint(beam.theta_deg, beam.phi_deg, kx, ky, complex(kz), is_propagating))

    k_required_max = None
    if len(beams) >= 2:
        k_required_max = required_transverse_spatial_frequency_per_m(
            tuple((b.theta_deg, b.phi_deg) for b in beams), frequency_hz,
        )
    return ObjectWaveSpectrum(
        beams=tuple(points), k_required_max_per_m=k_required_max,
        all_beams_propagating=all(p.is_propagating for p in points),
    )
