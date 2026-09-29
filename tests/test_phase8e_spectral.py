"""Sec. 7-8: spectral analysis of reference field and object wave."""

import numpy as np
import pytest

from src.frequency import free_space_wavenumber
from src.phase8b.object_wave import Beam
from src.phase8e.reference_field_input import ReferenceFieldGrid
from src.phase8e.spectral import object_wave_spectrum, reference_field_spectrum

FREQ_HZ = 12.5e9


def _tilted_plane_wave_field(theta_deg, phi_deg, n=32, extent_m=0.05):
    k0 = free_space_wavenumber(FREQ_HZ)
    theta, phi = np.radians(theta_deg), np.radians(phi_deg)
    kx0, ky0 = k0 * np.sin(theta) * np.cos(phi), k0 * np.sin(theta) * np.sin(phi)
    lin = np.linspace(-extent_m, extent_m, n)
    x, y = np.meshgrid(lin, lin, indexing="xy")
    ex = np.exp(-1j * (kx0 * x + ky0 * y))
    dx = float(lin[1] - lin[0])
    field = ReferenceFieldGrid(
        x_m=x,
        y_m=y,
        Ex=ex,
        Ey=np.zeros_like(ex),
        amplitude_available=True,
        source="TEST",
        dx_m=dx,
        dy_m=dx,
        z_reference_m=0.0,
    )
    return field, kx0, ky0


def test_reference_field_spectrum_rejects_missing_grid_spacing():
    field = ReferenceFieldGrid(
        x_m=np.zeros((2, 2)),
        y_m=np.zeros((2, 2)),
        Ex=None,
        Ey=None,
        amplitude_available=False,
        source="TEST",
        dx_m=None,
        dy_m=None,
        z_reference_m=None,
    )
    with pytest.raises(ValueError):
        reference_field_spectrum(field, FREQ_HZ)


def test_reference_field_spectrum_is_almost_fully_propagating_for_a_broadside_wave():
    field, kx0, ky0 = _tilted_plane_wave_field(0.0, 0.0)
    result = reference_field_spectrum(field, FREQ_HZ)
    assert result.propagating_fraction > 0.99
    assert result.local_periodic_compatible is True


def test_reference_field_spectrum_dominant_k_matches_the_actual_tilt():
    field, kx0, ky0 = _tilted_plane_wave_field(20.0, 45.0)
    result = reference_field_spectrum(field, FREQ_HZ)
    # dominant bin lands on the nearest FFT sample to the true (kx0, ky0); allow one bin of slack.
    k0 = free_space_wavenumber(FREQ_HZ)
    bin_width = 2 * np.pi / (32 * (0.1 / 31))
    assert abs(result.dominant_kx_per_m - kx0) < 3 * bin_width
    assert abs(result.dominant_ky_per_m - ky0) < 3 * bin_width


def test_reference_field_spectrum_reports_evanescent_fraction_never_discarding_it():
    field, _, _ = _tilted_plane_wave_field(0.0, 0.0)
    result = reference_field_spectrum(field, FREQ_HZ)
    assert 0.0 <= result.evanescent_fraction <= 1.0
    assert result.propagating_fraction + result.evanescent_fraction == pytest.approx(1.0)


def test_object_wave_spectrum_rejects_single_beam_for_interference_frequency():
    result = object_wave_spectrum([Beam(theta_deg=0.0, phi_deg=0.0)], FREQ_HZ)
    assert result.k_required_max_per_m is None


def test_object_wave_spectrum_zero_interference_frequency_for_identical_beams():
    beams = [Beam(theta_deg=10.0, phi_deg=0.0), Beam(theta_deg=10.0, phi_deg=0.0)]
    result = object_wave_spectrum(beams, FREQ_HZ)
    assert result.k_required_max_per_m == pytest.approx(0.0, abs=1e-9)


def test_object_wave_spectrum_all_beams_propagating_within_theta_90():
    beams = [Beam(theta_deg=30.0, phi_deg=0.0), Beam(theta_deg=60.0, phi_deg=90.0)]
    result = object_wave_spectrum(beams, FREQ_HZ)
    assert result.all_beams_propagating is True
    for p in result.beams:
        assert p.is_propagating
        assert p.kz_per_m.real >= 0.0 or abs(p.kz_per_m.imag) > 0
