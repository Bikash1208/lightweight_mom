"""Sec. 6: reference-plane compatibility and angular-spectrum propagation."""

import numpy as np

from src.phase8e.reference_plane import check_reference_plane_compatibility, propagate_angular_spectrum

FREQ_HZ = 12.5e9


def test_zero_height_is_compatible():
    result = check_reference_plane_compatibility(0.0)
    assert result.status == "COMPATIBLE"


def test_nonzero_height_is_a_mismatch():
    result = check_reference_plane_compatibility(1.7e-3)
    assert result.status == "BRIDGE_REFERENCE_PLANE_MISMATCH"
    assert result.delta_z_m == 1.7e-3


def test_unknown_height_is_reported_as_unknown_not_assumed_compatible():
    result = check_reference_plane_compatibility(None)
    assert result.status == "UNKNOWN_REFERENCE_PLANE"


def test_propagate_angular_spectrum_round_trip_recovers_the_original_field():
    n = 16
    lin = np.linspace(-0.02, 0.02, n)
    x, y = np.meshgrid(lin, lin, indexing="xy")
    dx = float(lin[1] - lin[0])
    field = np.exp(-0.5 * ((x / 0.01) ** 2 + (y / 0.01) ** 2)) + 0j  # smooth, mostly-propagating field

    forward = propagate_angular_spectrum(field, dx, dx, FREQ_HZ, 0.005)
    back = propagate_angular_spectrum(forward, dx, dx, FREQ_HZ, -0.005)
    assert np.allclose(back, field, atol=1e-9)


def test_propagate_angular_spectrum_zero_distance_is_identity():
    n = 8
    lin = np.linspace(-0.01, 0.01, n)
    x, y = np.meshgrid(lin, lin, indexing="xy")
    dx = float(lin[1] - lin[0])
    field = np.exp(1j * 3.0 * x)
    result = propagate_angular_spectrum(field, dx, dx, FREQ_HZ, 0.0)
    assert np.allclose(result, field, atol=1e-9)
