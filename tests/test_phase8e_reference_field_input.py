"""Sec. 3-4: reference-field input construction."""

import numpy as np
import pytest

from src.phase8e.reference_field_input import (
    monopole_demo_reference_field,
    phase_only_reference_field,
    spherical_feed_reference_field,
)


def test_spherical_feed_field_is_reference_plane_compatible_by_construction():
    x, y = np.meshgrid(np.linspace(-0.02, 0.02, 5), np.linspace(-0.02, 0.02, 5), indexing="xy")
    field = spherical_feed_reference_field(x, y, feed_height_m=0.05, frequency_hz=12.5e9)
    assert field.z_reference_m == 0.0
    assert field.amplitude_available
    assert field.source == "SYNTHETIC_DEMO_SPHERICAL_FEED"


def test_spherical_feed_amplitude_decays_with_distance_from_boresight():
    x, y = np.meshgrid(np.linspace(-0.05, 0.05, 11), np.linspace(-0.05, 0.05, 11), indexing="xy")
    field = spherical_feed_reference_field(x, y, feed_height_m=0.05, frequency_hz=12.5e9)
    center_idx = 5, 5
    edge_idx = 0, 0
    assert abs(field.Ex[center_idx]) > abs(field.Ex[edge_idx])


def test_spherical_feed_rejects_nonpositive_feed_height():
    x, y = np.meshgrid(np.linspace(-0.01, 0.01, 3), np.linspace(-0.01, 0.01, 3), indexing="xy")
    with pytest.raises(ValueError):
        spherical_feed_reference_field(x, y, feed_height_m=0.0, frequency_hz=12.5e9)


def test_monopole_demo_field_requires_k_swg():
    x, y = np.meshgrid(np.linspace(-0.01, 0.01, 3), np.linspace(-0.01, 0.01, 3), indexing="xy")
    with pytest.raises(ValueError):
        monopole_demo_reference_field(x, y, k_swg=None)


def test_monopole_demo_field_z_reference_is_zero_surface_wave():
    x, y = np.meshgrid(np.linspace(-0.01, 0.01, 3), np.linspace(-0.01, 0.01, 3), indexing="xy")
    field = monopole_demo_reference_field(x, y, k_swg=300.0)
    assert field.z_reference_m == 0.0


def test_phase_only_field_marks_amplitude_unavailable_and_never_fabricates_ex():
    x, y = np.meshgrid(np.linspace(-0.01, 0.01, 3), np.linspace(-0.01, 0.01, 3), indexing="xy")
    phase = np.zeros_like(x)
    field = phase_only_reference_field(x, y, phase)
    assert field.amplitude_available is False
    assert field.Ex is None
    assert field.Ey is None
    assert np.array_equal(field.phase(), phase)


def test_phase_method_derives_from_ex_when_amplitude_available():
    x, y = np.meshgrid(np.linspace(-0.01, 0.01, 3), np.linspace(-0.01, 0.01, 3), indexing="xy")
    field = spherical_feed_reference_field(x, y, feed_height_m=0.05, frequency_hz=12.5e9)
    assert np.allclose(field.phase(), np.degrees(np.angle(field.Ex)))
