"""Phase 15 Part A: sector (quadrant) multiplexing partition -- `sector_index_map` must reduce to
equal quadrants for 4 beams 90deg apart (the brief's own worked example) and stay well-defined for
an uneven beam spacing."""

import numpy as np
import pytest

from src.hologram.sector_multiplexing import sector_boundaries_deg, sector_index_map
from src.phase8b.object_wave import Beam


def test_four_beams_90deg_apart_give_four_equal_quadrants():
    beams = [Beam(theta_deg=20.0, phi_deg=p) for p in (0.0, 90.0, 180.0, 270.0)]
    # one point at each beam's own azimuth (phi' = 0, 90, 180, 270 deg)
    x = np.array([1.0, 0.0, -1.0, 0.0])
    y = np.array([0.0, 1.0, 0.0, -1.0])
    idx = sector_index_map(x, y, (0.0, 0.0), beams)
    assert list(idx) == [0, 1, 2, 3]


def test_sector_boundaries_are_the_midpoints_between_beams():
    beams = [Beam(theta_deg=20.0, phi_deg=p) for p in (0.0, 90.0, 180.0, 270.0)]
    assert sector_boundaries_deg((0.0, 0.0), beams) == pytest.approx([45.0, 135.0, 225.0, 315.0])


def test_uneven_beam_spacing_still_partitions_the_full_circle():
    beams = [Beam(theta_deg=20.0, phi_deg=p) for p in (0.0, 30.0, 200.0)]
    phi_grid = np.linspace(0.0, 359.0, 360)
    x = np.cos(np.radians(phi_grid))
    y = np.sin(np.radians(phi_grid))
    idx = sector_index_map(x, y, (0.0, 0.0), beams)
    assert set(idx.tolist()) == {0, 1, 2}  # every beam owns at least one point, no 4th index appears


def test_sector_index_map_is_about_the_given_feed():
    beams = [Beam(theta_deg=20.0, phi_deg=0.0), Beam(theta_deg=20.0, phi_deg=180.0)]
    x, y = np.array([3.0]), np.array([0.0])
    idx_from_origin = sector_index_map(x, y, (0.0, 0.0), beams)
    idx_from_far_feed = sector_index_map(x, y, (10.0, 0.0), beams)  # point now "behind" this feed
    assert idx_from_origin[0] == 0
    assert idx_from_far_feed[0] == 1


def test_empty_beams_raises():
    with pytest.raises(ValueError):
        sector_index_map(np.array([0.0]), np.array([0.0]), (0.0, 0.0), [])
    with pytest.raises(ValueError):
        sector_boundaries_deg((0.0, 0.0), [])
