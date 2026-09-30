"""Phase 15 Part A: `feed_position_from_nan_mask` -- the feed (x0, y0) a reference field's own
NaN/singularity sits at, used as the default origin for polarization-compensation's rho_hat."""

import numpy as np

from src.fields.reference_field_grid import feed_position_from_nan_mask


def test_feed_position_is_the_centroid_of_the_nan_cells():
    x, y = np.meshgrid(np.linspace(-0.02, 0.02, 5), np.linspace(-0.02, 0.02, 5), indexing="xy")
    nan_mask = np.zeros_like(x, dtype=bool)
    nan_mask[2, 2] = True  # the single feed-point NaN, at the grid centre by construction
    x0, y0 = feed_position_from_nan_mask(x, y, nan_mask)
    assert x0 == 0.0 and y0 == 0.0


def test_feed_position_off_centre():
    x, y = np.meshgrid(np.linspace(0.0, 0.04, 5), np.linspace(0.0, 0.04, 5), indexing="xy")
    nan_mask = np.zeros_like(x, dtype=bool)
    nan_mask[1, 3] = True  # x=0.03, y=0.01
    x0, y0 = feed_position_from_nan_mask(x, y, nan_mask)
    assert x0 == 0.03 and y0 == 0.01


def test_no_nan_defaults_to_origin():
    x, y = np.meshgrid(np.linspace(-0.02, 0.02, 5), np.linspace(-0.02, 0.02, 5), indexing="xy")
    nan_mask = np.zeros_like(x, dtype=bool)
    x0, y0 = feed_position_from_nan_mask(x, y, nan_mask)
    assert (x0, y0) == (0.0, 0.0)
