"""Geometry checks for the 5.8 GHz meander screen (no RF solver needed)."""

import math

import numpy as np
import pytest

from src.dispersion.elliptical_cells import MeanderSlotPatch, VerticalLegMeanderSlotPatch


def _patch(**kw):
    return MeanderSlotPatch(4.2e-3, 2.4e-3, 2, 0.6e-3, 0.2e-3, **kw)


def test_slot_stays_inside_disc_and_rotates():
    g = _patch()
    axis = np.linspace(-2.1e-3, 2.1e-3, 421)
    x, y = np.meshgrid(axis, axis)
    metal = g.in_feature(x, y)
    area = metal.sum() * (axis[1] - axis[0]) ** 2
    assert area == pytest.approx(g.feature_area, rel=0.015)
    assert not g.in_feature(0.0, -0.8e-3)  # first slot run
    assert g.in_feature(0.0, -0.4e-3)  # copper finger
    rotated = _patch(orientation=math.pi / 2)
    # A few samples land exactly on a floating-point boundary after rotation.
    assert np.mean(rotated.in_feature(x, y) != g.in_feature(y, -x)) < 0.001


def test_removed_tooth_makes_isolated_island():
    g = _patch(removed_tooth=0, island_length=0.8e-3)
    assert g.in_feature(0.0, -0.4e-3)
    assert not g.in_feature(0.5e-3, -0.4e-3)
    assert not g.in_feature(0.0, -0.65e-3)
    assert not g.in_feature(0.0, -0.15e-3)


def test_fabrication_minima_enforced():
    with pytest.raises(ValueError, match="0.2 mm"):
        MeanderSlotPatch(4.2e-3, 2.4e-3, 2, 0.19e-3, 0.2e-3)
    with pytest.raises(ValueError, match="0.2 mm"):
        MeanderSlotPatch(4.2e-3, 2.4e-3, 2, 0.6e-3, 0.19e-3)
    with pytest.raises(ValueError, match="clear"):
        _patch(island_length=0.8e-3, removed_tooth=0, island_clearance=0.1e-3)


def test_proposed_seven_vertical_legs_and_continuous_cut():
    g = VerticalLegMeanderSlotPatch(11.2e-3, 7, 1.2e-3, 3.0e-3, 0.4e-3)
    for i in range(7):
        x = (i - 3) * 1.2e-3
        assert not g.in_feature(x, 0.0)
        assert not g.in_feature(x, -1.0e-3)
    for i in range(6):
        x = (i - 2.5) * 1.2e-3
        y = 1.3e-3 if i % 2 == 0 else -1.3e-3
        assert not g.in_feature(x, y)  # alternating bridge joins adjacent cuts
        assert g.in_feature(x, -y)  # copper finger at opposite end
    assert g.in_feature(5.5e-3, 0.0)  # rim intact
    assert g.in_feature(-5.5e-3, 0.0)
    axis = np.linspace(-5.6e-3, 5.6e-3, 1121)
    x, y = np.meshgrid(axis, axis)
    assert g.in_feature(x, y).sum() * (axis[1] - axis[0])**2 == pytest.approx(g.feature_area, rel=0.01)
