"""Phase 12A: ``ReferencePlanes.bottom_air_thickness=None`` (no bottom port, grounded stack)."""

import pytest

from src.excitation import ReferencePlanes


def test_existing_two_port_construction_unchanged():
    rp = ReferencePlanes(10e-3, 5e-3, 1.575e-3)
    assert rp.z_top == 10e-3
    assert rp.z_bottom == -(1.575e-3 + 5e-3)
    assert rp.grounded is False


def test_none_bottom_thickness_marks_grounded():
    rp = ReferencePlanes(10e-3, None, 1.575e-3)
    assert rp.grounded is True
    assert rp.bottom_air_thickness is None


def test_z_bottom_raises_when_grounded():
    rp = ReferencePlanes(10e-3, None, 1.575e-3)
    with pytest.raises(ValueError, match="bottom"):
        rp.z_bottom


def test_negative_bottom_thickness_still_rejected():
    with pytest.raises(ValueError):
        ReferencePlanes(10e-3, -1e-3, 1.575e-3)


def test_negative_top_thickness_still_rejected():
    with pytest.raises(ValueError):
        ReferencePlanes(-1e-3, 5e-3, 1.575e-3)
