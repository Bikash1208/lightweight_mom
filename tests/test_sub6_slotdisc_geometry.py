"""SlotDiscNotch (scripts/research/sub6_slotdisc_lib.py) must reproduce SlottedCircularPatch exactly without a notch; notch/bridge masks."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts/research"))
import sub6_slotdisc_lib as S  # noqa: E402
from src.dispersion.slotted_cells import SlottedCircularPatch  # noqa: E402


def test_no_notch_equals_slotted_circular_patch():
    g = S.SlotDiscNotch(3e-3, 4.8e-3, 1.6e-3)
    h = SlottedCircularPatch(6e-3, 4.8e-3, 1.6e-3, 0.0)
    x, y = np.meshgrid(np.linspace(-3.4e-3, 3.4e-3, 301), np.linspace(-3.4e-3, 3.4e-3, 301))
    assert np.array_equal(g.in_feature(x, y), h.in_feature(x, y))


def test_notch_and_bridge_masks():
    r, L, W, gc = 5.4e-3, 8.64e-3, 2.9e-3, 1.0e-3
    g = S.SlotDiscNotch(r, L, W, gc, bridge=False)
    b = S.SlotDiscNotch(r, L, W, gc, bridge=True)
    xr = 0.5 * (L / 2 + r)
    assert not g.in_feature(np.array([xr]), np.array([0.0]))[0]  # empty notch
    assert b.in_feature(np.array([xr]), np.array([0.0]))[0]  # bridge region is metal
    assert not b.in_feature(np.array([L / 2 + 0.01e-3]), np.array([0.0]))[0]  # notch between slot end and bridge stays open


def test_off_bar50_impedance():
    z = S.z_off_bar50(5.8e9)
    assert abs(z - (6.7 - 160.9j)) < 1.0  # (0.15 pF || 5 kohm) + 0.6 nH at 5.8 GHz
