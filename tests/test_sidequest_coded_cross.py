"""Crossed two-diode cell (x and y diodes): the four states differ and the cell is symmetric under x <-> y."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "research"))
import sidequest_coded_cross_cell as X  # noqa: E402
import sidequest_gate_cell as gc  # noqa: E402
from src.hologram.coded_cell_eval import sheet_tensor  # noqa: E402


@pytest.fixture(scope="module")
def ys():
    return {(a, b): sheet_tensor(X.cross_cell(a, b), gc.F0, 32) for a in ("OFF", "ON") for b in ("OFF", "ON")}


def test_each_diode_switches_its_own_axis(ys):
    assert (ys["ON", "OFF"][0, 0] - ys["OFF", "OFF"][0, 0]).imag > 1e-3  # x diode raises Y_xx by > 1 mS
    assert abs(ys["ON", "OFF"][1, 1] - ys["OFF", "OFF"][1, 1]) < 1e-4  # ... and leaves Y_yy alone
    assert (ys["OFF", "ON"][1, 1] - ys["OFF", "OFF"][1, 1]).imag > 1e-3


def test_x_y_transposition_symmetry_and_reciprocity(ys):
    assert np.allclose(ys["ON", "OFF"], ys["OFF", "ON"][::-1, ::-1], rtol=1e-6, atol=1e-9)
    for y in ys.values():
        assert abs(y[0, 1] - y[1, 0]) < 1e-9 * abs(y).max()
