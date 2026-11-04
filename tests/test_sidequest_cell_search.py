"""Cell-search builders: notch superpixel geometry and the loaded Bloch tensor (diode state must move the on-axis admittance)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "research"))
import sidequest_cell_search as C  # noqa: E402


def test_notch_superpixel_geometry_and_loads():
    c = C.superpixel_notch("ro3010", "ON", "OFF", wd=0.4)
    assert len(c.loads) == 2 and {ld.direction for ld in c.loads} == {"x", "y"}
    x = (np.arange(64) + 0.5) * 0.1e-3 - 3.2e-3
    gx, gy = np.meshgrid(x, x, indexing="ij")
    m = c.is_metal(gx, gy)
    assert 0.75 < m.mean() < 0.95  # four patches with 0.2 mm gaps
    assert not c.is_metal(np.array([0.0]), np.array([2.0e-3]))[0]  # central gap away from the diode window is open


def test_diode_state_moves_the_on_axis_admittance_rt5880_disc():
    from src.frequency import free_space_wavenumber

    f0 = C.STACKS["rt5880"][1]
    k0 = free_space_wavenumber(f0)
    off = C.bloch_loaded(C.cell("rt5880", "OFF", "OFF", diam=6.0, ls=2.8, pos=0.8), f0, 32, 1.08 * k0)
    on = C.bloch_loaded(C.cell("rt5880", "ON", "OFF", diam=6.0, ls=2.8, pos=0.8), f0, 32, 1.08 * k0)
    assert (on[0, 0] - off[0, 0]).imag < -2e-3  # ON lowers Y_xx by > 2 mS
    assert abs(on[1, 1] - off[1, 1]) < 5e-4  # and leaves Y_yy alone
