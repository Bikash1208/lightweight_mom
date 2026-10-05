"""Phase 15 Part E: the paper's 26 GHz cell pair helpers (beta modes, empty-state geometry).

The pair values are the patch-plane reflections from examples/make_paper_board_26ghz.py (26 GHz,
2.57 mm period); the beta checks use the analytic bare-slab TM0 root, so they run without a MoM solve."""

import numpy as np
import pytest

from src.app.routers.hmsa import _two_state_geometry_to_cell_m
from src.dispersion.analytic_tm0 import tm0_beta
from src.dispersion.sheet_model import find_beta_from_gamma
from src.frequency import free_space_wavenumber
from src.hologram.two_state_cells import TwoStateCell, TwoStateChoice, mixed_surface_mode, paper_surface_wave_beta
from src.materials import Substrate

FREQ = 26e9
SUB = Substrate("user", 2.2, 0.0009, 1.575e-3)


def _bare_slab_gamma(frequency_hz, eps_r, thickness_m):
    eta0 = 376.730313668
    z = 1j * eta0 / np.sqrt(eps_r) * np.tan(free_space_wavenumber(frequency_hz) * np.sqrt(eps_r) * thickness_m)
    return complex((z - eta0) / (z + eta0))


def test_paper_beta_is_mishra_eq1():
    assert paper_surface_wave_beta(FREQ, 2.2) == pytest.approx(free_space_wavenumber(FREQ) * np.sqrt(1.6), rel=1e-12)
    assert paper_surface_wave_beta(FREQ, 2.2) == pytest.approx(689.27, abs=0.05)


def test_mixed_surface_of_identical_states_is_the_single_state_beta():
    gamma = _bare_slab_gamma(FREQ, 2.2, 1.575e-3)
    state = TwoStateCell("bare", "bare", gamma, {"kind": "empty"})
    choice = TwoStateChoice(state_0=state, state_1=state, achieved_separation_deg=0.0)
    mixed = mixed_surface_mode(choice, SUB, FREQ, period_m=2.57e-3)
    single = find_beta_from_gamma(gamma, SUB, FREQ, period_m=2.57e-3)
    assert mixed.beta_rad_per_m == pytest.approx(single.beta_rad_per_m, rel=1e-6)
    assert mixed.beta_rad_per_m == pytest.approx(tm0_beta(FREQ, 2.2, 1.575e-3), rel=1e-4)


def test_empty_state_maps_to_no_patch_cell():
    assert _two_state_geometry_to_cell_m({"kind": "empty"}) == (0.0, 0.0)
