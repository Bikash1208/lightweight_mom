"""Phase 20 R3 item 3: LP-3 cell and the optimisation helpers (scripts/research/phase20_r3_opt*.py; the optimiser itself lives outside the solver core)."""

import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

from src.basis_functions import MeshSettings
from src.excitation import ReferencePlanes
from src.ris.cells import OFF, ON, Lp3Cell, madp_kit, smp1345_kit
from src.sparameters import compute_s_parameters_oblique

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "research"))
import phase20_r3_opt as opt  # noqa: E402


def test_lattice_fit_recognises_an_ideal_ladder_in_any_order_and_offset():
    ph = np.array([30.0, 120.0, 210.0, 300.0])
    for perm in ([0, 1, 2, 3], [2, 0, 3, 1]):
        dev, _ = opt.lattice_fit(ph[perm] + 17.0)
        assert dev == pytest.approx(0.0, abs=1e-9)


def test_lattice_fit_reports_the_deviation_of_a_distorted_ladder():
    dev, _ = opt.lattice_fit([0.0, 90.0, 180.0, 270.0 + 24.0])
    assert 8.0 < dev < 24.0  # best offset spreads the 24 deg error over the four states


def test_step_error_keeps_the_nominal_assignment_and_counts_reversed_ladders():
    perm = (0, 1, 2, 3)
    assert opt.step_error([0, 90, 180, 270], perm) == pytest.approx(0.0, abs=1e-9)
    assert opt.step_error([0, 90, 200, 270], perm) == pytest.approx(20.0, abs=1e-9)
    assert opt.step_error([0, -90, -180, -270], perm) == pytest.approx(
        0.0, abs=1e-9
    )  # mirrored ladder is the same ladder


def test_corner_kits_smp_has_ct_and_ls_corners_madp_an_assumed_ls_offset():
    smp = opt.corner_kits(smp1345_kit())
    assert set(smp) == {"nominal", "C_T-20%", "C_T+20%", "L_s-20%", "L_s+20%"}
    assert pytest.approx(1.2 * smp["nominal"].off.C_T) == smp["C_T+20%"].off.C_T
    assert smp["L_s+20%"].on.L_s == pytest.approx(1.2 * smp["nominal"].on.L_s)
    madp = opt.corner_kits(madp_kit())
    assert set(madp) == {"nominal", "C_T-20%", "C_T+20%", "L_s+0.1nH(assumed)"}
    assert madp["nominal"].off.L_s == 0.0 and madp["L_s+0.1nH(assumed)"].off.L_s == pytest.approx(0.1e-9)


def test_lp3_cell_has_three_diodes_eight_codes_and_fits_the_cell():
    c = Lp3Cell()
    assert c.n_diodes == 3 and len(c.states()) == 8 and len(c.diode_sites()) == 3
    uc = c.unit_cell((ON, OFF, ON), madp_kit())
    assert len(uc.loads) == 3
    xs = sorted(s[0] for s in c.diode_sites())
    assert xs[0] < 0 < xs[1] < xs[2]  # d1 on the -x side, the chain on +x with d3 behind d2


def test_lp3_dead_segment_behind_an_open_diode_changes_little_but_a_closed_chain_does_more():
    c = Lp3Cell()
    pl = ReferencePlanes(0.0, None, 1.575e-3)

    def tm(state):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            s = compute_s_parameters_oblique(
                c.unit_cell(state, madp_kit()), 10e9, MeshSettings(56, 56, 29), pl, 0.0, 0.0
            )
        return s.s11[1, 1]

    base = tm((OFF, OFF, OFF))
    floating = tm((OFF, OFF, ON))  # segment 2 connected to an unconnected segment 1: nearly the same as all-OFF
    chained = tm((OFF, ON, ON))
    assert abs(floating - base) < abs(chained - base)
    assert abs(tm((OFF, OFF, OFF))) <= 1.0 + 1e-6
