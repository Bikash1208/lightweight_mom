"""Side quest: physics checks of the gate-strip line model (scripts/research/sidequest_gate_cell.py)."""
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "research"))
warnings.simplefilter("ignore")
import sidequest_gate_cell as gc  # noqa: E402

YREF = 1j * 7.82e-3


def test_identical_region_is_transparent():
    r = gc.strip_transmission(YREF, YREF, 3)
    assert abs(r["S21_db"]) < 0.05 and r["S11_db"] < -60  # tiny loss only (tan d)


def test_lossless_gate_conserves_power():
    for b in (-200e-3, -20e-3, 2e-3, 20e-3):
        r = gc.strip_transmission(YREF, 1j * b, 3)
        p = 10 ** (r["S21_db"] / 10) + 10 ** (r["S11_db"] / 10)
        assert 0.99 < p <= 1.0 + 1e-6, (b, p)  # deficit = dielectric loss (tan d 0.0009) only


def test_more_inductive_sheet_blocks_more():
    s = [gc.strip_transmission(YREF, -1j * b * 1e-3, 3)["S21_db"] for b in (20, 100, 500)]
    assert s[0] > s[1] > s[2]


def test_beta_matches_sheet_model_for_reactive_sheet():
    from src.dispersion.sheet_model import find_beta_from_sheet_admittance

    ref = find_beta_from_sheet_admittance(YREF, gc.SUB, gc.F0).beta_rad_per_m
    assert abs(gc.beta_complex(YREF).real - ref) / ref < 1e-6


def test_cell_dimensions_must_be_on_the_mesh_grid():
    import pytest

    with pytest.raises(ValueError):
        gc.make_dogbone(1.23e-3, 0.8e-3, 6.0e-3, 0.4e-3, "ON", lt=5.6e-3)


def test_conductive_sheet_has_damped_quasi_tem_root():
    roots = gc.roots_all((100.1 + 68.0j) * 1e-3)
    assert any(1.40 < r.real / gc.K0 < 1.55 and 5.0 < -r.imag < 15.0 for r in roots)  # sheet over ground: PPW-like mode, alpha ~ 8.6 Np/m


def test_transition_row_matches_the_pass_state():
    import sidequest_gate_match as gm

    yg = 0.0009e-3 + 2.5e-3j
    s11_bare, _ = gm.cascade(YREF, [(yg, 3)])
    tr = [(1j * 5.64e-3, 1)]
    s11, s21 = gm.cascade(YREF, tr + [(yg, 3)] + tr)
    assert -20 * np.log10(abs(s11_bare)) < 10 and -20 * np.log10(abs(s11)) > 15 and -20 * np.log10(abs(s21)) < 0.5
