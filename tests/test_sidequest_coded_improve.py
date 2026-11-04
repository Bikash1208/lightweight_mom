"""Side quest round 11: efficiency-improvement model (I1-I7) sanity: library mapping, option design, power bookkeeping, launcher baseline."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "research"))
import sidequest_coded_improve as I  # noqa: E402

BD = I.Board()
VS, PAIRS = I.bit_sets()
AMP = np.sqrt(2 * I.A0 / BD.rho) * np.exp(-0.5 * (I.A0 + I.AA0) * (BD.rho - I.st.RHO_KEEP))


def test_library_maps_the_baseline_cell_to_plus_minus_one():
    lib = I.Lib()
    v = lib.sw[(2.8, 1.6)]
    assert abs(v["m_off"] - 1.0) < 1e-9 and abs(v["m_on"] + 1.0) < 1e-9 and abs(v["g_on"] - 1.0) < 1e-9
    assert min(lib.static.values()) < -1.0 < 1.0  # the solid disc / shorted slit reaches beyond the diode ON value


def test_bit_design_values_and_diode_count():
    m, nd, loss = I.design(BD, AMP, VS, PAIRS, [0.0, 0.0], 0.6)
    assert set(np.unique(m[:, :2])) <= {-1.0, 1.0} and not m[:, 2].any()
    static_el = (np.abs(m[0, :2] - m[1, :2]) < 1e-12)  # (2 elements, n)
    assert int((~static_el).sum()) == int(nd.sum())
    assert np.all(loss >= 0) and loss.max() <= 1.0 + 1e-12


def test_more_penalty_means_fewer_diodes():
    n = [I.design(BD, AMP, VS, PAIRS, [0.0, 0.0], lam)[1].sum() for lam in (0.0, 0.3, 1.0)]
    assert n[0] > n[1] > n[2]


def test_power_bookkeeping_is_bounded():
    m, nd, loss = I.design(BD, AMP, VS, PAIRS, [0.0, 0.0], 0.6)
    amp, conv, edge, alpha, aa = I.amplitude_profile(BD, m, loss)
    assert 0.0 < conv < 1.0 and 0.0 < edge < 1.0 and conv + edge < 1.0
    assert abs(alpha.mean() - I.A0) < 1e-9  # bit cells: |m| = 1 on both axes -> alpha = A0 for every cell


def test_launcher_baseline_matches_the_model_value():
    from src.dispersion.pin_launch import Launcher, PinLaunch

    s = PinLaunch(10.54e9, 2.2, 1.575e-3, Launcher(top_above_surface_m=6e-3, pin_radius_m=0.305e-3), tan_delta=0.0009).split()
    assert abs(s["launch_fraction"] - 0.627) < 0.01


def test_e1_tensor_options_radiate_the_same_power_in_every_direction():
    tens = I.e1_options(1.3)
    mx, my, mxy, nd, ld, ft = I.option_arrays(np.array([0.0]), [], tens)
    assert ft.sum() == len(tens) and np.all(nd[ft] == 0)
    rx, ry = np.cos(np.linspace(0, 2 * np.pi, 7)), np.sin(np.linspace(0, 2 * np.pi, 7))
    for i in np.flatnonzero(ft)[1:]:
        p = (mx[i, 0] * rx + mxy[i, 0] * ry) ** 2 + (mxy[i, 0] * rx + my[i, 0] * ry) ** 2  # |Delta Y rho_hat|^2
        assert np.allclose(p, p[0]) and abs(p[0] - (mx[i, 0] ** 2 + mxy[i, 0] ** 2)) < 1e-12  # trace-free tensor: reflection, norm preserved
