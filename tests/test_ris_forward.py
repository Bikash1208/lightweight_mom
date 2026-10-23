"""Phase 20 R5-2: plane-wave forward model of a coded reflecting RIS -- analytic and cross-checks."""

import math

import numpy as np
import pytest

from src.hologram.coded_synthesis import textbook_quantization_loss_db
from src.ris import forward as fw
from src.ris.calibrated_library import CalibratedLibrary

F = 10.5
LAM = fw.C0 / (F * 1e9)
P = 15e-3


def _ideal(rows, cols, phase, f=F, theta_i=0.0, phi_i=0.0):
    return fw.aperture_from_amplitudes(np.exp(1j * phase), np.zeros((rows, cols)), f, theta_i, phi_i, P)


@pytest.fixture(scope="module")
def cal():
    return CalibratedLibrary.load()


def test_uniform_aperture_peak_gain_and_analytic_pattern():
    ap = _ideal(20, 20, np.zeros((20, 20)))
    g0, _ = fw.pattern_cut(ap, 0.0, np.array([0.0]))
    assert g0[0] == pytest.approx(4 * math.pi * ap.area_m2 / LAM**2, rel=1e-12)
    th = np.radians(np.linspace(0.5, 80.0, 40))
    g, _ = fw.pattern_cut(ap, 0.0, np.degrees(th))
    u = np.sin(th)
    k0 = 2 * math.pi / LAM
    af = np.abs(np.sin(20 * k0 * u * P / 2) / (20 * np.sin(k0 * u * P / 2)))
    ef = np.sinc(k0 * u * P / 2 / math.pi) * np.sinc(0.0)  # uy = 0
    expected = g0[0] * (af * ef * np.cos(th)) ** 2
    assert np.allclose(g, expected, rtol=1e-9, atol=1e-12)


def test_ideal_gradient_steers_to_the_analytic_angle_with_full_aperture_gain():
    n, ts = 20, 30.0
    x, y = fw.cell_positions(n, n, P)
    k0 = 2 * math.pi / LAM
    ap = _ideal(n, n, -k0 * math.sin(math.radians(ts)) * x[None, :] * np.ones((n, 1)))
    m = fw.metrics(ap, obliquity="none", target=(ts, 0.0))
    assert m.pointing_error_deg < 0.1
    ef = (np.sinc(k0 * math.sin(math.radians(ts)) * P / 2 / math.pi)) ** 2  # cell element factor, power
    assert 10 ** (m.peak_gain_dbi / 10) == pytest.approx(4 * math.pi * ap.area_m2 / LAM**2 * ef, rel=2e-3)
    assert m.cross_peak_db < -100  # no opposite-handed power for an ideal cell


def test_lossless_aperture_radiates_the_incident_power_with_the_cos_obliquity():
    n = 20
    ap = _ideal(n, n, np.zeros((n, n)))
    assert fw.metrics(ap).radiated_over_incident == pytest.approx(1.0, abs=0.05)
    assert fw.metrics(ap).aperture_efficiency == pytest.approx(1.0, abs=0.05)
    x, _ = fw.cell_positions(n, n, P)
    k0 = 2 * math.pi / LAM
    steered = _ideal(n, n, -k0 * math.sin(math.radians(30.0)) * x[None, :] * np.ones((n, 1)))
    assert (
        0.6 < fw.metrics(steered).radiated_over_incident < 1.02
    )  # cos(theta) x cell element factor^2 ~ 0.87 x 0.81, never above the incident power


def test_two_bit_quantisation_loss_matches_theory():
    rows, cols = 12, 80
    x, _ = fw.cell_positions(rows, cols, P)
    k0 = 2 * math.pi / LAM
    ts = math.asin(0.2713)  # irrational-ish period (many distinct phase errors)
    psi = -k0 * math.sin(ts) * x[None, :] * np.ones((rows, 1))
    q = np.round(psi / (math.pi / 2)) * (math.pi / 2)
    cont = fw.metrics(_ideal(rows, cols, psi), obliquity="none")
    quant = fw.metrics(_ideal(rows, cols, q), obliquity="none")
    loss = cont.peak_gain_dbi - quant.peak_gain_dbi
    assert loss == pytest.approx(textbook_quantization_loss_db(2), abs=0.15)


def test_plain_two_bit_gradient_has_the_known_quantisation_lobe():
    rows, cols = 8, 64
    x, _ = fw.cell_positions(rows, cols, P)
    k0 = 2 * math.pi / LAM
    psi = -k0 * math.sin(math.radians(20.0)) * x[None, :] * np.ones((rows, 1))
    q = np.round(psi / (math.pi / 2)) * (math.pi / 2)
    m = fw.metrics(_ideal(rows, cols, q), obliquity="none")
    # a staircase of 4 levels: strongest quantisation lobe about 20 log10(1/3)... within a few dB of -9.5 dB
    assert -13.0 < m.sll_db < -7.0


def test_energy_bookkeeping_per_state(cal):
    for theta in (0.0, 30.0):
        for f in (10.0, 10.5, 11.0):
            for s in (1, 2, 3, 4):
                c = cal.response(f, theta, 0.0, s).c
                assert (
                    abs(c[0, 0]) ** 2 + abs(c[1, 0]) ** 2 <= 1.0 + 0.05
                )  # HFSS-calibrated: passive within the calibration scatter


def test_periodic_gradient_orders_equal_the_isolated_cell_prediction(cal):
    """The R3 prediction b_m = (1/N) sum_site r_site exp(+j 2 pi m x/N) is the large-array limit of the forward model at the order direction."""
    code = [1, 2, 3, 4]
    nx, reps = 4, 11
    states = np.tile(np.array(code * reps), (3, 1))
    ap = fw.aperture_from_states(cal, states, F, 0.0, 0.0, P)
    cs = {s: cal.response(F, 0.0, 0.0, s).c for s in (1, 2, 3, 4)}
    k0 = 2 * math.pi / LAM
    for m in (-1, 0, 1):
        pred = {
            w: sum(cs[code[i]][w, 0] * np.exp(2j * math.pi * m * ((i - (nx - 1) / 2) / nx)) for i in range(nx)) / nx
            for w in (0, 1)
        }
        ux = m * LAM / (nx * P)
        ec, ex = fw._fields(ap, np.array([ux]), np.array([0.0]), "none")
        norm = ap.area_m2 * np.sinc(0.5 * k0 * ux * P / math.pi)
        assert ec[0] / norm == pytest.approx(pred[0], abs=1e-9)
        assert ex[0] / norm == pytest.approx(pred[1], abs=1e-9)


def test_pattern_symmetry_under_phi_plus_180_for_a_symmetric_coding():
    rng = np.random.default_rng(3)
    half = rng.uniform(0, 2 * np.pi, (6, 6))
    ph = np.block([[half, half[:, ::-1]], [half[::-1], half[::-1, ::-1]]])  # symmetric under x -> -x, y -> -y
    ap = _ideal(12, 12, ph)
    th = np.linspace(1, 70, 25)
    g1, _ = fw.pattern_cut(ap, 20.0, th)
    g2, _ = fw.pattern_cut(ap, 200.0, th)
    assert np.allclose(g1, g2, rtol=1e-9)


def test_refusals_and_option_checks():
    ap = _ideal(4, 4, np.zeros((4, 4)))
    with pytest.raises(ValueError):
        fw._fields(ap, np.array([0.0]), np.array([0.0]), "bogus")
    e, _ = fw._fields(ap, np.array([1.2]), np.array([0.0]), "cos")
    assert np.isnan(e[0])  # invisible region
