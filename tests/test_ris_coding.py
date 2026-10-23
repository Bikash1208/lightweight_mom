"""Phase 20 R5-3: coding synthesis -- required phase, state choice on the real responses, determinism, bandwidth helper."""

import math

import numpy as np
import pytest

from src.ris import coding as cd
from src.ris import forward as fw
from src.ris.calibrated_library import CalibratedLibrary

F = 10.5


@pytest.fixture(scope="module")
def cal():
    return CalibratedLibrary.load()


def test_required_phase_is_the_steering_gradient():
    psi = cd.required_phase(4, 6, 15e-3, F, 30.0, 0.0)
    k0 = 2 * math.pi * F * 1e9 / fw.C0
    assert psi[0, 1] - psi[0, 0] == pytest.approx(-k0 * math.sin(math.radians(30.0)) * 15e-3)
    assert np.allclose(psi[1:] - psi[:-1], 0.0)  # phi_s = 0: no variation along y
    psi90 = cd.required_phase(4, 6, 15e-3, F, 30.0, 90.0)
    assert np.allclose(psi90[:, 1:] - psi90[:, :-1], 0.0, atol=1e-12)


def test_normal_incidence_broadside_uses_one_state_and_is_the_uniform_aperture(cal):
    c = cd.design(cal, "plain", F, 0.0, 0.0)
    assert len(np.unique(c.states)) == 1
    m, _ = cd.evaluate(cal, c)
    assert m.pointing_error_deg < 0.2


def test_nearest_state_maximises_the_coherent_sum_over_all_choices(cal):
    table, _ = cd._state_table(cal, F, 0.0, 0.0)
    psi = cd.required_phase(6, 6, 15e-3, F, 25.0, 0.0)
    st = cd.nearest_states(table, psi)
    rng = np.random.default_rng(0)

    def coherent(codes):  # best-offset coherent sum toward the target
        return abs(np.sum(table[codes - 1] * np.exp(-1j * psi)))

    assert coherent(st) >= max(coherent(rng.integers(1, 5, st.shape)) for _ in range(200))
    assert coherent(st) <= st.size  # |C| <= 1 per cell


def test_designs_are_deterministic_and_dither_depends_on_the_seed(cal):
    a = cd.design(cal, "dither", F, 30.0, 0.0, seed=1)
    b = cd.design(cal, "dither", F, 30.0, 0.0, seed=1)
    c = cd.design(cal, "dither", F, 30.0, 0.0, seed=2)
    assert np.array_equal(a.states, b.states) and not np.array_equal(a.states, c.states)
    assert a.seed == 1 and "LOCAL_RESPONSE_AT_INCIDENT_ANGLE_ASSUMED" in a.notes


def test_real_coding_beam_points_near_the_target_and_cos_scan_loss(cal):
    c = cd.design(cal, "plain", F, 30.0, 0.0)
    m, gt = cd.evaluate(cal, c)
    assert m.pointing_error_deg < 3.0
    ref = cd.ideal_references(30.0, 0.0, F)["ideal_continuous"]
    assert gt < ref.peak_gain_dbi  # the real 2-bit coding with diode loss is below the ideal continuous aperture
    assert m.peak_gain_dbi > ref.peak_gain_dbi - 3.5


def test_refinement_never_lowers_its_own_objective(cal):
    c = cd.design(cal, "refined", F, 50.0, 0.0)
    assert c.log["objective_after"] >= c.log["objective_before"] - 1e-9
    again = cd.design(cal, "refined", F, 50.0, 0.0)
    assert np.array_equal(c.states, again.states)


def test_gain_bandwidth_interpolates_and_flags_open_ends():
    f = np.arange(10.0, 11.01, 0.1)
    g = 30.0 - 4.0 * ((f - 10.5) / 0.5) ** 2
    bw = cd.gain_bandwidth(f, g, drops=(1.0, 3.0))
    assert bw[1.0]["f_lo"] == pytest.approx(10.25, abs=0.02) and bw[1.0]["f_hi"] == pytest.approx(10.75, abs=0.02)
    assert not bw[1.0]["open_low"]
    assert bw[3.0]["open_low"] is False and bw[3.0]["width"] == pytest.approx(0.866, abs=0.03)
    flat = cd.gain_bandwidth(f, np.full_like(f, 20.0))
    assert flat[1.0]["open_low"] and flat[1.0]["open_high"]


def test_grating_onset_matches_the_recorded_limits():
    assert cd.grating_onset_deg(10.5) == pytest.approx(64.7, abs=0.2)  # 60 deg is grating-free at 10.5 GHz
    assert cd.grating_onset_deg(11.0) < 60.0
    assert cd.grating_onset_deg(10.71) == pytest.approx(60.0, abs=0.3)
