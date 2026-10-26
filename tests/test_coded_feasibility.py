"""Phase 28 item 2: quantisers on known inputs, the anchor to the Phase 26 single-hexagon prediction, contrast / bias-budget bookkeeping."""

import json
from pathlib import Path

import numpy as np
import pytest

from src.dispersion.leakage_model import leakage_alpha_first_order, leakage_alpha_first_order_tensor
from src.hologram import coded_feasibility as cf
from src.hologram import tensor_cp_design as t26
from src.hologram.tensor_leaky_cp import SupercellTable
from src.ris.diodes import bias_current_budget

ROOT = Path(__file__).resolve().parents[1]


# ---------------------------------------------------------------------------------------------- quantisers
def test_uniform_levels_are_truly_2_to_the_bits_levels_including_the_end_points():
    assert list(cf.uniform_levels(1)) == [-1.0, 1.0]
    np.testing.assert_allclose(cf.uniform_levels(2), [-1.0, -1.0 / 3.0, 1.0 / 3.0, 1.0])
    _, q = cf.quantise_nearest(np.linspace(-1, 1, 10001), cf.uniform_levels(2))
    assert len(np.unique(q)) == 4  # a dense ramp uses exactly four values, not three


@pytest.mark.parametrize(
    "value, expected_level",
    [
        (0.0, 1.0 / 3.0),
        (-1e-9, -1.0 / 3.0),
        (0.2, 1.0 / 3.0),
        (0.65, 1.0 / 3.0),
        (0.67, 1.0),
        (-0.67, -1.0),
        (-0.5, -1.0 / 3.0),
        (5.0, 1.0),
        (-5.0, -1.0),
    ],
)
def test_nearest_level_on_known_values_2bit(value, expected_level):
    idx, q = cf.quantise_nearest(value, cf.uniform_levels(2))
    assert q == pytest.approx(expected_level)


def test_nearest_level_1bit_is_the_sign_with_ties_to_the_lower_level():
    _, q = cf.quantise_nearest(np.array([-0.3, 0.3, 0.0]), cf.uniform_levels(1))
    assert list(q) == [-1.0, 1.0, -1.0]


@pytest.mark.parametrize(
    "angle, bits, expected_index",
    [
        (0.1, 2, 0),
        (2 * np.pi - 0.1, 2, 0),  # wraps to level 0, not to the 3 pi/2 level
        (-0.1, 2, 0),
        (np.pi / 4, 2, 1),  # exact half-way goes to the higher level
        (np.pi / 4 - 1e-9, 2, 0),
        (3.5, 2, 2),  # pi (3.14) is nearer than 3 pi/2 (4.71)
        (2 * np.pi - np.pi / 4 + 1e-9, 2, 0),
        (2 * np.pi - np.pi / 4 - 1e-9, 2, 3),
        (np.pi - 0.1, 1, 1),
        (np.pi + 0.1, 1, 1),
        (3 * np.pi / 2 + 0.1, 1, 0),  # 1-bit: levels 0 and pi, 3 pi/2 + eps is nearer 2 pi = 0
    ],
)
def test_circular_quantiser_wraps_correctly(angle, bits, expected_index):
    idx, q = cf.quantise_circular(angle, 2**bits, 2 * np.pi)
    assert int(idx) == expected_index
    assert 0.0 <= float(q) < 2 * np.pi


def test_circular_quantiser_error_is_at_most_half_a_step_on_the_circle():
    a = np.linspace(-10, 10, 20001)
    for bits in (1, 2, 3):
        _, q = cf.quantise_circular(a, 2**bits, np.pi)  # the slot-axis period is pi
        err = np.abs(((a - q + np.pi / 2) % np.pi) - np.pi / 2)
        assert err.max() <= np.pi / 2**bits / 2 + 1e-12


def test_bounded_db_reports_a_bound_not_a_clamped_number():
    assert cf.bounded_db(12.3)["value"] == pytest.approx(12.3)
    over = cf.bounded_db(240.0)
    assert over["value"] is None and over["text"] == "> 40 dB"
    assert cf.bounded_db(float("inf"))["value"] is None


# ---------------------------------------------------------------------------------------------- geometry
def test_hexagon_and_lattice_are_the_phase26_ones():
    np.testing.assert_allclose(cf.hexagon_m(0.150), t26.hexagon())
    x, y = cf.lattice(t26.hexagon(), 6.4e-3)
    x0, y0 = t26.lattice_cells(t26.hexagon())
    assert x.size == 1430 and np.array_equal(x, x0) and np.array_equal(y, y0)


# ---------------------------------------------------------------------------------------------- the anchor
@pytest.fixture(scope="module")
def anchor():
    table = SupercellTable.from_files(t26.TABLE_GLOB)
    poly = cf.hexagon_m(0.150)
    x, y, d, march = cf.design_cells(table, poly, 6.4e-3, 30.0, 225.0)
    return table, x, y, d, march


def test_unquantised_cp_currents_equal_the_phase26_currents(anchor):
    _, _, _, d, _ = anchor
    _, axis = cf.cp_states(d, None)
    jx, jy = cf.cp_currents(d, axis)
    np.testing.assert_allclose(jx, d["jx"], rtol=0, atol=1e-12 * np.abs(d["jx"]).max())
    np.testing.assert_allclose(jy, d["jy"], rtol=0, atol=1e-12 * np.abs(d["jy"]).max())


def test_anchor_no_quantisation_reproduces_the_phase26_single_hexagon_prediction(anchor):
    """Tolerances: D_ap and SLL 1e-6 dB, AR / XPD 1e-6 dB, eta_conv 1e-12 (same code path, same grid: the refinement is reported separately)."""
    _, x, y, d, march = anchor
    ref = json.loads((ROOT / "results/phase26/tensor_cp_hexagon.json").read_text(encoding="utf-8"))["selected"]
    m26 = ref["metrics_ap"]
    got = cf.analyse(d["jx"], d["jy"], x, y, cf.green_function(), 6.4e-3, "cp", (30.0, 225.0))
    assert march["eta_conv"] == pytest.approx(ref["eta_conv"], abs=1e-12)
    assert got["d_ap_coarse_grid_dbi"] == pytest.approx(m26["directivity_ap_rhcp_dbi"], abs=1e-6)
    assert cf.launch_gain_dbi(got["d_ap_coarse_grid_dbi"], march["eta_conv"]) == pytest.approx(
        ref["gain_rhcp_dbi_launch_bare"], abs=1e-6
    )
    assert got["ar_at_coarse_peak_db"] == pytest.approx(m26["ar_at_peak_db"], abs=1e-6)
    assert got["xpd_at_coarse_peak_db"] == pytest.approx(m26["xpd_at_peak_db"], abs=1e-6)
    assert got["sll_coarse_grid_db"] == pytest.approx(m26["rhcp_sll_db"], abs=1e-6)
    assert got["rhcp_fraction"] == pytest.approx(m26["rhcp_fraction"], abs=1e-9)
    # the 0.25 deg refinement can only raise the peak (by less than the coarse-grid scalloping)
    assert 0.0 <= got["d_ap_dbi"] - got["d_ap_coarse_grid_dbi"] < 0.3
    assert got["pointing_error_deg"] < 1.0


# ---------------------------------------------------------------------------------------------- quantised cases behave as they should
def test_cp_quantisation_never_beats_the_continuous_design_and_2bit_beats_1bit(anchor):
    _, x, y, d, _ = anchor
    g = cf.green_function()
    res = {}
    for bits in (None, 2, 1):
        _, axis = cf.cp_states(d, bits)
        jx, jy = cf.cp_currents(d, axis)
        res[bits] = cf.analyse(jx, jy, x, y, g, 6.4e-3, "cp", (30.0, 225.0))
    assert res[None]["d_ap_dbi"] >= res[2]["d_ap_dbi"] >= res[1]["d_ap_dbi"] - 1e-9
    assert res[2]["d_ap_dbi"] > res[1]["d_ap_dbi"] + 0.1
    # polarisation purity degrades with quantisation (the continuous design is the cleanest) and an unbounded value is never a number
    assert res[None]["xpd_at_peak"]["value"] > res[1]["xpd_at_peak"]["value"]
    assert res[None]["xpd_at_peak"]["value"] > 10.0


def test_lp_levels_used_and_quantisation_ordering(anchor):
    _, x, y, d, _ = anchor
    m = cf.lp_modulation(d, x, y, 30.0, 225.0)
    assert np.all(np.abs(m) <= 1.0 + 1e-12)
    for bits in (1, 2):
        idx, mq = cf.lp_states(m, bits)
        assert len(np.unique(idx)) == 2**bits and len(np.unique(mq)) == 2**bits
    g = cf.green_function()
    dd = {}
    for bits in (None, 2, 1):
        _, mq = cf.lp_states(m, bits)
        jx, jy = cf.lp_currents(d, mq)
        dd[bits] = cf.analyse(jx, jy, x, y, g, 6.4e-3, "lp", (30.0, 225.0))["d_ap_dbi"]
    assert dd[None] >= dd[2] - 1e-9 and dd[2] >= dd[1] - 1e-9


def test_first_harmonic_gain_of_the_levels():
    assert cf.first_harmonic_gain(cf.uniform_levels(1)) == pytest.approx(4 / np.pi, abs=1e-3)
    g2 = cf.first_harmonic_gain(cf.uniform_levels(2))
    assert 1.0 < g2 < 4 / np.pi


# ---------------------------------------------------------------------------------------------- contrast
def test_contrast_inversion_reproduces_the_requested_alpha_and_is_quadratic():
    c = cf.contrast_from_alpha(2.0, 1.06, 30.0, bits=1)
    beta = c["beta_over_k0_mean"] * cf.free_space_wavenumber(cf.F)
    ybar = 1j * c["y_mean_ms"] * 1e-3
    th = np.radians(30.0)
    a_t = leakage_alpha_first_order_tensor(beta, ybar, c["tensor"]["delta_ms"] * 1e-3, th, cf.SUB, cf.F)
    assert a_t == pytest.approx(2.0, rel=1e-9)
    db = c["scalar"]["d_b_full_swing_ms"] * 1e-3
    c1 = db * c["scalar"]["first_harmonic_gain_g"] / 4.0
    a_s = leakage_alpha_first_order(beta, ybar, c1, th, np.radians(c["scalar"]["worst_zeta_deg"]), cf.SUB, cf.F)
    assert a_s == pytest.approx(2.0, rel=1e-9)
    c4 = cf.contrast_from_alpha(8.0, 1.06, 30.0, bits=1)
    assert c4["tensor"]["delta_ms"] == pytest.approx(2 * c["tensor"]["delta_ms"], rel=1e-9)  # alpha x4 -> delta x2
    assert c["beta_over_k0_mean"] == pytest.approx(1.06, abs=1e-6)
    assert c["scalar"]["dx_switch_ohm"] > 0 and c["tensor"]["dx_aniso_ohm"] > 0


# ---------------------------------------------------------------------------------------------- bias budget
def test_bias_budget_uses_one_dc_power_method_and_matches_the_repo_diode_budget():
    b = cf.bias_budget(1430, 2, 0.5, i_forward_a=10e-3, v_rail_v=3.3)
    assert b["n_diodes"] == 2860 and b["independent_bias_lines"] == 2860
    _, p_repo = bias_current_budget(
        2860, 0.5, i_forward=10e-3, v_forward=3.3
    )  # same formula with the rail instead of V_F
    assert b["dc_power_expected_w"] == pytest.approx(p_repo)
    assert b["dc_power_all_on_w"] == pytest.approx(2860 * 10e-3 * 3.3)
    assert b["shift_register"]["outputs_per_chip_8"]["chips"] == 358  # ceil(2860 / 8)
    assert b["bias_line_sharing"]["cells_per_line_4"] == 715


def test_state_bit_stats():
    assert cf.state_bit_stats(np.array([0, 1, 2, 3]), 2) == pytest.approx(0.5)
    assert cf.state_bit_stats(np.array([3, 3]), 2) == 1.0
    assert cf.state_bit_stats(np.array([0, 0, 0]), 1) == 0.0


def test_capped_table_limits_alpha_and_nothing_else():
    table = SupercellTable.from_files(t26.TABLE_GLOB)
    cap = cf.capped_table(table, 1.5)
    assert (
        np.nanmax(cap.alpha) <= 1.5
        and np.array_equal(cap.beta, table.beta, equal_nan=True)
        and cap.slots == table.slots
    )


def test_tensor_principal_axes_differ_in_reactance_and_beta_and_lpx_current_is_along_x(anchor):
    c = cf.contrast_from_alpha(2.0, 1.06, 30.0, bits=1)
    assert (
        c["tensor"]["dx_aniso_ohm"] > 10.0 and c["tensor"]["d_beta_over_k0_principal_axes"] > 1e-3
    )  # a susceptance (reactive) anisotropy, not a conductance
    _, _, _, d, _ = anchor
    jx, jy = cf.lpx_currents(d, np.ones(d["amp"].size))
    assert np.all(jy == 0) and np.allclose(np.abs(jx), d["amp"])


def test_lp_scalar_radial_current_has_the_radial_direction(anchor):
    _, _, _, d, _ = anchor
    jx, jy = cf.lp_currents(d, np.ones(d["amp"].size))
    np.testing.assert_allclose(
        jy * np.cos(d["phi"]), jx * np.sin(d["phi"]), atol=1e-12 * np.abs(jx).max()
    )  # J parallel to (cos phi, sin phi)
