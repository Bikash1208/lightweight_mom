"""Phase wrapping / unwrapping: mathematical properties, branch cut, and the sampling limit."""

import numpy as np
import pytest

from src.phase import (
    max_wrapped_step_deg,
    principal_phase_deg,
    unwrap_deg,
    wrap_deg,
    wrapped_difference_deg,
    wrapped_steps_deg,
)

RNG = np.random.default_rng(20240921)


# ----------------------------------------------------------------------------- wrapping


@pytest.mark.parametrize(
    "raw, expected",
    [
        (0.0, 0.0),
        (-0.0, 0.0),
        (90.0, 90.0),
        (-90.0, -90.0),
        (179.0, 179.0),
        (-179.0, -179.0),
        (180.0, 180.0),
        (-180.0, 180.0),
        (181.0, -179.0),
        (-181.0, 179.0),
        (360.0, 0.0),
        (-360.0, 0.0),
        (540.0, 180.0),
        (-540.0, 180.0),
        (720.5, 0.5),
        (-3600.25, -0.25),
    ],
)
def test_wrap_known_values(raw, expected):
    assert wrap_deg(raw) == pytest.approx(expected, abs=1e-12)


def test_wrap_lands_in_half_open_interval_and_differs_by_multiple_of_360():
    x = RNG.uniform(-1e5, 1e5, 20000)
    w = wrap_deg(x)
    assert np.all(w > -180.0) and np.all(w <= 180.0)
    k = (x - w) / 360.0
    assert np.allclose(k, np.round(k), atol=1e-8)


def test_wrap_is_idempotent_and_periodic():
    x = RNG.uniform(-2000, 2000, 5000)
    assert np.array_equal(wrap_deg(wrap_deg(x)), wrap_deg(x))
    assert np.allclose(wrap_deg(x + 360.0), wrap_deg(x), atol=1e-9)


def test_wrap_scalar_in_scalar_out_array_in_array_out():
    assert isinstance(wrap_deg(200.0), float)
    out = wrap_deg(np.array([[190.0, -190.0]]))
    assert out.shape == (1, 2) and np.allclose(out, [[-170.0, 170.0]])


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
def test_wrap_rejects_non_finite(bad):
    with pytest.raises(ValueError):
        wrap_deg(bad)
    with pytest.raises(ValueError):
        wrap_deg([1.0, bad])


def test_branch_cut_negative_real_axis_is_plus_180_for_both_signed_zeros():
    """arg(-1 + 0j) = +pi but arg(-1 - 0j) = -pi in numpy: the database must hold one spelling.

    The negative zero is built with `complex(-1.0, -0.0)`: the literal `-1 - 0j` yields +0.0 or -0.0 in the
    imaginary part depending on the Python version's int/complex arithmetic and constant folding."""
    neg_zero = complex(-1.0, -0.0)
    assert np.degrees(np.angle(neg_zero)) == -180.0  # documents the numpy behaviour
    assert principal_phase_deg(complex(-1.0, 0.0)) == 180.0
    assert principal_phase_deg(neg_zero) == 180.0
    assert principal_phase_deg(1e-3 * neg_zero) == 180.0


def test_principal_phase_of_axes_and_quadrants():
    z = np.array([1, 1j, -1, -1j, 1 + 1j, -1 + 1j, -1 - 1j, 1 - 1j], dtype=complex)
    assert np.allclose(principal_phase_deg(z), [0, 90, 180, -90, 45, 135, -135, -45])


def test_principal_phase_matches_exp_of_phase_and_ignores_magnitude():
    phi = RNG.uniform(-179.0, 179.0, 200)
    z = RNG.uniform(0.1, 3.0, 200) * np.exp(1j * np.radians(phi))
    assert np.allclose(principal_phase_deg(z), phi, atol=1e-10)


def test_principal_phase_rejects_non_finite():
    with pytest.raises(ValueError):
        principal_phase_deg(complex(np.nan, 0.0))


def test_wrapped_difference_takes_the_short_way_round_the_cut():
    assert wrapped_difference_deg(179.0, -179.0) == pytest.approx(-2.0)
    assert wrapped_difference_deg(-179.0, 179.0) == pytest.approx(2.0)
    assert wrapped_difference_deg(10.0, 350.0) == pytest.approx(20.0)
    assert wrapped_difference_deg(45.0, 45.0) == 0.0


def test_wrapped_difference_antisymmetric_except_at_exactly_180():
    a, b = RNG.uniform(-180, 180, 500), RNG.uniform(-180, 180, 500)
    assert np.allclose(wrapped_difference_deg(a, b), -wrapped_difference_deg(b, a), atol=1e-9)


# --------------------------------------------------------------------------- unwrapping


@pytest.mark.parametrize("step", [-170.0, -100.0, -30.0, 0.0, 45.0, 100.0, 170.0])
def test_unwrap_recovers_linear_ramp_when_step_below_180(step):
    true = 25.0 + step * np.arange(40)
    assert np.allclose(unwrap_deg(wrap_deg(true)), true, atol=1e-9)


@pytest.mark.parametrize("step", [190.0, 250.0, -300.0])
def test_unwrap_aliases_when_step_exceeds_180_documented_limit(step):
    """Above 180 deg per sample the answer is *wrong* and looks fine: the only guard is the
    step diagnostic, which cannot see the true step either."""
    true = step * np.arange(10)
    got = unwrap_deg(wrap_deg(true))
    assert not np.allclose(got, true)
    expected_step = wrap_deg(step)  # what the samples actually look like
    assert np.allclose(np.diff(got), expected_step)


def test_unwrap_smooth_function_recovers_continuous_phase_anchored_at_first_sample():
    t = np.linspace(0.0, 1.0, 400)
    true = -2000.0 * t + 300.0 * np.sin(6.0 * t)  # max step ~ 5 deg
    got = unwrap_deg(wrap_deg(true))
    assert np.allclose(got - got[0], true - true[0], atol=1e-8)
    assert got[0] == wrap_deg(true[0])


def test_unwrap_is_continuous_and_consistent_with_wrapping():
    t = np.linspace(0.0, 1.0, 300)
    true = 1500.0 * t**2
    got = unwrap_deg(wrap_deg(true))
    assert np.max(np.abs(np.diff(got))) < 180.0
    diff = (got - wrap_deg(true)) / 360.0  # a whole number of turns at every sample
    assert np.allclose(diff, np.round(diff), atol=1e-8)


def test_unwrap_leaves_already_continuous_input_unchanged():
    x = np.cumsum(RNG.uniform(-50, 50, 100))
    assert np.allclose(unwrap_deg(x), x, atol=1e-9)


def test_unwrap_axis_argument_operates_per_row_or_column():
    t = np.arange(30)
    rows = np.stack([-70.0 * t, 40.0 * t + 10.0, 0.0 * t])
    wrapped = wrap_deg(rows)
    assert np.allclose(unwrap_deg(wrapped, axis=1), rows - rows[:, :1] + wrapped[:, :1], atol=1e-9)
    assert np.allclose(unwrap_deg(wrapped.T, axis=0), unwrap_deg(wrapped, axis=1).T)


def test_unwrap_short_inputs_and_bad_arguments():
    assert unwrap_deg([12.0]).tolist() == [12.0]
    assert unwrap_deg(np.array([])).size == 0
    with pytest.raises(ValueError):
        unwrap_deg(5.0)
    with pytest.raises(ValueError):
        unwrap_deg([1.0, np.nan])


def test_unwrap_step_of_exactly_180_is_resolved_to_plus_180():
    assert np.allclose(unwrap_deg([0.0, 180.0, 0.0]), [0.0, 180.0, 360.0])
    assert np.allclose(unwrap_deg([0.0, -180.0]), [0.0, 180.0])


def test_wrapped_steps_and_max_step():
    p = np.array([170.0, -170.0, -130.0, 100.0])
    assert np.allclose(wrapped_steps_deg(p), [20.0, 40.0, -130.0])
    assert max_wrapped_step_deg(p) == pytest.approx(130.0)
    assert max_wrapped_step_deg([5.0]) == 0.0
    assert max_wrapped_step_deg(np.zeros((3, 1)), axis=1) == 0.0


def test_delay_is_negative_phase_and_unwraps_monotonically_decreasing():
    """exp(+jwt) convention: a pure delay tau gives S = exp(-j 2 pi f tau); unwrapped phase falls."""
    f = np.linspace(1e9, 2e9, 201)
    tau = 1.5e-9
    z = np.exp(-2j * np.pi * f * tau)
    u = unwrap_deg(principal_phase_deg(z))
    assert np.all(np.diff(u) < 0.0)
    assert u[-1] - u[0] == pytest.approx(-360.0 * tau * (f[-1] - f[0]), abs=1e-6)
