"""Phase 23 B: the cell-wise far-field decomposition equals the full vector pattern, and a one-cell
change updates the field by exactly that cell's contribution."""

import numpy as np
import pytest

from src.green_function import LayeredGreenFunction
from src.materials import Substrate
from src.pattern.cellwise_pattern import (
    cell_contributions,
    combined_field,
    theta_phi_direction_set,
    uv_direction_set,
    verify_against_vector_pattern,
)
from src.pattern.vector_aperture_pattern import _vector_far_field, aperture_current_from_scalar, copol_component

F = 10.54e9
P = 0.0064


@pytest.fixture(scope="module")
def green():
    return LayeredGreenFunction(Substrate("user", 2.2, 0.0009, 1.575e-3), F, backing="pec")


def _aperture(n=9, seed=3):
    rng = np.random.default_rng(seed)
    x = (np.arange(n) - (n - 1) / 2) * P
    xx, yy = np.meshgrid(x, x, indexing="xy")
    s = np.exp(1j * rng.uniform(-np.pi, np.pi, xx.shape)) * rng.uniform(0.3, 1.0, xx.shape)
    return xx, yy, s


@pytest.mark.parametrize("pol", ["rhcp", "lhcp", "hp", "vp", "slant_p45"])
def test_cellwise_field_matches_the_vector_pattern(green, pol):
    xx, yy, s = _aperture()
    feed = (0.4 * P, -0.7 * P)  # off-centre feed, not on a cell centre
    dirs = uv_direction_set(
        du_far=0.1, beams_deg=[(30.0, 45.0)], du_near=0.05, near_radius_deg=15.0, patch_radius_deg=3.0, patch_step=0.02
    )
    cw, direct = verify_against_vector_pattern(xx, yy, s, feed, dirs, F, green, P, P, pol)
    assert cw.shape == direct.shape == (dirs.size,)
    assert np.max(np.abs(cw - direct)) <= 1e-10 * np.max(np.abs(direct))


def test_theta_component_and_default_dtype(green):
    xx, yy, s = _aperture(6)
    dirs, _ = theta_phi_direction_set(10.0, 30.0)
    cw, direct = verify_against_vector_pattern(xx, yy, s, (0.0, 0.0), dirs, F, green, P, P, "theta")
    assert np.max(np.abs(cw - direct)) <= 1e-10 * np.max(np.abs(direct))


def test_single_cell_update_is_the_cells_own_contribution(green):
    xx, yy, s = _aperture(7)
    feed = (0.0, 0.0)
    dirs = uv_direction_set(du_far=0.15)
    c = cell_contributions(xx.ravel(), yy.ravel(), feed, dirs, F, green, P, P, "rhcp")
    field = s.ravel() @ c
    s2 = s.copy()
    s2[2, 5] = -s2[2, 5] * np.exp(1j * 0.7)  # change one cell
    jx, jy = aperture_current_from_scalar(xx, yy, s2, feed_xy_m=feed)
    et, ep = _vector_far_field(dirs.theta_deg, dirs.phi_deg, xx, yy, jx, jy, F, green, P, P)
    new_direct = copol_component(et, ep, "rhcp")
    updated = field + (s2[2, 5] - s[2, 5]) * c[2 * 7 + 5]  # O(D) update
    assert np.max(np.abs(updated - new_direct)) <= 1e-10 * np.max(np.abs(new_direct))


def test_complex64_storage_is_close(green):
    xx, yy, s = _aperture(7)
    dirs = uv_direction_set(du_far=0.15)
    c128 = cell_contributions(xx.ravel(), yy.ravel(), (0, 0), dirs, F, green, P, P, "rhcp")
    c64 = cell_contributions(xx.ravel(), yy.ravel(), (0, 0), dirs, F, green, P, P, "rhcp", dtype=np.complex64)
    assert c64.dtype == np.complex64
    assert np.max(np.abs(c64 - c128)) <= 1e-6 * np.max(np.abs(c128))


def test_combined_excitation_is_the_weighted_sum_of_port_patterns(green):
    xx, yy, s1 = _aperture(5, seed=1)
    _, _, s2 = _aperture(5, seed=2)
    dirs = uv_direction_set(du_far=0.2)
    f1 = verify_against_vector_pattern(xx, yy, s1, (0, 0), dirs, F, green, P, P, "rhcp")[1]
    f2 = verify_against_vector_pattern(xx + 0.1, yy, s2, (0.1, 0), dirs, F, green, P, P, "rhcp")[1]
    w = [1.0, 0.8 * np.exp(1j * np.radians(70.0))]
    both = combined_field([f1, f2], w)
    # one aperture holding both ports' cells with the excitation folded into the cell values
    x_all = np.concatenate([xx.ravel(), (xx + 0.1).ravel()])
    y_all = np.concatenate([yy.ravel(), yy.ravel()])
    c1 = cell_contributions(xx.ravel(), yy.ravel(), (0, 0), dirs, F, green, P, P, "rhcp")
    c2 = cell_contributions((xx + 0.1).ravel(), yy.ravel(), (0.1, 0), dirs, F, green, P, P, "rhcp")
    assert x_all.size == c1.shape[0] + c2.shape[0]
    total = (w[0] * s1.ravel()) @ c1 + (w[1] * s2.ravel()) @ c2
    assert np.max(np.abs(total - both)) <= 1e-12 * np.max(np.abs(both))


def test_direction_set_structure():
    beams = [(30.0, 45.0), (30.0, 225.0)]
    d = uv_direction_set(
        du_far=0.05, beams_deg=beams, du_near=0.02, near_radius_deg=15.0, patch_radius_deg=3.0, patch_step=0.01
    )
    assert np.all(np.hypot(d.u, d.v) <= 1.0 + 1e-9)
    assert set(np.unique(d.tag)) == {-1, 0, 1}
    # no direction twice
    key = np.round(d.u * 1e6).astype(np.int64) * 10_000_000 + np.round(d.v * 1e6).astype(np.int64)
    assert np.unique(key).size == d.size
    # patch points are within 3 deg of their beam; the lattice is denser there than far away
    for k, (th, ph) in enumerate(beams):
        sel = d.tag == k
        t0, p0 = np.radians(th), np.radians(ph)
        t, p = np.radians(d.theta_deg[sel]), np.radians(d.phi_deg[sel])
        ang = np.degrees(np.arccos(np.clip(np.cos(t) * np.cos(t0) + np.sin(t) * np.sin(t0) * np.cos(p - p0), -1, 1)))
        assert ang.max() <= 3.0 + 1e-6 and sel.sum() > 50
    assert uv_direction_set(du_far=0.05).size == pytest.approx(np.pi / 0.05**2, rel=0.05)
    held, shape = theta_phi_direction_set()
    assert shape == (181, 360) and held.size == 181 * 360
