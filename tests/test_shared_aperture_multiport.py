"""Phase 29 design point B: shared-aperture 7-port model (MODEL_ONLY): known-answer checks of the tensor fit, the cell current, the pattern-domain operator and the scan helpers."""

import numpy as np
import pytest

from src.hologram import seven_hex_cluster as c
from src.hologram import shared_aperture_multiport as sa
from src.hologram import tensor_cp_design as t
from src.hologram.tensor_leaky_cp import SupercellTable


def test_port_layout_centre_plus_ring_opposite_to_the_beam():
    p = sa.port_layout(30.0)
    assert p[0] == (0.0, 0.0, 0.0, 0.0) and len(p) == 7
    for x, y, th, ph in p[1:]:
        assert th == 30.0 and np.hypot(x, y) == pytest.approx(30.0)
        assert (np.degrees(np.arctan2(y, x)) - (ph + 180.0) + 180.0) % 360.0 - 180.0 == pytest.approx(0.0, abs=1e-9)
    q = sa.port_layout(30.0, ring_side="same")
    assert np.degrees(np.arctan2(q[2][1], q[2][0])) % 360.0 == pytest.approx(60.0)


def test_cells_keep_out_every_monopole_and_stay_in_the_hexagon():
    ports = sa.port_layout(24.0)
    x, y = sa.aperture_cells(60.0, ports)
    assert x.size > 100
    for px, py, _, _ in ports:
        assert np.hypot(x - px, y - py).min() >= sa.KEEPOUT_MM - 1e-9
    assert t.inside(sa.hexagon_m(60.0) * 1e3, x, y).all()


def test_tensor_current_is_a_reflection_of_the_incident_vector():
    rng = np.random.default_rng(3)
    s = rng.normal(size=(5, 2)) + 1j * rng.normal(size=(5, 2))
    axis = rng.uniform(0, np.pi, 5)
    kap = rng.uniform(1, 3, 5)
    j = sa.tensor_current(s, kap, axis)
    assert np.allclose(np.linalg.norm(j, axis=1), kap * np.linalg.norm(s, axis=1))  # |J| = kappa |s|
    assert np.allclose(
        sa.tensor_current(s, kap, np.zeros(5)), kap[:, None] * s * np.array([1, -1])
    )  # axis 0: (sx, -sy)
    # an incident LINEAR direction psi is reflected to 2a - psi
    psi = np.array([0.3])
    sl = np.stack([np.cos(psi), np.sin(psi)], axis=-1).astype(complex)
    jl = sa.tensor_current(sl, np.ones(1), np.array([0.7]))
    assert np.arctan2(jl[0, 1].real, jl[0, 0].real) == pytest.approx(2 * 0.7 - 0.3)


def test_fit_tensor_recovers_a_known_tensor_and_weights_select_the_port():
    rng = np.random.default_rng(5)
    n = 6
    s = rng.normal(size=(4, n, 2)) + 1j * rng.normal(size=(4, n, 2))
    p0, q0 = rng.normal(size=n), rng.normal(size=n)
    t_ = np.stack([p0 * s[..., 0] + q0 * s[..., 1], q0 * s[..., 0] - p0 * s[..., 1]], axis=-1)
    p, q, _ = sa.fit_tensor(s, t_)
    assert np.allclose(p, p0) and np.allclose(q, q0)
    # two conflicting targets: weight 1 on port 0 only -> port 0's tensor
    t2 = t_.copy()
    t2[1] *= -1.0
    om = np.zeros((4, n))
    om[0] = 1.0
    p, q, _ = sa.fit_tensor(s, t2, om)
    assert np.allclose(p, p0) and np.allclose(q, q0)


def test_pattern_domain_operator_equals_the_phase26_pattern_routine():
    green = c.green_function()
    rng = np.random.default_rng(7)
    ports = sa.port_layout(20.0)
    x, y = (a * 1e-3 for a in sa.aperture_cells(40.0, ports))
    s = rng.normal(size=(x.size, 2)) + 1j * rng.normal(size=(x.size, 2))
    p, q = rng.normal(size=x.size), rng.normal(size=x.size)
    setup = sa._ff_setup(green, x, y, theta_step=15.0, phi_step=45.0)
    m = sa._ff_operator(setup, s) @ np.concatenate([p, q])
    jx, jy = p * s[:, 0] + q * s[:, 1], q * s[:, 0] - p * s[:, 1]
    et, ep = t._vector_far_field(setup["th"], setup["ph"], x, y, jx, jy, t.F, green, t.PITCH, t.PITCH)
    nd = setup["th"].size
    assert np.allclose(m[:nd], sa.rhcp_field(et, ep))
    from src.pattern.vector_aperture_pattern import copol_component

    assert np.allclose(m[nd:], copol_component(et, ep, "lhcp"))


def test_mrc_weights_are_unit_norm_with_the_requested_number_of_active_ports():
    theta, phi = np.arange(0.0, 90.1, 3.0), np.arange(0.0, 360.0, 6.0)
    rng = np.random.default_rng(1)
    fields = [
        (rng.normal(size=(theta.size, phi.size)) + 0j, rng.normal(size=(theta.size, phi.size)) + 0j) for _ in range(7)
    ]
    w = sa.mrc_weights(fields, theta, phi, (30.0, 60.0), n_active=3)
    assert np.linalg.norm(w) == pytest.approx(1.0) and np.count_nonzero(w) == 3
    # maximum-ratio: no other unit-norm weights on the same ports give more RHCP field at the direction
    r = np.array([sa.rhcp_field(*f)[10, 10] for f in fields])
    assert abs(w @ r) == pytest.approx(np.linalg.norm(r[np.nonzero(w)[0]]))


def test_scan_summary_pass_logic():
    rows = [
        {"theta_deg": th, "peak_gain_dbi": 8.0 - dip, "pointing_error_deg": err, "sll_db": -10.0}
        for th, dip, err in ((0.0, 0.0, 1.0), (15.0, 2.5, 2.9), (45.0, 6.0, 20.0))
    ]
    s = sa.scan_summary(rows)
    assert s["pass_inner"] and s["band_30_45"]["max_pointing_error_deg"] == 20.0
    rows[1]["pointing_error_deg"] = 3.5
    assert not sa.scan_summary(rows)["pass_inner"]


@pytest.mark.slow
def test_small_shared_design_runs_and_ports_superpose():
    table = SupercellTable.from_files(t.TABLE_GLOB)
    d = sa.design_shared(table, 45.0, 18.0, 1000.0, 20000.0, gamma=4.0, n_iter=1, n_psi=120)
    assert len(d["J"]) == 7 and all(0.0 < e < 1.0 for e in d["eta_conv"])
    theta, phi = np.arange(0.0, 90.1, 5.0), np.arange(0.0, 360.0, 10.0)
    f = sa.port_fields(d, c.green_function(), theta, phi)
    for j, fj in enumerate(f):  # each port's field is normalised to LAUNCH * eta_conv of 1 W
        assert t.prad(*fj, theta, phi) == pytest.approx(t.LAUNCH * d["eta_conv"][j], rel=1e-9)
    w = np.array([1, 0, 0.5j, 0, 0, 0, 0], dtype=complex)
    et, ep = sa.combine(f, w)
    assert np.allclose(et, f[0][0] + 0.5j * f[2][0]) and np.allclose(ep, f[0][1] + 0.5j * f[2][1])
