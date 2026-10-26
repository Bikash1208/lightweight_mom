"""Phase 20 R4: state-library assembly (src/ris/state_library.py) on synthetic data with known answers."""

import numpy as np
import pytest

from src.ris.polarization import cartesian_reflection, circular_reflection, jones_to_te_tm
from src.ris.state_library import (
    assemble,
    correction_map,
    fill_normal_incidence,
    interpolate_correction,
    structure_frame_circular,
    te_tm_from_structure_circular,
)

F = np.linspace(8.0, 12.0, 5)
TH = np.array([0.0, 15.0, 30.0, 45.0])
PHI = np.arange(0.0, 360.0, 45.0)


def _jones(k):
    a = np.radians(45.0 * k)
    return 0.9 * np.array([[np.cos(2 * a), np.sin(2 * a)], [np.sin(2 * a), -np.cos(2 * a)]]) * np.exp(1j * 0.3)


def _grid_with_only_phi0_at_normal_incidence():
    s = np.zeros((len(F), len(TH), len(PHI), 4, 2, 2), dtype=complex)
    rng = np.random.default_rng(1)
    for it, th in enumerate(TH):
        for ip, phi in enumerate(PHI):
            for k in range(4):
                if th == 0.0 and phi != 0.0:
                    continue  # to be filled by the exact rotation
                s[:, it, ip, k] = jones_to_te_tm(_jones(k), np.radians(phi)) * (1 + 0.01 * th) + 0 * rng.normal()
    return s


def test_normal_incidence_azimuths_follow_from_the_exact_frame_rotation():
    s = fill_normal_incidence(_grid_with_only_phi0_at_normal_incidence(), TH, PHI)
    for ip, phi in enumerate(PHI):
        for k in range(4):
            assert np.allclose(s[2, 0, ip, k], jones_to_te_tm(_jones(k), np.radians(phi)), atol=1e-12)
            assert np.allclose(cartesian_reflection(s[2, 0, ip, k], np.radians(phi)), _jones(k), atol=1e-12)


def test_structure_frame_roundtrip_and_pb_phase_in_the_fixed_frame():
    s = fill_normal_incidence(_grid_with_only_phi0_at_normal_incidence(), TH, PHI)
    c = structure_frame_circular(s, PHI)
    assert np.allclose(te_tm_from_structure_circular(c, PHI), s, atol=1e-12)
    # at normal incidence the structure-frame circular matrix does not depend on the azimuth at all
    assert np.allclose(c[:, 0, :, :], c[:, 0:1, 0:1, :][:, 0, :, :], atol=1e-12)
    # and the handedness-preserving term of the state turned by alpha carries exp(+j 2 alpha) (R3 convention)
    r = c[0, 0, 0, :, 0, 0]
    steps = np.angle(r[1:] / r[:-1], deg=True) % 360.0
    assert np.allclose(steps, 90.0, atol=1e-9)


def test_correction_map_applies_a_known_constant_shift_everywhere_and_reports_the_mesh_step():
    s = fill_normal_incidence(_grid_with_only_phi0_at_normal_incidence(), TH, PHI)
    cf = np.linspace(8.0, 12.0, 9)
    pts = [(0.0, 0.0)] + [(t, p) for t in (15.0, 30.0, 45.0) for p in (0.0, 90.0)]
    shift = (0.02 + 0.03j) * np.ones((4, 2, 2))
    c60 = np.zeros((len(cf), len(pts), 4, 2, 2), dtype=complex)
    c120 = c60 + shift
    c90 = c60 + 0.5 * shift
    lib = assemble(F, TH, PHI, s, cf, pts, c60, c90, c120)
    from src.ris.state_library import R2_RICHARDSON

    expect = (1 + 0.5 * R2_RICHARDSON) * shift  # C_inf(p=2) - C60 with C120 = shift, C90 = shift / 2, C60 = 0
    assert np.allclose(lib["correction"], expect, atol=1e-12)
    assert np.allclose(lib["c_corr"], lib["c_raw"] + expect, atol=1e-12)
    from src.ris.state_library import R1_RICHARDSON

    assert np.allclose(
        lib["uncertainty"], (R1_RICHARDSON - R2_RICHARDSON) * 0.5 * abs(shift), atol=1e-12
    )  # |C_inf(p=1) - C_inf(p=2)|, zero phi spread
    assert np.allclose(lib["s_corr"], te_tm_from_structure_circular(lib["c_corr"], PHI), atol=1e-12)


def test_cubic_interpolation_follows_a_smooth_cubic_in_f_and_theta():
    cf = np.linspace(8.0, 12.0, 17)
    tn = np.array([0.0, 15.0, 30.0, 45.0])
    ff, tt = np.meshgrid(cf, tn, indexing="ij")
    base = (1 + 0.1 * ff - 0.02 * ff**2 + 0.001 * ff**3 + 0.01 * tt - 1e-4 * tt**2 + 1e-6 * tt**3).astype(complex)
    vals = np.broadcast_to(base[:, :, None, None, None], (*base.shape, 4, 2, 2)).copy()
    fg, tg = np.linspace(8.0, 12.0, 41), np.arange(0.0, 46.0, 5.0)
    out = interpolate_correction(cf, tn, vals, fg, tg)
    FF, TT = np.meshgrid(fg, tg, indexing="ij")
    expect = 1 + 0.1 * FF - 0.02 * FF**2 + 0.001 * FF**3 + 0.01 * TT - 1e-4 * TT**2 + 1e-6 * TT**3
    # scipy's 'cubic' tensor-product scheme is a smooth cubic interpolant, not exact on every cubic: 1.4e-5 here
    assert np.allclose(out[:, :, 0, 0, 0].real, expect, atol=1e-4)


def test_correction_map_phi_spread_is_the_disagreement_between_azimuths():
    cf = np.linspace(8.0, 12.0, 5)
    pts = [(15.0, 0.0), (15.0, 90.0)]
    c60 = np.zeros((5, 2, 4, 2, 2), dtype=complex)
    c120 = np.zeros_like(c60)
    c120[:, 1] = 0.1  # the second azimuth sees a correction of 0.1, the first none
    _, d, spread, _ = correction_map(cf, pts, c60, c60, c120)
    from src.ris.state_library import R2_RICHARDSON

    assert np.allclose(d, 0.05 * (1 + R2_RICHARDSON)) and np.allclose(spread, 0.1 * (1 + R2_RICHARDSON))


def test_pb_circular_check_of_the_helper_matrix():
    assert circular_reflection(-np.eye(2))[1, 0] == pytest.approx(-1.0)


def _synthetic_library():
    from src.ris.state_library import StateLibrary

    f, th, phi = np.linspace(8, 12, 41), np.arange(0.0, 46, 5.0), np.arange(0.0, 360, 15.0)
    FF, TT, PP = np.meshgrid(f, th, phi, indexing="ij")
    c = np.zeros((*FF.shape, 4, 2, 2), dtype=complex)
    for k in range(4):
        base = np.exp(1j * (np.radians(90 * k) + 0.4 * (FF - 10) + 0.002 * TT + 0.2 * np.cos(np.radians(PP))))
        c[:, :, :, k, 0, 0] = 0.9 * base
        c[:, :, :, k, 1, 0] = 0.15 * np.exp(1j * np.radians(2 * PP)) * (1 + 0.01 * TT)
        c[:, :, :, k, 0, 1] = 0.05
        c[:, :, :, k, 1, 1] = 0.02
    arrays = {
        "f_ghz": f,
        "theta_deg": th,
        "phi_deg": phi,
        "state": np.arange(1, 5),
        "c_corr": c,
        "uncertainty": np.abs(c) * 0.01,
    }
    return StateLibrary(arrays), (
        lambda ff, tt, pp, k: (
            0.9 * np.exp(1j * (np.radians(90 * k) + 0.4 * (ff - 10) + 0.002 * tt + 0.2 * np.cos(np.radians(pp))))
        )
    )


def test_interpolator_follows_a_smooth_field_and_is_periodic_in_phi():
    lib, truth = _synthetic_library()
    rng = np.random.default_rng(5)
    f, th, ph = rng.uniform(8, 12, 40), rng.uniform(0, 45, 40), rng.uniform(0, 360, 40)
    for k in range(4):
        c = lib.structure_circular(f, th, ph, k + 1)
        assert np.max(np.abs(c[:, 0, 0] - truth(f, th, ph, k))) < 2e-3
    a = lib.structure_circular(10.0, 20.0, 5.0, 2)
    b = lib.structure_circular(10.0, 20.0, 365.0, 2)
    assert np.allclose(a, b, atol=1e-12)
    assert np.allclose(
        lib.structure_circular(10.0, 20.0, 359.0, 1), lib.structure_circular(10.0, 20.0, -1.0, 1), atol=1e-12
    )  # across the seam


def test_interpolator_refuses_outside_the_grid_with_the_named_state():
    from src.ris.state_library import OutOfDomain

    lib, _ = _synthetic_library()
    with pytest.raises(OutOfDomain, match="NOT_AVAILABLE_FROM_CURRENT_DATA"):
        lib.structure_circular(12.5, 10.0, 0.0, 1)
    with pytest.raises(OutOfDomain):
        lib.structure_circular(10.0, 50.0, 0.0, 1)
    assert lib.label(45.0) == "EXTRAPOLATED_BEYOND_R3_TARGETS" and lib.label(30.0) == "within_R3_targets"


def test_interpolated_te_tm_matches_the_inverse_transform_of_the_interpolated_circular():
    lib, _ = _synthetic_library()
    c = lib.structure_circular(9.3, 12.0, 77.0, 3)[0]
    from src.ris.polarization import from_structure_frame, te_tm_from_circular

    assert np.allclose(
        lib.te_tm(9.3, 12.0, 77.0, 3)[0], te_tm_from_circular(from_structure_frame(c, np.radians(77.0))), atol=1e-14
    )


def test_interpolator_seam_boundaries_arrays_and_uncertainty():
    lib, truth = _synthetic_library()
    # a dense sweep across the phi seam (345..360 and 0..15) against the periodic truth: a wrong pad width or wrap would show here
    ph = np.concatenate([np.linspace(340.0, 360.0, 41), np.linspace(0.0, 20.0, 41)])
    f, th = np.full_like(ph, 9.37), np.full_like(ph, 21.3)
    c = lib.structure_circular(f, th, ph, 3)
    assert np.max(np.abs(c[:, 0, 0] - truth(f, th, ph, 2))) < 2e-3
    # 1e-9 slack at the grid edge is clipped onto the grid, never a plain scipy error
    edge = lib.structure_circular(
        8.0 - 5e-10, 45.0 + 5e-10, 100.0, 1
    )  # (the 12 GHz / 45 deg corner has a propagating order and is refused)
    assert np.all(np.isfinite(edge))
    # array f with scalar phi, and f / phi grids of different shapes: te_tm uses the same broadcast and point order as structure_circular
    f2 = np.linspace(8.5, 11.5, 7)[:, None]
    phi2 = np.array([10.0, 100.0, 250.0])
    ts = lib.te_tm(f2, 12.0, phi2, 2)
    assert ts.shape == (21, 2, 2)
    cs = lib.structure_circular(f2, 12.0, phi2, 2)
    from src.ris.polarization import from_structure_frame, te_tm_from_circular

    ff, pp = np.broadcast_arrays(f2, phi2)
    k = 8
    assert np.allclose(ts[k], te_tm_from_circular(from_structure_frame(cs[k], np.radians(pp.ravel()[k]))), atol=1e-13)
    # the real-valued uncertainty table goes through the single-part path
    u = lib.uncertainty(9.0, 10.0, 77.0, 1)
    assert u.shape == (1, 2, 2) and np.all(u >= 0)


def test_higher_order_margin_matches_the_brute_force_rule_and_the_known_onset():
    from src.ris.state_library import higher_order_margin

    # onset at phi = 0: f = c / (P (1 + sin theta)) = 11.7076 GHz at theta = 45 deg for P = 15 mm
    assert higher_order_margin(11.70, 45.0, 0.0) > 1.0 > higher_order_margin(11.72, 45.0, 0.0)
    assert higher_order_margin(10.0, 30.0, 0.0) > 1.0  # the design point has only the specular order
    rng = np.random.default_rng(2)
    f, th, ph = rng.uniform(8, 12, 500), rng.uniform(0, 45, 500), rng.uniform(0, 360, 500)
    brute = np.full(500, np.inf)
    a = 299.792458 / (f * 15.0)
    for m in range(-3, 4):
        for n in range(-3, 4):
            if (m, n) != (0, 0):
                brute = np.minimum(
                    brute,
                    np.hypot(
                        np.sin(np.radians(th)) * np.cos(np.radians(ph)) + a * m,
                        np.sin(np.radians(th)) * np.sin(np.radians(ph)) + a * n,
                    ),
                )
    assert np.all((higher_order_margin(f, th, ph) < 1.0) == (brute < 1.0))  # the first orders are enough on this grid


def test_interpolator_refuses_propagating_orders_and_flags_the_guard_band():
    from src.ris.state_library import OutOfDomain

    lib, _ = _synthetic_library()
    with pytest.raises(OutOfDomain, match="higher diffraction order"):
        lib.structure_circular(11.9, 45.0, 0.0, 1)
    assert lib.guard(11.65, 45.0, 0.0) and not lib.guard(10.0, 20.0, 0.0)
    lib.structure_circular(11.65, 45.0, 0.0, 1)  # below onset: still answered (flagged by guard())
