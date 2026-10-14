"""Phase 26: launch fraction of a pin on a grounded slab -- known-answer limits, power balance, the top-hat
disk against an independent far-field integral, convergence."""

import numpy as np
import pytest
from scipy.special import j1

from src.constants import ETA0
from src.dispersion.analytic_tm0 import tm0_beta
from src.dispersion.pin_launch import Launcher, PinLaunch, _current_profile
from src.frequency import free_space_wavenumber

F = 10.54e9
K0 = free_space_wavenumber(F)
H = 1.575e-3


def _direct_power(pl, top, a, h):
    """Radiated power of the same current system (pin + disk + PEC-ground images) from the free-space far
    field -- independent of the transmission-line machinery. eps_r = 1 only."""
    th = np.linspace(1e-4, np.pi / 2, 601)
    zg = -h
    zz = np.linspace(zg, top, 1501)
    cur, _, _ = _current_profile(zz, pl.launcher, F, 1.0, h, zg)
    itop = cur[-1]
    rho = np.linspace(0.0, a, 300) if a > 0 else None
    out = []
    for t in th:
        kz = K0 * np.cos(t)
        v = -np.sin(t) * (
            np.trapezoid(cur * np.exp(1j * kz * zz), zz) + np.trapezoid(cur * np.exp(1j * kz * (2 * zg - zz)), zz)
        )
        d = 0.0
        if a > 0:
            fj = np.trapezoid(itop * (1 - rho**2 / a**2) * j1(K0 * np.sin(t) * rho), rho)
            d = 1j * np.cos(t) * fj * (np.exp(1j * kz * top) - np.exp(1j * kz * (2 * zg - top)))
        out.append(v + d)
    a_th = np.array(out)
    return ETA0 * K0**2 / (32 * np.pi**2) * np.trapezoid(abs(a_th) ** 2 * np.sin(th), th) * 2 * np.pi


def test_short_dipole_over_ground_known_answer():
    length, h = 0.5e-3, 0.2e-3
    pl = PinLaunch(F, 1.0, h, Launcher(top_above_surface_m=length - h, pin_radius_m=0.1e-3), n_seg=12)
    moment = 0.5 * length  # sinusoidal current of an electrically tiny pin is a triangle
    assert pl.space_power() == pytest.approx(ETA0 * K0**2 * moment**2 / (6 * np.pi), rel=0.01)


@pytest.mark.parametrize("top,a", [(2.2e-3, 0.0), (2.2e-3, 1.5e-3), (4e-3, 3e-3)])
def test_pin_and_top_hat_match_the_free_space_far_field(top, a):
    pl = PinLaunch(F, 1.0, H, Launcher(top_above_surface_m=top, disk_radius_m=a), n_seg=30)
    assert pl.space_power() == pytest.approx(_direct_power(pl, top, a, H), rel=2e-3)


def test_pole_is_the_bare_tm0_constant_and_power_balances():
    pl = PinLaunch(F, 2.2, H, Launcher(top_above_surface_m=1.2e-3, disk_radius_m=1e-3), n_seg=16)
    assert pl.beta() == pytest.approx(tm0_beta(F, 2.2, H), rel=1e-9)
    r = pl.split()
    total = r["p_space_w_per_a2"] + r["p_sw_w_per_a2"]
    assert r["p_sw_w_per_a2"] > 0.0 and 0.0 < r["launch_fraction"] < 1.0
    assert pl.power_balance(tan_delta=1e-3, n_per_k0=2000) == pytest.approx(total, rel=0.02)


def test_thin_slab_limit_launches_nothing_and_no_slab_power_is_recovered():
    thin = 1e-3 * H
    pl = PinLaunch(F, 2.2, thin, Launcher(top_above_surface_m=2e-3), n_seg=12)
    ref = PinLaunch(F, 1.0, thin, Launcher(top_above_surface_m=2e-3), n_seg=12)
    r = pl.split()
    assert r["launch_fraction"] < 1e-3
    assert r["p_space_w_per_a2"] == pytest.approx(ref.space_power(), rel=0.01)


def test_launch_fraction_grows_with_slab_thickness():
    fr = [
        PinLaunch(F, 2.2, h, Launcher(top_above_surface_m=1.2e-3), n_seg=12).split()["launch_fraction"]
        for h in (0.1e-3, 0.4e-3, 0.8e-3, 1.575e-3)
    ]
    assert all(b > a for a, b in zip(fr, fr[1:], strict=False))


def test_converged_in_segments_and_angles():
    base = Launcher(top_above_surface_m=3e-3, disk_radius_m=1e-3)
    a = PinLaunch(F, 2.2, H, base, n_seg=12).split()["launch_fraction"]
    b = PinLaunch(F, 2.2, H, base, n_seg=48).split()["launch_fraction"]
    assert a == pytest.approx(b, rel=0.01)
    pl = PinLaunch(F, 2.2, H, base, n_seg=24)
    assert pl.space_power(n_theta=200) == pytest.approx(pl.space_power(n_theta=800), rel=2e-3)


def test_zero_radius_disk_is_no_disk():
    a = PinLaunch(F, 2.2, H, Launcher(top_above_surface_m=2e-3, disk_radius_m=0.0), n_seg=12).split()
    b = PinLaunch(F, 2.2, H, Launcher(top_above_surface_m=2e-3, disk_radius_m=1e-9), n_seg=12).split()
    assert a["launch_fraction"] == pytest.approx(b["launch_fraction"], rel=1e-4)


def test_loaded_slab_pole_is_the_sheet_model_beta_and_zero_sheet_is_unchanged():
    from src.dispersion.sheet_model import find_beta_from_sheet_admittance
    from src.materials import Substrate

    y = 1j * 4e-3
    launcher = Launcher(top_above_surface_m=1.2e-3)
    loaded = PinLaunch(F, 2.2, H, launcher, n_seg=16, y_sheet_s=y)
    ref = find_beta_from_sheet_admittance(y, Substrate("s", 2.2, 0.0, H), F).beta_rad_per_m
    assert loaded.beta() == pytest.approx(ref, rel=1e-8)
    bare = PinLaunch(F, 2.2, H, launcher, n_seg=16).split()
    zero = PinLaunch(F, 2.2, H, launcher, n_seg=16, y_sheet_s=0.0).split()
    assert bare == zero
    r = loaded.split()
    total = r["p_space_w_per_a2"] + r["p_sw_w_per_a2"]
    assert loaded.power_balance(tan_delta=1e-3, n_per_k0=2000) == pytest.approx(total, rel=0.03)
    assert r["beta_over_k0"] > bare["beta_over_k0"]
