"""Phase 26 C: rotating-slot CP hologram design march and currents."""

import numpy as np
import pytest

from src.hologram.leaky_forward import polygon_ray_distance
from src.hologram.tensor_leaky_cp import SupercellTable, cell_design, march_rays, modulation_k

K0 = 2 * np.pi * 10.54e9 / 299792458.0


def _table(alpha_scale=1.0):
    k = np.array([0.556, 0.741, 0.889, 1.111, 1.481])
    a = np.array([[0.4, 0.5, 0.6, 0.5, 0.7], [1.5, 1.8, 1.7, 1.6, 2.0], [3.0, 3.0, 3.0, 2.5, 4.2]]) * alpha_scale
    b = np.full_like(a, 1.06)
    return SupercellTable(((4.4, 1.2), (4.6, 1.6), (4.8, 2.0)), k, a, b, 6.0, 6.4)


def _hex(r=0.150):
    return np.array([[r * np.cos(np.radians(90 + 60 * i)), r * np.sin(np.radians(90 + 60 * i))] for i in range(6)])


def test_modulation_k_limits():
    th = np.radians(30.0)
    assert modulation_k(1.06, 0.0, th, 0.0) == pytest.approx(0.56)
    assert modulation_k(1.06, np.pi, th, 0.0) == pytest.approx(1.56)
    assert modulation_k(1.06, np.pi / 2, th, 0.0) == pytest.approx(np.hypot(1.06, 0.5))


def test_power_conservation_radiated_plus_edge_equals_launched():
    m = march_rays(_table(), _hex(), (0.0, 0.0), K0, np.radians(30), np.radians(225), 0.6, n_psi=72)
    # radiated per radian = int 2 alpha W d rho ; edge remainder = W_edge
    rad = np.sum(2.0 * m["alpha"] * m["w"], axis=1) * m["drho"]
    # first-order explicit march: agreement to the step size
    assert np.allclose(rad + m["w_edge"], 1.0, atol=0.01)


def test_conversion_is_capped_by_the_table_and_grows_with_available_alpha():
    kw = dict(vertices_m=_hex(), feed_xy_m=(0.0, 0.0), k0=K0, theta0_rad=np.radians(30), phi0_rad=np.radians(225), eta_design=0.9, n_psi=36)
    e1 = march_rays(_table(1.0), **kw)["eta_conv"]
    e2 = march_rays(_table(2.0), **kw)["eta_conv"]
    assert 0.0 < e1 < e2 < 1.0
    # no deeper slot than the table's: even a 0.9 request cannot exceed 1 - exp(-2 a_max rho_edge) per ray
    a_max = np.nanmax(_table(1.0).alpha)
    assert e1 <= np.mean(1 - np.exp(-2 * a_max * polygon_ray_distance(_hex(), (0, 0), np.linspace(-np.pi, np.pi, 36, endpoint=False)))) + 1e-6


def test_circular_decomposition_one_component_in_phase_the_other_evanescent():
    # chi = phi_ref - kt.r: J.(1,-j)/2 e^{j chi} has phase -kt.r; the other has wavevector 2 beta rho - kt
    th, ph = np.radians(30), np.radians(225)
    rng = np.random.default_rng(0)
    r = rng.uniform(0.02, 0.14, 50)
    a = rng.uniform(-np.pi, np.pi, 50)
    x, y = r * np.cos(a), r * np.sin(a)
    phi_ref = 1.06 * K0 * r
    kdotr = K0 * np.sin(th) * (x * np.cos(ph) + y * np.sin(ph))
    chi = phi_ref - kdotr
    j = np.exp(-1j * phi_ref)
    jx, jy = j * np.cos(chi), j * np.sin(chi)
    comp_a = (jx - 1j * jy) / 2.0  # coefficient of (1, -j)... projection onto (1, +j)/sqrt2 conj
    comp_a = 0.5 * (jx + 1j * jy)  # = e^{-j(phi_ref)} e^{j chi}/... magnitude 1/2 with phase chi - phi_ref
    assert np.allclose(np.angle(comp_a * np.exp(1j * kdotr)), 0.0, atol=1e-9)
    assert np.allclose(np.abs(comp_a), 0.5)


def test_cell_design_gives_beam_phase_and_unit_slot_axis_range():
    t = _table()
    m = march_rays(t, _hex(), (0.0, 0.0), K0, np.radians(30), np.radians(225), 0.7, n_psi=72)
    xs = np.array([0.03, -0.05, 0.0, 0.07])
    ys = np.array([0.0, 0.04, -0.06, 0.05])
    d = cell_design(m, t, xs, ys, (0.0, 0.0), K0, np.radians(30), np.radians(225))
    assert np.all((d["slot_angle_rad"] >= 0) & (d["slot_angle_rad"] < np.pi))
    assert set(d["slot_L_mm"]) <= {4.4, 4.6, 4.8}
    assert np.all(d["amp"] >= 0)
    # |J|^2 = amp^2 only for a circularly-phased unit direction vector: |cos|^2 + |sin|^2 = 1 for real chi
    assert np.allclose(np.abs(d["jx"]) ** 2 + np.abs(d["jy"]) ** 2, d["amp"] ** 2, rtol=1e-9)
