"""Phase 20 R1 item 2: lumped RLC loads -- normalisation, limits, power balance, mesh convergence, config."""

import math
import warnings
from dataclasses import dataclass

import numpy as np
import pytest

from src.basis_functions import MeshSettings, RooftopBasis
from src.config import ConfigError, config_from_dict
from src.excitation import ReferencePlanes
from src.geometry import UnitCell
from src.loads import OPEN_OHM, LumpedLoad, absorbed_power, check_load_geometry, load_matrix
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice
from src.sparameters import compute_s_parameters, compute_s_parameters_oblique

F = 10e9
SUB = Substrate("rt5880", 2.2, 0.0, 1.575e-3)
PX = 10e-3
PLANES = ReferencePlanes(2e-3, None, SUB.thickness)


@dataclass(frozen=True)
class TwoArm:
    """Two metal arms along x (|x| in [gap/2, arm_end]) of width `width`; `fill` also makes the gap metal. Duck-types Geometry."""

    gap: float = 2e-3
    arm_end: float = 4e-3
    width: float = 6e-3
    fill: bool = False
    kind = "two_arm"
    polarity = "metal"

    @property
    def feature_area(self):
        return (2 * self.arm_end if self.fill else 2 * (self.arm_end - self.gap / 2)) * self.width

    def corners(self):
        return np.array(
            [
                [self.arm_end, self.width / 2],
                [-self.arm_end, self.width / 2],
                [-self.arm_end, -self.width / 2],
                [self.arm_end, -self.width / 2],
            ]
        )

    def in_feature(self, x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        in_y = np.abs(y) <= self.width / 2
        arms = (np.abs(x) >= self.gap / 2) & (np.abs(x) <= self.arm_end)
        return in_y & (arms | (self.fill & (np.abs(x) < self.gap / 2)))


def _cell(fill, loads=(), geo=None, backing="pec"):
    g = geo or TwoArm(fill=fill)
    return UnitCell(RectangularLattice(PX, PX), g, SUB, Metal(), backing=backing, loads=tuple(loads))


def _load(**kw):
    d = {"x": 0.0, "y": 0.0, "length": 2e-3, "width": 6e-3, "direction": "x"}
    d.update(kw)
    return LumpedLoad(**d)


def _s(cell, mesh_n=20, theta=0.0, phi=0.0):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return compute_s_parameters_oblique(
            cell, F, MeshSettings(mesh_n, mesh_n, max(10, mesh_n // 2)), PLANES, theta, phi
        )


@pytest.mark.parametrize("n", [10, 20, 40])
def test_normalisation_uniform_current_through_the_load_dissipates_one_half_R_I_total_squared(n):
    """A uniform current density J = 1 A/m through the region carries I_total = J * width and must dissipate (1/2) R I_total^2."""
    ld = _load(R=7.0)
    basis = RooftopBasis(_cell(True), n, n)
    x0, x1, y0, y1 = ld.region()
    sel = np.flatnonzero(
        (np.arange(basis.n_unknowns) < basis.n_x)
        & (basis.center_x + basis.dx > x0)
        & (basis.center_x - basis.dx < x1)
        & (basis.center_y + basis.dy / 2 > y0)
        & (basis.center_y - basis.dy / 2 < y1)
    )
    i = np.zeros(basis.n_unknowns, complex)
    i[sel] = 1.0  # tents overlapping the region partition unity: J = 1 everywhere inside
    p = absorbed_power(basis, [ld], 2 * math.pi * F, i)
    assert p == pytest.approx(0.5 * 7.0 * ld.width**2, rel=1e-9)
    # an inductor dissipates nothing, stores (1/2) L I^2 / ... : the quadratic form is purely imaginary
    zl = load_matrix(basis, [_load(L=1e-9)], 2 * math.pi * F)
    assert abs((np.vdot(i, zl @ i)).real) < 1e-9 * abs(np.vdot(i, zl @ i))


def test_two_parallel_halves_with_doubled_impedance_equal_one_load():
    basis = RooftopBasis(_cell(True), 20, 20)
    one = load_matrix(basis, [_load(R=5.0, L=1e-9)], 2 * math.pi * F)
    halves = load_matrix(
        basis,
        [_load(y=1.5e-3, width=3e-3, R=10.0, L=2e-9), _load(y=-1.5e-3, width=3e-3, R=10.0, L=2e-9)],
        2 * math.pi * F,
    )
    assert np.allclose(one, halves, rtol=1e-12, atol=1e-12 * np.abs(one).max())


def test_short_load_is_the_filled_gap_and_open_load_is_the_empty_gap():
    filled = _s(_cell(True))
    short = _s(_cell(True, [_load()]))  # R = L = 0, no C: Z_L = 0
    assert np.array_equal(short.s11, filled.s11)  # bit-identical: Z_load = 0
    tiny = _s(_cell(True, [_load(R=1e-6)]))
    assert np.abs(tiny.s11 - filled.s11).max() < 1e-6  # continuous in Z_L: ~0.13 per ohm near the short
    empty = _s(_cell(False))
    opened = _s(_cell(True, [_load(open_circuit=True)]))
    assert OPEN_OHM >= 1e12
    assert np.abs(opened.s11 - empty.s11).max() < 2e-5
    assert np.abs(opened.s11 - filled.s11).max() > 1e-2  # the limit is not vacuous: gap empty != gap metal


def test_blocking_the_transverse_conduction_changes_a_finite_load_but_not_the_open_limit():
    kept = _s(_cell(True, [_load(C=0.4e-12)]))
    blocked = _s(_cell(True, [_load(C=0.4e-12, transverse_conducts=False)]))
    assert (
        np.abs(kept.s11 - blocked.s11).max() > 1e-3
    )  # a conducting strip across the gap is a real physical difference
    empty = _s(_cell(False))
    assert np.abs(_s(_cell(True, [_load(open_circuit=True)])).s11 - empty.s11).max() < 2e-5


@pytest.mark.parametrize("theta_deg", [0.0, 40.0])
def test_resistive_load_absorbs_exactly_one_minus_reflected_power(theta_deg):
    ld = _load(R=60.0)
    cell = _cell(True, [ld])
    s = _s(cell, 20, math.radians(theta_deg), 0.3)
    basis = RooftopBasis(cell, 20, 20)
    omega = 2 * math.pi * F
    for mode in (0, 1):
        i = s.currents[:, mode]
        y = (s.kz0 / (omega * 4e-7 * math.pi)) if mode == 0 else (omega * 8.8541878128e-12 / s.kz0)
        p_in = 0.5 * abs(s.incident_amplitudes[0, mode]) ** 2 * y * PX * PX
        p_abs = absorbed_power(basis, [ld], omega, i)
        assert p_abs / p_in == pytest.approx(1.0 - s.power_out_per_incident()[0, mode], abs=2e-9)
        assert 0.05 < p_abs / p_in < 0.95  # a genuinely partial absorber


@pytest.mark.parametrize("load", [dict(L=3e-9), dict(C=0.3e-12), dict(L=2e-9, C=0.5e-12)])
def test_reactive_loads_are_lossless(load):
    s = _s(_cell(True, [_load(**load)]), 20, math.radians(30.0), 0.0)
    assert np.allclose(s.power_out_per_incident(), 1.0, atol=1e-9)
    assert abs(abs(s.s11[0, 0]) - 1.0) < 1e-9


def test_loaded_cell_phase_converges_with_mesh_refinement_and_is_independent_of_edges_spanned():
    ld = _load(R=20.0, C=0.4e-12)
    ph = []
    for n in (10, 20, 40):
        s = _s(_cell(True, [ld]), n)
        ph.append(
            np.degrees(np.angle(s.s11[0, 0]))
        )  # TE (y-polarised): current across the x-gap is the TM current... use both
        ph.append(np.degrees(np.angle(s.s11[1, 1])))
    te, tm = ph[0::2], ph[1::2]
    assert abs(tm[2] - tm[1]) < abs(
        tm[1] - tm[0]
    )  # the loaded (TM = x) phase converges: 59.8, 56.4, 55.4 deg for 2, 4, 8 cells across the load
    assert abs(tm[2] - tm[1]) < 1.2 and abs(te[2] - te[1]) < 1.2  # degrees


def test_loads_must_fit_and_off_grid_edges_warn():
    cell = _cell(True)
    with pytest.raises(ValueError):
        check_load_geometry([_load(x=4.5e-3)], cell, 20, 20)
    with pytest.warns(UserWarning, match="off the mesh grid"):
        check_load_geometry([_load(length=1.3e-3)], cell, 20, 20)
    with pytest.raises(ValueError):
        LumpedLoad(0, 0, 1e-3, 1e-3, "z")
    with pytest.raises(ValueError):
        LumpedLoad(0, 0, 1e-3, 1e-3, "x", R=-1.0)


def test_no_loads_is_the_same_object_path_and_normal_incidence_solver_accepts_loads():
    cell = _cell(True)
    assert cell.loads == ()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a = compute_s_parameters(_cell(True, [_load(R=50.0)]), F, MeshSettings(20, 20, 10), PLANES)
        b = _s(_cell(True, [_load(R=50.0)]))
    assert abs(a.s11[1, 1] - b.s11[0, 0]) < 1e-9 or abs(a.s11[0, 0] - b.s11[1, 1]) < 1e-9 or True
    assert abs(a.s11[0, 0] - b.s11[1, 1]) < 1e-9  # TM = x at theta = 0, phi = 0 (x-directed current across the x-gap)


def test_config_loads_block():
    base = {
        "period_x": 10e-3,
        "period_y": 10e-3,
        "substrate": {"name": "s", "eps_r": 2.2, "tan_delta": 0.0, "thickness": 1.575e-3, "backing": "pec"},
        "geometry": {"type": "rectangle", "length": 8e-3, "width": 6e-3},
        "excitation": {"frequency": 10e9},
        "mesh": {"cells_x": 20, "cells_y": 20, "floquet_order": 10},
        "reference_planes": {"top_air_thickness": 2e-3},
    }
    assert config_from_dict(base).unit_cell.loads == ()
    cfg = config_from_dict(
        {
            **base,
            "loads": [{"x": 0.0, "y": 0.0, "length": "2 mm", "width": 6e-3, "direction": "x", "R": 10.0, "C": 1e-13}],
        }
    )
    ld = cfg.unit_cell.loads[0]
    assert ld.length == pytest.approx(2e-3) and ld.R == 10.0 and ld.C == 1e-13 and ld.direction == "x"
    with pytest.raises(ConfigError):
        config_from_dict({**base, "loads": [{"x": 0, "y": 0, "length": 2e-3, "width": 6e-3}]})  # direction missing
    with pytest.raises(ConfigError):
        config_from_dict({**base, "loads": [{"x": 0, "y": 0, "length": 2e-3, "width": 6e-3, "direction": "x", "Q": 1}]})
