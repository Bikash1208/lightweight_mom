import math

import numpy as np
import pytest

from src.geometry import Geometry, UnitCell
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice

PX, PY = 8e-3, 6e-3
LAT = RectangularLattice(PX, PY)


def test_reciprocal_lattice():
    assert LAT.reciprocal_x == pytest.approx(2 * math.pi / PX, rel=1e-15)
    assert LAT.reciprocal_y == pytest.approx(2 * math.pi / PY, rel=1e-15)
    assert LAT.cell_area == pytest.approx(PX * PY, rel=1e-15)
    # reciprocal_x * Px = 2*pi : one full phase cycle per period
    assert LAT.reciprocal_x * PX == pytest.approx(2 * math.pi, rel=1e-15)


@pytest.mark.parametrize("bad", [0.0, -8e-3, float("nan"), float("inf"), "8e-3", None, True])
def test_lattice_rejects_bad_periods(bad):
    with pytest.raises(ValueError):
        RectangularLattice(bad, PY)
    with pytest.raises(ValueError):
        RectangularLattice(PX, bad)


def test_floquet_wavenumber_spacing_and_symmetry_at_normal_incidence():
    m = np.arange(-4, 5)
    kx, ky = LAT.transverse_wavenumbers(0.0, 0.0, m, m)
    assert np.allclose(np.diff(kx), LAT.reciprocal_x, rtol=1e-14)
    assert np.allclose(np.diff(ky), LAT.reciprocal_y, rtol=1e-14)
    assert np.allclose(kx, -kx[::-1], atol=1e-9)  # kx_{-m} = -kx_m when kx0 = 0
    assert kx[4] == 0.0 and ky[4] == 0.0  # (0,0) mode is the incident direction


def test_floquet_wavenumbers_oblique_offset():
    kx0, ky0 = 123.4, -56.7
    kx, ky = LAT.transverse_wavenumbers(kx0, ky0, 2, -3)
    assert kx == pytest.approx(kx0 + 2 * 2 * math.pi / PX, rel=1e-14)
    assert ky == pytest.approx(ky0 - 3 * 2 * math.pi / PY, rel=1e-14)


def test_floquet_indices_must_be_integers():
    with pytest.raises(ValueError):
        LAT.transverse_wavenumbers(0.0, 0.0, 0.5, 0)
    with pytest.raises(ValueError):
        LAT.transverse_wavenumbers(0.0, 0.0, 0, np.array([0.0, 1.0]))


@pytest.mark.parametrize("kx0, ky0", [(0.0, 0.0), (137.0, -42.0), (-511.0, 300.0)])
def test_floquet_harmonics_are_quasi_periodic(kx0, ky0):
    """psi_mn(x+Px, y) = psi_mn(x, y) * exp(-j kx0 Px) for every harmonic (m, n)."""
    rng = np.random.default_rng(1)
    x = rng.uniform(-PX, PX, 20)
    y = rng.uniform(-PY, PY, 20)
    phase_x, phase_y = LAT.bloch_phase(kx0, ky0)
    for m in range(-3, 4):
        for n in range(-3, 4):
            kx, ky = LAT.transverse_wavenumbers(kx0, ky0, m, n)

            def psi(xx, yy):
                return np.exp(-1j * (kx * xx + ky * yy))

            assert np.allclose(psi(x + PX, y), psi(x, y) * phase_x, atol=1e-12)
            assert np.allclose(psi(x, y + PY), psi(x, y) * phase_y, atol=1e-12)


def test_bloch_phase_unit_modulus_and_normal_incidence():
    px, py = LAT.bloch_phase(0.0, 0.0)
    assert px == 1.0 and py == 1.0
    px, py = LAT.bloch_phase(211.0, -97.0)
    assert abs(px) == pytest.approx(1.0, rel=1e-14)
    assert abs(py) == pytest.approx(1.0, rel=1e-14)
    assert px == pytest.approx(np.exp(-1j * 211.0 * PX), rel=1e-14)


def test_floquet_harmonics_orthogonal_over_the_cell():
    """(1/Px) * integral of psi_m * conj(psi_m') dx = delta_mm' (uniform sampling is exact)."""
    kx0 = 173.0
    n_pts = 64
    x = (np.arange(n_pts) + 0.5) / n_pts * PX
    for m in range(-3, 4):
        for mp in range(-3, 4):
            kxm, _ = LAT.transverse_wavenumbers(kx0, 0.0, m, 0)
            kxp, _ = LAT.transverse_wavenumbers(kx0, 0.0, mp, 0)
            inner = np.mean(np.exp(-1j * kxm * x) * np.conj(np.exp(-1j * kxp * x)))
            assert inner == pytest.approx(1.0 if m == mp else 0.0, abs=1e-12)


def test_mode_grid():
    m, n = LAT.mode_grid(2)
    assert m.shape == n.shape == (5, 5)
    assert m[0, 0] == -2 and m[-1, 0] == 2 and n[0, 0] == -2 and n[0, -1] == 2
    assert np.issubdtype(m.dtype, np.integer)
    m0, n0 = LAT.mode_grid(0)
    assert m0.shape == (1, 1) and m0[0, 0] == 0 and n0[0, 0] == 0
    for bad in (-1, 1.5, "2"):
        with pytest.raises(ValueError):
            LAT.mode_grid(bad)


@pytest.mark.parametrize("kind", ["rectangle", "slot"])
def test_pattern_is_exactly_periodic(kind):
    sub = Substrate("t", 2.2, 0.0009, 1.575e-3)
    geo = Geometry(kind, 4.5e-3, 1.5e-3, center_x=0.4e-3, center_y=-0.3e-3, rotation=math.radians(25))
    cell = UnitCell(LAT, geo, sub, Metal())
    rng = np.random.default_rng(7)
    x = rng.uniform(-2 * PX, 2 * PX, 5000)
    y = rng.uniform(-2 * PY, 2 * PY, 5000)
    base = cell.is_metal(x, y)
    for i, j in [(1, 0), (0, 1), (-3, 2), (5, -4)]:
        assert np.array_equal(cell.is_metal(x + i * PX, y + j * PY), base)
