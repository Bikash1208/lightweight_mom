"""Phase 30 S1: direct MoM leaky-mode solve of a rotating-axis supercell (`src/dispersion/leaky_supercell_mom.py`).

Mathematical / physical properties, not "it runs": the FFT-filled complex-analytic Z equals the production assembler; Z is symmetric
(reciprocity) and analytic in kx0 and f; a supercell of identical cells has the single cell's mode (no leakage, lossless -> real);
leakage needs the modulation (rotating slots) and a radiating space harmonic; loss only adds attenuation.
"""

import math

import numpy as np

from src.dispersion import leaky_supercell_mom as L
from src.green_function import LayeredGreenFunction
from src.materials import Substrate
from src.mom_matrix import assemble_impedance_matrix

H = 0.635e-3
GEOM = {"major": 2.6e-3, "slot_length": 1.6e-3, "slot_width": 0.5e-3}  # gap 0.6 mm: resolved at 12 cells / 3.2 mm


def _sub(tand=0.0):
    return Substrate("ro", 11.2, tand, H)


def _peak_mode(basis, sub, kx0, f_lo, f_hi, order_x, order_y, n=61):
    """Real-frequency mode of a bound (lossless) problem: largest probe peak, refined by the eigenvalue secant."""
    v, w = L.probe_vectors(basis.n_unknowns)
    fs = np.linspace(f_lo, f_hi, n)
    g = np.array(
        [abs(L.probe_function(L.assemble_z_complex(basis, sub, f * 1e9, kx0, 0.0, order_x, order_y), v, w)) for f in fs]
    )
    i = int(np.argmax(g[1:-1])) + 1

    def zf(p):
        return L.assemble_z_complex(basis, sub, p * 1e9, kx0, 0.0, order_x, order_y)

    p, _, _ = L.find_mode_eig(zf, fs[i], fs[i] + 0.002)
    return p


def test_fft_filled_z_equals_the_production_assembler():
    sub = _sub(0.0022)
    _, basis = L.build_supercell(L.s1_angles(3), 8, sub)
    f, kx0, ky0, order = 10.0e9, 37.0, 21.0, 9
    z = L.assemble_z_complex(basis, sub, f, kx0, ky0, order, order)
    z0 = assemble_impedance_matrix(basis, LayeredGreenFunction(sub, f, "pec"), kx0, ky0, order)
    assert np.abs(z - z0).max() / np.abs(z0).max() < 1e-12


def test_z_is_reciprocal_and_analytic_in_kx0_and_f():
    sub = _sub(0.0022)
    _, basis = L.build_supercell(L.s1_angles(3), 6, sub)
    o = (20, 6)
    f, kx0 = 10.0e9, 50.0 - 8.0j
    z = L.assemble_z_complex(basis, sub, f, kx0, 0.0, *o)
    zm = L.assemble_z_complex(basis, sub, f, -kx0, 0.0, *o)
    assert (
        np.abs(z.T - zm).max() / np.abs(z).max() < 1e-12
    )  # reciprocity: Z(k)^T = Z(-k), also at complex k (Z is symmetric only at k = 0)
    z0 = L.assemble_z_complex(basis, sub, f, 0.0, 0.0, *o)
    assert np.abs(z0 - z0.T).max() / np.abs(z0).max() < 1e-12
    # Cauchy-Riemann: dZ/d(Re p) = -j dZ/d(Im p) for p = kx0 and for p = f
    d = 1e-3
    zr = (
        L.assemble_z_complex(basis, sub, f, kx0 + d, 0.0, *o) - L.assemble_z_complex(basis, sub, f, kx0 - d, 0.0, *o)
    ) / (2 * d)
    zi = (
        L.assemble_z_complex(basis, sub, f, kx0 + 1j * d, 0.0, *o)
        - L.assemble_z_complex(basis, sub, f, kx0 - 1j * d, 0.0, *o)
    ) / (2 * d)
    assert np.abs(zr + 1j * zi).max() / np.abs(zr).max() < 1e-6
    df = 1e5
    fr = (
        L.assemble_z_complex(basis, sub, f + df, kx0, 0.0, *o) - L.assemble_z_complex(basis, sub, f - df, kx0, 0.0, *o)
    ) / (2 * df)
    fi = (
        L.assemble_z_complex(basis, sub, f + 1j * df, kx0, 0.0, *o)
        - L.assemble_z_complex(basis, sub, f - 1j * df, kx0, 0.0, *o)
    ) / (2 * df)
    assert np.abs(fr + 1j * fi).max() / np.abs(fr).max() < 1e-6


def test_real_k_is_unchanged_by_the_outgoing_branch_choice():
    # for real kt and a real frequency the outgoing/decaying choice coincides with the production Im(kz) <= 0 branch
    from src.periodic_bc import normal_wavenumber

    k = 300.0
    kt = np.array([100.0, 500.0, 299.0, 301.0])
    ours = L._kz_air(k * k - kt**2)
    theirs = normal_wavenumber(k, kt, 0.0)
    assert np.allclose(ours, theirs, rtol=1e-12, atol=0)


def test_supercell_of_identical_cells_has_the_single_cell_mode_and_does_not_leak():
    sub = _sub(0.0)
    kx0 = math.radians(52.72) / 3.2e-3
    ang0 = (0.0,)
    _, b1 = L.build_supercell(ang0, 12, sub, **GEOM)
    f1 = _peak_mode(b1, sub, kx0, 8.0, 14.0, *L.orders_for(b1, 1, 2))
    _, b3 = L.build_supercell(ang0 * 3, 12, sub, **GEOM)
    o3 = L.orders_for(b3, 3, 2)

    def zf(p):
        return L.assemble_z_complex(b3, sub, p * 1e9, kx0, 0.0, *o3)

    f3, _, _ = L.find_mode_eig(zf, f1.real + 0.003, f1.real - 0.003)
    assert (
        abs(f3 - f1) / f1.real < 1e-3
    )  # same Bloch wave, three copies of the cell (Floquet truncation differs slightly)
    assert abs(f3.imag) / f3.real < 1e-6  # lossless, bound (no radiating harmonic): real frequency


def test_leakage_needs_the_modulation_and_loss_only_adds_attenuation():
    # N = 5 rotating slots: 2*pi/L = 393 rad/m puts the (-1) harmonic of a beta ~ 280 wave inside the light cone (radiates)
    n, m = 5, 12
    lossless, lossy = _sub(0.0), _sub(0.01)
    kx0 = 0.3
    res = {}
    for name, angles in (("rot", L.s1_angles(n)), ("same", (0.0,) * n)):
        _, b = L.build_supercell(angles, m, lossless, **GEOM)
        o = L.orders_for(b, n, 2)
        f0 = _peak_mode(b, lossless, kx0, 9.0, 14.0, *o, n=51)
        for tag, sub in (("0", lossless), ("d", lossy)):
            _, bb = L.build_supercell(angles, m, sub, **GEOM)

            def zf(p, bb=bb, sub=sub, o=o):
                return L.assemble_z_complex(bb, sub, p * 1e9, kx0, 0.0, *o)

            p, _, _ = L.find_mode_eig(zf, f0.real + 0.002 + 0.001j, f0.real - 0.002 + 0.0005j)
            res[name + tag] = p
    q = {k: v.real / (2.0 * v.imag) for k, v in res.items()}
    assert q["rot0"] > 0 and q["rot0"] < 1e5  # a finite radiative Q
    assert (
        q["same0"] > 1e3 * q["rot0"] or res["same0"].imag < 1e-9 * res["same0"].real
    )  # no modulation -> (numerically) no leakage
    # tan(delta) adds attenuation: 1/Q_tot > 1/Q_rad, and the extra loss is (approximately) additive
    assert 1.0 / q["rotd"] > 1.0 / q["rot0"]
    assert 1.0 / q["samed"] > 1.0 / q["same0"] + 1e-6
