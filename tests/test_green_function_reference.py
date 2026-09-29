"""Green's function vs. an independent full-Maxwell plane-wave solver (reference_layered.py)."""

import numpy as np
import pytest

from reference_layered import solve_current_sheet
from src.green_function import LayeredGreenFunction
from src.materials import Substrate

F = 10e9
H = 1.575e-3

SUBSTRATES = {
    "lossless_RT5880": (2.2, 0.0),
    "lossy_RT5880": (2.2, 0.0009),
    "lossy_high_eps": (10.2, 0.0023),
    "very_lossy": (4.4, 0.05),
}
CURRENTS = [(1.0, 0.0), (0.0, 1.0), (0.3 + 0.2j, -0.7j)]


def make(name, thickness=H):
    eps_r, tand = SUBSTRATES[name]
    sub = Substrate(name, eps_r, tand, thickness)
    return sub, LayeredGreenFunction(sub, F)


# kt/k0 chosen to cover: kt=0, radiating (<1), between air and substrate cutoff
# (1 < kt/k0 < sqrt(eps_r)), and strongly evanescent everywhere.
KT_OVER_K0 = [0.0, 0.2, 0.7, 0.99, 1.1, 1.4, 1.6, 2.0, 5.0, 20.0]
ANGLES = [0.0, 0.4, 1.3, -2.0]


@pytest.mark.parametrize("name", SUBSTRATES)
def test_dyad_matches_independent_field_matching(name):
    sub, g = make(name)
    for kt_r in KT_OVER_K0:
        for ang in ANGLES:
            kx, ky = kt_r * g.k0 * np.cos(ang), kt_r * g.k0 * np.sin(ang)
            G = g.spectral_dyad(kx, ky)
            for J in CURRENTS:
                ref = solve_current_sheet(kx, ky, J, F, sub.eps_complex, sub.thickness)["Et"]
                mine = G @ np.array(J)
                assert mine == pytest.approx(ref, rel=1e-9, abs=1e-9 * np.linalg.norm(ref))


@pytest.mark.parametrize("thickness", [1e-5, 0.4e-3, 1.575e-3, 6e-3, 0.03])
@pytest.mark.parametrize("name", ["lossless_RT5880", "lossy_RT5880"])
def test_dyad_matches_reference_for_various_thickness(name, thickness):
    sub, g = make(name, thickness)
    for kt_r in (0.0, 0.5, 1.3, 3.0):
        kx, ky = kt_r * g.k0 * 0.6, kt_r * g.k0 * 0.8
        ref = solve_current_sheet(kx, ky, (1.0, 0.5j), F, sub.eps_complex, thickness)["Et"]
        assert g.spectral_dyad(kx, ky) @ np.array([1.0, 0.5j]) == pytest.approx(ref, rel=1e-8)


# ------------------------------------------------------- energy conservation

@pytest.mark.parametrize("name", SUBSTRATES)
def test_energy_balance_source_equals_radiated_plus_dissipated(name):
    """P_source = -1/2 Re(J* . E_t)  must equal  P_up + P_down + P_dissipated."""
    sub, g = make(name)
    for kt_r in KT_OVER_K0:
        for J in CURRENTS:
            kx, ky = kt_r * g.k0 * 0.6, kt_r * g.k0 * 0.8
            r = solve_current_sheet(kx, ky, J, F, sub.eps_complex, H)
            # power from the *closed-form* dyad, not from the reference solver's own E_t
            p_from_dyad = -0.5 * np.real(np.conj(J) @ (g.spectral_dyad(kx, ky) @ np.array(J)))
            assert p_from_dyad == pytest.approx(r["P_source"], rel=1e-8, abs=1e-12)
            total = r["P_up"] + r["P_down"] + r["P_diss"]
            assert total == pytest.approx(r["P_source"], rel=1e-5, abs=1e-9)


@pytest.mark.parametrize("name", ["lossless_RT5880", "lossy_RT5880"])
def test_radiation_only_for_kt_below_k0_and_outgoing(name):
    """Air-side power flow is outward (>0 up and down) for kt < k0 and zero above."""
    sub, g = make(name)
    for kt_r in (0.1, 0.5, 0.95):
        r = solve_current_sheet(kt_r * g.k0, 0.0, (0.0, 1.0), F, sub.eps_complex, H)
        assert r["P_up"] > 0.0 and r["P_down"] > 0.0
    for kt_r in (1.2, 3.0):
        r = solve_current_sheet(kt_r * g.k0, 0.0, (0.0, 1.0), F, sub.eps_complex, H)
        assert abs(r["P_up"]) < 1e-12 and abs(r["P_down"]) < 1e-12


def test_lossless_has_no_dissipation_and_lossy_has():
    sub, g = make("lossless_RT5880")
    r = solve_current_sheet(0.3 * g.k0, 0.0, (1.0, 0.0), F, sub.eps_complex, H)
    assert r["P_diss"] == 0.0
    assert r["P_source"] == pytest.approx(r["P_up"] + r["P_down"], rel=1e-9)
    sub, g = make("lossy_RT5880")
    r = solve_current_sheet(0.3 * g.k0, 0.0, (1.0, 0.0), F, sub.eps_complex, H)
    assert r["P_diss"] > 0.0


def test_reference_solver_detects_wrong_kz_branch():
    """Sanity check on the test itself: with kz -> -kz the radiated power is not outgoing."""
    sub, g = make("lossless_RT5880")
    bad = solve_current_sheet(0.3 * g.k0, 0.0, (1.0, 0.0), F, sub.eps_complex, H, branch_sign=-1.0)
    assert bad["P_up"] < 0.0 or bad["P_down"] < 0.0
    assert bad["Et"] != pytest.approx(g.spectral_dyad(0.3 * g.k0, 0.0) @ np.array([1.0, 0.0]), rel=1e-3)
