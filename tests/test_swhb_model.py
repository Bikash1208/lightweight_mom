"""Physics checks of the SWHB board model and the mushroom-branch cell model (sidequest Route B).  Fast, no MoM solves."""
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "research"))
import sidequest_cellB_model as M  # noqa: E402
import sidequest_swhb as S  # noqa: E402

ENV = S.Env(3.38, 0.0027, 1.524e-3)


def test_single_cell_far_field_power_matches_element_q():
    """The incoherent element power Q |J~|^2 equals the hemisphere integral of the coherent far field of one cell (same Green function)."""
    ysh = 3e-3j
    th = (np.arange(180) + 0.5) * (math.pi / 2) / 180
    ph = (np.arange(72) + 0.5) * 2 * math.pi / 72
    xy = np.array([[0.0, 0.0]])
    j = np.array([[1.0, 0.0]])  # J along x, 1 A/m, cell 1 mm so that sinc ~ 1
    p = 1e-3
    eth, eph = S.farfield(ENV, xy, j, ysh, ysh, th, ph, p)
    prad, _ = S.radiated_power(eth, eph, th, ph)
    q = S.element_q(ENV, ysh, ysh, nth=180, nph=72)
    jt2 = (1.0 * p * p) ** 2
    assert abs(prad / (q * jt2) - 1.0) < 0.03


def test_unmodulated_lossless_sheet_does_not_radiate():
    """Equal cells: no modulation current, the SW power is only reduced by the dielectric loss; no radiated power."""
    xy = np.array([[x, y] for x in np.arange(-0.05, 0.0501, 6.4e-3) for y in np.arange(-0.05, 0.0501, 6.4e-3) if 0.01 <= math.hypot(x, y) <= 0.05])
    y = np.full(len(xy), -3e-3j)
    r = S.run_state(ENV, xy, y, lambda kr: 0j, (0j, 0j), 6.4e-3, 1.0)
    assert r["ok"]
    assert r["P_rad_unit"] < 1e-12
    expected_edge = math.exp(-2 * r["alpha_loss"] * (0.05 - 0.01))
    assert abs(r["P_edge"] - expected_edge) < 0.1


def test_power_closure_modulated_board():
    """Radiated + edge + loss accounts for the launched power within the march / kappa iteration tolerance."""
    rng = np.random.default_rng(1)
    xy = np.array([[x, y] for x in np.arange(-0.06, 0.0601, 6.4e-3) for y in np.arange(-0.06, 0.0601, 6.4e-3) if 0.01 <= math.hypot(x, y) <= 0.06])
    y = np.where(rng.random(len(xy)) < 0.5, 3e-3j, -4e-3j)
    r = S.run_state(ENV, xy, y, lambda kr: 0j, (0j, 0j), 6.4e-3, 1.0)
    assert r["ok"]
    loss = 1.0 - r["P_rad_unit"] - r["P_edge"]
    assert -0.05 < loss < 0.4  # nonnegative within tolerance, bounded by the loss estimate
    assert r["P_rad_unit"] > 1e-3  # a random 1-bit pattern on a weakly bound mode radiates little but not zero


def test_branch_limits():
    """ON branch -> 1/(j w mu0 h) for R -> 0, OFF branch -> j w C for large series reactance; contrast grows when L_ser shrinks."""
    h = 1.524e-3
    yon = M.branch("ON", h, "MA4AGP907_measured", dR=1e-6)
    assert abs(yon - 1 / (1j * M.W0 * M.MU0 * h)) < 1e-9
    yoff = M.branch("OFF", h, "MA4AGP907_measured")
    c = 28.16e-15
    assert abs(yoff.imag - M.W0 * c / (1 - M.W0 ** 2 * M.MU0 * h * c)) < 1e-6
    d0 = abs(M.branch("ON", h, "MA4AGP907_measured", L_ser=0) - M.branch("OFF", h, "MA4AGP907_measured", L_ser=0))
    d1 = abs(M.branch("ON", h, "MA4AGP907_measured", L_ser=2e-9) - M.branch("OFF", h, "MA4AGP907_measured", L_ser=2e-9))
    assert d0 > d1


def test_tm0_bare_slab_beta_matches_transverse_resonance():
    """Bare grounded slab (zero sheet): the root of the transverse resonance sits between the light line and the slab line."""
    r = ENV.beta_alpha(lambda kr: 0j)
    assert r is not None and 1.0 < r[0] < math.sqrt(3.38)
