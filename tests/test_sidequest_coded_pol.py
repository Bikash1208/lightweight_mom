"""Side quest round 8: polarisation-switched K = 2 coded hologram (rotated 1-diode cells vs 2 diodes per cell), model level."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "research"))
import sidequest_coded_pol as P  # noqa: E402


def test_cell_constants_are_physical():
    assert 1.0 < P.BETA_K0 < 1.1 and 0.2 < P.A0 < 2.0 and P.AABS0 >= 0.0
    assert abs(P.DELTA) * 1e3 > 3.0  # mS: the rotated disc cell is a weaker modulator than the x-only split patch (7.6 mS)


def test_zero_penalty_design_picks_the_larger_contribution_and_is_deterministic():
    B = P.PolBoard(P.CHIS, 0.8)
    d1, d2 = P.design(B, "R1", [0.0, 0.0], 0.0), P.design(B, "R1", [0.0, 0.0], 0.0)
    assert all(np.array_equal(d1[k], d2[k]) for k in d1)
    r = np.real(B.c)  # (2, nchi, n) unit-sign contributions at eps = 0
    ii = np.arange(B.n)
    chosen = d1["sa"][0] * r[0, d1["idx"][0], ii] + d1["sb"][0] * r[1, d1["idx"][0], ii]
    assert np.all(chosen >= (np.abs(r[0]) + np.abs(r[1])).max(axis=0) - 1e-12)


def test_penalty_trades_diodes_for_static_cells():
    B = P.PolBoard(P.CHIS, 0.8)
    fr = [P.summarize(P.design(B, "R1", [0.0, 0.0], lam))["static_frac"] for lam in (0.0, 0.5, 1.2)]
    assert fr[0] < fr[1] < fr[2]


def test_two_diode_route_uses_x_and_y_elements_only():
    B = P.PolBoard(np.radians([0.0, 90.0]), 1.0)
    d = P.design(B, "R2", [0.0, 0.0], 0.0)
    assert d["idx"].shape == (2, B.n) and set(d["idx"].ravel()) == {0, 1}
    s = P.summarize(d)
    assert s["diodes"] <= 2 * B.n and s["lines_with_inverters"] <= 1 and s["lines_no_inverters"] <= 2
